"use strict";
/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
var toscad;
(function (toscad) {
    // easy to use constant for later reference.
    const resolveLatticeKind = (lattice) => (lattice ?? '').toLowerCase() === 'square' ? 'square' : 'honeycomb';
    const LATTICE_CONFIG = {
        honeycomb: {
            basesPerTurn: 10.5,
            phases: {
                0: [0.0, 3.5, 7.0, 10.5],
                // 1: [1.75, 5.25, 8.75]
                // 1: [2.77, 6.27, 9.77]
                // numbers found by using 215/145 asymmetry in major/minor groove and allowing void phase to break ties.
                1: [0.7292, 4.229, 7.7292]
            },
            voidPhase: 5.25,
            tieEpsilon: 0.01
        },
        square: {
            basesPerTurn: 32 / 3,
            phases: {
                // Symmetric and therefore simple.
                0: [0.0, 8 / 3, 16 / 3, 8.0, 32 / 3],
                1: [4 / 3, 4.0, 20 / 3, 28 / 3]
            },
            voidPhase: null,
            tieEpsilon: 0.01
        }
    };
    function findHelixID(targetId, helices) {
        for (let i = 0; i < helices.length; i++) {
            const helix = helices[i];
            for (const nt of helix) {
                if (nt.id === targetId)
                    return i;
            }
        }
        return null;
    }
    toscad.findHelixID = findHelixID;
    // helper function
    // Finds the local angle distribution of a given helix relative to its neighbors.
    // Uses crossover offsets defined by LATTICE_CONFIG to determine ideal angle buckets. 
    function getAngleHelix(grid, helices, helixId, lattice) {
        void helices;
        const result = new Map();
        const latticeType = resolveLatticeKind(lattice);
        const latticeConfig = LATTICE_CONFIG[latticeType];
        const nearestPhase = (value, parity) => {
            const candidates = latticeConfig.phases[parity] ?? [];
            if (!candidates.length)
                return value;
            return candidates.reduce((prev, curr) => {
                const distPrev = Math.abs(prev - value);
                const distCurr = Math.abs(curr - value);
                if (parity === 1
                    && latticeConfig.voidPhase !== null
                    && Math.abs(distPrev - distCurr) < latticeConfig.tieEpsilon) {
                    if (prev === latticeConfig.voidPhase && curr !== latticeConfig.voidPhase)
                        return curr;
                    if (curr === latticeConfig.voidPhase && prev !== latticeConfig.voidPhase)
                        return prev;
                }
                return distCurr < distPrev ? curr : prev;
            }, candidates[0]);
        };
        const hubCrossovers = [];
        for (const crossover of crossoverNts(grid)) {
            if (crossover.fromHelix === helixId) {
                const hubMark = grid.get(crossover.fromNt.id);
                if (!hubMark)
                    continue;
                hubCrossovers.push({
                    adj_helix: crossover.toHelix,
                    offset: crossover.fromOffset,
                    direction: hubMark.direction
                });
            }
            else if (crossover.toHelix === helixId) {
                const hubMark = grid.get(crossover.toNt.id);
                if (!hubMark)
                    continue;
                hubCrossovers.push({
                    adj_helix: crossover.fromHelix,
                    offset: crossover.toOffset,
                    direction: hubMark.direction
                });
            }
        }
        if (!hubCrossovers.length)
            return result;
        const groupedByNeighbor = new Map();
        for (const crossover of hubCrossovers) {
            if (!groupedByNeighbor.has(crossover.adj_helix)) {
                groupedByNeighbor.set(crossover.adj_helix, []);
            }
            groupedByNeighbor.get(crossover.adj_helix).push(crossover);
        }
        const neighbors = Array.from(groupedByNeighbor.keys()).sort((a, b) => a - b);
        if (!neighbors.length)
            return result;
        const pairTallies = new Map();
        for (let i = 0; i < neighbors.length; i++) {
            const neighborA = neighbors[i];
            const groupA = groupedByNeighbor.get(neighborA) ?? [];
            for (let j = i + 1; j < neighbors.length; j++) {
                const neighborB = neighbors[j];
                const groupB = groupedByNeighbor.get(neighborB) ?? [];
                const pairKey = `${neighborA}|${neighborB}`;
                const bucket = new Map();
                for (const crossoverA of groupA) {
                    for (const crossoverB of groupB) {
                        const offsetA = crossoverA.offset;
                        const dirA = crossoverA.direction;
                        const offsetB = crossoverB.offset;
                        const dirB = crossoverB.direction;
                        // Here's the math i spent so long figuring out:
                        const rawX = offsetB - offsetA;
                        const y = dirA === dirB ? 0 : 1;
                        const phase = ((rawX % latticeConfig.basesPerTurn) + latticeConfig.basesPerTurn) % latticeConfig.basesPerTurn;
                        const idealPhase = nearestPhase(phase, y);
                        const angleRaw = (360 / latticeConfig.basesPerTurn) * idealPhase + (y * 215);
                        const relativeAngle = (Math.round(angleRaw % 360) + 360) % 360;
                        bucket.set(relativeAngle, (bucket.get(relativeAngle) ?? 0) + 1);
                    }
                }
                pairTallies.set(pairKey, bucket);
            }
        }
        const pairConsensus = new Map();
        for (const [pairKey, bucket] of pairTallies.entries()) {
            if (!bucket.size)
                continue;
            let modeAngle = 0;
            let modeCount = -1;
            for (const [angle, count] of bucket.entries()) {
                if (count > modeCount || (count === modeCount && angle < modeAngle)) {
                    modeAngle = angle;
                    modeCount = count;
                }
            }
            pairConsensus.set(pairKey, modeAngle);
        }
        const baseReference = neighbors[0];
        for (const neighbor of neighbors) {
            const angle = neighbor === baseReference
                ? 0
                : (pairConsensus.get(`${baseReference}|${neighbor}`) ?? 0);
            result.set(neighbor, {
                helixId,
                adj_helix: neighbor,
                angle
            });
        }
        return result;
    }
    // Helper function to collect all backbone crossovers with their helix and offset info.
    function crossoverNts(grid) {
        const allNtIds = new Set();
        for (const [ntId] of grid.entries())
            allNtIds.add(ntId);
        const visited = new Set();
        const crossovers = [];
        for (const [ntId] of grid.entries()) {
            if (visited.has(ntId))
                continue;
            const startNt = elements.get(ntId);
            if (!startNt || !(startNt instanceof Nucleotide))
                continue;
            // Find 5' end
            let fivePrime = startNt;
            const walkBack = new Set();
            walkBack.add(fivePrime.id);
            while (true) {
                const prev = fivePrime.n5;
                if (!prev || !(prev instanceof Nucleotide))
                    break;
                if (!allNtIds.has(prev.id))
                    break;
                if (walkBack.has(prev.id))
                    break;
                walkBack.add(prev.id);
                fivePrime = prev;
            }
            // Walk 5' -> 3' and record backbone helix transitions
            let curr = fivePrime;
            const walkForward = new Set();
            let prevNt = null;
            let prevMark = null;
            while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                if (walkForward.has(curr.id))
                    break;
                walkForward.add(curr.id);
                visited.add(curr.id);
                const mark = grid.get(curr.id);
                if (mark) {
                    if (prevNt && prevMark && prevMark.helixId !== mark.helixId) {
                        crossovers.push({
                            fromHelix: prevMark.helixId,
                            toHelix: mark.helixId,
                            fromOffset: prevMark.offset,
                            toOffset: mark.offset,
                            fromNt: prevNt,
                            toNt: curr
                        });
                    }
                    prevNt = curr;
                    prevMark = mark;
                }
                else {
                    prevNt = null;
                    prevMark = null;
                }
                const n3ref = curr.n3;
                curr = (n3ref && n3ref instanceof Nucleotide) ? n3ref : null;
            }
        }
        return crossovers;
    }
    toscad.crossoverNts = crossoverNts;
    // Auto-detect whether a structure was built on a honeycomb or square lattice.
    //
    // Per-helix, sort crossover-endpoint nucleotides by offset, then for every
    // consecutive pair compute the gap and test:
    //     fitsHC = (gap mod 28) ∈ {6, 13, 20, 27}
    //     fitsSQ = (gap mod 32) ∈ {7, 15, 23, 31}
    // Exactly one match → +1 vote for that lattice. Matches both → ambiguous
    // (skip). Matches neither → reject. Whichever lattice has more votes wins;
    // ties (including zero votes) fall back to honeycomb.
    //
    // Binder helices (passed from setGrid) attach outside the main lattice, so
    // any crossover with either endpoint on a binder helix is excluded.
    function detectLatticeKind(grid, binderHelices = []) {
        const binderSet = new Set(binderHelices);
        // helixId -> ntId -> point. Map-by-ntId dedupes nucleotides that appear
        // as both the "to" of one crossover and the "from" of another.
        const pointsByHelix = new Map();
        const addPoint = (helixId, ntId, offset, direction) => {
            let inner = pointsByHelix.get(helixId);
            if (!inner) {
                inner = new Map();
                pointsByHelix.set(helixId, inner);
            }
            inner.set(ntId, { offset, direction });
        };
        let skippedBinderCrossovers = 0;
        for (const crossover of crossoverNts(grid)) {
            // Exclude crossovers that touch a binder helix on either side.
            if (binderSet.has(crossover.fromHelix) || binderSet.has(crossover.toHelix)) {
                skippedBinderCrossovers++;
                continue;
            }
            const fromMark = grid.get(crossover.fromNt.id);
            const toMark = grid.get(crossover.toNt.id);
            if (fromMark) {
                addPoint(crossover.fromHelix, crossover.fromNt.id, crossover.fromOffset, fromMark.direction);
            }
            if (toMark) {
                addPoint(crossover.toHelix, crossover.toNt.id, crossover.toOffset, toMark.direction);
            }
        }
        const HC_PERIOD = 28;
        const HC_OFFSETS = {
            0: new Set([6, 13, 20, 27]),
            1: new Set([1, 4, 8])
        };
        const SQ_PERIOD = 32;
        const SQ_OFFSETS = {
            0: new Set([7, 15, 23, 31]),
            1: new Set([3, 11, 19, 27])
        };
        let honeycombVotes = 0;
        let squareVotes = 0;
        let totalGaps = 0;
        // Per-gap-size breakdown of *accepted* gaps only. Map<gap, {hc, sq}>.
        const acceptedByGap = new Map();
        const bumpAccepted = (gap, kind) => {
            let entry = acceptedByGap.get(gap);
            if (!entry) {
                entry = { hc: 0, sq: 0 };
                acceptedByGap.set(gap, entry);
            }
            entry[kind]++;
        };
        console.log(`[detectLatticeKind] starting; helices with crossover endpoints=${pointsByHelix.size}, ` +
            `binder helices=${binderHelices.length} (${binderHelices.join(',') || 'none'}), ` +
            `binder-touching crossovers skipped=${skippedBinderCrossovers}`);
        for (const [helixId, inner] of pointsByHelix.entries()) {
            if (inner.size < 2) {
                console.log(`[detectLatticeKind]   helix=${helixId} has ${inner.size} crossover endpoint(s) — skipping`);
                continue;
            }
            const sorted = Array.from(inner.values()).sort((a, b) => a.offset - b.offset);
            const summary = sorted.map(p => `${p.offset}${p.direction === 'forward' ? 'f' : 'b'}`).join(', ');
            console.log(`[detectLatticeKind]   helix=${helixId} endpoints (offset+dir, sorted): ${summary}`);
            for (let i = 1; i < sorted.length; i++) {
                const a = sorted[i - 1];
                const b = sorted[i];
                const rawX = b.offset - a.offset;
                if (rawX === 0)
                    continue;
                const parity = a.direction === b.direction ? 0 : 1;
                // TODO: gap=1 skipped as a proxy for "same crossover junction"; add a real same-crossover grouping fn.
                // For parity=1 rawX=1 is a valid signal (parity-1 honeycomb list includes 1), so only skip for parity=0.
                if (parity === 0 && rawX === 1)
                    continue;
                totalGaps++;
                const hcResid = ((rawX % HC_PERIOD) + HC_PERIOD) % HC_PERIOD;
                const sqResid = ((rawX % SQ_PERIOD) + SQ_PERIOD) % SQ_PERIOD;
                const fitsHC = HC_OFFSETS[parity].has(hcResid);
                const fitsSQ = SQ_OFFSETS[parity].has(sqResid);
                let verdict;
                if (fitsHC && fitsSQ) {
                    verdict = 'ambiguous';
                }
                else if (fitsHC) {
                    verdict = 'HC';
                    honeycombVotes++;
                    bumpAccepted(rawX, 'hc');
                }
                else if (fitsSQ) {
                    verdict = 'SQ';
                    squareVotes++;
                    bumpAccepted(rawX, 'sq');
                }
                else {
                    verdict = 'rejected';
                }
                console.log(`[detectLatticeKind]     helix=${helixId} gap=${rawX} (offsets ${a.offset}->${b.offset}) ` +
                    `parity=${parity} | gap%28=${hcResid} (HC ${fitsHC ? 'yes' : 'no'}) ` +
                    `gap%32=${sqResid} (SQ ${fitsSQ ? 'yes' : 'no'}) -> ${verdict}`);
            }
        }
        console.log(`[detectLatticeKind] total gaps inspected (post gap=0,1 skip): ${totalGaps}`);
        const acceptedSizes = Array.from(acceptedByGap.keys()).sort((a, b) => a - b);
        if (acceptedSizes.length === 0) {
            console.log(`[detectLatticeKind] accepted gap-size breakdown: none`);
        }
        else {
            console.log(`[detectLatticeKind] accepted gap-size breakdown:`);
            for (const size of acceptedSizes) {
                const { hc, sq } = acceptedByGap.get(size);
                const parts = [];
                if (hc)
                    parts.push(`HC=${hc}`);
                if (sq)
                    parts.push(`SQ=${sq}`);
                console.log(`[detectLatticeKind]   gap=${size}: ${parts.join(', ')}`);
            }
        }
        const detected = squareVotes > honeycombVotes ? 'square' : 'honeycomb';
        console.log(`[detectLatticeKind] totals: HC=${honeycombVotes}, SQ=${squareVotes} -> ${detected}`);
        return detected;
    }
    toscad.detectLatticeKind = detectLatticeKind;
    // Now run getAngleHelix for every helix to get a full network map.
    // Note that every helix is has RELATIVE angles to its neighbors, they are not globally aligned to anything yet.
    // The global alignment is done in calculateGridPositions. 
    function getAngles(grid, helices, lattice = 'honeycomb') {
        const networkMap = new Map();
        const helixIds = new Set();
        for (const [, mark] of grid.entries()) {
            helixIds.add(mark.helixId);
        }
        const sortedhids = Array.from(helixIds).sort((a, b) => a - b);
        for (const currentHID of sortedhids) {
            const helixAngles = getAngleHelix(grid, helices, currentHID, lattice);
            const angleMap = new Map();
            for (const [adjHelixId, angleInfo] of helixAngles.entries()) {
                angleMap.set(adjHelixId, angleInfo.angle);
            }
            networkMap.set(currentHID, angleMap);
        }
        return networkMap;
    }
    toscad.getAngles = getAngles;
    function tempGlobalPos(networkMap, grid, lattice = 'honeycomb') {
        const latticeType = resolveLatticeKind(lattice);
        const ANGLES = latticeType === 'square' ? [0, 90, 180, 270] : [0, 120, 240];
        const HONEYCOMB_STEP_BY_PARITY = {
            even: {
                0: { dCol: 1, dRow: 0 },
                120: { dCol: -1, dRow: 0 },
                240: { dCol: 0, dRow: -1 }
            },
            odd: {
                0: { dCol: 1, dRow: 0 },
                120: { dCol: 0, dRow: 1 },
                240: { dCol: -1, dRow: 0 }
            }
        };
        const SQUARE_STEPS = {
            0: { dCol: 1, dRow: 0 },
            90: { dCol: 0, dRow: 1 },
            180: { dCol: -1, dRow: 0 },
            270: { dCol: 0, dRow: -1 }
        };
        // Helper functions
        const normalizeAngle = (angle) => ((angle % 360) + 360) % 360;
        // Returns the smallest angular distance between two angles (in degrees)
        const angleDistance = (a, b) => {
            const diff = Math.abs(normalizeAngle(a) - normalizeAngle(b));
            return Math.min(diff, 360 - diff);
        };
        // checks the parity, based on row
        const parityAt = (col, row) => (((col + row) & 1) === 0 ? 'even' : 'odd');
        const getStep = (col, row, angle) => {
            if (latticeType === 'square')
                return SQUARE_STEPS[angle];
            return HONEYCOMB_STEP_BY_PARITY[parityAt(col, row)][angle];
        };
        const snapToLatticeAngle = (angle) => {
            let best = ANGLES[0];
            let bestDist = Number.POSITIVE_INFINITY;
            for (const cand of ANGLES) {
                const d = angleDistance(angle, cand);
                if (d < bestDist || (d === bestDist && cand < best)) {
                    best = cand;
                    bestDist = d;
                }
            }
            return best;
        };
        const latticeAngleFromDelta = (col, row, dCol, dRow) => {
            for (const a of ANGLES) {
                const s = getStep(col, row, a);
                if (s.dCol === dCol && s.dRow === dRow)
                    return a;
            }
            return null;
        };
        // Get number of crossovers
        const connectionCounts = getConnectionCounts(grid);
        const getWeight = (a, b) => connectionCounts.get(a)?.get(b) ?? 1;
        // Collect every helix referenced by the network map
        const allHelixIds = new Set();
        allHelixIds.add(0);
        for (const [from, row] of networkMap.entries()) {
            allHelixIds.add(from);
            for (const to of row.keys())
                allHelixIds.add(to);
        }
        // Derive a helix's global orientation from its parent's.
        const deriveChildOrientation = (parentId, childId, parentCol, parentRow, childCol, childRow, parentLocalAngle, snappedForwardAngle) => {
            const dColBack = parentCol - childCol;
            const dRowBack = parentRow - childRow;
            const backGlobal = latticeAngleFromDelta(childCol, childRow, dColBack, dRowBack);
            const desiredBack = backGlobal !== null
                ? backGlobal
                : normalizeAngle(snappedForwardAngle + 180);
            const childLocalBack = networkMap.get(childId)?.get(parentId);
            if (typeof childLocalBack === 'number') {
                return normalizeAngle(desiredBack - childLocalBack);
            }
            const inferredBackLocal = normalizeAngle(parentLocalAngle + 180);
            return normalizeAngle(desiredBack - inferredBackLocal);
        };
        // Initialize these. In this case, "canonical" means BFS-derived. Thus, non-canonical mean that the helix is not where the MST placed it.
        const canonicalCoord = new Map();
        const canonicalOrient = new Map();
        // if there are 2 totally disconnected structures, they will be placed at (ROOT_SEPARATION,0). It is just a safegaurd.
        // The first root is always placed at (0,0)
        const ROOT_SEPARATION = 1000;
        const sortedHelixIds = [...allHelixIds].sort((a, b) => a - b);
        let nextRootCol = ROOT_SEPARATION;
        // BFS from every root
        for (const rootId of sortedHelixIds) {
            if (canonicalCoord.has(rootId))
                continue;
            const rootCol = rootId === 0 ? 0 : nextRootCol;
            canonicalCoord.set(rootId, [rootCol, 0]);
            canonicalOrient.set(rootId, 0);
            if (rootId !== 0)
                nextRootCol += ROOT_SEPARATION;
            const queue = [rootId];
            // BFS manual instead of Array.shift() for performance
            let qi = 0;
            while (qi < queue.length) {
                const parentId = queue[qi++];
                const parentCoord = canonicalCoord.get(parentId);
                const parentOrient = canonicalOrient.get(parentId);
                const [pCol, pRow] = parentCoord;
                const edges = networkMap.get(parentId);
                if (!edges || edges.size === 0)
                    continue;
                // Sort candidates by crossover weight, then by ascending id.
                const candidates = [...edges.entries()]
                    .filter(([nid]) => nid !== parentId && !canonicalCoord.has(nid))
                    .map(([nid, localAngle]) => ({
                    nid,
                    localAngle: normalizeAngle(localAngle),
                    weight: getWeight(parentId, nid)
                }))
                    .sort((a, b) => b.weight - a.weight || a.nid - b.nid);
                for (const c of candidates) {
                    if (canonicalCoord.has(c.nid))
                        continue;
                    const predictedGlobal = normalizeAngle(c.localAngle + parentOrient);
                    const snapped = snapToLatticeAngle(predictedGlobal);
                    const step = getStep(pCol, pRow, snapped);
                    const childCol = pCol + step.dCol;
                    const childRow = pRow + step.dRow;
                    canonicalCoord.set(c.nid, [childCol, childRow]);
                    canonicalOrient.set(c.nid, deriveChildOrientation(parentId, c.nid, pCol, pRow, childCol, childRow, c.localAngle, snapped));
                    queue.push(c.nid);
                }
            }
        }
        // record every edge's prediction for the far endpoint
        // So far, we only have the tree edges (i.e., canonical edges). But sometimes, there are non-tree edges, caused by 2 helices giving different angles for the same helix.
        const predictions = new Map();
        const push = (childId, pred) => {
            let arr = predictions.get(childId);
            if (!arr) {
                arr = [];
                predictions.set(childId, arr);
            }
            arr.push(pred);
        };
        for (const [parentId, edges] of networkMap.entries()) {
            const pCoord = canonicalCoord.get(parentId);
            const pOrient = canonicalOrient.get(parentId);
            if (!pCoord || pOrient === undefined)
                continue;
            const [pCol, pRow] = pCoord;
            for (const [childId, localAngle] of edges.entries()) {
                if (childId === parentId)
                    continue;
                const cCoord = canonicalCoord.get(childId);
                const cOrient = canonicalOrient.get(childId);
                if (!cCoord || cOrient === undefined)
                    continue;
                const localAngleN = normalizeAngle(localAngle);
                const predictedGlobal = normalizeAngle(localAngleN + pOrient);
                const snapped = snapToLatticeAngle(predictedGlobal);
                const step = getStep(pCol, pRow, snapped);
                const predCol = pCol + step.dCol;
                const predRow = pRow + step.dRow;
                const [cCol, cRow] = cCoord;
                push(childId, {
                    coord: [predCol, predRow],
                    globalOrientation: cOrient,
                    viaParent: parentId,
                    parentCoord: [pCol, pRow],
                    parentOrientation: pOrient,
                    edgeLocalAngle: localAngleN,
                    edgeWeight: getWeight(parentId, childId),
                    isCanonical: predCol === cCol && predRow === cRow
                });
            }
        }
        // Synthetic root predictions for helices with no incoming edge (backup)
        for (const helixId of allHelixIds) {
            if (predictions.has(helixId))
                continue;
            const coord = canonicalCoord.get(helixId);
            const orient = canonicalOrient.get(helixId);
            if (!coord || orient === undefined)
                continue;
            push(helixId, {
                coord: [coord[0], coord[1]],
                globalOrientation: orient,
                viaParent: -1,
                parentCoord: [coord[0], coord[1]],
                parentOrientation: 0,
                edgeLocalAngle: 0,
                edgeWeight: 0,
                isCanonical: true
            });
        }
        return predictions;
    }
    toscad.tempGlobalPos = tempGlobalPos;
    // For every helix, count how many backbone crossovers it shares with each neighbor. Useful to remove bad combinations (such as end-only, which happens when a helix is broken in 2 pieces)
    function getConnectionCounts(grid) {
        const counts = new Map();
        const bump = (a, b) => {
            if (!counts.has(a))
                counts.set(a, new Map());
            const inner = counts.get(a);
            inner.set(b, (inner.get(b) ?? 0) + 1);
        };
        for (const crossover of crossoverNts(grid)) {
            if (crossover.fromHelix === crossover.toHelix)
                continue;
            bump(crossover.fromHelix, crossover.toHelix);
            bump(crossover.toHelix, crossover.fromHelix);
        }
        return counts;
    }
    toscad.getConnectionCounts = getConnectionCounts;
    // Cosine similarity between the crossover-count vectors of two helices.
    //
    //   cos(A, B) = Σ_i (w_{A,i} · w_{B,i})  /  ( ||w_A|| · ||w_B|| )
    //
    // where w_{X,i} = crossovers(X, i) and i ranges over "neighbor helices"
    // (any helix ≠ A, ≠ B). A and B themselves are excluded so the score
    // reflects overlap in EXTERNAL connectivity — how similarly A and B are
    // wired to the rest of the graph — rather than being inflated by any
    // direct A↔B crossovers.
    //
    // Returns 0 if either vector has zero norm (i.e., the helix has no
    // strong external neighbors). Result is in [0, 1] since all weights
    // are non-negative.
    function helixPairCosine(A, B, connectionCounts) {
        const nbA = connectionCounts.get(A);
        const nbB = connectionCounts.get(B);
        if (!nbA || !nbB)
            return 0;
        // Neighbor universe = (A's neighbors ∪ B's neighbors) \ {A, B}.
        const neighbors = new Set();
        for (const n of nbA.keys())
            if (n !== A && n !== B)
                neighbors.add(n);
        for (const n of nbB.keys())
            if (n !== A && n !== B)
                neighbors.add(n);
        let dot = 0;
        let sumSqA = 0;
        let sumSqB = 0;
        for (const i of neighbors) {
            const wA = nbA.get(i) ?? 0;
            const wB = nbB.get(i) ?? 0;
            dot += wA * wB;
            sumSqA += wA * wA;
            sumSqB += wB * wB;
        }
        if (sumSqA === 0 || sumSqB === 0)
            return 0;
        const raw = dot / (Math.sqrt(sumSqA) * Math.sqrt(sumSqB));
        return Math.max(0, Math.min(1, raw));
    }
    toscad.helixPairCosine = helixPairCosine;
    // helper function to check for angle collisions in the map.
    function angleCollisions(networkMap) {
        const overlappingHelices = [];
        for (const [helixId, angleMap] of networkMap.entries()) {
            if (!angleMap || angleMap.size === 0)
                continue;
            const anglesToNeighbors = new Map();
            for (const [adjHelix, angle] of angleMap.entries()) {
                if (!anglesToNeighbors.has(angle))
                    anglesToNeighbors.set(angle, []);
                anglesToNeighbors.get(angle).push(adjHelix);
            }
            const localConflicts = [];
            for (const [angle, collidedHelices] of anglesToNeighbors.entries()) {
                const uniqueCollided = Array.from(new Set(collidedHelices)).sort((a, b) => a - b);
                if (uniqueCollided.length < 2)
                    continue;
                localConflicts.push({
                    angle,
                    colliding_adj_helices: uniqueCollided
                });
            }
            if (localConflicts.length > 0) {
                localConflicts.sort((a, b) => a.angle - b.angle);
                overlappingHelices.push({
                    helixId,
                    conflicts: localConflicts
                });
            }
        }
        overlappingHelices.sort((a, b) => a.helixId - b.helixId);
        return overlappingHelices;
    }
    toscad.angleCollisions = angleCollisions;
    // helper function. Checks if 2 helices are mutually disjoint (through offsets in the grid). 
    // typically used AFTER aligning the grid otherwise this is nonsense.
    function disjoint(helices, h1, h2, grid) {
        const buildSignedOffsets = (helixId) => {
            const signedOffsets = new Set();
            const nts = helices[helixId] ?? [];
            for (const nt of nts) {
                const mark = grid.get(nt.id);
                if (!mark || mark.helixId !== helixId)
                    continue;
                const x = mark.offset;
                const y = mark.direction === 'backward' ? 1 : 0;
                signedOffsets.add((2 * x) + y);
            }
            return signedOffsets;
        };
        const set1 = buildSignedOffsets(h1);
        const set2 = buildSignedOffsets(h2);
        if (set1.size === 0 || set2.size === 0)
            return false;
        const smaller = set1.size <= set2.size ? set1 : set2;
        const larger = set1.size <= set2.size ? set2 : set1;
        for (const offs of smaller) {
            if (larger.has(offs))
                return false;
        }
        return true;
    }
    toscad.disjoint = disjoint;
    // Variant of disjoint() that pretends each helix's grid marks have been shifted by shifts.get(id).
    // Doesn't mutate anything; used to preview whether a proposed reorder would collide with a
    // still-unshifted helix (so we know whether to expand the reorder set).
    function disjointWithShifts(helices, h1, h2, grid, shifts) {
        const build = (helixId) => {
            const s = new Set();
            const delta = shifts.get(helixId) ?? 0;
            const nts = helices[helixId] ?? [];
            for (const nt of nts) {
                const mark = grid.get(nt.id);
                if (!mark || mark.helixId !== helixId)
                    continue;
                const x = mark.offset + delta;
                const y = mark.direction === 'backward' ? 1 : 0;
                s.add(2 * x + y);
            }
            return s;
        };
        const s1 = build(h1);
        const s2 = build(h2);
        if (s1.size === 0 || s2.size === 0)
            return false;
        const smaller = s1.size <= s2.size ? s1 : s2;
        const larger = s1.size <= s2.size ? s2 : s1;
        for (const off of smaller) {
            if (larger.has(off))
                return false;
        }
        return true;
    }
    // Bitmask DP over orderings of the helices we need to reorder. For each subset S of the
    // reorder set, f[S] is the minimum external-crossover misalignment achievable by placing
    // exactly the helices in S end-to-end starting at offset 0, in some order. The offset each
    // helix ends up at is determined by the sum of spans of the helices before it, which is a
    // set-property (not an ordering property), so 2^N states with N transitions each captures
    // every permutation without redundancy.
    //
    // "External" here means "not currently in the reorder set" — includes helices in the combine
    // set that we chose to leave fixed AND helices completely outside the combine set. Their
    // grid offsets are treated as anchors we want to keep alignment with.
    //
    // Returns a shift-per-helix map, or null if the DP cannot run (missing footprint info).
    function runReorderDP(helices, reorderIds, grid) {
        const N = reorderIds.length;
        if (N < 2)
            return new Map();
        // Per-helix footprint on the unsigned offset axis.
        const minOff = new Array(N);
        const span = new Array(N);
        for (let i = 0; i < N; i++) {
            let mn = Infinity;
            let mx = -Infinity;
            const nts = helices[reorderIds[i]] ?? [];
            for (const nt of nts) {
                const mark = grid.get(nt.id);
                if (!mark || mark.helixId !== reorderIds[i])
                    continue;
                if (mark.offset < mn)
                    mn = mark.offset;
                if (mark.offset > mx)
                    mx = mark.offset;
            }
            if (mn === Infinity)
                return null;
            minOff[i] = mn;
            span[i] = mx - mn + 1;
        }
        const externalCX = reorderIds.map(() => []);
        const idToLocal = new Map();
        reorderIds.forEach((id, idx) => idToLocal.set(id, idx));
        for (const cx of crossoverNts(grid)) {
            const fromLocal = idToLocal.get(cx.fromHelix);
            const toLocal = idToLocal.get(cx.toHelix);
            // Both endpoints inside reorder set → internal, ignore.
            if (fromLocal !== undefined && toLocal !== undefined)
                continue;
            // Neither → irrelevant.
            if (fromLocal === undefined && toLocal === undefined)
                continue;
            if (fromLocal !== undefined) {
                externalCX[fromLocal].push({ offInternal: cx.fromOffset, offExternal: cx.toOffset });
            }
            if (toLocal !== undefined) {
                externalCX[toLocal].push({ offInternal: cx.toOffset, offExternal: cx.fromOffset });
            }
        }
        // For helix i placed at position p, an external crossover contributes:
        //   | offExternal - (offInternal + p - minOff[i]) |
        // = | (offExternal + minOff[i] - offInternal) - p |
        // Precompute targetP_k = offExternal + minOff[i] - offInternal per crossover.
        const targets = reorderIds.map((_, i) => externalCX[i].map(cx => cx.offExternal + minOff[i] - cx.offInternal));
        const costOf = (i, p) => {
            let total = 0;
            const arr = targets[i];
            for (let k = 0; k < arr.length; k++) {
                total += Math.abs(arr[k] - p);
            }
            return total;
        };
        // Sum of spans per subset. Grows by one bit at a time.
        const size = 1 << N;
        const sumSpan = new Int32Array(size);
        for (let S = 1; S < size; S++) {
            const low = S & -S;
            const li = Math.log2(low) | 0;
            sumSpan[S] = sumSpan[S ^ low] + span[li];
        }
        // f[S] = best score reaching subset S. parent[S] = the local helix index that was
        // placed last to achieve f[S]; used to recover the winning ordering.
        const f = new Float64Array(size);
        const parent = new Int32Array(size);
        f[0] = 0;
        parent[0] = -1;
        for (let S = 1; S < size; S++) {
            let best = Infinity;
            let bestBit = -1;
            let bits = S;
            while (bits !== 0) {
                const low = bits & -bits;
                const i = Math.log2(low) | 0;
                bits ^= low;
                const pred = S ^ (1 << i);
                const pos = sumSpan[pred];
                const c = f[pred] + costOf(i, pos);
                if (c < best) {
                    best = c;
                    bestBit = i;
                }
            }
            f[S] = best;
            parent[S] = bestBit;
        }
        // Walk parent[] backwards to reconstruct the ordering (last helix → first helix).
        const orderRev = [];
        let S = size - 1;
        while (S !== 0) {
            const i = parent[S];
            if (i < 0)
                break;
            orderRev.push(i);
            S = S ^ (1 << i);
        }
        const order = orderRev.reverse();
        // Assign final positions end-to-end starting at 0.
        const position = new Array(N).fill(0);
        let running = 0;
        for (const i of order) {
            position[i] = running;
            running += span[i];
        }
        // Global translation T: the whole block can slide along the axis. Given each helix's
        // shift Δ_i(T) = position[i] + T - minOff[i], each external crossover's misalignment is
        //   | offExternal - offInternal - (position[i] - minOff[i]) - T |
        // Set r_k = offExternal - offInternal - (position[i] - minOff[i]). T* = median(r_k)
        // minimizes Σ |r_k - T|.
        const residuals = [];
        for (let i = 0; i < N; i++) {
            const base = position[i] - minOff[i];
            for (const cx of externalCX[i]) {
                residuals.push(cx.offExternal - cx.offInternal - base);
            }
        }
        let T = 0;
        if (residuals.length > 0) {
            residuals.sort((a, b) => a - b);
            T = residuals[Math.floor(residuals.length / 2)];
        }
        const shifts = new Map();
        for (let i = 0; i < N; i++) {
            const delta = position[i] + T - minOff[i];
            if (delta !== 0)
                shifts.set(reorderIds[i], delta);
        }
        return shifts;
    }
    // Public entry point. Given the full combine set, decides which subset of helices actually
    // needs to shift to make the group mutually disjoint on the grid, runs the DP over that
    // subset, and returns per-helix shift deltas. Helices in the combine set that are already
    // disjoint from everyone (before AND after the DP) get no shift.
    //
    // MAX_REORDER caps the bitmask DP at 2^20 = ~1M states. Larger overlap groups return an
    // empty map, which lets the caller fall back to its normal disjointness error path.
    function computeCombineShifts(helices, combineIds, grid) {
        const empty = new Map();
        if (!Array.isArray(combineIds) || combineIds.length < 2)
            return empty;
        const uniqueIds = [...new Set(combineIds)].filter(id => id >= 0 && id < helices.length && Array.isArray(helices[id]) && helices[id].length > 0).sort((a, b) => a - b);
        if (uniqueIds.length < 2)
            return empty;
        // Seed the reorder set with helices that already collide with someone else in the group.
        const mustReorder = new Set();
        for (let i = 0; i < uniqueIds.length; i++) {
            for (let j = i + 1; j < uniqueIds.length; j++) {
                if (!disjoint(helices, uniqueIds[i], uniqueIds[j], grid)) {
                    mustReorder.add(uniqueIds[i]);
                    mustReorder.add(uniqueIds[j]);
                }
            }
        }
        if (mustReorder.size < 2)
            return empty;
        const MAX_REORDER = 20;
        if (mustReorder.size > MAX_REORDER) {
            console.warn(`[computeCombineShifts] Overlap group of ${mustReorder.size} helices exceeds cap of ${MAX_REORDER}; skipping auto-shift.`);
            return empty;
        }
        // Iterate: run DP, then check if any leave-alone helix now collides with the reordered
        // block. If yes, pull it into the reorder set and re-run. Terminates in ≤ |uniqueIds|
        // iterations because mustReorder only grows.
        let shifts = empty;
        while (true) {
            const arr = [...mustReorder];
            const dp = runReorderDP(helices, arr, grid);
            if (!dp)
                return empty;
            shifts = dp;
            let expanded = false;
            for (const id of uniqueIds) {
                if (mustReorder.has(id))
                    continue;
                for (const oid of mustReorder) {
                    if (!disjointWithShifts(helices, oid, id, grid, shifts)) {
                        mustReorder.add(id);
                        expanded = true;
                        break;
                    }
                }
                if (expanded)
                    break;
            }
            if (!expanded)
                break;
            if (mustReorder.size > MAX_REORDER) {
                console.warn(`[computeCombineShifts] Reorder set grew past cap of ${MAX_REORDER}; skipping auto-shift.`);
                return empty;
            }
        }
        return shifts;
    }
    toscad.computeCombineShifts = computeCombineShifts;
    // Applies a shift map to `grid` in-place. Every GridMark whose helixId is in `shifts` has
    // its `offset` incremented by the mapped delta. Used both by the combine-time apply and by
    // undo (with negated deltas). Returns the number of marks touched, for logging.
    function applyCombineShifts(grid, shifts) {
        if (!(grid instanceof Map) || shifts.size === 0)
            return 0;
        let touched = 0;
        for (const mark of grid.values()) {
            const delta = shifts.get(mark.helixId);
            if (delta !== undefined && delta !== 0) {
                mark.offset += delta;
                touched++;
            }
        }
        return touched;
    }
    toscad.applyCombineShifts = applyCombineShifts;
    // Chebyshev-ring expanding search for the closest (col,row) cell not present in `occupied`.
    // Rings d = 1, 2, 3, ... are traversed in a stable (dCol, dRow) order so the result is
    // deterministic. Lattice-agnostic: the visual scadnano grid places one helix per cell
    // regardless of honeycomb/square adjacency, so plain Chebyshev distance is sufficient.
    function findNearestOpenPos(anchor, occupied, maxRadius = 10000) {
        const keyOf = (c, r) => `${c},${r}`;
        const [ac, ar] = anchor;
        if (!occupied.has(keyOf(ac, ar)))
            return [ac, ar];
        for (let d = 1; d <= maxRadius; d++) {
            for (let dc = -d; dc <= d; dc++) {
                for (let dr = -d; dr <= d; dr++) {
                    if (Math.max(Math.abs(dc), Math.abs(dr)) !== d)
                        continue;
                    const c = ac + dc;
                    const r = ar + dr;
                    if (!occupied.has(keyOf(c, r)))
                        return [c, r];
                }
            }
        }
        // Fallback (unreachable in practice): walk along +col until an empty column is found.
        let c = ac + maxRadius + 1;
        while (occupied.has(keyOf(c, ar)))
            c += 1;
        return [c, ar];
    }
    // Rejects "side-by-side" helix merges by measuring how much of each helix's
    // axis segment is covered by the other's shadow. Endpoints of each helix are
    // picked as the two farthest-apart nucleotides (2-pass diameter approximation),
    // using the basepair midpoint if the pair is in the same helix and the
    // nucleotide's own position otherwise. Each endpoint of one helix is projected
    // onto the other's axis, the resulting interval is clipped to the segment
    // [0, L], and the clipped fraction of L is the coverage. Returns true when
    // max(covA, covB) >= threshold — i.e., the helices lie on top of each other
    // rather than meeting end-to-end, and the merge should be rejected.
    function axisShadowOverlap(helixA, helixB, threshold = 0.3) {
        if (!Array.isArray(helixA) || !Array.isArray(helixB))
            return false;
        if (helixA.length < 2 || helixB.length < 2)
            return false;
        const endpointsOf = (helix) => {
            const ids = new Set();
            for (const nt of helix)
                ids.add(nt.id);
            const pointOf = (nt) => {
                const own = nt.getPos();
                if (nt.pair && ids.has(nt.pair.id)) {
                    return own.clone().add(nt.pair.getPos()).multiplyScalar(0.5);
                }
                return own.clone();
            };
            // Farthest from seed, then farthest from that — approximate diameter in O(n).
            const seedPos = helix[0].getPos();
            let farA = helix[0];
            let farADist = -1;
            for (const nt of helix) {
                const d = seedPos.distanceToSquared(nt.getPos());
                if (d > farADist) {
                    farADist = d;
                    farA = nt;
                }
            }
            const farAPos = farA.getPos();
            let farB = farA;
            let farBDist = -1;
            for (const nt of helix) {
                const d = farAPos.distanceToSquared(nt.getPos());
                if (d > farBDist) {
                    farBDist = d;
                    farB = nt;
                }
            }
            if (farA === farB)
                return null;
            return [pointOf(farA), pointOf(farB)];
        };
        const epA = endpointsOf(helixA);
        const epB = endpointsOf(helixB);
        if (!epA || !epB)
            return false;
        const [P0, P1] = epA;
        const [Q0, Q1] = epB;
        // Coverage of segment S0→S1 by the shadow of X0/X1 on its axis, clipped to [0, L].
        const coverage = (S0, S1, X0, X1) => {
            const axis = S1.clone().sub(S0);
            const L = axis.length();
            if (L < 1e-9)
                return 0;
            axis.divideScalar(L);
            const t0 = X0.clone().sub(S0).dot(axis);
            const t1 = X1.clone().sub(S0).dot(axis);
            const lo = Math.max(0, Math.min(t0, t1));
            const hi = Math.min(L, Math.max(t0, t1));
            return Math.max(0, hi - lo) / L;
        };
        const covA = coverage(P0, P1, Q0, Q1);
        const covB = coverage(Q0, Q1, P0, P1);
        console.log(`[axisShadowOverlap] covA=${covA.toFixed(3)} covB=${covB.toFixed(3)}`);
        return Math.max(covA, covB) >= threshold;
    }
    toscad.axisShadowOverlap = axisShadowOverlap;
    // Splits `helixId` into two helices along the boundary defined by `nucleotides`:
    //   - The nucleotides passed in are moved into a brand-new helix appended to `helices[]`.
    //   - The remaining nucleotides stay in the original helix (helixApos = original position).
    //   - GridMap marks for the moved nucleotides get their `helixId` remapped to the new helix,
    //     while their `offset` and `direction` are left untouched — so the physical layout of
    //     every base (its column and strand direction) is preserved across the split.
    //   - The new helix is placed in the visual grid at the closest available empty cell to
    //     the original helix's position (Chebyshev-ring search via findNearestOpenPos).
    //
    // Returns the two positions plus references to the same (now mutated) helices, GridMap, and
    // helixPos. Caller is expected to trigger downstream refreshes (editor node list, angle
    // recalc, etc.) themselves.
    //
    // Returns null on invalid input:
    //   - helixId out of range,
    //   - no nucleotides passed,
    //   - no helixPos entry for helixId,
    //   - none of the supplied nucleotides actually belong to helixId,
    //   - all of helixId's nucleotides are being moved (nothing left to keep).
    function splitHelix(grid, helixPos, helixId, helices, nucleotides) {
        if (helixId < 0 || helixId >= helices.length) {
            console.warn(`[splitHelix] helixId=${helixId} out of range (helices.length=${helices.length})`);
            return null;
        }
        if (!Array.isArray(nucleotides) || nucleotides.length === 0) {
            console.warn('[splitHelix] no nucleotides supplied to split off');
            return null;
        }
        const anchor = helixPos.get(helixId);
        if (!anchor) {
            console.warn(`[splitHelix] no helixPos entry for helixId=${helixId}`);
            return null;
        }
        // Build lookup of ids to move.
        const moveIds = new Set();
        for (const nt of nucleotides) {
            if (nt instanceof Nucleotide)
                moveIds.add(nt.id);
        }
        // Partition the current helix's nucleotides. A nucleotide is only moved when both
        //   (a) it's in the caller's move list, AND
        //   (b) its GridMark still claims membership in helixId.
        // (b) protects against stale selections from a previous split/combine.
        const currentHelixNts = helices[helixId] ?? [];
        const keepList = [];
        const moveList = [];
        for (const nt of currentHelixNts) {
            if (!moveIds.has(nt.id)) {
                keepList.push(nt);
                continue;
            }
            const mark = grid.get(nt.id);
            if (!mark || mark.helixId !== helixId) {
                keepList.push(nt);
                continue;
            }
            moveList.push(nt);
        }
        if (moveList.length === 0) {
            console.warn(`[splitHelix] none of the supplied nucleotides belong to helix ${helixId}`);
            return null;
        }
        if (keepList.length === 0) {
            console.warn(`[splitHelix] all nucleotides of helix ${helixId} are being moved; nothing to keep`);
            return null;
        }
        // Apply the structural split. The new helix takes the next free slot at the end of the
        // helices[] array.
        const newHelixId = helices.length;
        helices[helixId] = keepList;
        helices.push(moveList);
        // Remap GridMap membership. Offsets and directions stay put per spec.
        for (const nt of moveList) {
            const mark = grid.get(nt.id);
            if (mark)
                mark.helixId = newHelixId;
        }
        // Place the new helix in the closest empty visual-grid cell to the original.
        const helixApos = [anchor[0], anchor[1]];
        const occupiedKeys = new Set();
        for (const pos of helixPos.values()) {
            occupiedKeys.add(`${pos[0]},${pos[1]}`);
        }
        const helixBpos = findNearestOpenPos(helixApos, occupiedKeys);
        helixPos.set(newHelixId, helixBpos);
        console.log(`[splitHelix] split helix ${helixId} (${currentHelixNts.length} nts) -> ` +
            `kept ${keepList.length} nts at [${helixApos[0]},${helixApos[1]}], ` +
            `moved ${moveList.length} nts into new helix ${newHelixId} at [${helixBpos[0]},${helixBpos[1]}]`);
        return {
            keptHelixId: helixId,
            newHelixId,
            helixApos,
            helixBpos,
            helices,
            grid,
            helixPos
        };
    }
    toscad.splitHelix = splitHelix;
    // anglecomb3 — three phases per iteration:
    //
    //   Phase 1 — hash pairs: MERGE UNCONDITIONALLY.
    //             No gates. No cos threshold. If hashAxisOverlap nominated
    //             the pair (its own end-to-end cylinder + antiparallel-axis
    //             tests already passed upstream), we trust it. All hash
    //             pairs are pIdx-exclusive by construction, so ordering
    //             doesn't matter — we retranslate pIdx→helixId per pair.
    //             cos / axisDot / disjoint status are logged for visibility
    //             but do NOT gate anything.
    //
    //   Phase 2 — anglecomb2 strict: bucket by (viaParent, exactLocalAngle),
    //             gates 1+2+3, strict-parent unanimity. One merge per pass,
    //             loops to fixed-point.
    //
    //   Phase 3 — QUEUE-DRAIN of anglecomb2 relaxed passers:
    //             candidates = { bucket pair + gates 1+2+3 + ≥1 strong
    //                            mutual witness (gate 4b) }
    //             Sort by cos desc, axisDot desc (tiebreak).
    //             Drain greedy: each helix participates in at most one
    //             merge per drain (any subsequent entry touching an already
    //             consumed helix is dropped and gets a fresh chance in the
    //             next outer iteration).
    //
    // Outer iterations repeat until zero merges occur. Merge mechanic is
    // anglecomb2's light-touch: helixId remap on the grid + splice on the
    // helices array. Phase 1 can produce non-disjoint merges (nucleotides
    // at overlapping offsets in the merged helix); a warning is logged in
    // that case.
    function anglecomb3(grid, helices, lattice = 'honeycomb', partials, hashMergePairs, angleMap = getAngles(grid, helices, lattice)) {
        const AXIS_DOT_THRESHOLD = 0.935;
        const MAX_ITERATIONS = 200;
        let networkMap = angleMap;
        const mergedPairs = [];
        // ── shared helpers ──────────────────────────────────────────────
        // (kept local; anglecomb2 has near-identical inline versions, but
        // duplicating rather than sharing avoids coupling anglecomb2's evolution.)
        // Light merge — matches anglecomb2's mergeHelixInto exactly.
        const mergeHelixInto = (keepHelix, mergedHelix) => {
            const keepNts = helices[keepHelix] ?? [];
            const mergedNts = helices[mergedHelix] ?? [];
            const keepNtIds = keepNts.map(n => n.id);
            const mergedNtIds = mergedNts.map(n => n.id);
            const remapHelixId = (helixId) => {
                if (helixId === mergedHelix)
                    return keepHelix;
                if (helixId > mergedHelix)
                    return helixId - 1;
                return helixId;
            };
            for (const [, mark] of grid.entries()) {
                mark.helixId = remapHelixId(mark.helixId);
            }
            helices[keepHelix] = keepNts.concat(mergedNts);
            helices.splice(mergedHelix, 1);
            return { keepNtIds, mergedNtIds };
        };
        // Helix 3D axis — same as anglecomb2's helper.
        const helixAxis = (helixId) => {
            const nts = helices[helixId] ?? [];
            let minOffset = Number.POSITIVE_INFINITY;
            let maxOffset = Number.NEGATIVE_INFINITY;
            let minNt = null;
            let maxNt = null;
            for (const nt of nts) {
                if (!(nt instanceof Nucleotide))
                    continue;
                const mark = grid.get(nt.id);
                if (!mark || mark.helixId !== helixId)
                    continue;
                if (mark.offset < minOffset) {
                    minOffset = mark.offset;
                    minNt = nt;
                }
                if (mark.offset > maxOffset) {
                    maxOffset = mark.offset;
                    maxNt = nt;
                }
            }
            if (!minNt || !maxNt || minNt === maxNt)
                return null;
            const dir = minNt.getPos().clone().sub(maxNt.getPos());
            const len = dir.length();
            if (!isFinite(len) || len === 0)
                return null;
            return dir.divideScalar(len);
        };
        const pIdxToHelixId = (pIdx) => {
            const p = partials?.[pIdx];
            if (!Array.isArray(p) || p.length === 0)
                return -1;
            const mark = grid.get(p[0].id);
            return (mark && typeof mark.helixId === 'number') ? mark.helixId : -1;
        };
        // Gate 2 (axis collinearity, |dot| >= threshold) and Gate 3
        // (axis-shadow non-overlap). Returned axisDot is |dot| — used for
        // tiebreaking in phase 4.
        const gate23 = (A, B) => {
            const axisA = helixAxis(A);
            const axisB = helixAxis(B);
            if (!axisA || !axisB)
                return { pass: false, axisDot: 0, reason: 'missing-axis' };
            const axisDot = Math.abs(axisA.dot(axisB));
            if (axisDot < AXIS_DOT_THRESHOLD)
                return { pass: false, axisDot, reason: 'axisDot-below-threshold' };
            if (axisShadowOverlap(helices[A], helices[B]))
                return { pass: false, axisDot, reason: 'shadow-overlap' };
            return { pass: true, axisDot, reason: 'ok' };
        };
        // (viaParent, exactLocalAngle) buckets from preds — same signal as
        // anglecomb2's collision detector.
        const collectBuckets = (preds) => {
            const bucket = new Map();
            for (const [hid, plist] of preds.entries()) {
                for (const p of plist) {
                    if (p.viaParent < 0)
                        continue;
                    const key = `${p.viaParent}|${p.edgeLocalAngle}`;
                    let b = bucket.get(key);
                    if (!b) {
                        b = { parent: p.viaParent, localAngle: p.edgeLocalAngle, coord: p.coord, hids: new Set() };
                        bucket.set(key, b);
                    }
                    b.hids.add(hid);
                }
            }
            return [...bucket.values()];
        };
        const inAnyBucket = (A, B, buckets) => {
            for (const b of buckets) {
                if (b.hids.has(A) && b.hids.has(B))
                    return true;
            }
            return false;
        };
        // Gate 4a — strict-parent unanimity for a child at a specific cell.
        // Reused from anglecomb2 verbatim.
        const strictParentsAgree = (childHid, cell, preds, cc) => {
            let sawStrong = false;
            for (const p of preds.get(childHid) ?? []) {
                if (p.viaParent < 0)
                    continue;
                const w = cc.get(p.viaParent)?.get(childHid) ?? 0;
                if (w < 2)
                    continue;
                sawStrong = true;
                if (p.coord[0] !== cell[0] || p.coord[1] !== cell[1])
                    return false;
            }
            return sawStrong;
        };
        // Gate 4b — strong mutual witnesses (≥2 crossovers to BOTH A and B).
        // Reused from anglecomb2.
        const computeStrongMutuals = (A, B, preds, cc) => {
            const parentsA = new Set();
            for (const p of preds.get(A) ?? [])
                if (p.viaParent >= 0)
                    parentsA.add(p.viaParent);
            const out = [];
            for (const n of parentsA) {
                if (n === A || n === B)
                    continue;
                let isBParent = false;
                for (const p of preds.get(B) ?? []) {
                    if (p.viaParent === n) {
                        isBParent = true;
                        break;
                    }
                }
                if (!isBParent)
                    continue;
                const nToA = cc.get(n)?.get(A) ?? 0;
                const nToB = cc.get(n)?.get(B) ?? 0;
                if (nToA >= 2 && nToB >= 2)
                    out.push(n);
            }
            return out.sort((a, b) => a - b);
        };
        const minNtIdOf = (hid) => {
            const nts = helices[hid] ?? [];
            let m = Number.POSITIVE_INFINITY;
            for (const nt of nts)
                if (nt && nt.id < m)
                    m = nt.id;
            return isFinite(m) ? m : -1;
        };
        const haveHashPairs = Array.isArray(hashMergePairs) && hashMergePairs.length > 0;
        const havePartials = Array.isArray(partials) && partials.length > 0;
        if (!haveHashPairs || !havePartials) {
            console.log(`[anglecomb3] no partials / no hash pairs — hash-driven phases will be skipped`);
        }
        // ── outer iteration ────────────────────────────────────────────
        let iteration = 0;
        let anyMergeThisIteration = true;
        while (anyMergeThisIteration && iteration < MAX_ITERATIONS) {
            iteration++;
            anyMergeThisIteration = false;
            // ── PHASE 1: hash pairs — UNCONDITIONAL merge ──────────────
            // hashAxisOverlap already ran its own end-to-end tests (cylinder
            // overlap + antiparallel free-side axes) upstream. If it says
            // yes, we merge. cos / axisDot / gate outcomes are logged for
            // visibility but do NOT gate anything.
            if (haveHashPairs && havePartials) {
                let phase1Merged = true;
                while (phase1Merged) {
                    phase1Merged = false;
                    const cc = getConnectionCounts(grid);
                    for (const mp of hashMergePairs) {
                        const hA = pIdxToHelixId(mp.a);
                        const hB = pIdxToHelixId(mp.b);
                        if (hA < 0 || hB < 0 || hA === hB)
                            continue;
                        if (!helices[hA] || !helices[hB])
                            continue;
                        const cos = helixPairCosine(hA, hB, cc);
                        const g23 = gate23(hA, hB);
                        const wasDisjoint = disjoint(helices, hA, hB, grid);
                        const shadowOverlap = helices[hA] && helices[hB]
                            ? axisShadowOverlap(helices[hA], helices[hB])
                            : false;
                        const keep = Math.min(hA, hB);
                        const merged = Math.max(hA, hB);
                        const ids = mergeHelixInto(keep, merged);
                        mergedPairs.push({
                            keepHelix: keep, mergedHelix: merged,
                            keepNtIds: ids.keepNtIds, mergedNtIds: ids.mergedNtIds,
                            phase: 'hash', cos, axisDot: g23.axisDot,
                            wasDisjoint, shadowOverlap,
                        });
                        console.log(`[anglecomb3] iter=${iteration} phase=hash merged ${merged}→${keep} ` +
                            `cos=${cos.toFixed(3)} axisDot=${g23.axisDot.toFixed(3)} ` +
                            `disjoint=${wasDisjoint} shadowOverlap=${shadowOverlap}` +
                            (!g23.pass ? ` (gate23 would have rejected: ${g23.reason})` : ``));
                        if (!wasDisjoint) {
                            console.warn(`[anglecomb3] ⚠️ non-disjoint hash merge ${merged}→${keep} — ` +
                                `merged helix will have overlapping offsets on its number line.`);
                        }
                        if (shadowOverlap) {
                            console.warn(`[anglecomb3] ⚠️ hash merge ${merged}→${keep} shadow-overlaps — ` +
                                `axisShadowOverlap would have flagged this as side-by-side. Trusting hash.`);
                        }
                        networkMap = getAngles(grid, helices, lattice);
                        phase1Merged = true;
                        anyMergeThisIteration = true;
                        break; // restart phase-1 loop with fresh cc / helixIds
                    }
                }
            }
            // ── PHASE 2: anglecomb2 strict (one merge per pass, to fixed-point) ──
            {
                let phase2Merged = true;
                while (phase2Merged) {
                    phase2Merged = false;
                    const preds = tempGlobalPos(networkMap, grid, lattice);
                    const cc = getConnectionCounts(grid);
                    const buckets = collectBuckets(preds);
                    outer2: for (const bucket of buckets) {
                        if (bucket.hids.size < 2)
                            continue;
                        const sorted = [...bucket.hids].sort((a, b) => a - b);
                        for (let i = 0; i < sorted.length; i++) {
                            const A = sorted[i];
                            for (let j = i + 1; j < sorted.length; j++) {
                                const B = sorted[j];
                                if (!disjoint(helices, A, B, grid))
                                    continue;
                                const g23 = gate23(A, B);
                                if (!g23.pass)
                                    continue;
                                if (!strictParentsAgree(A, bucket.coord, preds, cc))
                                    continue;
                                if (!strictParentsAgree(B, bucket.coord, preds, cc))
                                    continue;
                                const cos = helixPairCosine(A, B, cc);
                                const keep = Math.min(A, B);
                                const merged = Math.max(A, B);
                                const ids = mergeHelixInto(keep, merged);
                                mergedPairs.push({
                                    keepHelix: keep, mergedHelix: merged,
                                    keepNtIds: ids.keepNtIds, mergedNtIds: ids.mergedNtIds,
                                    phase: 'strict', cos, axisDot: g23.axisDot,
                                    source: bucket.parent, cell: bucket.coord,
                                });
                                console.log(`[anglecomb3] iter=${iteration} phase=strict merged ${merged}→${keep} ` +
                                    `at (${bucket.coord[0]},${bucket.coord[1]}) source=${bucket.parent} ` +
                                    `cos=${cos.toFixed(3)} axisDot=${g23.axisDot.toFixed(3)}`);
                                networkMap = getAngles(grid, helices, lattice);
                                phase2Merged = true;
                                anyMergeThisIteration = true;
                                break outer2;
                            }
                        }
                    }
                }
            }
            // ── PHASE 3: queue-drain of relaxed anglecomb2 passers ──────
            {
                const preds = tempGlobalPos(networkMap, grid, lattice);
                const cc = getConnectionCounts(grid);
                const buckets = collectBuckets(preds);
                // dedup key: unordered pair. Prefer the higher (cos, axisDot).
                const queueMap = new Map();
                const enqueue = (e) => {
                    const key = `${Math.min(e.A, e.B)}|${Math.max(e.A, e.B)}`;
                    const existing = queueMap.get(key);
                    if (!existing ||
                        e.cos > existing.cos ||
                        (e.cos === existing.cos && e.axisDot > existing.axisDot)) {
                        queueMap.set(key, e);
                    }
                };
                // Relaxed anglecomb2 — bucket pairs passing gates 1+2+3+4b
                for (const bucket of buckets) {
                    if (bucket.hids.size < 2)
                        continue;
                    const sorted = [...bucket.hids].sort((a, b) => a - b);
                    for (let i = 0; i < sorted.length; i++) {
                        const A = sorted[i];
                        for (let j = i + 1; j < sorted.length; j++) {
                            const B = sorted[j];
                            if (!disjoint(helices, A, B, grid))
                                continue;
                            const g23 = gate23(A, B);
                            if (!g23.pass)
                                continue;
                            const witnesses = computeStrongMutuals(A, B, preds, cc);
                            if (witnesses.length === 0)
                                continue;
                            const cos = helixPairCosine(A, B, cc);
                            const minNtA = minNtIdOf(A);
                            const minNtB = minNtIdOf(B);
                            if (minNtA < 0 || minNtB < 0)
                                continue;
                            enqueue({
                                A, B, minNtA, minNtB,
                                cos, axisDot: g23.axisDot,
                                cell: bucket.coord, source: bucket.parent, witnesses,
                            });
                        }
                    }
                }
                // Sort: cos desc, axisDot desc.
                const queue = [...queueMap.values()].sort((x, y) => y.cos - x.cos || y.axisDot - x.axisDot);
                if (queue.length > 0) {
                    console.log(`[anglecomb3] iter=${iteration} phase=queue built ${queue.length} candidate(s) ` +
                        `(top cos=${queue[0].cos.toFixed(3)}, bottom cos=${queue[queue.length - 1].cos.toFixed(3)})`);
                }
                // Drain greedy. Track consumed helices by their pre-drain
                // representative nt id (immutable across merges).
                const consumedNts = new Set();
                let queueMergedAny = false;
                for (const entry of queue) {
                    if (consumedNts.has(entry.minNtA) || consumedNts.has(entry.minNtB)) {
                        console.log(`[anglecomb3] iter=${iteration} phase=queue-relaxed skipped (${entry.A},${entry.B}) ` +
                            `cos=${entry.cos.toFixed(3)} — helix consumed earlier in drain`);
                        continue;
                    }
                    // Resolve to CURRENT helixIds via nt-id → grid mark.
                    const currA = grid.get(entry.minNtA)?.helixId ?? -1;
                    const currB = grid.get(entry.minNtB)?.helixId ?? -1;
                    if (currA < 0 || currB < 0 || currA === currB)
                        continue;
                    if (!helices[currA] || !helices[currB])
                        continue;
                    const keep = Math.min(currA, currB);
                    const merged = Math.max(currA, currB);
                    const ids = mergeHelixInto(keep, merged);
                    mergedPairs.push({
                        keepHelix: keep, mergedHelix: merged,
                        keepNtIds: ids.keepNtIds, mergedNtIds: ids.mergedNtIds,
                        phase: 'queue-relaxed',
                        cos: entry.cos,
                        axisDot: entry.axisDot,
                        source: entry.source,
                        cell: entry.cell,
                        witnesses: entry.witnesses,
                    });
                    console.log(`[anglecomb3] iter=${iteration} phase=queue-relaxed merged ${merged}→${keep} ` +
                        `cos=${entry.cos.toFixed(3)} axisDot=${entry.axisDot.toFixed(3)}` +
                        (entry.cell ? ` cell=(${entry.cell[0]},${entry.cell[1]})` : '') +
                        (entry.source !== undefined ? ` source=${entry.source}` : '') +
                        (entry.witnesses && entry.witnesses.length > 0
                            ? ` witnesses=[${entry.witnesses.join(',')}]` : ''));
                    consumedNts.add(entry.minNtA);
                    consumedNts.add(entry.minNtB);
                    queueMergedAny = true;
                }
                if (queueMergedAny) {
                    networkMap = getAngles(grid, helices, lattice);
                    anyMergeThisIteration = true;
                }
            }
            console.log(`[anglecomb3] iter=${iteration} complete — ` +
                (anyMergeThisIteration ? 'had merges, continuing' : 'no merges, done'));
        }
        if (iteration >= MAX_ITERATIONS) {
            console.warn(`[anglecomb3] hit MAX_ITERATIONS=${MAX_ITERATIONS}`);
        }
        return { networkMap, mergedPairs };
    }
    toscad.anglecomb3 = anglecomb3;
    // anglecorr2 — global-position variant of anglecorr.
    //
    // Consumes the multimap from tempGlobalPos. Detection is identical to
    // anglecomb2 (bucket by (viaParent, exact local angle)), but this
    // function handles the OTHER branch of the shared signal: when the
    // colliding pair is NOT disjoint on the offset axis, they physically
    // overlap and can't coexist at the same lattice slot. One must move.
    //
    // Loser identification is unchanged from current anglecorr: source's
    // baseHelix (0° neighbor) is inviolate; otherwise, consensus mass
    // (Σ crossover counts to other neighbors) picks the winner. Rotation
    // steps 120° (honeycomb) or 90° (square) through source's row until an
    // unused slot is found.
    //
    // Two rounds:
    //   Strict: skip ties in the mass tiebreak. Only clear winners rotate.
    //   Relaxed: on tie, higher-id loses (matches current anglecorr).
    //
    // Mutates networkMap in place and returns it. Does NOT re-run getAngles
    // between passes — that would wipe the rotations, since getAngles
    // derives from the grid and rotations don't touch the grid.
    function anglecorr2(grid, helices, lattice = 'honeycomb', angleMap = getAngles(grid, helices, lattice)) {
        const MAX_PASSES = 200;
        const latticeType = resolveLatticeKind(lattice);
        const correctionStep = latticeType === 'square' ? 90 : 120;
        const maxLatticeAngles = latticeType === 'square' ? 4 : 3;
        const normalizeAngle = (a) => ((a % 360) + 360) % 360;
        let networkMap = angleMap;
        const correctedPairs = [];
        for (const round of ['strict', 'relaxed']) {
            let pass = 0;
            while (pass++ < MAX_PASSES) {
                const preds = tempGlobalPos(networkMap, grid, lattice);
                const connectionCounts = getConnectionCounts(grid);
                // Same bucket detection as anglecomb2: (viaParent, exact
                // local angle). All predictions with the same source and
                // same local angle land in one bucket.
                const bucket = new Map();
                for (const [hid, plist] of preds.entries()) {
                    for (const p of plist) {
                        if (p.viaParent < 0)
                            continue;
                        const key = `${p.viaParent}|${p.edgeLocalAngle}`;
                        let b = bucket.get(key);
                        if (!b) {
                            b = { parent: p.viaParent, localAngle: p.edgeLocalAngle, hids: new Set() };
                            bucket.set(key, b);
                        }
                        b.hids.add(hid);
                    }
                }
                const candidates = [...bucket.values()].filter(b => b.hids.size >= 2);
                if (candidates.length === 0) {
                    console.log(`[anglecorr2] round=${round} pass=${pass} no conflicts`);
                    break;
                }
                let correctedThisPass = false;
                outer: for (const { parent: source, localAngle: conflictAngle, hids } of candidates) {
                    const sourceMap = networkMap.get(source);
                    if (!sourceMap)
                        continue;
                    // baseHelix = source's neighbor at local angle 0. If a
                    // conflict pair contains it, it stays put (its 0° is
                    // the anchor of source's local frame).
                    const baseCandidates = [];
                    for (const [nid, ang] of sourceMap.entries()) {
                        if (normalizeAngle(ang) === 0)
                            baseCandidates.push(nid);
                    }
                    const baseHelix = baseCandidates.length > 0 ? Math.min(...baseCandidates) : null;
                    const sorted = [...hids].sort((a, b) => a - b);
                    for (let i = 0; i < sorted.length; i++) {
                        const A = sorted[i];
                        for (let j = i + 1; j < sorted.length; j++) {
                            const B = sorted[j];
                            // Disjoint → anglecomb2's territory, not ours.
                            if (disjoint(helices, A, B, grid))
                                continue;
                            // Source can hold at most maxLatticeAngles
                            // distinct-angle neighbors. If it's already
                            // full, there's no slot to rotate into.
                            if (sourceMap.size > maxLatticeAngles) {
                                console.warn(`[anglecorr2] source=${source} has ${sourceMap.size} connections ` +
                                    `and non-disjoint conflict at angle=${conflictAngle}°; no free ` +
                                    `${latticeType} slots available`);
                                continue;
                            }
                            // ── Identify the loser (the helix to rotate) ──
                            let adjustedHelix;
                            let decisionLog;
                            if (A === baseHelix) {
                                adjustedHelix = B;
                                decisionLog = `baseHelix=${A} inviolate`;
                            }
                            else if (B === baseHelix) {
                                adjustedHelix = A;
                                decisionLog = `baseHelix=${B} inviolate`;
                            }
                            else {
                                // Consensus mass: sum of crossover counts to each
                                // helix's OTHER neighbors (excluding source and
                                // the competing candidate). A frayed fragment
                                // scores near-zero; a well-supported structural
                                // body scores high.
                                const massOf = (candidateId) => {
                                    const nMap = networkMap.get(candidateId);
                                    if (!nMap)
                                        return { mass: 0, voters: 0 };
                                    const wRow = connectionCounts.get(candidateId);
                                    let mass = 0;
                                    let voters = 0;
                                    for (const nid of nMap.keys()) {
                                        if (nid === source)
                                            continue;
                                        if (nid === A || nid === B)
                                            continue;
                                        const w = wRow?.get(nid) ?? 0;
                                        if (w > 0) {
                                            mass += w;
                                            voters++;
                                        }
                                    }
                                    return { mass, voters };
                                };
                                const { mass: massA, voters: votersA } = massOf(A);
                                const { mass: massB, voters: votersB } = massOf(B);
                                if (massA === massB) {
                                    if (round === 'strict')
                                        continue; // strict: no tiebreak
                                    // Relaxed tiebreak: structural weight =
                                    // (num helices) × (total crossovers) =
                                    // voters × mass. When mass ties, this
                                    // rewards the more-distributed candidate
                                    // over the more-concentrated one.
                                    const swA = massA * votersA;
                                    const swB = massB * votersB;
                                    if (swA !== swB) {
                                        adjustedHelix = swA > swB ? B : A;
                                        decisionLog = `tied mass=${massA}; structural weight A=${swA} vs B=${swB}`;
                                    }
                                    else {
                                        adjustedHelix = Math.max(A, B);
                                        decisionLog = `tied mass=${massA} and structural weight=${swA}; higher id loses`;
                                    }
                                }
                                else if (massA > massB) {
                                    adjustedHelix = B;
                                    decisionLog = `A wins mass=${massA}/${votersA} over B mass=${massB}/${votersB}`;
                                }
                                else {
                                    adjustedHelix = A;
                                    decisionLog = `B wins mass=${massB}/${votersB} over A mass=${massA}/${votersA}`;
                                }
                            }
                            // ── Rotate loser's angle into an unused slot ──
                            const oldAngle = normalizeAngle(sourceMap.get(adjustedHelix) ?? 0);
                            const usedAngles = new Set();
                            for (const [nid, ang] of sourceMap.entries()) {
                                if (nid === adjustedHelix)
                                    continue;
                                usedAngles.add(normalizeAngle(ang));
                            }
                            let newAngle = oldAngle;
                            let foundOpen = false;
                            for (let attempt = 1; attempt < maxLatticeAngles; attempt++) {
                                const cand = normalizeAngle(oldAngle + attempt * correctionStep);
                                if (!usedAngles.has(cand)) {
                                    newAngle = cand;
                                    foundOpen = true;
                                    break;
                                }
                            }
                            if (!foundOpen) {
                                console.warn(`[anglecorr2] source=${source} could not place corrected ` +
                                    `angle for helix=${adjustedHelix}; all ${correctionStep}° ` +
                                    `slots occupied`);
                                continue;
                            }
                            sourceMap.set(adjustedHelix, newAngle);
                            correctedPairs.push({
                                source, adjustedHelix, oldAngle, newAngle,
                                conflictLocalAngle: conflictAngle, round, decision: decisionLog
                            });
                            console.log(`[anglecorr2] round=${round} source=${source} rotated helix ${adjustedHelix}: ` +
                                `${oldAngle}° → ${newAngle}° (conflict=${conflictAngle}°, ${decisionLog})`);
                            correctedThisPass = true;
                            break outer;
                        }
                    }
                }
                if (!correctedThisPass) {
                    console.log(`[anglecorr2] round=${round} pass=${pass} ${candidates.length} conflict(s) ` +
                        `but none passable`);
                    break;
                }
                // Deliberately no getAngles refresh here. Rotations mutate
                // networkMap in place; getAngles would rebuild it from the
                // grid and erase our work.
            }
        }
        return { networkMap, correctedPairs };
    }
    toscad.anglecorr2 = anglecorr2;
    // ── Shared helper: collect all crossover shift observations ──────
    // For each pair of helices connected by backbone crossovers, returns
    // the list of observed shifts (offsetA - offsetB for each crossover
    // from A→B). Used by both alignGridPrim and alignGridDP.
    function collectShiftObservations(grid) {
        const allNtIds = new Set();
        for (const [ntId] of grid.entries())
            allNtIds.add(ntId);
        const visited = new Set();
        const helixIds = new Set();
        // shifts[a][b] = array of (offsetA - offsetB) values
        const shifts = new Map();
        const ensurePair = (a, b) => {
            if (!shifts.has(a))
                shifts.set(a, new Map());
            if (!shifts.get(a).has(b))
                shifts.get(a).set(b, []);
            return shifts.get(a).get(b);
        };
        for (const [ntId] of grid.entries()) {
            if (visited.has(ntId))
                continue;
            const startNt = elements.get(ntId);
            if (!startNt || !(startNt instanceof Nucleotide))
                continue;
            // Find 5' end
            let fivePrime = startNt;
            const walkBack = new Set();
            walkBack.add(fivePrime.id);
            while (true) {
                const prev = fivePrime.n5;
                if (!prev || !(prev instanceof Nucleotide))
                    break;
                if (!allNtIds.has(prev.id))
                    break;
                if (walkBack.has(prev.id))
                    break;
                walkBack.add(prev.id);
                fivePrime = prev;
            }
            // Walk 5'→3'
            let curr = fivePrime;
            const walkForward = new Set();
            let prevMark = null;
            while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                if (walkForward.has(curr.id))
                    break;
                walkForward.add(curr.id);
                visited.add(curr.id);
                const mark = grid.get(curr.id);
                if (mark) {
                    helixIds.add(mark.helixId);
                    if (prevMark && prevMark.helixId !== mark.helixId) {
                        // Crossover: shift = offsetFrom - offsetTo
                        // If we add this value to helix "to", the crossover aligns.
                        const shiftVal = prevMark.offset - mark.offset;
                        ensurePair(prevMark.helixId, mark.helixId).push(shiftVal);
                        // Reverse: offsetTo - offsetFrom = -shiftVal
                        ensurePair(mark.helixId, prevMark.helixId).push(-shiftVal);
                    }
                    prevMark = mark;
                }
                else {
                    prevMark = null;
                }
                const n3ref = curr.n3;
                curr = (n3ref && n3ref instanceof Nucleotide) ? n3ref : null;
            }
        }
        return { shifts, helixIds };
    }
    toscad.collectShiftObservations = collectShiftObservations;
    // ── Helper: compute median of a sorted-or-unsorted number array ─
    function median(arr) {
        if (arr.length === 0)
            return 0;
        const sorted = arr.slice().sort((a, b) => a - b);
        const mid = Math.floor(sorted.length / 2);
        return sorted.length % 2 === 0
            ? Math.round((sorted[mid - 1] + sorted[mid]) / 2)
            : sorted[mid];
    }
    function buildCellIndex(grid) {
        const idx = new Map();
        for (const [, m] of grid.entries()) {
            let hm = idx.get(m.helixId);
            if (!hm) {
                hm = new Map();
                idx.set(m.helixId, hm);
            }
            const key = `${m.direction}|${m.offset}`;
            if (!hm.has(key))
                hm.set(key, new Set());
            hm.get(key).add(m.helixId);
        }
        return idx;
    }
    // A cell is "free for selfHid" when nothing occupies it, or only marks of
    // selfHid occupy it (a helix never conflicts with itself — its own marks
    // are presumed internally consistent before any move).
    function cellFree(idx, hid, selfHid, dir, off) {
        const hm = idx.get(hid);
        if (!hm)
            return true;
        const s = hm.get(`${dir}|${off}`);
        if (!s)
            return true;
        if (selfHid === hid)
            return true;
        for (const other of s)
            if (other !== selfHid)
                return false;
        return true;
    }
    // Would shifting every mark of helix hid by delta cause a collision?
    function collidesAt(grid, idx, hid, delta) {
        if (delta === 0)
            return false;
        for (const [, m] of grid.entries()) {
            if (m.helixId !== hid)
                continue;
            if (!cellFree(idx, hid, hid, m.direction, m.offset + delta))
                return true;
        }
        return false;
    }
    // Pick the best non-colliding delta for helix hid.
    // Tries the ideal shift first, then walks outward ±1, ±2, … up to
    // maxRadius (±1 keeps crossover parity so even ±2 crossovers still align).
    // Returns null when no collision-free placement exists within the radius —
    // the caller must then keep the helix's previous position (i.e. not move).
    function pickCollisionFreeDelta(grid, hid, idealDelta, maxRadius = 24) {
        const idx = buildCellIndex(grid);
        if (!collidesAt(grid, idx, hid, idealDelta))
            return idealDelta;
        for (let radius = 1; radius <= maxRadius; radius++) {
            for (const sign of [1, -1]) {
                const cand = idealDelta + sign * radius;
                if (!collidesAt(grid, idx, hid, cand))
                    return cand;
            }
        }
        return null;
    }
    // Apply a single helix's shift and then re-normalize all offsets so the
    // global minimum sits at 0 (keeps offsets non-negative even when the
    // minimum was produced by a rejected move's neighbors).
    function applySingleShift(grid, hid, delta) {
        if (delta === 0)
            return;
        for (const [, m] of grid.entries())
            if (m.helixId === hid)
                m.offset += delta;
        let globalMin = Infinity;
        for (const [, m] of grid.entries())
            if (m.offset < globalMin)
                globalMin = m.offset;
        if (globalMin !== 0 && globalMin !== Infinity) {
            for (const [, m] of grid.entries())
                m.offset -= globalMin;
        }
    }
    // ── Helper: apply shift to every nt on a helix, then normalize ──
    function applyHelixShifts(grid, shiftMap) {
        for (const [, mark] of grid.entries()) {
            const s = shiftMap.get(mark.helixId);
            if (s !== undefined && s !== 0) {
                mark.offset += s;
            }
        }
        // Normalize: find global minimum, shift everything so min = 0
        let globalMin = Infinity;
        for (const [, mark] of grid.entries()) {
            if (mark.offset < globalMin)
                globalMin = mark.offset;
        }
        if (globalMin !== 0 && globalMin !== Infinity) {
            for (const [, mark] of grid.entries()) {
                mark.offset -= globalMin;
            }
        }
    }
    /**
     * alignGridPrim — MST-based offset alignment.
     *
     * Builds a complete graph where nodes = helices and edge weight =
     * number of crossover observations between them. Prim's MST selects
     * the most reliable (most-observed) edges. Walks the MST from helix 0
     * to compute cumulative shifts using the median of observed shifts
     * per edge (robust to outlier crossovers).
     *
     * If binderHelices is provided, those helices are excluded from the
     * global MST alignment. Instead, each strand segment on a binder
     * helix is aligned individually to the crossover offset of its parent.
     */
    function alignGridPrim(grid, binderHelices) {
        const { shifts, helixIds } = collectShiftObservations(grid);
        const countGridConflicts = (currentGrid) => {
            const checkMap = new Map();
            let conflicts = 0;
            for (const [ntId, pos] of currentGrid.entries()) {
                if (!checkMap.has(pos.helixId)) {
                    checkMap.set(pos.helixId, {
                        forward: new Map(),
                        backward: new Map()
                    });
                }
                const strandMap = checkMap.get(pos.helixId)[pos.direction];
                if (strandMap.has(pos.offset) && strandMap.get(pos.offset) !== ntId) {
                    conflicts++;
                }
                else {
                    strandMap.set(pos.offset, ntId);
                }
            }
            return conflicts;
        };
        // ── Build weighted edge list for Prim's ─────────────────────────
        // weight = number of crossover observations (higher = more reliable)
        const helixList = Array.from(helixIds).sort((a, b) => a - b);
        if (helixList.length <= 1) {
            console.log('[alignGridPrim] Only 0-1 helices, nothing to align.');
            return;
        }
        // Prim's MST starting from helix 0
        const inMST = new Set();
        // mstEdges: parent → child with median shift
        const mstParent = new Map();
        // Priority: pick the edge with the highest weight (most observations)
        inMST.add(0);
        while (inMST.size < helixList.length) {
            let bestNeighbor = -1;
            let bestFrom = -1;
            let bestWeight = 0;
            for (const inNode of inMST) {
                const neighbors = shifts.get(inNode);
                if (!neighbors)
                    continue;
                for (const [neighbor, observations] of neighbors.entries()) {
                    if (inMST.has(neighbor))
                        continue;
                    if (observations.length > bestWeight) {
                        bestWeight = observations.length;
                        bestNeighbor = neighbor;
                        bestFrom = inNode;
                    }
                }
            }
            if (bestNeighbor === -1) {
                // Disconnected graph — pick an unvisited helix, anchor it
                for (const hId of helixList) {
                    if (!inMST.has(hId)) {
                        inMST.add(hId);
                        // No parent (disconnected), shift = 0 relative to itself
                        break;
                    }
                }
                continue;
            }
            // The median shift from bestFrom → bestNeighbor
            const observations = shifts.get(bestFrom).get(bestNeighbor);
            const medianShift = median(observations);
            mstParent.set(bestNeighbor, { parent: bestFrom, shift: medianShift });
            inMST.add(bestNeighbor);
        }
        // ── Walk MST from helix 0 to compute cumulative shifts ──────────
        const cumulativeShift = new Map();
        cumulativeShift.set(0, 0); // anchor
        // BFS order: process nodes so parent's cumulative shift is known
        const bfsQueue = [0];
        let qi = 0;
        // Build children adjacency from mstParent
        const children = new Map();
        for (const [child, { parent }] of mstParent.entries()) {
            if (!children.has(parent))
                children.set(parent, []);
            children.get(parent).push(child);
        }
        while (qi < bfsQueue.length) {
            const node = bfsQueue[qi++];
            const nodeShift = cumulativeShift.get(node) ?? 0;
            const kids = children.get(node) ?? [];
            for (const child of kids) {
                const edge = mstParent.get(child);
                // edge.shift = offsetParent - offsetChild at crossover
                // To align child with parent: child += edge.shift + parentCumulativeShift
                // Actually: cumulativeShift[child] = cumulativeShift[parent] + edge.shift
                cumulativeShift.set(child, nodeShift + edge.shift);
                bfsQueue.push(child);
            }
        }
        // Handle disconnected helices (not in MST tree from 0)
        for (const hId of helixList) {
            if (!cumulativeShift.has(hId)) {
                cumulativeShift.set(hId, 0);
            }
        }
        // ── Apply shifts (non-binder helices) ─────────────────────────
        // Sequential, overlap-aware placement: each helix is moved to the
        // collision-free delta closest to its MST-ideal cumulative shift.
        // Applying one helix at a time lets later helices see (and avoid)
        // the cells occupied by earlier ones — including pseudo-helices of a
        // split merged helix, which would otherwise fold back onto the same
        // helixId with colliding offsets. A helix with no valid placement
        // within the search radius keeps its previous position.
        {
            const orderedHelices = helixList
                .filter(h => h !== 0 && (cumulativeShift.get(h) ?? 0) !== 0)
                .sort((a, b) => {
                const da = Math.abs(cumulativeShift.get(a) ?? 0);
                const db = Math.abs(cumulativeShift.get(b) ?? 0);
                return db !== da ? db - da : a - b;
            });
            for (const hid of orderedHelices) {
                const ideal = cumulativeShift.get(hid) ?? 0;
                const delta = pickCollisionFreeDelta(grid, hid, ideal);
                if (delta === null) {
                    console.warn(`[alignGridPrim] No collision-free placement for helix ${hid} (ideal shift ${ideal}); keeping current position.`);
                    continue;
                }
                if (delta !== ideal) {
                    console.log(`[alignGridPrim] Helix ${hid}: ideal shift ${ideal} would overlap; placed at ${delta} instead.`);
                }
                applySingleShift(grid, hid, delta);
            }
        }
        // ── Binder correction: align each strand segment individually ───
        const binderSet = new Set(binderHelices ?? []);
        if (binderSet.size > 0) {
            console.log(`[alignGridPrim] Aligning binder helices: [${Array.from(binderSet).sort((a, b) => a - b).join(', ')}]`);
            // Walk all strands to find crossovers INTO binder helices.
            // For each strand segment on a binder helix, compute the
            // per-segment shift from its crossover parent.
            const allNtIds = new Set();
            for (const [ntId] of grid.entries())
                allNtIds.add(ntId);
            const binderVisited = new Set();
            for (const [ntId] of grid.entries()) {
                if (binderVisited.has(ntId))
                    continue;
                const startNt = elements.get(ntId);
                if (!startNt || !(startNt instanceof Nucleotide))
                    continue;
                // Find 5' end
                let fivePrime = startNt;
                const walkBack = new Set();
                walkBack.add(fivePrime.id);
                while (true) {
                    const prev = fivePrime.n5;
                    if (!prev || !(prev instanceof Nucleotide))
                        break;
                    if (!allNtIds.has(prev.id))
                        break;
                    if (walkBack.has(prev.id))
                        break;
                    walkBack.add(prev.id);
                    fivePrime = prev;
                }
                const runs = [];
                let currentRun = null;
                let curr = fivePrime;
                const walkFwd = new Set();
                while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                    if (walkFwd.has(curr.id))
                        break;
                    walkFwd.add(curr.id);
                    binderVisited.add(curr.id);
                    const mark = grid.get(curr.id);
                    if (mark) {
                        if (currentRun && currentRun.helixId === mark.helixId) {
                            currentRun.ntIds.push(curr.id);
                        }
                        else {
                            currentRun = { helixId: mark.helixId, ntIds: [curr.id] };
                            runs.push(currentRun);
                        }
                    }
                    else {
                        currentRun = null;
                    }
                    const n3ref = curr.n3;
                    curr = (n3ref && n3ref instanceof Nucleotide) ? n3ref : null;
                }
                // For each run on a binder helix, find the crossover offset
                // from its adjacent non-binder run and align.
                for (let i = 0; i < runs.length; i++) {
                    const run = runs[i];
                    if (!binderSet.has(run.helixId))
                        continue;
                    // Look for the adjacent non-binder run to get the
                    // crossover offset. Check the run before and after.
                    let parentOffset = null;
                    let binderCrossoverOffset = null;
                    // Check previous run (crossover INTO binder)
                    if (i > 0 && !binderSet.has(runs[i - 1].helixId)) {
                        const prevRun = runs[i - 1];
                        const lastNtId = prevRun.ntIds[prevRun.ntIds.length - 1];
                        const lastMark = grid.get(lastNtId);
                        const firstMark = grid.get(run.ntIds[0]);
                        if (lastMark && firstMark) {
                            parentOffset = lastMark.offset;
                            binderCrossoverOffset = firstMark.offset;
                        }
                    }
                    // Check next run (crossover OUT of binder) if no prev
                    if (parentOffset === null && i < runs.length - 1 && !binderSet.has(runs[i + 1].helixId)) {
                        const nextRun = runs[i + 1];
                        const firstNtId = nextRun.ntIds[0];
                        const firstMark = grid.get(firstNtId);
                        const lastMark = grid.get(run.ntIds[run.ntIds.length - 1]);
                        if (firstMark && lastMark) {
                            parentOffset = firstMark.offset;
                            binderCrossoverOffset = lastMark.offset;
                        }
                    }
                    if (parentOffset !== null && binderCrossoverOffset !== null) {
                        const segmentShift = parentOffset - binderCrossoverOffset;
                        if (segmentShift !== 0) {
                            for (const ntId of run.ntIds) {
                                const mark = grid.get(ntId);
                                if (mark)
                                    mark.offset += segmentShift;
                            }
                        }
                    }
                }
            }
        }
        // ── Final normalization ──────────────────────────────────────────
        let globalMin2 = Infinity;
        for (const [, mark] of grid.entries()) {
            if (mark.offset < globalMin2)
                globalMin2 = mark.offset;
        }
        if (globalMin2 !== 0 && globalMin2 !== Infinity) {
            for (const [, mark] of grid.entries()) {
                mark.offset -= globalMin2;
            }
        }
        // ── Binder post-pass (after alignGridPrim) ─────────────────────
        const binderPostSet = new Set(binderHelices ?? []);
        if (binderPostSet.size > 0) {
            const binderList = Array.from(binderPostSet).sort((a, b) => a - b);
            console.log(`[alignGridPrim] Binder helices noted: [${binderList.join(', ')}]`);
            const allNtIds = new Set();
            for (const [ntId] of grid.entries())
                allNtIds.add(ntId);
            const visited = new Set();
            const binderRuns = [];
            const ntToRun = new Map();
            // Rebuild strand runs and keep only runs on binder helices.
            for (const [ntId] of grid.entries()) {
                if (visited.has(ntId))
                    continue;
                const startNt = elements.get(ntId);
                if (!startNt || !(startNt instanceof Nucleotide))
                    continue;
                let fivePrime = startNt;
                const walkBack = new Set();
                walkBack.add(fivePrime.id);
                while (true) {
                    const prev = fivePrime.n5;
                    if (!prev || !(prev instanceof Nucleotide))
                        break;
                    if (!allNtIds.has(prev.id))
                        break;
                    if (walkBack.has(prev.id))
                        break;
                    walkBack.add(prev.id);
                    fivePrime = prev;
                }
                let curr = fivePrime;
                const walkFwd = new Set();
                let currentRunNtIds = [];
                let currentHelix = null;
                while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                    if (walkFwd.has(curr.id))
                        break;
                    walkFwd.add(curr.id);
                    visited.add(curr.id);
                    const mark = grid.get(curr.id);
                    if (!mark) {
                        if (currentHelix !== null && binderPostSet.has(currentHelix) && currentRunNtIds.length > 0) {
                            const runIdx = binderRuns.length;
                            binderRuns.push({ helixId: currentHelix, ntIds: currentRunNtIds });
                            currentRunNtIds.forEach((id) => ntToRun.set(id, runIdx));
                        }
                        currentRunNtIds = [];
                        currentHelix = null;
                    }
                    else if (currentHelix === mark.helixId || currentHelix === null) {
                        currentHelix = mark.helixId;
                        currentRunNtIds.push(curr.id);
                    }
                    else {
                        if (binderPostSet.has(currentHelix) && currentRunNtIds.length > 0) {
                            const runIdx = binderRuns.length;
                            binderRuns.push({ helixId: currentHelix, ntIds: currentRunNtIds });
                            currentRunNtIds.forEach((id) => ntToRun.set(id, runIdx));
                        }
                        currentHelix = mark.helixId;
                        currentRunNtIds = [curr.id];
                    }
                    const n3ref = curr.n3;
                    curr = (n3ref && n3ref instanceof Nucleotide) ? n3ref : null;
                }
                if (currentHelix !== null && binderPostSet.has(currentHelix) && currentRunNtIds.length > 0) {
                    const runIdx = binderRuns.length;
                    binderRuns.push({ helixId: currentHelix, ntIds: currentRunNtIds });
                    currentRunNtIds.forEach((id) => ntToRun.set(id, runIdx));
                }
            }
            // Force all binder substrands to forward direction.
            for (const run of binderRuns) {
                for (const ntId of run.ntIds) {
                    const mark = grid.get(ntId);
                    if (mark)
                        mark.direction = 'forward';
                }
            }
            const runCenter = (run) => {
                let sum = 0;
                let count = 0;
                for (const ntId of run.ntIds) {
                    const mark = grid.get(ntId);
                    if (!mark)
                        continue;
                    sum += mark.offset;
                    count++;
                }
                return count > 0 ? sum / count : 0;
            };
            const shiftRunBy = (run, delta) => {
                if (delta === 0)
                    return;
                for (const ntId of run.ntIds) {
                    const mark = grid.get(ntId);
                    if (mark)
                        mark.offset += delta;
                }
            };
            // Keep pushing the front run by +5 until no binder run overlaps remain.
            let guard = 0;
            while (guard++ < 2000) {
                const overlapBuckets = new Map();
                for (const run of binderRuns) {
                    for (const ntId of run.ntIds) {
                        const mark = grid.get(ntId);
                        if (!mark || !binderPostSet.has(mark.helixId))
                            continue;
                        const key = `${mark.helixId}|${mark.offset}`;
                        if (!overlapBuckets.has(key))
                            overlapBuckets.set(key, new Set());
                        const runIdx = ntToRun.get(ntId);
                        if (runIdx !== undefined)
                            overlapBuckets.get(key).add(runIdx);
                    }
                }
                let moved = false;
                for (const [, runSet] of overlapBuckets.entries()) {
                    if (runSet.size <= 1)
                        continue;
                    let frontRun = null;
                    let frontCenter = -Infinity;
                    for (const runIdx of runSet.values()) {
                        const run = binderRuns[runIdx];
                        if (!run)
                            continue;
                        const c = runCenter(run);
                        if (c > frontCenter) {
                            frontCenter = c;
                            frontRun = run;
                        }
                    }
                    if (frontRun) {
                        shiftRunBy(frontRun, 5);
                        moved = true;
                        break;
                    }
                }
                if (!moved)
                    break;
            }
        }
        console.log(`[alignGridPrim] Aligned ${helixList.length} helices. Shifts:`, Object.fromEntries(Array.from(cumulativeShift.entries()).sort((a, b) => a[0] - b[0])));
        toscad.validateGrid(grid);
        const conflictsAfter = countGridConflicts(grid);
        if (conflictsAfter > 0) {
            throw new Error(`[alignGridPrim] Overlaps remain after post-pass: ${conflictsAfter}`);
        }
        return { shifts: cumulativeShift };
    }
    toscad.alignGridPrim = alignGridPrim;
    /**
     * alignMergedGroups — merged-helix group alignment pass.
     *
     * Companion to alignGridPrim, and the only place that knows about merge
     * provenance. setGrid is responsible purely for *placement*: it copies
     * each merged helix's origin-groups into the grid with their pre-merge
     * (separate-helix) offsets, which is collision-free but unaligned. This
     * function performs the alignment.
     *
     * Each merged helix is temporarily split into pseudo-helices — one per
     * origin-group — and alignGridPrim is run over the whole grid, so each
     * group is anchored by its own crossovers to the full (already-placed)
     * lattice. The groups were mutually disjoint on the offset axis at merge
     * time (disjointness is a merge precondition), so once aligned they can be
     * folded back under the real helixId without overlap. Afterwards,
     * crossovers into the merged helix observe ~0 shift, so the pipeline's own
     * alignGridPrim pass leaves the merged helix in place instead of shifting
     * it as a unit and misaligning one of the groups.
     *
     * Must be called immediately after setGrid and BEFORE directionAlign2:
     * directionAlign2 flips helices and rewrites offsets, which changes what
     * collectShiftObservations sees.
     *
     * No-op unless at least one merged helix has 2+ origin-groups, so it is
     * safe to call unconditionally on every pipeline iteration.
     */
    function alignMergedGroups(grid, mergedGroups, binderHelices) {
        if (!mergedGroups || mergedGroups.size === 0)
            return;
        let maxHelixId = -1;
        for (const [, m] of grid) {
            if (m.helixId > maxHelixId)
                maxHelixId = m.helixId;
        }
        // pseudo helixId -> real (merged) helixId
        const pseudoToReal = new Map();
        let nextPseudo = maxHelixId + 1;
        mergedGroups.forEach((groups, helixId) => {
            if (groups.length < 2)
                return;
            // Group 0 keeps the real helixId; the rest become pseudo-helices.
            for (let gi = 1; gi < groups.length; gi++) {
                const pseudo = nextPseudo++;
                pseudoToReal.set(pseudo, helixId);
                for (const ntId of groups[gi]) {
                    const m = grid.get(ntId);
                    if (m)
                        m.helixId = pseudo;
                }
            }
        });
        // Every merged helix has a single origin-group: nothing to align, and
        // no pseudo-helices were created, so leave the grid untouched.
        if (pseudoToReal.size === 0)
            return;
        alignGridPrim(grid, binderHelices);
        for (const [, m] of grid) {
            const real = pseudoToReal.get(m.helixId);
            if (real !== undefined)
                m.helixId = real;
        }
        // Inter-group conflict resolution: NEVER distort a helix.
        // An origin-group is moved ONLY as a whole unit, and ONLY
        // its offsets change — direction/orientation is never
        // touched. For every group that overlaps others on the
        // same helix, compute the minimum-|delta| shift (searching
        // BOTH +/- directions) that leaves it disjoint from every
        // other group on that helix, then apply it in a single
        // step. Searching both directions matters: alignGridPrim
        // can shift a group either way relative to its siblings,
        // and the "correct" resolution is whichever direction has
        // the group already close to disjoint. A pure +delta nudge
        // can either explode (walking a small group through a
        // large one) or leave the crash-later state where the
        // required delta exceeds the sweep budget.
        // Shift every mark of one origin-group by delta (offset only).
        const shiftGroup = (helixId, gi, delta) => {
            if (delta === 0)
                return;
            const group = mergedGroups.get(helixId)?.[gi];
            if (!group)
                return;
            for (const ntId of group) {
                const m = grid.get(ntId);
                if (m && m.helixId === helixId)
                    m.offset += delta;
            }
        };
        // For (helixId, gi), find the min-|delta| shift such that
        // no mark of gi collides (same direction + same offset)
        // with any mark of ANY OTHER group on the same helix.
        // Returns 0 if already disjoint, or null if no delta
        // within the search bound resolves the conflict.
        const findCleanShift = (helixId, gi) => {
            const groups = mergedGroups.get(helixId);
            if (!groups)
                return 0;
            const myGroup = groups[gi];
            if (!myGroup || myGroup.length === 0)
                return 0;
            const mySet = new Set(myGroup);
            // My occupied cells, indexed by direction.
            const myCells = new Map();
            let myMin = Infinity, myMax = -Infinity;
            for (const ntId of myGroup) {
                const m = grid.get(ntId);
                if (!m || m.helixId !== helixId)
                    continue;
                let s = myCells.get(m.direction);
                if (!s) {
                    s = new Set();
                    myCells.set(m.direction, s);
                }
                s.add(m.offset);
                if (m.offset < myMin)
                    myMin = m.offset;
                if (m.offset > myMax)
                    myMax = m.offset;
            }
            if (myCells.size === 0)
                return 0;
            // Cells occupied by OTHER groups on the same helix.
            const otherCells = new Map();
            let otherMin = Infinity, otherMax = -Infinity;
            for (const [ntId, m] of grid) {
                if (m.helixId !== helixId)
                    continue;
                if (mySet.has(ntId))
                    continue;
                let s = otherCells.get(m.direction);
                if (!s) {
                    s = new Set();
                    otherCells.set(m.direction, s);
                }
                s.add(m.offset);
                if (m.offset < otherMin)
                    otherMin = m.offset;
                if (m.offset > otherMax)
                    otherMax = m.offset;
            }
            if (otherCells.size === 0)
                return 0;
            const conflictsAt = (delta) => {
                for (const [dir, myOffs] of myCells) {
                    const others = otherCells.get(dir);
                    if (!others)
                        continue;
                    for (const off of myOffs) {
                        if (others.has(off + delta))
                            return true;
                    }
                }
                return false;
            };
            if (!conflictsAt(0))
                return 0;
            // Bound: any collision requires (my_off + delta) to
            // land on an other-off, so |delta| ≤ (other-span +
            // my-span). Add a small pad so we can step JUST past
            // the far edge in either direction.
            const mySpan = (myMax - myMin) || 0;
            const otherSpan = (otherMax - otherMin) || 0;
            const maxSearch = Math.max(mySpan + otherSpan + 8, 32);
            for (let mag = 1; mag <= maxSearch; mag++) {
                if (!conflictsAt(mag))
                    return mag;
                if (!conflictsAt(-mag))
                    return -mag;
            }
            return null;
        };
        // Iterate: resolve one group at a time, using the current
        // grid state (so gi=2 sees where gi=1 landed). Each
        // successful shift makes that group globally disjoint on
        // its helix, so the outer loop terminates when a full
        // pass moves nothing. Guard cap protects against
        // pathological chains between helices.
        let movedAny = false;
        let unresolvable = 0;
        let sweepGuard = 0;
        while (sweepGuard++ < 200) {
            let didMove = false;
            mergedGroups.forEach((groups, helixId) => {
                for (let gi = 1; gi < groups.length; gi++) {
                    const delta = findCleanShift(helixId, gi);
                    if (delta === null) {
                        unresolvable++;
                        continue;
                    }
                    if (delta === 0)
                        continue;
                    shiftGroup(helixId, gi, delta);
                    didMove = true;
                    movedAny = true;
                }
            });
            if (!didMove)
                break;
        }
        if (movedAny) {
            console.warn(`[alignMergedGroups] Merged-helix inter-group overlap resolved by min-|delta| shifts` +
                ` (sweeps=${sweepGuard - 1}${unresolvable > 0 ? `, unresolved=${unresolvable}` : ''}).`);
        }
        if (unresolvable > 0) {
            console.error(`[alignMergedGroups] ${unresolvable} merged-helix group(s) had no collision-free placement` +
                ` within the search bound. Grid will fail validation downstream.`);
        }
    }
    toscad.alignMergedGroups = alignMergedGroups;
    /**
     * Calculates absolute grid coordinates from local helix-to-helix angles.
     *
     * Phase 1: Build a strict lattice spanning tree from helix 0 using
     * weighted BFS (top-3 children only at each node).
     *
     * Phase 2: Place deferred/artefact helices in nearest open coordinates
     * around their parent once the phase-1 core is locked.
     */
    function calculateGlobalPositions(networkMap, crossoverWeights, options, lattice = 'honeycomb') {
        const latticeType = resolveLatticeKind(lattice);
        const HONEYCOMB_ANGLES = [0, 120, 240];
        const SQUARE_ANGLES = [0, 90, 180, 270];
        const ANGLES = latticeType === 'square'
            ? SQUARE_ANGLES.slice()
            : HONEYCOMB_ANGLES.slice();
        const maxChildrenPerNode = ANGLES.length;
        const HONEYCOMB_STEP_BY_PARITY = {
            even: {
                0: { dCol: 1, dRow: 0 },
                120: { dCol: -1, dRow: 0 },
                240: { dCol: 0, dRow: -1 }
            },
            odd: {
                0: { dCol: 1, dRow: 0 },
                120: { dCol: 0, dRow: 1 },
                240: { dCol: -1, dRow: 0 }
            }
        };
        const SQUARE_STEPS = {
            0: { dCol: 1, dRow: 0 },
            90: { dCol: 0, dRow: 1 },
            180: { dCol: -1, dRow: 0 },
            270: { dCol: 0, dRow: -1 }
        };
        const normalizeAngle = (angle) => ((angle % 360) + 360) % 360;
        const angleDistance = (a, b) => {
            const diff = Math.abs(normalizeAngle(a) - normalizeAngle(b));
            return Math.min(diff, 360 - diff);
        };
        const parityAt = (coord) => (((coord.col + coord.row) & 1) === 0 ? 'even' : 'odd');
        const getStep = (coord, angle) => {
            if (latticeType === 'square') {
                return SQUARE_STEPS[angle];
            }
            return HONEYCOMB_STEP_BY_PARITY[parityAt(coord)][angle];
        };
        const snapToLatticeAngle = (angle) => {
            let best = ANGLES[0];
            let bestDist = Number.POSITIVE_INFINITY;
            for (const candidate of ANGLES) {
                const dist = angleDistance(angle, candidate);
                if (dist < bestDist || (dist === bestDist && candidate < best)) {
                    best = candidate;
                    bestDist = dist;
                }
            }
            return best;
        };
        const latticeAngleFromDelta = (coord, dCol, dRow) => {
            for (const angle of ANGLES) {
                const delta = getStep(coord, angle);
                if (delta.dCol === dCol && delta.dRow === dRow)
                    return angle;
            }
            return null;
        };
        const weightMap = new Map();
        const ensureWeightRow = (helixId) => {
            if (!weightMap.has(helixId))
                weightMap.set(helixId, new Map());
            return weightMap.get(helixId);
        };
        const setWeight = (from, to, weight) => {
            ensureWeightRow(from).set(to, weight);
        };
        if (crossoverWeights instanceof Map) {
            for (const [from, row] of crossoverWeights.entries()) {
                for (const [to, weight] of row.entries()) {
                    setWeight(from, to, weight);
                }
            }
        }
        else if (Array.isArray(crossoverWeights)) {
            for (const edge of crossoverWeights) {
                if ('fromHelix' in edge && 'toHelix' in edge) {
                    setWeight(edge.fromHelix, edge.toHelix, edge.weight);
                }
                else {
                    setWeight(edge.from, edge.to, edge.weight);
                }
            }
        }
        const getWeight = (a, b) => {
            const direct = weightMap.get(a)?.get(b);
            if (direct !== undefined)
                return direct;
            const reverse = weightMap.get(b)?.get(a);
            if (reverse !== undefined)
                return reverse;
            return 1;
        };
        const allHelixIds = new Set();
        allHelixIds.add(0);
        for (const [from, row] of networkMap.entries()) {
            allHelixIds.add(from);
            for (const to of row.keys())
                allHelixIds.add(to);
        }
        for (const [from, row] of weightMap.entries()) {
            allHelixIds.add(from);
            for (const to of row.keys())
                allHelixIds.add(to);
        }
        const positions = new Map();
        const occupied = new Map();
        const globalRotationOffsets = new Map();
        const deferredQueue = [];
        const deferredSeen = new Set();
        const keyOf = (coord) => `${coord.col},${coord.row}`;
        const enqueueDeferred = (parentId, helixId, inheritedOffset) => {
            const key = `${parentId}|${helixId}`;
            if (deferredSeen.has(key) || positions.has(helixId))
                return;
            deferredSeen.add(key);
            deferredQueue.push({ parentId, helixId, inheritedOffset });
        };
        const maxSearchRadius = Math.max(1, options?.maxSearchRadius ?? 256);
        const latticeNeighbors = (coord) => {
            return ANGLES.map((angle) => {
                const step = getStep(coord, angle);
                return {
                    col: coord.col + step.dCol,
                    row: coord.row + step.dRow
                };
            });
        };
        const findNearestOpen = (anchor) => {
            if (!occupied.has(keyOf(anchor)))
                return { col: anchor.col, row: anchor.row };
            const visited = new Set();
            const queue = [{ coord: anchor, dist: 0 }];
            visited.add(keyOf(anchor));
            let qIdx = 0;
            while (qIdx < queue.length) {
                const { coord, dist } = queue[qIdx++];
                if (dist >= maxSearchRadius)
                    continue;
                const neighbors = latticeNeighbors(coord).sort((a, b) => {
                    const da = Math.abs(a.col - anchor.col) + Math.abs(a.row - anchor.row);
                    const db = Math.abs(b.col - anchor.col) + Math.abs(b.row - anchor.row);
                    if (da !== db)
                        return da - db;
                    if (a.col !== b.col)
                        return a.col - b.col;
                    return a.row - b.row;
                });
                for (const next of neighbors) {
                    const nextKey = keyOf(next);
                    if (visited.has(nextKey))
                        continue;
                    visited.add(nextKey);
                    if (!occupied.has(nextKey))
                        return next;
                    queue.push({ coord: next, dist: dist + 1 });
                }
            }
            const fallback = { col: anchor.col + maxSearchRadius + 1, row: anchor.row };
            while (occupied.has(keyOf(fallback))) {
                fallback.col += 1;
            }
            return fallback;
        };
        const computeChildOffset = (parentId, childId, parentCoord, childCoord, parentLocalAngle, snappedForwardAngle) => {
            const dColBack = parentCoord.col - childCoord.col;
            const dRowBack = parentCoord.row - childCoord.row;
            const backGlobal = latticeAngleFromDelta(childCoord, dColBack, dRowBack);
            const desiredBackGlobalAngle = backGlobal !== null
                ? backGlobal
                : normalizeAngle(snappedForwardAngle + 180);
            const childLocalBack = networkMap.get(childId)?.get(parentId);
            if (typeof childLocalBack === 'number') {
                return normalizeAngle(desiredBackGlobalAngle - childLocalBack);
            }
            const inferredBackLocal = normalizeAngle(parentLocalAngle + 180);
            return normalizeAngle(desiredBackGlobalAngle - inferredBackLocal);
        };
        const processNode = (node, targetQueue) => {
            const parentCoord = positions.get(node.helixId);
            if (!parentCoord)
                return;
            const localEdges = networkMap.get(node.helixId);
            if (!localEdges || localEdges.size === 0)
                return;
            const candidates = Array.from(localEdges.entries())
                .filter(([neighborId]) => neighborId !== node.helixId && !positions.has(neighborId))
                .map(([neighborId, localAngle]) => ({
                neighborId,
                localAngle: normalizeAngle(localAngle),
                weight: getWeight(node.helixId, neighborId)
            }))
                .sort((a, b) => {
                if (b.weight !== a.weight)
                    return b.weight - a.weight;
                return a.neighborId - b.neighborId;
            });
            const selected = candidates.slice(0, maxChildrenPerNode);
            const overflow = candidates.slice(maxChildrenPerNode);
            for (const item of overflow) {
                enqueueDeferred(node.helixId, item.neighborId, node.offset);
            }
            for (const item of selected) {
                if (positions.has(item.neighborId))
                    continue;
                const predictedGlobal = normalizeAngle(item.localAngle + node.offset);
                const snappedAngle = snapToLatticeAngle(predictedGlobal);
                const step = getStep(parentCoord, snappedAngle);
                const childCoord = {
                    col: parentCoord.col + step.dCol,
                    row: parentCoord.row + step.dRow
                };
                const cellKey = keyOf(childCoord);
                const occupant = occupied.get(cellKey);
                if (occupant !== undefined && occupant !== item.neighborId) {
                    enqueueDeferred(node.helixId, item.neighborId, node.offset);
                    continue;
                }
                positions.set(item.neighborId, childCoord);
                occupied.set(cellKey, item.neighborId);
                const childOffset = computeChildOffset(node.helixId, item.neighborId, parentCoord, childCoord, item.localAngle, snappedAngle);
                globalRotationOffsets.set(item.neighborId, childOffset);
                targetQueue.push({ helixId: item.neighborId, offset: childOffset });
            }
        };
        // Phase 1: strict weighted BFS spanning tree from helix 0.
        positions.set(0, { col: 0, row: 0 });
        occupied.set('0,0', 0);
        globalRotationOffsets.set(0, 0);
        const mainQueue = [{ helixId: 0, offset: 0 }];
        let mainIdx = 0;
        while (mainIdx < mainQueue.length) {
            processNode(mainQueue[mainIdx++], mainQueue);
        }
        // Phase 2: place deferred artefacts nearest to their parent.
        const runDeferredSubBfs = options?.runDeferredSubBfs ?? false;
        const deferredSubQueue = [];
        let deferredIdx = 0;
        while (deferredIdx < deferredQueue.length) {
            const item = deferredQueue[deferredIdx++];
            if (positions.has(item.helixId))
                continue;
            const parentCoord = positions.get(item.parentId) ?? positions.get(0) ?? { col: 0, row: 0 };
            const coord = findNearestOpen(parentCoord);
            positions.set(item.helixId, coord);
            occupied.set(keyOf(coord), item.helixId);
            const inherited = globalRotationOffsets.get(item.parentId);
            const chosenOffset = inherited !== undefined ? inherited : item.inheritedOffset;
            globalRotationOffsets.set(item.helixId, chosenOffset);
            if (runDeferredSubBfs) {
                deferredSubQueue.push({ helixId: item.helixId, offset: chosenOffset });
            }
        }
        if (runDeferredSubBfs) {
            let subIdx = 0;
            while (subIdx < deferredSubQueue.length) {
                processNode(deferredSubQueue[subIdx++], deferredSubQueue);
            }
        }
        // Ensure every helix in the input graph gets a coordinate.
        const rootCoord = positions.get(0) ?? { col: 0, row: 0 };
        const sortedHelixIds = Array.from(allHelixIds).sort((a, b) => a - b);
        for (const helixId of sortedHelixIds) {
            if (positions.has(helixId))
                continue;
            const coord = findNearestOpen(rootCoord);
            positions.set(helixId, coord);
            occupied.set(keyOf(coord), helixId);
        }
        const result = new Map();
        for (const [helixId, coord] of positions.entries()) {
            result.set(helixId, [coord.col, coord.row]);
        }
        return result;
    }
    toscad.calculateGlobalPositions = calculateGlobalPositions;
    const RENUMBER_JUMP_PENALTY = 1000;
    function applyHelixRenumber(helices, grid, helixPos, remap) {
        const n = helices.length;
        const newHelices = new Array(n);
        for (let oldId = 0; oldId < n; oldId++) {
            const newId = remap.get(oldId);
            const slot = (newId !== undefined && newId >= 0 && newId < n) ? newId : oldId;
            newHelices[slot] = helices[oldId];
        }
        for (const mark of grid.values()) {
            const newId = remap.get(mark.helixId);
            if (newId !== undefined)
                mark.helixId = newId;
        }
        const newHelixPos = new Map();
        for (const [oldId, pos] of helixPos.entries()) {
            const newId = remap.get(oldId);
            newHelixPos.set(newId !== undefined ? newId : oldId, pos);
        }
        return { helices: newHelices, helixPos: newHelixPos };
    }
    toscad.applyHelixRenumber = applyHelixRenumber;
    // ─────────────────────────────────────────────────────────────────────
    //  renumberHelicesGNN — Greedy Nearest-Neighbor renumber
    //
    //  Same input/output contract as renumberHelices, but the path is built
    //  by walking from helix 0 and at each step extending to the cheapest
    //  unvisited helix under the same edgeCost (jump penalty + Euclidean
    //  grid distance). No MST, no DFS preorder, no 2-Opt — the linearization
    //  step itself respects locality, so consecutive IDs prefer crossover
    //  neighbors and otherwise fall back to nearest physical neighbor.
    //
    //  Trade-off vs MST+DFS+2Opt: a greedy tour can leave one or two long
    //  "tail" hops when the last unvisited node sits far from the current
    //  endpoint. In exchange we avoid the mid-path teleports that come from
    //  popping out of MST subtrees.
    // ─────────────────────────────────────────────────────────────────────
    function renumberHelicesGNN(grid, helixPos, lattice = 'honeycomb', binderHelices) {
        // ── 1. Adjacency + cost helpers (identical to renumberHelices) ──
        const { crossovers, helixIds } = toscad.collectCrossovers(grid);
        for (const id of helixPos.keys())
            helixIds.add(id);
        const nodes = Array.from(helixIds).sort((a, b) => a - b);
        const n = nodes.length;
        const idx = new Map();
        nodes.forEach((id, i) => idx.set(id, i));
        const isSquare = (lattice ?? '').toLowerCase() === 'square';
        const world = nodes.map(id => {
            const pos = helixPos.get(id) ?? [0, 0];
            const v = isSquare
                ? scadnano.squareToWorld(pos[0], pos[1])
                : scadnano.honeycombToWorld(pos[0], pos[1]);
            return [v.x, v.y];
        });
        // Anchor selection: lowest row (y), tie-broken by lowest col (x).
        // This pins helix 0 to the top-left corner of the layout for both
        // lattices. Because helix 0 also seeds calculateGlobalPositions'
        // Phase-1 BFS on the NEXT iteration, we must not let it be a binder:
        // binders don't participate in the offset MST or lattice-detection
        // passes, so a binder anchor would misroot the whole layout. When
        // the geometric top-left IS a binder, walk the same (row, col)
        // ordering to the next candidate that isn't. If every helix is a
        // binder (degenerate — no core structure), fall back to the pure
        // geometric top-left.
        const binderNodeSet = new Set();
        if (Array.isArray(binderHelices)) {
            for (const hid of binderHelices) {
                const i = idx.get(hid);
                if (i !== undefined)
                    binderNodeSet.add(i);
            }
        }
        const posOrder = Array.from({ length: n }, (_, i) => i).sort((a, b) => {
            const pa = helixPos.get(nodes[a]) ?? [0, 0];
            const pb = helixPos.get(nodes[b]) ?? [0, 0];
            if (pa[1] !== pb[1])
                return pa[1] - pb[1];
            return pa[0] - pb[0];
        });
        let startNode = posOrder[0] ?? 0;
        let skippedBinderAnchors = 0;
        if (binderNodeSet.size > 0) {
            let chosen = null;
            for (const cand of posOrder) {
                if (!binderNodeSet.has(cand)) {
                    chosen = cand;
                    break;
                }
                skippedBinderAnchors++;
            }
            if (chosen !== null) {
                startNode = chosen;
            }
            else {
                // Every helix is a binder — nothing to pick from. Log so the
                // caller can see something is off with the input.
                console.warn(`[renumberHelicesGNN] every helix is flagged as a binder; ` +
                    `falling back to geometric top-left (helixId=${nodes[startNode]}).`);
            }
        }
        const startPos = helixPos.get(nodes[startNode]) ?? [0, 0];
        console.log(`[renumberHelicesGNN] anchor node idx=${startNode} ` +
            `(helixId=${nodes[startNode]}, col=${startPos[0]}, row=${startPos[1]})` +
            (skippedBinderAnchors > 0
                ? ` — skipped ${skippedBinderAnchors} binder(s) higher in top-left ordering`
                : ''));
        const dist = (a, b) => {
            const [ax, ay] = world[a];
            const [bx, by] = world[b];
            return Math.hypot(ax - bx, ay - by);
        };
        const adj = Array.from({ length: n }, () => new Array(n).fill(false));
        for (const [from, inner] of crossovers.entries()) {
            const i = idx.get(from);
            if (i === undefined)
                continue;
            for (const [to, counts] of inner.entries()) {
                const j = idx.get(to);
                if (j === undefined)
                    continue;
                if ((counts.sameWalk + counts.diffWalk) > 0) {
                    adj[i][j] = true;
                    adj[j][i] = true;
                }
            }
        }
        const isJump = (a, b) => !adj[a][b];
        const edgeCost = (a, b) => (isJump(a, b) ? RENUMBER_JUMP_PENALTY : 0) + dist(a, b);
        if (n <= 1) {
            const order = nodes.slice();
            const remap = new Map();
            order.forEach((id, i) => remap.set(id, i));
            console.log(`[renumberHelicesGNN] n=${n}, nothing to renumber.`);
            return {
                order, remap,
                stats: {
                    n, components: n, lowerBoundJumps: 0,
                    initialJumps: 0, initialDistance: 0,
                    afterPhase1Jumps: 0, finalJumps: 0, finalDistance: 0
                }
            };
        }
        // ── 2. Connected components (lower bound on jumps, log only) ────
        const compId = new Array(n).fill(-1);
        let nComps = 0;
        for (let s = 0; s < n; s++) {
            if (compId[s] !== -1)
                continue;
            const stack = [s];
            compId[s] = nComps;
            while (stack.length) {
                const u = stack.pop();
                for (let v = 0; v < n; v++) {
                    if (compId[v] !== -1)
                        continue;
                    if (adj[u][v]) {
                        compId[v] = nComps;
                        stack.push(v);
                    }
                }
            }
            nComps++;
        }
        const lowerBoundJumps = Math.max(0, nComps - 1);
        // ── 3. Greedy nearest-neighbor walk from anchor node ────────────
        const visited = new Array(n).fill(false);
        const path = [startNode];
        visited[startNode] = true;
        let current = startNode;
        for (let step = 1; step < n; step++) {
            let bestV = -1;
            let bestCost = Infinity;
            for (let v = 0; v < n; v++) {
                if (visited[v])
                    continue;
                const c = edgeCost(current, v);
                if (c < bestCost) {
                    bestCost = c;
                    bestV = v;
                }
            }
            if (bestV === -1)
                break; // disconnected universe; shouldn't happen
            visited[bestV] = true;
            path.push(bestV);
            current = bestV;
        }
        // ── 4. Path stats ───────────────────────────────────────────────
        const pathStats = (p) => {
            let jumps = 0, distance = 0;
            for (let i = 0; i + 1 < p.length; i++) {
                if (isJump(p[i], p[i + 1]))
                    jumps++;
                distance += dist(p[i], p[i + 1]);
            }
            return { jumps, distance };
        };
        const final = pathStats(path);
        // ── 5. Build remap and report ──────────────────────────────────
        const order = path.map(i => nodes[i]);
        const remap = new Map();
        order.forEach((id, i) => remap.set(id, i));
        const stats = {
            n,
            components: nComps,
            lowerBoundJumps,
            initialJumps: final.jumps,
            initialDistance: final.distance,
            afterPhase1Jumps: final.jumps,
            finalJumps: final.jumps,
            finalDistance: final.distance,
        };
        console.log(`[renumberHelicesGNN] n=${n}, components=${nComps}, ` +
            `jumps lower bound=${lowerBoundJumps} | ` +
            `final=[${final.jumps} jumps, dist=${final.distance.toFixed(2)}]`);
        return { order, remap, stats };
    }
    toscad.renumberHelicesGNN = renumberHelicesGNN;
})(toscad || (toscad = {}));
