"use strict";
/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
/*
    We use "Col" and "Row" to refer to the 2D grid coordinates of a helix in the lattice. Refer to the honeycomb and square lattice diagrams in
    https://github.com/UC-Davis-molecular-computing/scadnano.
    Note that (col, row) == (x,y). The following codebase uses (col, row) simply to follow the scadnano convention.

    Ideas that are not yet implemented:
    - Use major/minor groove asymmetries as a gate for helices, or deterministically labelling helix directions.
*/
var toscad;
(function (toscad) {
    // easy to use constant for later reference.
    toscad.resolveLatticeKind = (lattice) => (lattice ?? '').toLowerCase() === 'square' ? 'square' : 'honeycomb';
    const LATTICE_CONFIG = {
        honeycomb: {
            basesPerTurn: 10.5,
            phases: {
                0: [0.0, 3.5, 7.0, 10.5],
                // numbers found by using 215/145 asymmetry in major/minor groove and allowing void phase to break ties.
                // Example: 360/10.5 * (0.7292) + 215 = 240º.
                1: {
                    145: [0.7292, 4.229, 7.7292],
                    215: [2.7708, 6.2708, 9.7708]
                }
            },
            voidPhase: 5.25,
            tieEpsilon: 0.01
        },
        square: {
            basesPerTurn: 32 / 3,
            phases: {
                0: [0.0, 8 / 3, 16 / 3, 8.0, 32 / 3],
                1: {
                    145: [4.2963, 6.9629, 9.6296, 1.6296],
                    215: [6.3704, 9.0370, 11.7037, 3.7037]
                }
            },
            voidPhase: null,
            tieEpsilon: 0.01
        }
    };
    // helper function
    // Finds the local angle distribution of a given helix relative to its neighbors.
    // TODO: Potentially allow it to output more than 1 angle per neighbor...
    function getAngleHelix(grid, helices, helixId, lattice) {
        void helices;
        const result = new Map();
        const latticeType = toscad.resolveLatticeKind(lattice);
        const latticeConfig = LATTICE_CONFIG[latticeType];
        const nearestPhase = (value, parity) => {
            const phaseEntry = latticeConfig.phases[parity];
            let candidates;
            if (parity === 0) {
                candidates = phaseEntry.map(v => ({ value: v, bucket: 215 }));
            }
            else {
                const split = phaseEntry;
                candidates = [
                    ...split[145].map(v => ({ value: v, bucket: 145 })),
                    ...split[215].map(v => ({ value: v, bucket: 215 }))
                ];
            }
            return candidates.reduce((prev, curr) => {
                const distPrev = Math.abs(prev.value - value);
                const distCurr = Math.abs(curr.value - value);
                if (parity === 1
                    && latticeConfig.voidPhase !== null
                    && Math.abs(distPrev - distCurr) < latticeConfig.tieEpsilon) {
                    if (prev.value === latticeConfig.voidPhase && curr.value !== latticeConfig.voidPhase)
                        return curr;
                    if (curr.value === latticeConfig.voidPhase && prev.value !== latticeConfig.voidPhase)
                        return prev;
                }
                return distCurr < distPrev ? curr : prev;
            }, candidates[0]);
        };
        // collect crossover informations, where it happens and which direction -> direction tells us whether crossover is going to be at current strand or it's pair.
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
        // local map for crossovers between current helix to adjacent helix. This is used to find the most common angle between two helices.
        const groupedByNeighbor = new Map();
        for (const crossover of hubCrossovers) {
            if (!groupedByNeighbor.has(crossover.adj_helix)) {
                groupedByNeighbor.set(crossover.adj_helix, []);
            }
            groupedByNeighbor.get(crossover.adj_helix).push(crossover);
        }
        const neighbors = Array.from(groupedByNeighbor.keys()).sort((a, b) => a - b);
        // should always have at least one neighbor... this is just a sanity check.
        if (!neighbors.length)
            return result;
        const pairTallies = new Map();
        // iterate over all unique pairs to find the most common angle between them.
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
                        const matchedPhase = nearestPhase(phase, y);
                        // `bucket` labels WHERE the groove sits (145 or 215), so the
                        // correction that rotates it back onto a lattice direction is
                        // the complement, 360 - bucket (145 <-> 215). Adding the label
                        // itself double-counts and lands every antiparallel (y=1) edge
                        // exactly one lattice step off (90 deg square / 120 deg honeycomb).
                        const grooveCorrection = 360 - matchedPhase.bucket;
                        const angleRaw = (360 / latticeConfig.basesPerTurn) * matchedPhase.value + (y * grooveCorrection);
                        const relativeAngle = (Math.round(angleRaw % 360) + 360) % 360;
                        // bucket that contains all angles, such as 90, 120, 215, etc. It is lattice-agnostic.
                        bucket.set(relativeAngle, (bucket.get(relativeAngle) ?? 0) + 1);
                    }
                }
                pairTallies.set(pairKey, bucket);
            }
        }
        // the mode angle for each pair of neighbors.
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
    // TODO: Theres a todo inside, go look for it.
    function detectLatticeKind(grid, binderHelices = []) {
        // add binders into a set so we can ignore them.
        const binderSet = new Set(binderHelices);
        // helixId -> ntId -> point. Map-by-ntId dedupes nucleotides that appear as both the "to" of one crossover and the "from" of another.
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
        for (const [, inner] of pointsByHelix.entries()) {
            if (inner.size < 2) {
                continue;
            }
            const sorted = Array.from(inner.values()).sort((a, b) => a.offset - b.offset);
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
            }
        }
        const acceptedSizes = Array.from(acceptedByGap.keys()).sort((a, b) => a - b);
        for (const size of acceptedSizes) {
            const { hc, sq } = acceptedByGap.get(size);
            const parts = [];
            if (hc)
                parts.push(`HC=${hc}`);
            if (sq)
                parts.push(`SQ=${sq}`);
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
    toscad.findNearestOpenPos = findNearestOpenPos;
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
        let RENUMBER_JUMP_PENALTY = 100; // Tuned for typical honeycomb/square spacing
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
