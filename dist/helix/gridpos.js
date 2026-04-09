/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
var toscad;
(function (toscad) {
    const resolveLatticeKind = (lattice) => (lattice ?? '').toLowerCase() === 'square' ? 'square' : 'honeycomb';
    const LATTICE_CONFIG = {
        honeycomb: {
            basesPerTurn: 10.5,
            phases: {
                0: [0.0, 3.5, 7.0, 10.5],
                // 1: [1.75, 5.25, 8.75]
                // 1: [2.77, 6.27, 9.77]
                1: [0.7292, 4.229, 7.7292]
            },
            voidPhase: 5.25,
            tieEpsilon: 0.01
        },
        square: {
            basesPerTurn: 32 / 3,
            phases: {
                0: [0.0, 8 / 3, 16 / 3, 8.0, 32 / 3],
                1: [4 / 3, 4.0, 20 / 3, 28 / 3]
            },
            voidPhase: null,
            tieEpsilon: 0.01
        }
    };
    // helper function
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
        for (const crossover of toscad.crossoverNts(grid)) {
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
    function anglecomb(grid, helices, lattice = 'honeycomb', angleMap = getAngles(grid, helices, lattice)) {
        let networkMap = angleMap;
        const mergedPairs = [];
        const mergeHelixInto = (keepHelix, mergedHelix) => {
            for (const [, mark] of grid.entries()) {
                if (mark.helixId === mergedHelix) {
                    mark.helixId = keepHelix;
                }
            }
            if (helices[mergedHelix] && helices[mergedHelix].length > 0) {
                if (!helices[keepHelix])
                    helices[keepHelix] = [];
                helices[keepHelix].push(...helices[mergedHelix]);
                helices[mergedHelix] = [];
            }
        };
        let pass = 0;
        while (pass++ < 200) {
            const collisionReports = angleCollisions(networkMap);
            if (collisionReports.length === 0) {
                console.log(`[anglecomb] pass=${pass} no collisions remain`);
                break;
            }
            let mergedInThisPass = false;
            outer: for (const report of collisionReports) {
                const sourceHelix = report.helixId;
                for (const conflict of report.conflicts) {
                    const collidedHelices = conflict.colliding_adj_helices
                        .filter((helixId) => helixId !== sourceHelix)
                        .sort((a, b) => a - b);
                    if (collidedHelices.length < 2)
                        continue;
                    for (let i = 0; i < collidedHelices.length; i++) {
                        const helixA = collidedHelices[i];
                        for (let j = i + 1; j < collidedHelices.length; j++) {
                            const helixB = collidedHelices[j];
                            if (!disjoint(helices, helixA, helixB, grid))
                                continue;
                            const keepHelix = Math.min(helixA, helixB);
                            const mergedHelix = Math.max(helixA, helixB);
                            mergeHelixInto(keepHelix, mergedHelix);
                            mergedPairs.push({
                                sourceHelix,
                                keepHelix,
                                mergedHelix,
                                angle: conflict.angle
                            });
                            console.log(`[anglecomb] combined helix ${mergedHelix} into ${keepHelix} in helices[][] and grid ` +
                                `(source=${sourceHelix}, angle=${conflict.angle})`);
                            mergedInThisPass = true;
                            break outer;
                        }
                    }
                }
            }
            if (!mergedInThisPass) {
                console.log(`[anglecomb] pass=${pass} collisions remain but no disjoint collision pairs could be merged`);
                break;
            }
            networkMap = getAngles(grid, helices, lattice);
        }
        networkMap = getAngles(grid, helices, lattice);
        return { networkMap, mergedPairs };
    }
    toscad.anglecomb = anglecomb;
    function anglecorr(grid, helices, lattice = 'honeycomb', angleMap = getAngles(grid, helices, lattice)) {
        let networkMap = angleMap;
        const correctedPairs = [];
        const normalizeAngle = (angle) => ((angle % 360) + 360) % 360;
        let pass = 0;
        while (pass++ < 200) {
            const collisionReports = angleCollisions(networkMap);
            if (collisionReports.length === 0) {
                console.log(`[anglecorr] pass=${pass} no collisions remain`);
                break;
            }
            let correctedInThisPass = false;
            outer: for (const report of collisionReports) {
                const sourceHelix = report.helixId;
                const sourceMap = networkMap.get(sourceHelix);
                if (!sourceMap || sourceMap.size === 0)
                    continue;
                const baseHelices = Array.from(sourceMap.entries())
                    .filter(([, angle]) => normalizeAngle(angle) === 0)
                    .map(([adjHelix]) => adjHelix);
                const baseHelix = baseHelices.length > 0 ? Math.min(...baseHelices) : null;
                for (const conflict of report.conflicts) {
                    const collidedHelices = conflict.colliding_adj_helices
                        .filter((helixId) => helixId !== sourceHelix)
                        .sort((a, b) => a - b);
                    if (collidedHelices.length < 2)
                        continue;
                    for (let i = 0; i < collidedHelices.length; i++) {
                        const helixA = collidedHelices[i];
                        for (let j = i + 1; j < collidedHelices.length; j++) {
                            const helixB = collidedHelices[j];
                            if (disjoint(helices, helixA, helixB, grid))
                                continue;
                            const connectionCount = sourceMap.size;
                            if (connectionCount >= 4) {
                                console.warn(`[anglecorr] source=${sourceHelix} has ${connectionCount} connections and non-disjoint ` +
                                    `collision at angle=${conflict.angle}; no empty spots available`);
                                continue;
                            }
                            const pairDescending = [helixA, helixB].sort((a, b) => b - a);
                            let adjustedHelix = pairDescending.find((hId) => hId !== baseHelix) ?? pairDescending[0];
                            const oldAngleValue = sourceMap.get(adjustedHelix);
                            if (oldAngleValue === undefined)
                                continue;
                            const oldAngle = normalizeAngle(oldAngleValue);
                            const usedAngles = new Set();
                            for (const [adjHelix, angle] of sourceMap.entries()) {
                                if (adjHelix === adjustedHelix)
                                    continue;
                                usedAngles.add(normalizeAngle(angle));
                            }
                            let newAngle = normalizeAngle(oldAngle + 120);
                            if (usedAngles.has(newAngle)) {
                                newAngle = normalizeAngle(newAngle + 120);
                            }
                            if (usedAngles.has(newAngle)) {
                                console.warn(`[anglecorr] source=${sourceHelix} could not place corrected angle for helix=${adjustedHelix}; ` +
                                    `all 120-degree alternatives occupied`);
                                continue;
                            }
                            sourceMap.set(adjustedHelix, newAngle);
                            correctedPairs.push({
                                sourceHelix,
                                adjustedHelix,
                                oldAngle,
                                newAngle,
                                conflictAngle: conflict.angle
                            });
                            console.log(`[anglecorr] corrected source=${sourceHelix}: helix ${adjustedHelix} angle ${oldAngle} -> ${newAngle} ` +
                                `(conflict angle=${conflict.angle})`);
                            correctedInThisPass = true;
                            break outer;
                        }
                    }
                }
            }
            if (!correctedInThisPass) {
                console.log(`[anglecorr] pass=${pass} collisions remain but no non-disjoint collision could be corrected`);
                break;
            }
        }
        return { networkMap, correctedPairs };
    }
    toscad.anglecorr = anglecorr;
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
        applyHelixShifts(grid, cumulativeShift);
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
    function computeHelixPcaAxis(helix) {
        if (!helix || helix.length < 2)
            return null;
        const points = [];
        for (const nt of helix) {
            if (!nt)
                continue;
            const pos = nt.getPos();
            if (!pos)
                continue;
            points.push(pos.clone());
        }
        if (points.length < 2)
            return null;
        const centroid = new THREE.Vector3();
        for (const p of points)
            centroid.add(p);
        centroid.multiplyScalar(1 / points.length);
        const centered = [];
        for (const p of points) {
            const rel = p.clone().sub(centroid);
            if (rel.lengthSq() > 1e-12)
                centered.push(rel);
        }
        if (centered.length < 2)
            return null;
        const cov = {
            xx: 0, xy: 0, xz: 0,
            yy: 0, yz: 0,
            zz: 0,
        };
        for (const r of centered) {
            cov.xx += r.x * r.x;
            cov.xy += r.x * r.y;
            cov.xz += r.x * r.z;
            cov.yy += r.y * r.y;
            cov.yz += r.y * r.z;
            cov.zz += r.z * r.z;
        }
        let axis = centered[0].clone().normalize();
        if (axis.lengthSq() <= 1e-12)
            return null;
        for (let i = 0; i < 16; i++) {
            const next = new THREE.Vector3(cov.xx * axis.x + cov.xy * axis.y + cov.xz * axis.z, cov.xy * axis.x + cov.yy * axis.y + cov.yz * axis.z, cov.xz * axis.x + cov.yz * axis.y + cov.zz * axis.z);
            if (next.lengthSq() <= 1e-12)
                break;
            next.normalize();
            if (next.dot(axis) < 0)
                next.multiplyScalar(-1);
            axis = next;
        }
        return axis.lengthSq() > 1e-12 ? axis.normalize() : null;
    }
    function computeHelixCentroid(helix) {
        if (!helix || helix.length === 0)
            return null;
        const sum = new THREE.Vector3();
        let count = 0;
        for (const nt of helix) {
            if (!nt)
                continue;
            const pos = nt.getPos();
            if (!pos)
                continue;
            sum.add(pos);
            count++;
        }
        return count > 0 ? sum.divideScalar(count) : null;
    }
    /**
     * Returns true when helix A and helix B are colinear — i.e., the
     * centre of helix B lies within `maxOffAxisDist` (simulation length
     * units) of the infinite line defined by helix A's centroid + its
     * PCA axis.  The check is symmetric: we test B against A *and*
     * A against B and require both to pass.
     *
     * Typical inter-helix spacing in DNA origami is ~2–3 oxDNA units,
     * so stacked (colinear) helices should have a rejection distance
     * close to 0; side-by-side helices will have ~2–3 units.
     */
    function areHelicesColinear(helixAId, helixBId, helices, pcaAxisCache, centroidCache, maxOffAxisDist = 10) {
        if (!helices)
            return true;
        const getAxis = (hId) => {
            if (!pcaAxisCache.has(hId))
                pcaAxisCache.set(hId, computeHelixPcaAxis(helices[hId]));
            return pcaAxisCache.get(hId) ?? null;
        };
        const getCentroid = (hId) => {
            if (!centroidCache.has(hId))
                centroidCache.set(hId, computeHelixCentroid(helices[hId]));
            return centroidCache.get(hId) ?? null;
        };
        const axisA = getAxis(helixAId);
        const centroidA = getCentroid(helixAId);
        const centroidB = getCentroid(helixBId);
        if (!axisA || !centroidA || !centroidB)
            return true; // no data — allow
        // Perpendicular (rejection) distance from centroidB to the line
        // through centroidA along axisA:
        //   d = |(centroidB - centroidA) × axisA|
        const diff = centroidB.clone().sub(centroidA);
        const rejection = diff.clone().cross(axisA).length();
        return rejection <= maxOffAxisDist;
    }
    function areHelixPcaAxesCompatible(helixAId, helixBId, helices, pcaAxisCache, maxAngleDeg = 45) {
        if (!helices)
            return true;
        const getAxis = (helixId) => {
            if (!pcaAxisCache.has(helixId)) {
                pcaAxisCache.set(helixId, computeHelixPcaAxis(helices[helixId]));
            }
            return pcaAxisCache.get(helixId) ?? null;
        };
        const axisA = getAxis(helixAId);
        const axisB = getAxis(helixBId);
        if (!axisA || !axisB)
            return false;
        const cosThreshold = Math.cos(maxAngleDeg * Math.PI / 180);
        const cosine = Math.min(1, Math.max(-1, Math.abs(axisA.dot(axisB))));
        return cosine >= cosThreshold;
    }
    // 2 helices are candidates for combining when one of their mutual connection has >3 helices connected to it.
    // They get combined only if they have no offset overalp and their PCA axes are compatible (roughly parallel, allowing for some angle).
    function combinedHelices(maxOffsetDist, grid, helices, binderHelices) {
        const retiredHelices = new Set();
        const touchedHelices = new Set();
        const combinedPairs = [];
        const binderSet = new Set(binderHelices ?? []);
        const pcaAxisCache = new Map();
        const centroidCache = new Map();
        const maxOffs = maxOffsetDist;
        const buildAdjacency = (crossovers) => {
            const adjacency = new Map();
            const ensure = (hId) => {
                if (!adjacency.has(hId))
                    adjacency.set(hId, new Set());
                return adjacency.get(hId);
            };
            for (const [from, toMap] of crossovers.entries()) {
                const fromSet = ensure(from);
                for (const [to, stats] of toMap.entries()) {
                    const totalConnections = stats.sameWalk + stats.diffWalk;
                    if (totalConnections <= 0)
                        continue;
                    fromSet.add(to);
                    ensure(to).add(from);
                }
            }
            return adjacency;
        };
        const buildOffsetSets = () => {
            const offsetSets = new Map();
            for (const [, mark] of grid.entries()) {
                if (!offsetSets.has(mark.helixId))
                    offsetSets.set(mark.helixId, new Set());
                offsetSets.get(mark.helixId).add(mark.offset);
            }
            return offsetSets;
        };
        // stupid??
        const offsetsDisjoint = (a, b) => {
            if (!a || !b || a.size === 0 || b.size === 0)
                return false;
            const smaller = a.size <= b.size ? a : b;
            const larger = a.size <= b.size ? b : a;
            for (const offs of smaller) {
                if (larger.has(offs))
                    return false;
            }
            return true;
        };
        const distanceSquared = (a, b, positions) => {
            const pa = positions.get(a);
            const pb = positions.get(b);
            if (!pa || !pb)
                return Number.POSITIVE_INFINITY;
            const dx = pb.x - pa.x;
            const dy = pb.y - pa.y;
            return dx * dx + dy * dy;
        };
        const mergeHelixInto = (keep, merged) => {
            for (const [, mark] of grid.entries()) {
                if (mark.helixId === merged) {
                    mark.helixId = keep;
                }
            }
            if (helices && helices[merged] && helices[merged].length > 0) {
                if (!helices[keep])
                    helices[keep] = [];
                helices[keep].push(...helices[merged]);
                helices[merged] = [];
            }
            // Recalculate the PCA axis and centroid for the merged helix after combining
            retiredHelices.add(merged);
            touchedHelices.add(keep);
            touchedHelices.add(merged);
            pcaAxisCache.delete(keep);
            pcaAxisCache.delete(merged);
            centroidCache.delete(keep);
            centroidCache.delete(merged);
            combinedPairs.push({ keep, merged });
        };
        const positions = helices ? getRelativePositions(helices) : new Map();
        let iteration = 0;
        while (iteration++ < 100) {
            const { crossovers, helixIds } = toscad.collectCrossovers(grid);
            const adjacency = buildAdjacency(crossovers);
            const offsetSets = buildOffsetSets();
            const areOffsetsDisjoint = (a, b) => helices ? disjoint(helices, a, b, grid) : offsetsDisjoint(offsetSets.get(a), offsetSets.get(b));
            for (const hId of helixIds) {
                if (!adjacency.has(hId))
                    adjacency.set(hId, new Set());
            }
            const hubs = Array.from(adjacency.entries())
                .filter(([helixId, neighbors]) => neighbors.size > 3 && !retiredHelices.has(helixId) && !binderSet.has(helixId))
                .map(([helixId]) => helixId)
                .sort((a, b) => a - b);
            if (hubs.length === 0)
                break;
            let mergedInThisIteration = false;
            for (const hub of hubs) {
                if (retiredHelices.has(hub))
                    continue;
                if (binderSet.has(hub))
                    continue;
                const neighbors = Array.from(adjacency.get(hub) ?? [])
                    .filter(n => !retiredHelices.has(n) && n !== hub && !binderSet.has(n));
                if (neighbors.length === 0)
                    continue;
                // Step 1: try combining one neighbor directly into the overloaded hub.
                const compatibleWithHub = neighbors.filter(n => areOffsetsDisjoint(hub, n)
                    && areHelixPcaAxesCompatible(hub, n, helices, pcaAxisCache, 45)
                    && areHelicesColinear(hub, n, helices, pcaAxisCache, centroidCache, maxOffs));
                if (compatibleWithHub.length > 0) {
                    let bestNeighbor = compatibleWithHub[0];
                    let bestDist = distanceSquared(hub, bestNeighbor, positions);
                    for (const candidate of compatibleWithHub) {
                        const distSq = distanceSquared(hub, candidate, positions);
                        if (distSq < bestDist) {
                            bestDist = distSq;
                            bestNeighbor = candidate;
                        }
                    }
                    // Refresh hub's offset snapshot BEFORE the merge so Step 2 sees the full combined range.
                    const neighborOffs = offsetSets.get(bestNeighbor);
                    if (neighborOffs) {
                        if (!offsetSets.has(hub))
                            offsetSets.set(hub, new Set());
                        for (const off of neighborOffs)
                            offsetSets.get(hub).add(off);
                    }
                    mergeHelixInto(hub, bestNeighbor);
                    mergedInThisIteration = true;
                }
                // Step 2: try combining overloaded hub neighbors with each other.
                const remainingNeighbors = neighbors
                    .filter(n => !retiredHelices.has(n));
                while (remainingNeighbors.length >= 2) {
                    let bestPair = null;
                    let bestPairDist = Number.POSITIVE_INFINITY;
                    for (let i = 0; i < remainingNeighbors.length; i++) {
                        for (let j = i + 1; j < remainingNeighbors.length; j++) {
                            const a = remainingNeighbors[i];
                            const b = remainingNeighbors[j];
                            if (!areOffsetsDisjoint(a, b))
                                continue;
                            if (!areHelixPcaAxesCompatible(a, b, helices, pcaAxisCache, 45))
                                continue;
                            if (!areHelicesColinear(a, b, helices, pcaAxisCache, centroidCache, maxOffs))
                                continue;
                            const distSq = distanceSquared(a, b, positions);
                            if (distSq < bestPairDist) {
                                bestPairDist = distSq;
                                bestPair = [a, b];
                            }
                        }
                    }
                    if (!bestPair)
                        break;
                    const keep = Math.min(bestPair[0], bestPair[1]);
                    const merged = Math.max(bestPair[0], bestPair[1]);
                    // Refresh keep's offset snapshot so the next pair check in this loop sees the combined range.
                    const mergedOffs = offsetSets.get(merged);
                    if (mergedOffs) {
                        if (!offsetSets.has(keep))
                            offsetSets.set(keep, new Set());
                        for (const off of mergedOffs)
                            offsetSets.get(keep).add(off);
                    }
                    mergeHelixInto(keep, merged);
                    mergedInThisIteration = true;
                    const mergedIdx = remainingNeighbors.indexOf(merged);
                    if (mergedIdx >= 0)
                        remainingNeighbors.splice(mergedIdx, 1);
                }
            }
            if (!mergedInThisIteration)
                break;
        }
        const activeAfterMerge = new Set();
        for (const [, mark] of grid.entries()) {
            activeAfterMerge.add(mark.helixId);
        }
        const maxHelixIdAfterMerge = activeAfterMerge.size > 0
            ? Math.max(...Array.from(activeAfterMerge))
            : -1;
        const removedHelices = [];
        for (let hId = 0; hId <= maxHelixIdAfterMerge; hId++) {
            if (!activeAfterMerge.has(hId))
                removedHelices.push(hId);
        }
        const sortedActive = Array.from(activeAfterMerge).sort((a, b) => a - b);
        const helixIdRemap = new Map();
        sortedActive.forEach((oldId, newId) => helixIdRemap.set(oldId, newId));
        for (const [, mark] of grid.entries()) {
            const remapped = helixIdRemap.get(mark.helixId);
            if (remapped !== undefined) {
                mark.helixId = remapped;
            }
        }
        if (helices) {
            const remappedHelices = sortedActive.map((oldId) => helices[oldId] ?? []);
            helices.length = 0;
            helices.push(...remappedHelices);
        }
        if (combinedPairs.length > 0) {
            console.log(`[combinedHelices] Combined ${combinedPairs.length} helix pairs: ` +
                combinedPairs.map(p => `${p.merged}->${p.keep}`).join(', '));
            if (removedHelices.length > 0) {
                console.log(`[combinedHelices] Removed empty helices: [${removedHelices.join(', ')}]`);
            }
        }
        else {
            console.log('[combinedHelices] No helix combinations applied.');
        }
        return {
            combinedPairs,
            combinedHelices: Array.from(touchedHelices).sort((a, b) => a - b),
            removedHelices,
            helixIdRemap: Object.fromEntries(helixIdRemap.entries())
        };
    }
    toscad.combinedHelices = combinedHelices;
    // a much more hand-wavy function, combines helices if they don't overlap and they have the same connections.
    function trialComb(grid, helices, binderHelices) {
        const combinedPairs = [];
        const retiredHelices = new Set();
        const touchedHelices = new Set();
        const binderSet = new Set(binderHelices ?? []);
        const pcaAxisCache = new Map();
        const buildAdjacency = (crossovers) => {
            const adjacency = new Map();
            const ensure = (hId) => {
                if (!adjacency.has(hId))
                    adjacency.set(hId, new Set());
                return adjacency.get(hId);
            };
            for (const [from, toMap] of crossovers.entries()) {
                const fromSet = ensure(from);
                for (const [to, stats] of toMap.entries()) {
                    const totalConnections = stats.sameWalk + stats.diffWalk;
                    if (totalConnections <= 0)
                        continue;
                    fromSet.add(to);
                    ensure(to).add(from);
                }
            }
            return adjacency;
        };
        const buildOffsetSets = () => {
            const offsetSets = new Map();
            for (const [, mark] of grid.entries()) {
                if (!offsetSets.has(mark.helixId))
                    offsetSets.set(mark.helixId, new Set());
                offsetSets.get(mark.helixId).add(mark.offset);
            }
            return offsetSets;
        };
        const setEqual = (a, b) => {
            if (a.size !== b.size)
                return false;
            for (const v of a) {
                if (!b.has(v))
                    return false;
            }
            return true;
        };
        const offsetsDisjoint = (a, b) => {
            if (!a || !b || a.size === 0 || b.size === 0)
                return false;
            const smaller = a.size <= b.size ? a : b;
            const larger = a.size <= b.size ? b : a;
            for (const off of smaller) {
                if (larger.has(off))
                    return false;
            }
            return true;
        };
        let iteration = 0;
        while (iteration++ < 100) {
            const { crossovers, helixIds } = toscad.collectCrossovers(grid);
            const adjacency = buildAdjacency(crossovers);
            const offsetSets = buildOffsetSets();
            const areOffsetsDisjoint = (a, b) => helices ? disjoint(helices, a, b, grid) : offsetsDisjoint(offsetSets.get(a), offsetSets.get(b));
            for (const hId of helixIds) {
                if (!adjacency.has(hId))
                    adjacency.set(hId, new Set());
            }
            const candidates = Array.from(adjacency.keys()).sort((a, b) => a - b);
            let mergedInIteration = false;
            for (let i = 0; i < candidates.length; i++) {
                const a = candidates[i];
                if (retiredHelices.has(a))
                    continue;
                if (binderSet.has(a))
                    continue;
                for (let j = i + 1; j < candidates.length; j++) {
                    const b = candidates[j];
                    if (retiredHelices.has(b))
                        continue;
                    if (binderSet.has(b))
                        continue;
                    const aNeighbors = adjacency.get(a) ?? new Set();
                    const bNeighbors = adjacency.get(b) ?? new Set();
                    if (!setEqual(aNeighbors, bNeighbors))
                        continue;
                    if (!areOffsetsDisjoint(a, b))
                        continue;
                    if (!areHelixPcaAxesCompatible(a, b, helices, pcaAxisCache, 45))
                        continue;
                    const keep = Math.min(a, b);
                    const merged = Math.max(a, b);
                    // Refresh keep's offset snapshot before merging so later candidates in this pass
                    // see the full combined range and don't incorrectly pass the disjoint check.
                    const trialMergedOffs = offsetSets.get(merged);
                    if (trialMergedOffs) {
                        if (!offsetSets.has(keep))
                            offsetSets.set(keep, new Set());
                        for (const off of trialMergedOffs)
                            offsetSets.get(keep).add(off);
                    }
                    for (const [, mark] of grid.entries()) {
                        if (mark.helixId === merged) {
                            mark.helixId = keep;
                        }
                    }
                    if (helices && helices[merged] && helices[merged].length > 0) {
                        if (!helices[keep])
                            helices[keep] = [];
                        helices[keep].push(...helices[merged]);
                        helices[merged] = [];
                    }
                    retiredHelices.add(merged);
                    touchedHelices.add(keep);
                    touchedHelices.add(merged);
                    pcaAxisCache.delete(keep);
                    pcaAxisCache.delete(merged);
                    combinedPairs.push({ keep, merged });
                    mergedInIteration = true;
                    break;
                }
            }
            if (!mergedInIteration)
                break;
        }
        const activeAfterMerge = new Set();
        for (const [, mark] of grid.entries()) {
            activeAfterMerge.add(mark.helixId);
        }
        const maxHelixIdAfterMerge = activeAfterMerge.size > 0
            ? Math.max(...Array.from(activeAfterMerge))
            : -1;
        const removedHelices = [];
        for (let hId = 0; hId <= maxHelixIdAfterMerge; hId++) {
            if (!activeAfterMerge.has(hId))
                removedHelices.push(hId);
        }
        const sortedActive = Array.from(activeAfterMerge).sort((a, b) => a - b);
        const helixIdRemap = new Map();
        sortedActive.forEach((oldId, newId) => helixIdRemap.set(oldId, newId));
        for (const [, mark] of grid.entries()) {
            const remapped = helixIdRemap.get(mark.helixId);
            if (remapped !== undefined) {
                mark.helixId = remapped;
            }
        }
        if (helices) {
            const remappedHelices = sortedActive.map((oldId) => helices[oldId] ?? []);
            helices.length = 0;
            helices.push(...remappedHelices);
        }
        if (combinedPairs.length > 0) {
            console.log(`[trialComb] Combined ${combinedPairs.length} helix pairs: ` +
                combinedPairs.map(p => `${p.merged}->${p.keep}`).join(', '));
            if (removedHelices.length > 0) {
                console.log(`[trialComb] Removed empty helices: [${removedHelices.join(', ')}]`);
            }
        }
        else {
            console.log('[trialComb] No helix combinations applied.');
        }
        return {
            combinedPairs,
            combinedHelices: Array.from(touchedHelices).sort((a, b) => a - b),
            removedHelices,
            helixIdRemap: Object.fromEntries(helixIdRemap.entries())
        };
    }
    toscad.trialComb = trialComb;
    function subtractPositions(map, keyA, keyB) {
        const posA = map.get(keyA);
        const posB = map.get(keyB);
        if (!posA || !posB)
            return null;
        return {
            x: posA.x - posB.x,
            y: posA.y - posB.y
        };
    }
    /**
     * Calculates absolute grid coordinates from local helix-to-helix angles.
     *
     * Phase 1: Build a strict lattice spanning tree from helix 0 using
     * weighted BFS (top-3 children only at each node).
     *
     * Phase 2: Place deferred/artefact helices in nearest open coordinates
     * around their parent once the phase-1 core is locked.
     */
    function calculateGlobalPositions(networkMap, crossoverWeights, options) {
        const ANGLES = [0, 120, 240];
        const STEP_BY_PARITY = {
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
        const normalizeAngle = (angle) => ((angle % 360) + 360) % 360;
        const angleDistance = (a, b) => {
            const diff = Math.abs(normalizeAngle(a) - normalizeAngle(b));
            return Math.min(diff, 360 - diff);
        };
        const parityAt = (coord) => (((coord.col + coord.row) & 1) === 0 ? 'even' : 'odd');
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
            const parity = parityAt(coord);
            for (const angle of ANGLES) {
                const delta = STEP_BY_PARITY[parity][angle];
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
            const parity = parityAt(coord);
            return ANGLES.map((angle) => {
                const step = STEP_BY_PARITY[parity][angle];
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
            const selected = candidates.slice(0, 3);
            const overflow = candidates.slice(3);
            for (const item of overflow) {
                enqueueDeferred(node.helixId, item.neighborId, node.offset);
            }
            for (const item of selected) {
                if (positions.has(item.neighborId))
                    continue;
                const predictedGlobal = normalizeAngle(item.localAngle + node.offset);
                const snappedAngle = snapToLatticeAngle(predictedGlobal);
                const step = STEP_BY_PARITY[parityAt(parentCoord)][snappedAngle];
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
    function HelixPosByRelativeBfs(grid, helices) {
        const placed = new Map();
        const occupied = new Set();
        const placedByCoord = new Map();
        const orphaned = [];
        const orphanedIds = new Set();
        const relativePositions = getRelativePositions(helices);
        const { crossovers, helixIds } = toscad.collectCrossovers(grid);
        const helixCount = Math.max(helices.length, ...Array.from(grid.values()).map((mark) => mark.helixId + 1), ...Array.from(helixIds).map((helixId) => helixId + 1), 0);
        const positionKey = (coord) => `${coord.x},${coord.y}`;
        const isOpen = (coord) => !occupied.has(positionKey(coord));
        const getOccupant = (coord) => placedByCoord.get(positionKey(coord));
        const tryPlace = (helixId, coord) => {
            if (placed.has(helixId) || !isOpen(coord))
                return false;
            placed.set(helixId, coord);
            occupied.add(positionKey(coord));
            placedByCoord.set(positionKey(coord), helixId);
            return true;
        };
        const enqueueOrphan = (helixId, parentId) => {
            if (placed.has(helixId) || orphanedIds.has(helixId))
                return;
            orphaned.push({ helixId, parentId });
            orphanedIds.add(helixId);
        };
        const buildAdjacency = () => {
            const adjacency = new Map();
            const ensure = (helixId) => {
                if (!adjacency.has(helixId))
                    adjacency.set(helixId, new Set());
                return adjacency.get(helixId);
            };
            for (let helixId = 0; helixId < helixCount; helixId++)
                ensure(helixId);
            // Primary: backbone crossover connections
            for (const [from, toMap] of crossovers.entries()) {
                const row = ensure(from);
                for (const to of toMap.keys()) {
                    row.add(to);
                    ensure(to).add(from);
                }
            }
            // Secondary: pair-bridge connections.
            // Binder helices have NO backbone crossovers (n3/n5 never leave the
            // same helix), so collectCrossovers misses them entirely.  Instead
            // they attach to double-stranded helices via base-pair links.
            // Scan every nucleotide in the grid: if nt A is on helixX and its
            // pair is on helixY (different), add X<->Y to adjacency.
            for (const [ntId, mark] of grid.entries()) {
                const nt = elements.get(ntId);
                if (!nt || !(nt instanceof Nucleotide))
                    continue;
                const pairNt = nt.pair;
                if (!pairNt || !(pairNt instanceof Nucleotide))
                    continue;
                const pairMark = grid.get(pairNt.id);
                if (!pairMark)
                    continue;
                if (pairMark.helixId === mark.helixId)
                    continue;
                ensure(mark.helixId).add(pairMark.helixId);
                ensure(pairMark.helixId).add(mark.helixId);
            }
            return adjacency;
        };
        const adjacency = buildAdjacency();
        // Spiral outward from `origin` in true Manhattan order (distance 1, 2, 3, …)
        // so orphans are always placed as close to their parent as possible.
        const findNearestOpenAround = (origin) => {
            if (isOpen(origin))
                return origin;
            const maxDist = helixCount + 32;
            for (let dist = 1; dist <= maxDist; dist++) {
                // Enumerate all integer coords at exact Manhattan distance = dist
                // and pick the first open one, preferring low |y| then low |x|
                // (i.e. stay close to the same row as origin first).
                let best = null;
                let bestKey = Infinity;
                for (let dx = -dist; dx <= dist; dx++) {
                    const dyAbs = dist - Math.abs(dx);
                    for (const dy of dyAbs === 0 ? [0] : [-dyAbs, dyAbs]) {
                        const cand = { x: origin.x + dx, y: origin.y + dy };
                        if (!isOpen(cand))
                            continue;
                        // Sort key: |dy| first (favour same row), then |dx|
                        const key = Math.abs(dy) * (maxDist * 2 + 1) + Math.abs(dx);
                        if (key < bestKey) {
                            bestKey = key;
                            best = cand;
                        }
                    }
                }
                if (best)
                    return best;
            }
            return { x: origin.x + maxDist + 1, y: origin.y };
        };
        const chooseRoot = () => {
            for (let helixId = 0; helixId < helixCount; helixId++) {
                if ((adjacency.get(helixId)?.size ?? 0) === 3)
                    return helixId;
            }
            for (let helixId = 0; helixId < helixCount; helixId++) {
                if ((adjacency.get(helixId)?.size ?? 0) > 0)
                    return helixId;
            }
            return 0;
        };
        const getSeedCoordinate = (helixId) => {
            const degree = adjacency.get(helixId)?.size ?? 0;
            if (degree === 3)
                return { x: 0, y: 0 };
            // Do NOT use raw relativePositions as grid coords — those are in
            // simulation units (2-3 oxDNA units per helix spacing) and would
            // place the seed 10-30+ cells away from the origin.
            // Just find the nearest open cell to the origin.
            return findNearestOpenAround({ x: 0, y: 0 });
        };
        // Returns the nearest already-placed neighbor's grid coord, or falls
        // back to the given default.  Used to anchor orphan/straggler placement.
        const nearestPlacedNeighbor = (helixId, fallback) => {
            for (const nbId of adjacency.get(helixId) ?? []) {
                const nb = placed.get(nbId);
                if (nb)
                    return nb;
            }
            return fallback;
        };
        const placeNeighborGroup = (parentId, queue) => {
            const parentCoord = placed.get(parentId);
            if (!parentCoord)
                return;
            const neighbors = Array.from(adjacency.get(parentId) ?? [])
                .map((helixId) => ({
                helixId,
                delta: subtractPositions(relativePositions, helixId, parentId)
            }))
                .filter((entry) => !!entry.delta)
                .sort((a, b) => {
                const dyDelta = Math.abs(b.delta.y) - Math.abs(a.delta.y);
                if (dyDelta !== 0)
                    return dyDelta;
                const dxDelta = Math.abs(b.delta.x) - Math.abs(a.delta.x);
                if (dxDelta !== 0)
                    return dxDelta;
                return a.helixId - b.helixId;
            });
            if (neighbors.length === 0)
                return;
            const planned = new Set();
            const assignedTargets = new Map();
            // Only unplaced neighbors are candidates for new slots.
            const unplacedNeighbors = neighbors.filter((entry) => !placed.has(entry.helixId));
            // Assign Y only when at least one candidate is truly vertical,
            // i.e. |dy| > |dx|. Otherwise skip Y and keep placements on ±X.
            const verticalCandidates = unplacedNeighbors
                .filter((entry) => Math.abs(entry.delta.y) > Math.abs(entry.delta.x))
                .sort((a, b) => {
                const dyDelta = Math.abs(b.delta.y) - Math.abs(a.delta.y);
                if (dyDelta !== 0)
                    return dyDelta;
                const ratioA = Math.abs(a.delta.y) / (Math.abs(a.delta.x) + 1e-9);
                const ratioB = Math.abs(b.delta.y) / (Math.abs(b.delta.x) + 1e-9);
                if (ratioB !== ratioA)
                    return ratioB - ratioA;
                return a.helixId - b.helixId;
            });
            if (verticalCandidates.length > 0) {
                const yCandidate = verticalCandidates[0];
                // For helices with ≤3 connections the y direction is forced by grid parity:
                //   (x + y) odd  → +1 (e.g. even-x/odd-y, odd-x/even-y)
                //   (x + y) even → -1 (e.g. even-x/even-y, odd-x/odd-y)
                // Using & 1 instead of % 2 so negatives resolve correctly.
                // For ≥4 connections use the deltaY sign as before.
                const totalConnections = adjacency.get(parentId)?.size ?? 0;
                const yStep = totalConnections <= 3
                    ? (((parentCoord.x + parentCoord.y) & 1) === 1 ? 1 : -1)
                    : (yCandidate.delta.y < 0 ? 1 : -1);
                assignedTargets.set(yCandidate.helixId, { x: parentCoord.x, y: parentCoord.y + yStep });
                planned.add(yCandidate.helixId);
            }
            const remaining = unplacedNeighbors.filter((entry) => !planned.has(entry.helixId));
            const positive = remaining
                .filter((entry) => entry.delta.x >= 0)
                .sort((a, b) => b.delta.x - a.delta.x || a.helixId - b.helixId);
            const negative = remaining
                .filter((entry) => entry.delta.x < 0)
                .sort((a, b) => a.delta.x - b.delta.x || a.helixId - b.helixId);
            const assignXSlot = (entry, xStep) => {
                const target = { x: parentCoord.x + xStep, y: parentCoord.y };
                assignedTargets.set(entry.helixId, target);
                planned.add(entry.helixId);
            };
            if (positive.length > 0 && negative.length > 0) {
                assignXSlot(positive[0], 1);
                assignXSlot(negative[0], -1);
            }
            else if (positive.length > 0) {
                assignXSlot(positive[0], 1);
            }
            else if (negative.length > 0) {
                assignXSlot(negative[0], -1);
            }
            for (const entry of remaining) {
                if (!planned.has(entry.helixId))
                    enqueueOrphan(entry.helixId, parentId);
            }
            for (const [helixId, target] of assignedTargets.entries()) {
                const placedCoord = placed.get(helixId);
                if (placedCoord) {
                    continue;
                }
                const occupant = getOccupant(target);
                if (occupant !== undefined && occupant !== helixId) {
                    enqueueOrphan(helixId, parentId);
                    continue;
                }
                if (tryPlace(helixId, target)) {
                    queue.push(helixId);
                }
                else {
                    enqueueOrphan(helixId, parentId);
                }
            }
        };
        const root = chooseRoot();
        occupied.add('0,0');
        const rootCoord = getSeedCoordinate(root);
        if (rootCoord.x === 0 && rootCoord.y === 0) {
            occupied.delete('0,0');
        }
        tryPlace(root, rootCoord);
        const queue = [root];
        let qIdx = 0;
        while (qIdx < queue.length) {
            placeNeighborGroup(queue[qIdx++], queue);
        }
        for (const orphan of orphaned) {
            if (placed.has(orphan.helixId))
                continue;
            // Anchor to the nearest already-placed neighbor (any connection),
            // not just the first parent that orphaned this helix.  This prevents
            // a helix from landing far away because it was orphaned early by a
            // distant parent.
            const fallback = placed.get(orphan.parentId) ?? placed.get(root) ?? { x: 0, y: 0 };
            const anchor = nearestPlacedNeighbor(orphan.helixId, fallback);
            const openCoord = findNearestOpenAround(anchor);
            if (tryPlace(orphan.helixId, openCoord)) {
                queue.push(orphan.helixId);
            }
        }
        while (qIdx < queue.length) {
            placeNeighborGroup(queue[qIdx++], queue);
        }
        for (let helixId = 0; helixId < helixCount; helixId++) {
            if (placed.has(helixId))
                continue;
            // Anchor stragglers to a placed neighbor, not raw relativePositions
            // (which are in simulation units and would scatter them far away).
            const anchor = nearestPlacedNeighbor(helixId, { x: 0, y: 0 });
            tryPlace(helixId, findNearestOpenAround(anchor));
        }
        const result = new Map();
        for (let helixId = 0; helixId < helixCount; helixId++) {
            const coord = placed.get(helixId) ?? { x: helixId, y: 0 };
            result.set(helixId, [coord.x, coord.y]);
        }
        return result;
    }
    toscad.HelixPosByRelativeBfs = HelixPosByRelativeBfs;
    function HelixPos(grid, helices) {
        const HEX_AXIAL_DIRS = [
            { q: 1, r: 0 },
            { q: 1, r: -1 },
            { q: 0, r: -1 },
            { q: -1, r: 0 },
            { q: -1, r: 1 },
            { q: 0, r: 1 }
        ];
        const axialAdd = (a, b) => ({ q: a.q + b.q, r: a.r + b.r });
        const axialDistance = (a, b) => {
            const dq = a.q - b.q;
            const dr = a.r - b.r;
            const ds = (-a.q - a.r) - (-b.q - b.r);
            return Math.max(Math.abs(dq), Math.abs(dr), Math.abs(ds));
        };
        const axialKey = (a) => `${a.q},${a.r}`;
        const edgeKey = (a, b) => `${a}|${b}`;
        const axialRToOddQRow = (q, r) => r + ((q - (q & 1)) / 2);
        // yes
        const helixCOM = (helix) => {
            if (!helix || helix.length === 0)
                return null;
            const sum = new THREE.Vector3();
            let count = 0;
            for (const nt of helix) {
                if (!(nt instanceof Nucleotide))
                    continue;
                sum.add(nt.getPos());
                count++;
            }
            if (count === 0)
                return null;
            return sum.divideScalar(count);
        };
        // Creates a coordinate system based on 
        const estimateHexBasis3D = (allHelices) => {
            const centers = [];
            const axisSamples = [];
            // "PCA" is basically using endpoints 1 and 2 as a rough axis direction, and then averaging the centroids to find the origin.
            for (const helix of allHelices) {
                if (!helix || helix.length === 0)
                    continue;
                const center = helixCOM(helix);
                if (center)
                    centers.push(center);
                const ep = toscad.helixEndpoints(helix);
                if (!ep)
                    continue;
                const axis = ep.end2.getPos().clone().sub(ep.end1.getPos());
                if (axis.lengthSq() > 1e-8)
                    axisSamples.push(axis.normalize());
            }
            if (centers.length === 0)
                return null;
            // find the average center to use as the origin
            const origin = new THREE.Vector3();
            for (const c of centers)
                origin.add(c);
            origin.divideScalar(centers.length);
            // Average helix direction to get a rough axis. Not really for any rigorous reason, it works fine for this and thus, we will use it.
            let axisVec = new THREE.Vector3(0, 0, 1);
            if (axisSamples.length > 0) {
                axisVec.set(0, 0, 0);
                for (const s of axisSamples)
                    axisVec.add(s);
                if (axisVec.lengthSq() < 1e-8)
                    axisVec.set(0, 0, 1);
                else
                    axisVec.normalize();
            }
            // smart and concise way of finding the perpendicular direction of the spread...
            // removes the entire axisVec vector from the axes, so only the perpendicular spread remains.
            const projectedCenters = centers.map((c) => {
                const rel = c.clone().sub(origin);
                // normalized already, so no need to divide by axisVec.lengthSq()
                return rel.sub(axisVec.clone().multiplyScalar(rel.dot(axisVec)));
            });
            let e1 = new THREE.Vector3(1, 0, 0);
            let bestD2 = Infinity;
            for (let i = 0; i < projectedCenters.length; i++) {
                for (let j = i + 1; j < projectedCenters.length; j++) {
                    const d = projectedCenters[j].clone().sub(projectedCenters[i]);
                    const d2 = d.lengthSq();
                    if (d2 > 1e-8 && d2 < bestD2) {
                        bestD2 = d2;
                        e1 = d.normalize();
                    }
                }
            }
            let e2 = axisVec.clone().cross(e1);
            if (e2.lengthSq() < 1e-8) {
                const fallback = Math.abs(axisVec.x) < 0.9
                    ? new THREE.Vector3(1, 0, 0)
                    : new THREE.Vector3(0, 1, 0);
                e2 = axisVec.clone().cross(fallback);
            }
            e2.normalize();
            e1 = e2.clone().cross(axisVec).normalize();
            const nearestDistances = [];
            for (let i = 0; i < projectedCenters.length; i++) {
                let nearest = Infinity;
                for (let j = 0; j < projectedCenters.length; j++) {
                    if (i === j)
                        continue;
                    const d = projectedCenters[j].clone().sub(projectedCenters[i]).length();
                    if (d > 1e-6 && d < nearest)
                        nearest = d;
                }
                if (nearest < Infinity)
                    nearestDistances.push(nearest);
            }
            let spacing = 1;
            if (nearestDistances.length > 0) {
                nearestDistances.sort((a, b) => a - b);
                spacing = nearestDistances[Math.floor(nearestDistances.length / 2)] || 1;
                if (spacing <= 1e-6)
                    spacing = 1;
            }
            const qVec = e1.clone().multiplyScalar(spacing);
            const rVec = e1.clone().multiplyScalar(0.5 * spacing)
                .add(e2.clone().multiplyScalar((Math.sqrt(3) / 2) * spacing));
            return { origin, qVec, rVec, axisVec };
        };
        const vectorToAxialContinuous = (v, basis) => {
            const inPlane = v.clone().sub(basis.axisVec.clone().multiplyScalar(v.dot(basis.axisVec)));
            const aa = basis.qVec.dot(basis.qVec);
            const ab = basis.qVec.dot(basis.rVec);
            const bb = basis.rVec.dot(basis.rVec);
            const ap = basis.qVec.dot(inPlane);
            const bp = basis.rVec.dot(inPlane);
            const det = aa * bb - ab * ab;
            if (Math.abs(det) < 1e-10)
                return null;
            return {
                q: (ap * bb - bp * ab) / det,
                r: (bp * aa - ap * ab) / det
            };
        };
        const quantizeToHexDirection = (v, basis) => {
            let best = HEX_AXIAL_DIRS[0];
            let bestScore = -Infinity;
            const dirVec = v.clone().sub(basis.axisVec.clone().multiplyScalar(v.dot(basis.axisVec)));
            if (dirVec.lengthSq() < 1e-12)
                return best;
            dirVec.normalize();
            for (const d of HEX_AXIAL_DIRS) {
                const world = basis.qVec.clone().multiplyScalar(d.q).add(basis.rVec.clone().multiplyScalar(d.r)).normalize();
                const score = dirVec.dot(world);
                if (score > bestScore) {
                    bestScore = score;
                    best = d;
                }
            }
            return { q: best.q, r: best.r };
        };
        const collectHelixAdjacency = (currentGrid) => {
            const adjacency = new Map();
            const ensure = (a, b) => {
                if (!adjacency.has(a))
                    adjacency.set(a, new Map());
                const row = adjacency.get(a);
                row.set(b, (row.get(b) ?? 0) + 1);
            };
            const allNtIds = new Set();
            for (const [ntId] of currentGrid.entries())
                allNtIds.add(ntId);
            const visited = new Set();
            for (const [ntId] of currentGrid.entries()) {
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
                const walkForward = new Set();
                let prevMark = null;
                while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                    if (walkForward.has(curr.id))
                        break;
                    walkForward.add(curr.id);
                    visited.add(curr.id);
                    const mark = currentGrid.get(curr.id);
                    if (mark) {
                        if (prevMark && prevMark.helixId !== mark.helixId) {
                            ensure(prevMark.helixId, mark.helixId);
                            ensure(mark.helixId, prevMark.helixId);
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
            return adjacency;
        };
        const nearestOpenAround = (target, occupied, maxRadius = 16) => {
            if (!occupied.has(axialKey(target)))
                return target;
            for (let radius = 1; radius <= maxRadius; radius++) {
                let best = null;
                let bestDist = Infinity;
                for (let dq = -radius; dq <= radius; dq++) {
                    const rMin = Math.max(-radius, -dq - radius);
                    const rMax = Math.min(radius, -dq + radius);
                    for (let dr = rMin; dr <= rMax; dr++) {
                        const cand = { q: target.q + dq, r: target.r + dr };
                        const k = axialKey(cand);
                        if (occupied.has(k))
                            continue;
                        const d = axialDistance(cand, target);
                        if (d < bestDist) {
                            bestDist = d;
                            best = cand;
                        }
                    }
                }
                if (best)
                    return best;
            }
            return { q: target.q + maxRadius + 1, r: target.r };
        };
        const computeOddQGridPositions = (currentGrid, allHelices) => {
            const positions = new Map();
            const helixCount = Math.max(0, ...Array.from(currentGrid.values()).map((m) => m.helixId + 1));
            if (helixCount === 0)
                return positions;
            const basis = estimateHexBasis3D(allHelices);
            const adjacency = collectHelixAdjacency(currentGrid);
            const centers = new Map();
            const projectedAxial = new Map();
            for (let h = 0; h < helixCount; h++) {
                const helix = allHelices[h] ?? [];
                const c = helixCOM(helix);
                if (!c)
                    continue;
                centers.set(h, c);
                if (basis) {
                    const rel = c.clone().sub(basis.origin);
                    const uv = vectorToAxialContinuous(rel, basis);
                    if (uv)
                        projectedAxial.set(h, uv);
                }
            }
            const preferredDelta = new Map();
            if (basis) {
                for (let a = 0; a < helixCount; a++) {
                    const row = adjacency.get(a);
                    if (!row)
                        continue;
                    const ca = centers.get(a);
                    if (!ca)
                        continue;
                    for (const [b] of row.entries()) {
                        const cb = centers.get(b);
                        if (!cb)
                            continue;
                        const d = quantizeToHexDirection(cb.clone().sub(ca), basis);
                        preferredDelta.set(edgeKey(a, b), d);
                        preferredDelta.set(edgeKey(b, a), { q: -d.q, r: -d.r });
                    }
                }
            }
            const placed = new Map();
            const occupied = new Set();
            const place = (helixId, coord) => {
                const k = axialKey(coord);
                if (occupied.has(k))
                    return false;
                placed.set(helixId, coord);
                occupied.add(k);
                return true;
            };
            const allIds = Array.from({ length: helixCount }, (_, i) => i);
            const roots = [0, ...allIds.filter((h) => h !== 0)].filter((h, i, arr) => h >= 0 && arr.indexOf(h) === i);
            for (const h of allIds) {
                const degree = adjacency.get(h)?.size ?? 0;
                if (degree > 3) {
                    console.warn(`[HelixPos] Helix ${h} has ${degree} neighbors (>3). This may indicate a combinedHelices issue.`);
                }
            }
            for (const root of roots) {
                if (placed.has(root))
                    continue;
                const preferredRoot = projectedAxial.get(root)
                    ? { q: Math.round(projectedAxial.get(root).q), r: Math.round(projectedAxial.get(root).r) }
                    : { q: 0, r: 0 };
                const rootCoord = (root === 0 && !occupied.has(axialKey({ q: 0, r: 0 })))
                    ? { q: 0, r: 0 }
                    : nearestOpenAround(preferredRoot, occupied);
                place(root, rootCoord);
                const queue = [root];
                let qi = 0;
                while (qi < queue.length) {
                    const current = queue[qi++];
                    const currentPos = placed.get(current);
                    if (!currentPos)
                        continue;
                    const neighbors = Array.from((adjacency.get(current) ?? new Map()).entries())
                        .sort((a, b) => b[1] - a[1])
                        .map(([id]) => id);
                    for (const nb of neighbors) {
                        if (placed.has(nb))
                            continue;
                        const pref = preferredDelta.get(edgeKey(current, nb)) ?? HEX_AXIAL_DIRS[0];
                        const base = axialAdd(currentPos, pref);
                        const candidates = HEX_AXIAL_DIRS
                            .map((d) => axialAdd(currentPos, d))
                            .sort((a, b) => {
                            const score = (coord) => {
                                let s = 0;
                                if (coord.q === base.q && coord.r === base.r)
                                    s -= 5;
                                const nbRow = adjacency.get(nb);
                                if (nbRow) {
                                    for (const [p] of nbRow.entries()) {
                                        const placedP = placed.get(p);
                                        if (!placedP)
                                            continue;
                                        const pd = preferredDelta.get(edgeKey(p, nb));
                                        if (!pd)
                                            continue;
                                        const expected = axialAdd(placedP, pd);
                                        s += axialDistance(coord, expected) * 10;
                                    }
                                }
                                const proj = projectedAxial.get(nb);
                                if (proj) {
                                    const dq = coord.q - proj.q;
                                    const dr = coord.r - proj.r;
                                    s += dq * dq + dr * dr;
                                }
                                return s;
                            };
                            return score(a) - score(b);
                        });
                        let placedNow = false;
                        for (const cand of candidates) {
                            if (place(nb, cand)) {
                                placedNow = true;
                                queue.push(nb);
                                break;
                            }
                        }
                        if (!placedNow) {
                            const fallback = nearestOpenAround(base, occupied);
                            if (place(nb, fallback)) {
                                queue.push(nb);
                            }
                        }
                    }
                }
            }
            for (let h = 0; h < helixCount; h++) {
                if (placed.has(h))
                    continue;
                const proj = projectedAxial.get(h);
                const preferred = proj
                    ? { q: Math.round(proj.q), r: Math.round(proj.r) }
                    : { q: 0, r: 0 };
                place(h, nearestOpenAround(preferred, occupied));
            }
            for (let h = 0; h < helixCount; h++) {
                const a = placed.get(h) ?? { q: h, r: 0 };
                const y = axialRToOddQRow(a.q, a.r);
                positions.set(h, [a.q, y]);
            }
            return positions;
        };
        return computeOddQGridPositions(grid, helices);
    }
    toscad.HelixPos = HelixPos;
    /**
     * Finds all crossovers from helix1 to helix2, computes the COM of the 4 nucleotides
     * involved in each crossover, and averages them to return a single 3D vector
     * representing the relative connection from helix1 to helix2.
     */
    function getCrossoverVector(helix1, helix2, helices) {
        // Collect all nucleotides in helix1 into a Set for fast lookup
        // const h1Set = new Set(helices[helix1].map(n => n.id));
        const h2Set = new Set(helices[helix2].map(n => n.id));
        const crossoverVectors = [];
        // Scan all nucleotides in helix 1 to find connections to helix 2
        for (const n1 of helices[helix1]) {
            // Check 5' backbone connection
            if (n1.n5 && n1.n5 instanceof Nucleotide && h2Set.has(n1.n5.id)) {
                // We found a backbone step from helix1 to helix2!
                const n2 = n1.n5;
                const n1pair = n1.pair;
                const n2pair = n2.pair;
                // Ensure it's a true 4-way Holliday Junction crossover (both have pairs)
                if (n1pair && n2pair && n1pair instanceof Nucleotide && n2pair instanceof Nucleotide) {
                    // Center of Helix 1 at this slice
                    const c1 = new THREE.Vector3();
                    c1.addVectors(n1.getPos(), n1pair.getPos()).divideScalar(2);
                    // Center of Helix 2 at this slice
                    const c2 = new THREE.Vector3();
                    c2.addVectors(n2.getPos(), n2pair.getPos()).divideScalar(2);
                    // True vector pointing from Helix 1 core to Helix 2 core
                    const v = new THREE.Vector3().subVectors(c2, c1);
                    crossoverVectors.push(v);
                }
            }
            // Check 3' backbone connection
            if (n1.n3 && n1.n3 instanceof Nucleotide && h2Set.has(n1.n3.id)) {
                const n2 = n1.n3;
                const n1pair = n1.pair;
                const n2pair = n2.pair;
                if (n1pair && n2pair && n1pair instanceof Nucleotide && n2pair instanceof Nucleotide) {
                    // Center of Helix 1 at this slice
                    const c1 = new THREE.Vector3();
                    c1.addVectors(n1.getPos(), n1pair.getPos()).divideScalar(2);
                    // Center of Helix 2 at this slice
                    const c2 = new THREE.Vector3();
                    c2.addVectors(n2.getPos(), n2pair.getPos()).divideScalar(2);
                    // True vector pointing from Helix 1 core to Helix 2 core
                    const v = new THREE.Vector3().subVectors(c2, c1);
                    crossoverVectors.push(v);
                }
            }
            // same as running an average over all 4 nucleotides and then running an average over THOSE vectors
        }
        if (crossoverVectors.length === 0)
            return null;
        // Average all crossover displacement vectors into a final definitive step vector
        const avgVector = new THREE.Vector3(0, 0, 0);
        for (const v of crossoverVectors) {
            avgVector.add(v);
        }
        avgVector.divideScalar(crossoverVectors.length);
        return avgVector;
    }
    toscad.getCrossoverVector = getCrossoverVector;
    function getCrossoverVectorWithBinders(helix1, helix2, helices) {
        // Collect all nucleotides in helix1 into a Set for fast lookup
        // const h1Set = new Set(helices[helix1].map(n => n.id));
        const h2Set = new Set(helices[helix2].map(n => n.id));
        const crossoverVectors = [];
        // For binder-like regions, pair may be missing on one/both sides.
        // Use the nucleotide position itself as fallback so these links still
        // contribute to relative placement.
        const localCenter = (nt) => {
            const p = nt.pair;
            if (p && p instanceof Nucleotide) {
                return nt.getPos().clone().add(p.getPos()).multiplyScalar(0.5);
            }
            return nt.getPos().clone();
        };
        const pushVector = (n1, n2) => {
            const c1 = localCenter(n1);
            const c2 = localCenter(n2);
            crossoverVectors.push(new THREE.Vector3().subVectors(c2, c1));
        };
        // Scan all nucleotides in helix 1 to find connections to helix 2
        for (const n1 of helices[helix1]) {
            // Check 5' backbone connection
            if (n1.n5 && n1.n5 instanceof Nucleotide && h2Set.has(n1.n5.id)) {
                const n2 = n1.n5;
                pushVector(n1, n2);
            }
            // Check 3' backbone connection
            if (n1.n3 && n1.n3 instanceof Nucleotide && h2Set.has(n1.n3.id)) {
                const n2 = n1.n3;
                pushVector(n1, n2);
            }
            // Check direct pair bridge across helices. This captures binder-like
            // attachments where backbone crossover signatures are sparse/absent.
            if (n1.pair && n1.pair instanceof Nucleotide && h2Set.has(n1.pair.id)) {
                const n2 = n1.pair;
                pushVector(n1, n2);
            }
            // same as running an average over all 4 nucleotides and then running an average over THOSE vectors
        }
        if (crossoverVectors.length === 0)
            return null;
        // Average all crossover displacement vectors into a final definitive step vector
        const avgVector = new THREE.Vector3(0, 0, 0);
        for (const v of crossoverVectors) {
            avgVector.add(v);
        }
        avgVector.divideScalar(crossoverVectors.length);
        return avgVector;
    }
    toscad.getCrossoverVectorWithBinders = getCrossoverVectorWithBinders;
    /**
     * Traverses the connections starting from Helix 0, and plots every connected helix
     * onto a 2D coordinate plane locally aligned relative to Helix 0's axis.
     * Returns a Map of HelixId -> { x, y }
     */
    function getRelativePositions(helices) {
        const positions = new Map();
        const visited = new Set();
        const queue = [];
        // Find anchor (Helix 0)
        positions.set(0, { x: 0, y: 0 });
        visited.add(0);
        queue.push(0);
        // 1. Helix long-axis: the direction the helices run along (their "Z").
        const endPts = toscad.helixEndpoints(helices[0]);
        let longAxis = new THREE.Vector3(0, 0, 1);
        if (endPts && endPts.end1 && endPts.end2) {
            longAxis.subVectors(endPts.end1.getPos(), endPts.end2.getPos()).normalize();
        }
        // 2. Build a deterministic 2D basis anchored to world axes, NOT to the
        //    first crossover found.  This ensures the output (x, y) map directly
        //    to world (X, Y) whenever possible, so plotting pos.x → world-X and
        //    pos.y → world-Y is always correct rather than being rotated by an
        //    arbitrary angle depending on BFS order.
        //
        //    Strategy: project world-X onto the plane perpendicular to longAxis.
        //    If longAxis is nearly parallel to world-X, fall back to world-Y,
        //    then world-Z.
        const deriveU0 = (axis) => {
            const candidates = [
                new THREE.Vector3(1, 0, 0),
                new THREE.Vector3(0, 1, 0),
                new THREE.Vector3(0, 0, 1),
            ];
            for (const c of candidates) {
                const proj = c.clone().projectOnPlane(axis);
                if (proj.lengthSq() > 0.01)
                    return proj.normalize();
            }
            return new THREE.Vector3(1, 0, 0); // degenerate fallback
        };
        const u0 = deriveU0(longAxis);
        const u1 = new THREE.Vector3().crossVectors(longAxis, u0).normalize();
        // BFS
        while (queue.length > 0) {
            const curr = queue.shift();
            const currPos = positions.get(curr);
            for (let i = 0; i < helices.length; i++) {
                if (i === curr)
                    continue;
                // getCrossoverVector gives us the true 3D spatial step between core axes
                const vec = getCrossoverVectorWithBinders(curr, i, helices);
                if (vec) { // connection exists
                    // Assign position if we haven't placed this helix yet
                    if (!visited.has(i)) {
                        visited.add(i);
                        queue.push(i);
                        // Flatten the 3D step onto our nice new 2D paper (coordinate basis u0, u1)
                        const dx = vec.dot(u0);
                        const dy = vec.dot(u1);
                        // The new position is simply the parent's position + the flat 2D step!
                        positions.set(i, {
                            x: currPos.x + dx,
                            y: currPos.y + dy
                        });
                    }
                }
            }
        }
        return positions;
    }
    toscad.getRelativePositions = getRelativePositions;
})(toscad || (toscad = {}));
