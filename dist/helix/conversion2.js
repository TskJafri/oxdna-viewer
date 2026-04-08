/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
/*
Here's an easy way to use this code:

    const helices = await honda.findHelices(nucleotideElements, 3);
    const { grid, binderHelices } = toscad.setGrid(helices);
    toscad.directionAlign2(grid);
    toscad.alignGridPrim(grid, binderHelices);
    toscad.combinedHelices(15, grid, helices, binderHelices); // The 15 is arbitrary, just works well for now.
    const { crossovers } = toscad.collectCrossovers(grid);
    let helixPos = toscad.HelixPosByRelativeBfs(grid, helices);
    let gridType = 'honeycomb'; // or 'square'.
    const scadnano = toscad.buildScadnano2(grid, helices, gridType, helixPos);
*/
var toscad;
(function (toscad) {
    function helixEndpoints(helix) {
        // Find the two most distant endpoints in the helix using BFS.
        // first we remove duplicates
        // best hope is that there never should be. All of the duplicates must necessarily be removed by findhelix2.ts.
        const nodes = Array.from(new Map(helix.map((n) => [n.id, n])).values());
        if (!nodes.length)
            return null;
        if (nodes.length !== helix.length) {
            console.log("Holy shit the world is doomed");
            console.log("Just kidding, there are duplicates in the helix");
            console.log("Did you give 1 helix as input or all of them? This function only takes 1.");
            console.warn("Pay attention something went wrong");
        }
        const nodeById = new Map(nodes.map((n) => [n.id, n]));
        const neighbors = (nt) => {
            const list = [];
            const n5 = nt.n5;
            if (n5 instanceof Nucleotide && nodeById.has(n5.id))
                list.push(n5);
            const n3 = nt.n3;
            if (n3 instanceof Nucleotide && nodeById.has(n3.id))
                list.push(n3);
            const pair = nt.pair;
            if (pair instanceof Nucleotide && nodeById.has(pair.id))
                list.push(pair);
            return list;
        };
        // BFS algo mentioned earlier
        const bfs = (start) => {
            // initialize
            const q = [start];
            const dist = new Map();
            dist.set(start.id, 0);
            let far = start;
            for (let i = 0; i < q.length; i++) {
                const cur = q[i];
                const d = dist.get(cur.id);
                if (d === undefined)
                    continue;
                const farthestDist = dist.get(far.id);
                if (farthestDist === undefined || d > farthestDist)
                    far = cur;
                for (const nb of neighbors(cur)) {
                    if (!dist.has(nb.id)) {
                        dist.set(nb.id, d + 1);
                        q.push(nb);
                    }
                }
            }
            return { far, dist };
        };
        const visited = new Set();
        let best = { end1: nodes[0], end2: nodes[0], diameter: 0 };
        for (const n of nodes) {
            if (visited.has(n.id))
                continue;
            const first = bfs(n);
            first.dist.forEach((_, id) => visited.add(id));
            const second = bfs(first.far);
            let farId = first.far.id;
            let diameter = 0;
            second.dist.forEach((d, id) => {
                if (d > diameter) {
                    diameter = d;
                    farId = id;
                }
            });
            if (diameter > best.diameter) {
                const farNode = nodeById.get(farId);
                if (farNode) {
                    best = { end1: first.far, end2: farNode, diameter };
                }
            }
        }
        return best;
    }
    toscad.helixEndpoints = helixEndpoints;
    ;
    function showHelixEndpoints(helices) {
        // const helices = await honda.findHelices(elements, 2);
        const endpoints = helices.map((helix, i) => {
            const res = helixEndpoints(helix);
            return {
                helixIndex: i,
                endpointA: res?.end1?.id,
                endpointB: res?.end2?.id,
                diameter: res?.diameter
            };
        });
        console.log(endpoints);
        return endpoints;
    }
    toscad.showHelixEndpoints = showHelixEndpoints;
    ;
    function crossoverEndpointsHelix(grid, helix, helixId) {
        const nodes = Array.from(new Map(helix.map((nt) => [nt.id, nt])).values());
        if (!nodes.length)
            return null;
        const nodeById = new Map(nodes.map((nt) => [nt.id, nt]));
        const neighbors = (nt) => {
            const list = [];
            const n5 = nt.n5;
            if (n5 instanceof Nucleotide && nodeById.has(n5.id))
                list.push(n5);
            const n3 = nt.n3;
            if (n3 instanceof Nucleotide && nodeById.has(n3.id))
                list.push(n3);
            const pair = nt.pair;
            if (pair instanceof Nucleotide && nodeById.has(pair.id))
                list.push(pair);
            return list;
        };
        const bfsDistances = (start) => {
            const q = [start];
            const dist = new Map();
            dist.set(start.id, 0);
            for (let i = 0; i < q.length; i++) {
                const cur = q[i];
                const d = dist.get(cur.id);
                if (d === undefined)
                    continue;
                for (const nb of neighbors(cur)) {
                    if (!dist.has(nb.id)) {
                        dist.set(nb.id, d + 1);
                        q.push(nb);
                    }
                }
            }
            return dist;
        };
        const crossoverNtIdsByHelix = new Map();
        const ensureSet = (hId) => {
            if (!crossoverNtIdsByHelix.has(hId))
                crossoverNtIdsByHelix.set(hId, new Set());
            return crossoverNtIdsByHelix.get(hId);
        };
        for (const crossover of crossoverNts(grid)) {
            ensureSet(crossover.fromHelix).add(crossover.fromNt.id);
            ensureSet(crossover.toHelix).add(crossover.toNt.id);
        }
        const crossoverIds = crossoverNtIdsByHelix.get(helixId);
        if (!crossoverIds || crossoverIds.size === 0)
            return null;
        const helixEnds = helixEndpoints(nodes);
        if (!helixEnds)
            return null;
        const closestCrossoverFrom = (start) => {
            const dist = bfsDistances(start);
            let bestNt = null;
            let bestDist = Infinity;
            for (const ntId of crossoverIds) {
                const d = dist.get(ntId);
                if (d === undefined)
                    continue;
                if (d < bestDist) {
                    const nt = nodeById.get(ntId);
                    if (!nt)
                        continue;
                    bestDist = d;
                    bestNt = nt;
                }
            }
            return bestNt;
        };
        const end1 = closestCrossoverFrom(helixEnds.end1);
        const end2 = closestCrossoverFrom(helixEnds.end2);
        if (!end1 || !end2)
            return null;
        const diameter = bfsDistances(end1).get(end2.id) ?? 0;
        return { end1, end2, diameter };
    }
    toscad.crossoverEndpointsHelix = crossoverEndpointsHelix;
    function setGrid(helices) {
        const grid = new Map();
        // --- Helpers ---
        const mark = (nt, helixId, offset, dir) => {
            if (!grid.has(nt.id)) {
                grid.set(nt.id, { helixId, offset, direction: dir });
            }
        };
        const isInHelix = (set, nt) => !!nt && set.has(nt.id);
        const getPair = (set, nt) => (nt && nt.pair && set.has(nt.pair.id)) ? nt.pair : null;
        // note: tracePath does NOT include the stopAt nucleotide... 
        // returns a single-segment path in the given direction.
        const tracePath = (start, dir, set, maxSteps = -1, stopAt) => {
            const path = [];
            const visited = new Set();
            let curr = start;
            while (curr && isInHelix(set, curr)) {
                if (visited.has(curr.id))
                    break;
                if (stopAt && curr.id === stopAt.id)
                    break;
                visited.add(curr.id);
                path.push(curr);
                if (maxSteps !== -1 && path.length >= maxSteps)
                    break;
                curr = curr[dir];
            }
            return path;
        };
        // note: this one DOES include the stopAt nucleotide...
        const tracePathWithStop = (start, dir, set, stopAt) => {
            const path = [];
            const visited = new Set();
            let curr = start;
            let safety = 0;
            while (curr && isInHelix(set, curr) && safety++ < 500) {
                if (visited.has(curr.id))
                    break;
                visited.add(curr.id);
                if (curr.id === stopAt.id) {
                    path.push(curr);
                    return { path, reachedStop: true };
                }
                path.push(curr);
                curr = curr[dir];
            }
            return { path, reachedStop: false };
        };
        // same as tracePath but at the end, the path[] is reversed.
        const traceBack = (start, revDir, set, stopAtId, maxSteps) => {
            const path = [];
            let curr = start[revDir];
            let steps = 0;
            while (curr && isInHelix(set, curr) && steps++ < maxSteps) {
                if (curr.id === stopAtId)
                    break;
                path.push(curr);
                curr = curr[revDir];
            }
            return path.reverse();
        };
        const findNextPaired = (start, dirs, set) => {
            let fwdCurr = start.fwd;
            let bwdCurr = start.bwd;
            let steps = 0;
            const visitedFwd = new Set();
            const visitedBwd = new Set();
            if (fwdCurr)
                visitedFwd.add(fwdCurr.id);
            if (bwdCurr)
                visitedBwd.add(bwdCurr.id);
            while (steps < 200) {
                const nextFwd = fwdCurr ? fwdCurr[dirs.fwd] : null;
                const nextBwd = bwdCurr ? bwdCurr[dirs.bwd] : null;
                const validFwd = isInHelix(set, nextFwd) && nextFwd && !visitedFwd.has(nextFwd.id);
                const validBwd = isInHelix(set, nextBwd) && nextBwd && !visitedBwd.has(nextBwd.id);
                if (!validFwd && !validBwd)
                    return null;
                steps++;
                if (validFwd && nextFwd) {
                    fwdCurr = nextFwd;
                    visitedFwd.add(fwdCurr.id);
                }
                else {
                    fwdCurr = null;
                }
                if (validBwd && nextBwd) {
                    bwdCurr = nextBwd;
                    visitedBwd.add(bwdCurr.id);
                }
                else {
                    bwdCurr = null;
                }
                if (fwdCurr) {
                    const p = getPair(set, fwdCurr);
                    if (p && !(start.bwd && p.id === start.bwd.id))
                        return { anchor: fwdCurr, steps, source: 'fwd' };
                }
                if (bwdCurr) {
                    const p = getPair(set, bwdCurr);
                    if (p && !(start.fwd && p.id === start.fwd.id))
                        return { anchor: p, steps, source: 'bwd_pair' };
                }
            }
            return null;
        };
        // Track which helices are binder-only (no internal base-pairing)
        const binderHelices = [];
        // --- Main Loop ---
        helices.forEach((helix, helixId) => {
            if (!helix.length)
                return;
            const helixSet = new Set(helix.map(n => n.id));
            const endpoints = helixEndpoints(helix);
            // Detect binder helix: no nucleotide has a pair within the helix
            const hasPairInHelix = helix.some(n => n.pair && n.pair instanceof Nucleotide && helixSet.has(n.pair.id));
            if (!hasPairInHelix) {
                binderHelices.push(helixId);
            }
            let offset = 0; // Local offset for the main backbone
            // --- A-D. MAIN BACKBONE LOGIC ---
            if (endpoints) {
                // start "forward" from any endpoint. They will be oriented later. Our main priority is to generate a grid without overlap and sufficient details.
                const helixFwd = endpoints.end1;
                const helixFwdDir = (isInHelix(helixSet, helixFwd.n3) ? 'n3' : 'n5');
                const helixBwdDir = (helixFwdDir === 'n3' ? 'n5' : 'n3');
                const revFwdDir = helixFwdDir === 'n3' ? 'n5' : 'n3';
                const revBwdDir = helixBwdDir === 'n3' ? 'n5' : 'n3';
                // Find Head
                let firstAnchor = null;
                let firstAnchorPair = null;
                // if it has a pair, set it as a head otherwise find a new anchorpoint.
                // findNextPaired finds an anchorpoint which does have a valid pair within the helix.
                if (getPair(helixSet, helixFwd)) {
                    firstAnchor = helixFwd;
                }
                else {
                    const result = findNextPaired({ fwd: helixFwd, bwd: null }, { fwd: helixFwdDir, bwd: helixBwdDir }, helixSet);
                    if (result)
                        firstAnchor = result.anchor;
                }
                if (firstAnchor) {
                    // by definition anchor's pair exists.
                    firstAnchorPair = getPair(helixSet, firstAnchor);
                    const headFwd = tracePath(firstAnchor, revFwdDir, helixSet);
                    const headBwd = tracePath(firstAnchorPair, revBwdDir, helixSet);
                    const fwdHeadLen = headFwd.length - 1;
                    const bwdHeadLen = headBwd.length - 1;
                    const startOffset = Math.max(fwdHeadLen, bwdHeadLen);
                    // remember we won't mark the anchor and it's pair, those will be marked to the grid later.
                    [...headFwd].reverse().forEach((n, i) => {
                        if (i < headFwd.length - 1)
                            mark(n, helixId, (startOffset - fwdHeadLen) + i, 'forward');
                    });
                    [...headBwd].reverse().forEach((n, i) => {
                        if (i < headBwd.length - 1)
                            mark(n, helixId, (startOffset - bwdHeadLen) + i, 'backward');
                    });
                    offset = startOffset;
                }
                else {
                    firstAnchor = helixFwd;
                }
                // By here we have the correct offset for the first anchorpoint.
                // The Body
                let currFwd = firstAnchor;
                let currBwd = firstAnchorPair;
                if (currFwd)
                    mark(currFwd, helixId, offset, 'forward');
                if (currBwd)
                    mark(currBwd, helixId, offset, 'backward');
                while (currFwd) {
                    const nextStep = findNextPaired({ fwd: currFwd, bwd: currBwd }, { fwd: helixFwdDir, bwd: helixBwdDir }, helixSet);
                    if (!nextStep)
                        break;
                    const nextAnchor = nextStep.anchor;
                    const nextPair = getPair(helixSet, nextAnchor);
                    // probably doesn't need this check but can happen due to cross/double pairing?
                    if (nextAnchor.id === currFwd.id)
                        break;
                    // Gap Detect
                    const fwdTrace = tracePathWithStop(currFwd, helixFwdDir, helixSet, nextAnchor);
                    // exclude the end points (the pairs themselves)
                    let fwdTail = fwdTrace.path.slice(1);
                    if (fwdTrace.reachedStop)
                        fwdTail.pop();
                    let fwdHead = [];
                    if (!fwdTrace.reachedStop || fwdTail.length === 0) {
                        const potentialHead = traceBack(nextAnchor, revFwdDir, helixSet, currFwd.id, 10);
                        if (potentialHead.length > fwdTail.length) {
                            fwdHead = potentialHead;
                            fwdTail = [];
                        }
                    }
                    // this backward trace is actually very significant. 
                    // without it, cross-pairing becomes a real issue.
                    let bwdTail = [];
                    let bwdHead = [];
                    if (currBwd && nextPair) {
                        const bwdTrace = tracePathWithStop(currBwd, helixBwdDir, helixSet, nextPair);
                        bwdTail = bwdTrace.path.slice(1);
                        if (bwdTrace.reachedStop)
                            bwdTail.pop();
                        if (!bwdTrace.reachedStop || bwdTail.length === 0) {
                            const potentialHead = traceBack(nextPair, revBwdDir, helixSet, currBwd.id, 10);
                            if (potentialHead.length > bwdTail.length) {
                                bwdHead = potentialHead;
                                bwdTail = [];
                            }
                        }
                    }
                    // Fill Grid
                    const gapLength = Math.max(fwdTail.length + fwdHead.length, bwdTail.length + bwdHead.length);
                    fwdTail.forEach((n, i) => mark(n, helixId, offset + i + 1, 'forward'));
                    const fwdHeadStart = offset + 1 + (gapLength - fwdHead.length);
                    fwdHead.forEach((n, i) => mark(n, helixId, fwdHeadStart + i, 'forward'));
                    bwdTail.forEach((n, i) => mark(n, helixId, offset + i + 1, 'backward'));
                    const bwdHeadStart = offset + 1 + (gapLength - bwdHead.length);
                    bwdHead.forEach((n, i) => mark(n, helixId, bwdHeadStart + i, 'backward'));
                    offset += gapLength + 1;
                    mark(nextAnchor, helixId, offset, 'forward');
                    if (nextPair)
                        mark(nextPair, helixId, offset, 'backward');
                    currFwd = nextAnchor;
                    currBwd = nextPair;
                }
                // The Tail
                if (currFwd) {
                    const fwdTail = tracePath(currFwd, helixFwdDir, helixSet).slice(1);
                    fwdTail.forEach((n, i) => mark(n, helixId, offset + i + 1, 'forward'));
                }
                if (currBwd) {
                    const bwdTail = tracePath(currBwd, helixBwdDir, helixSet).slice(1);
                    bwdTail.forEach((n, i) => mark(n, helixId, offset + i + 1, 'backward'));
                }
            }
            // --- E. THE BINDER SWEEPER (Cleanest Version) ---
            // 1. Identify what is missing from the Grid
            const unvisited = helix.filter(n => !grid.has(n.id));
            if (unvisited.length > 0) {
                console.log(`[setGrid] For binders, processing ${unvisited.length} disconnected items on Helix ${helixId}`);
                // 2. Determine "True" Start Offset from Grid State
                let currentBinderOffset = 0;
                let maxFoundOffset = -1;
                for (const val of grid.values()) {
                    if (val.helixId === helixId) {
                        if (val.offset > maxFoundOffset)
                            maxFoundOffset = val.offset;
                    }
                }
                // If the grid has content, start after it. If empty, start at 0.
                if (maxFoundOffset > -1) {
                    currentBinderOffset = maxFoundOffset + 4; // Add visual buffer
                }
                // 3. Group and Grid
                const unvisitedSet = new Set(unvisited.map(n => n.id));
                const processedExtras = new Set();
                for (const node of unvisited) {
                    if (processedExtras.has(node.id))
                        continue;
                    // Trace Back to find segment start
                    let startOfSegment = node;
                    // Walk 5' until we hit something that is NOT in the unvisited set (either in grid or null)
                    let bwdSearch = node.n5;
                    while (bwdSearch && unvisitedSet.has(bwdSearch.id)) {
                        startOfSegment = bwdSearch;
                        bwdSearch = bwdSearch.n5;
                    }
                    // Trace Forward (3') to grab the full segment
                    let currSeg = startOfSegment;
                    while (currSeg && unvisitedSet.has(currSeg.id)) {
                        if (processedExtras.has(currSeg.id))
                            break;
                        processedExtras.add(currSeg.id);
                        mark(currSeg, helixId, currentBinderOffset++, 'forward');
                        currSeg = currSeg.n3;
                    }
                    // Add small gap between distinct binder segments
                    currentBinderOffset += 4;
                }
            }
        });
        return { grid, binderHelices };
    }
    toscad.setGrid = setGrid;
    ;
    // grid flipper. Great helper function for the final output.
    function gridFlip(grid, helixId) {
        const ranges = new Map();
        for (const [, mark] of grid.entries()) {
            if (mark.helixId !== helixId)
                continue;
            const range = ranges.get(mark.helixId);
            if (!range) {
                ranges.set(mark.helixId, { min: mark.offset, max: mark.offset });
                continue;
            }
            if (mark.offset < range.min)
                range.min = mark.offset;
            if (mark.offset > range.max)
                range.max = mark.offset;
        }
        for (const [, mark] of grid.entries()) {
            if (mark.helixId !== helixId)
                continue;
            const range = ranges.get(mark.helixId);
            if (!range)
                continue;
            mark.offset = range.max + range.min - mark.offset;
            mark.direction = mark.direction === 'forward' ? 'backward' : 'forward';
        }
    }
    toscad.gridFlip = gridFlip;
    function getScaffoldStrand() {
        let maxLen = 0;
        let scaffold = null;
        systems.forEach(s => {
            s.strands.forEach(strand => {
                if (strand.getLength() > maxLen) {
                    maxLen = strand.getLength();
                    scaffold = strand;
                }
            });
        });
        return scaffold;
    }
    function collectCrossovers(grid) {
        const allNtIds = new Set();
        for (const [ntId] of grid.entries())
            allNtIds.add(ntId);
        const visited = new Set();
        // crossovers[fromHelix][toHelix] = { sameWalk: n, diffWalk: n }
        //   sameWalk  = both runs have same offset trend → need flip
        //   diffWalk  = runs have opposite offset trend → already correct
        const crossovers = new Map();
        const helixIds = new Set();
        const ensureEntry = (from, to) => {
            if (!crossovers.has(from))
                crossovers.set(from, new Map());
            const inner = crossovers.get(from);
            if (!inner.has(to))
                inner.set(to, { sameWalk: 0, diffWalk: 0 });
            return inner.get(to);
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
            const runs = [];
            let currentRun = null;
            let curr = fivePrime;
            const walkForward = new Set();
            while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                if (walkForward.has(curr.id))
                    break;
                walkForward.add(curr.id);
                visited.add(curr.id);
                const mark = grid.get(curr.id);
                if (mark) {
                    helixIds.add(mark.helixId);
                    if (currentRun && currentRun.helixId === mark.helixId) {
                        currentRun.offsets.push(mark.offset);
                    }
                    else {
                        currentRun = { helixId: mark.helixId, offsets: [mark.offset] };
                        runs.push(currentRun);
                    }
                }
                else {
                    currentRun = null;
                }
                const n3ref = curr.n3;
                curr = (n3ref && n3ref instanceof Nucleotide) ? n3ref : null;
            }
            // Now examine consecutive runs for crossovers
            for (let i = 0; i < runs.length - 1; i++) {
                const runA = runs[i];
                const runB = runs[i + 1];
                if (runA.helixId === runB.helixId)
                    continue;
                // Determine offset trend for each run.
                // For runs with ≥2 nts, compare first and last offset.
                // For single-nt runs, skip (can't determine trend).
                if (runA.offsets.length < 2 && runB.offsets.length < 2)
                    continue;
                // Use the trend near the crossover point:
                // runA trend: compare second-to-last offset to last offset
                // runB trend: compare first offset to second offset
                let trendA = null;
                let trendB = null;
                if (runA.offsets.length >= 2) {
                    const last = runA.offsets[runA.offsets.length - 1];
                    const prev = runA.offsets[runA.offsets.length - 2];
                    trendA = last > prev ? 'inc' : 'dec';
                }
                if (runB.offsets.length >= 2) {
                    const first = runB.offsets[0];
                    const second = runB.offsets[1];
                    trendB = second > first ? 'inc' : 'dec';
                }
                // If we can't determine one side, skip this crossover
                if (!trendA && !trendB)
                    continue;
                // If only one side is known, we still can't compare — skip
                if (!trendA || !trendB)
                    continue;
                const entry = ensureEntry(runA.helixId, runB.helixId);
                const entryRev = ensureEntry(runB.helixId, runA.helixId);
                if (trendA === trendB) {
                    // Same walk direction on both helices → need to flip one
                    entry.sameWalk++;
                    entryRev.sameWalk++;
                }
                else {
                    // Opposite walk direction → already correct
                    entry.diffWalk++;
                    entryRev.diffWalk++;
                }
            }
        }
        return { crossovers, helixIds };
    }
    toscad.collectCrossovers = collectCrossovers;
    ;
    /**
     * directionAlign2 — propagating BFS helix orientation alignment.
     *
     * Uses actual offset trends (increasing vs decreasing along the 5'→3'
     * walk) to determine strand direction on each helix — NOT grid.direction
     * labels (which can be wrong).
     *
     * At each crossover between consecutive runs on different helices,
     * checks whether the offset trend is the same or alternates:
     *   same trend (both increasing or both decreasing) → need to flip one
     *   opposite trend → already correct
     *
     * BFS from helix 0 (anchor). Flip immediately, re-scan, proceed.
     */
    function directionAlign2(grid) {
        // ── Initial scan to discover all helices and crossover stats ────
        const { crossovers, helixIds } = collectCrossovers(grid);
        const anchored = new Set();
        const flippedHelices = [];
        let edgeCount = 0;
        for (const [fromHelix, neighbors] of crossovers.entries()) {
            for (const [toHelix, stats] of neighbors.entries()) {
                if (fromHelix >= toHelix)
                    continue;
                if (stats.sameWalk + stats.diffWalk <= 0)
                    continue;
                edgeCount++;
            }
        }
        const applyFlipToCrossoverStats = (helixId) => {
            const neighbors = crossovers.get(helixId);
            if (!neighbors)
                return;
            for (const [neighborHelix, stats] of neighbors.entries()) {
                const reverseStats = crossovers.get(neighborHelix)?.get(helixId);
                if (!reverseStats)
                    continue;
                const same = stats.sameWalk;
                stats.sameWalk = stats.diffWalk;
                stats.diffWalk = same;
                const reverseSame = reverseStats.sameWalk;
                reverseStats.sameWalk = reverseStats.diffWalk;
                reverseStats.diffWalk = reverseSame;
            }
        };
        // Anchor helix 0
        anchored.add(0);
        const queue = [0];
        let qIdx = 0;
        while (qIdx < queue.length) {
            const currHelix = queue[qIdx++];
            const neighbors = crossovers.get(currHelix);
            if (!neighbors)
                continue;
            for (const [neighborHelix] of neighbors.entries()) {
                if (anchored.has(neighborHelix))
                    continue;
                // Collect votes from ALL anchored helices to this neighbor
                let totalSameWalk = 0;
                let totalDiffWalk = 0;
                for (const anchoredHelix of anchored) {
                    const anchoredNeighbors = crossovers.get(anchoredHelix);
                    if (!anchoredNeighbors)
                        continue;
                    const s = anchoredNeighbors.get(neighborHelix);
                    if (!s)
                        continue;
                    totalSameWalk += s.sameWalk;
                    totalDiffWalk += s.diffWalk;
                }
                // sameWalk = both sides increase (or both decrease) → flip
                // diffWalk = they alternate → already correct
                const shouldFlip = totalSameWalk > totalDiffWalk;
                if (shouldFlip) {
                    gridFlip(grid, neighborHelix);
                    flippedHelices.push(neighborHelix);
                    applyFlipToCrossoverStats(neighborHelix);
                }
                anchored.add(neighborHelix);
                queue.push(neighborHelix);
            }
        }
        // Handle disconnected helices
        for (const hId of helixIds) {
            if (anchored.has(hId))
                continue;
            anchored.add(hId);
            const subQueue = [hId];
            let subIdx = 0;
            while (subIdx < subQueue.length) {
                const currHelix = subQueue[subIdx++];
                const neighbors = crossovers.get(currHelix);
                if (!neighbors)
                    continue;
                for (const [neighborHelix] of neighbors.entries()) {
                    if (anchored.has(neighborHelix))
                        continue;
                    let totalSameWalk = 0;
                    let totalDiffWalk = 0;
                    for (const anchoredHelix of anchored) {
                        const anchoredNeighbors = crossovers.get(anchoredHelix);
                        if (!anchoredNeighbors)
                            continue;
                        const s = anchoredNeighbors.get(neighborHelix);
                        if (!s)
                            continue;
                        totalSameWalk += s.sameWalk;
                        totalDiffWalk += s.diffWalk;
                    }
                    if (totalSameWalk > totalDiffWalk) {
                        gridFlip(grid, neighborHelix);
                        flippedHelices.push(neighborHelix);
                        applyFlipToCrossoverStats(neighborHelix);
                    }
                    anchored.add(neighborHelix);
                    subQueue.push(neighborHelix);
                }
            }
        }
        console.log(`[directionAlign2] Flipped ${flippedHelices.length} helices: [${flippedHelices.sort((a, b) => a - b).join(', ')}]`);
        return {
            flippedHelices: flippedHelices.sort((a, b) => a - b),
            edgeCount
        };
    }
    toscad.directionAlign2 = directionAlign2;
    // Helper function to collect all backbone crossovers with their helix and offset info.
    // Slightly lengthy but quite useful.
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
    /**
     * buildScadnano2 — topology-driven scadnano export.
     *
     * Algorithm:
     *  1. Discover every strand by walking backbone links (n3/n5).
     *     - Find 5' ends (degree-1: has n3 but no n5, or n5 not in elements).
     *     - Walk n3 to build the 5'→3' ordered nucleotide list.
     *     - Handle circular strands (no degree-1 node).
     *  2. For each strand, split into domains whenever the helixId changes.
     *  3. For each domain (contiguous run on one helix):
     *     - start = min(offsets in run)
     *     - end   = max(offsets in run) + 1   (scadnano exclusive end)
     *     - forward = (first 5' offset in run === min offset)
     *  4. Sequence is built in backbone-walk order (guaranteed 5'→3').
     */
    // TODO: needs more testing.
    function buildScadnano2(grid, helices, gridType, helixPositions) {
        // ── Scaffold detection ──────────────────────────────────────────
        const scaffoldStrand = getScaffoldStrand();
        const SCAFFOLD_COLOR = '#0066cc';
        const STAPLE_COLORS = ['#f74308', '#57bb00', '#000000'];
        // ── Helix metadata ──────────────────────────────────────────────
        const helixCount = helices.length || Math.max(0, ...Array.from(grid.values()).map(m => m.helixId + 1));
        const helixMaxOffsets = new Map();
        for (const [, mark] of grid.entries()) {
            const current = helixMaxOffsets.get(mark.helixId) ?? -1;
            if (mark.offset > current)
                helixMaxOffsets.set(mark.helixId, mark.offset);
        }
        const scadHelices = Array.from({ length: helixCount }, (_, i) => ({
            max_offset: (helixMaxOffsets.get(i) ?? 0) + 1,
            grid_position: helixPositions?.get(i) ?? [0, i]
        }));
        // ── Step 1: Discover all strands via backbone topology ──────────
        // Build a set of all nucleotide ids that exist in the grid so we
        // only emit nucleotides that were actually placed.
        const allNtIds = new Set();
        for (const [ntId] of grid.entries()) {
            allNtIds.add(ntId);
        }
        // Track which nucleotides have been assigned to a strand already.
        const visited = new Set();
        // We'll collect strand data here.
        const scadStrands = [];
        // Iterate over every nucleotide in the grid and discover strands.
        for (const [ntId] of grid.entries()) {
            if (visited.has(ntId))
                continue;
            const startNt = elements.get(ntId);
            if (!startNt || !(startNt instanceof Nucleotide))
                continue;
            // ── 1a. Find the 5' end of this strand ──────────────────────
            // Walk n5 until we can't anymore (the node with no n5, or
            // whose n5 is not in the grid, is the 5' end).
            let fivePrime = startNt;
            const walkBack = new Set();
            walkBack.add(fivePrime.id);
            while (true) {
                const prev = fivePrime.n5;
                if (!prev || !(prev instanceof Nucleotide))
                    break;
                if (!allNtIds.has(prev.id))
                    break; // not in grid
                if (walkBack.has(prev.id))
                    break; // circular — stop
                walkBack.add(prev.id);
                fivePrime = prev;
            }
            // Detect circular: if fivePrime still has a valid n5 that
            // we stopped on because of the visited guard, it's circular.
            const n5OfFive = fivePrime.n5;
            const isCircular = n5OfFive instanceof Nucleotide &&
                allNtIds.has(n5OfFive.id) &&
                walkBack.has(n5OfFive.id);
            // ── 1b. Walk n3 from 5' end to build ordered nt list ────────
            const orderedNts = [];
            let curr = fivePrime;
            const walkForward = new Set();
            while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                if (walkForward.has(curr.id))
                    break; // full circle
                walkForward.add(curr.id);
                visited.add(curr.id);
                orderedNts.push(curr);
                const n3ref = curr.n3;
                curr = (n3ref && n3ref instanceof Nucleotide) ? n3ref : null;
            }
            if (orderedNts.length === 0)
                continue;
            const runs = [];
            let currentRun = null;
            for (const nt of orderedNts) {
                const mark = grid.get(nt.id);
                if (!mark) {
                    currentRun = null;
                    continue;
                }
                if (currentRun &&
                    currentRun.helixId === mark.helixId &&
                    currentRun.direction === mark.direction) {
                    currentRun.nts.push(nt);
                }
                else {
                    // New helix or new direction → new run
                    currentRun = { helixId: mark.helixId, direction: mark.direction, nts: [nt] };
                    runs.push(currentRun);
                }
            }
            // ── Step 3: Convert runs into scadnano domains ──────────────
            // Within a single run, offsets must be contiguous for a valid
            // scadnano domain ([start, end) claims every position in that
            // range). If there are gaps, split into sub-runs so each
            // sub-run has perfectly contiguous offsets.
            let sequence = '';
            const domains = [];
            for (const run of runs) {
                // Collect (offset, base, nt) tuples in walk order
                const entries = [];
                for (const nt of run.nts) {
                    const mark = grid.get(nt.id);
                    entries.push({ offset: mark.offset, base: nt.type || 'N', nt });
                }
                // Determine overall walk direction for this run:
                // forward = 5' end is at the smaller offset
                const firstOff = entries[0].offset;
                const lastOff = entries[entries.length - 1].offset;
                const forward = firstOff <= lastOff; // increasing or single-nt
                // Split into contiguous sub-runs.
                // Walk entries in order; a sub-run breaks when the next
                // offset isn't exactly ±1 from the previous.
                const step = forward ? 1 : -1;
                const subRuns = [];
                let currentSub = [entries[0]];
                for (let i = 1; i < entries.length; i++) {
                    const prev = entries[i - 1].offset;
                    const curr = entries[i].offset;
                    if (curr === prev + step) {
                        currentSub.push(entries[i]);
                    }
                    else {
                        subRuns.push(currentSub);
                        currentSub = [entries[i]];
                    }
                }
                subRuns.push(currentSub);
                // Emit a domain for each contiguous sub-run
                for (const sub of subRuns) {
                    const minOff = Math.min(sub[0].offset, sub[sub.length - 1].offset);
                    const maxOff = Math.max(sub[0].offset, sub[sub.length - 1].offset);
                    // Sequence in 5'→3' walk order (already correct)
                    for (const e of sub) {
                        sequence += e.base;
                    }
                    domains.push({
                        helix: run.helixId,
                        forward,
                        start: minOff,
                        end: maxOff + 1 // exclusive end
                    });
                }
            }
            if (domains.length > 0) {
                // Determine if this strand is the scaffold
                const isScaffold = scaffoldStrand !== null &&
                    fivePrime.strand === scaffoldStrand;
                const color = isScaffold
                    ? SCAFFOLD_COLOR
                    : STAPLE_COLORS[Math.floor(Math.random() * STAPLE_COLORS.length)];
                const strandObj = { color, sequence, domains };
                if (isScaffold) {
                    strandObj.is_scaffold = true;
                }
                if (isCircular) {
                    strandObj.circular = true;
                }
                scadStrands.push(strandObj);
            }
        }
        return {
            version: '0.20.1',
            grid: gridType,
            helices: scadHelices,
            strands: scadStrands
        };
    }
    toscad.buildScadnano2 = buildScadnano2;
    ;
    // Confirms whether every offset -> direction is unique. 
    function validateGrid(grid) {
        // Structure: Map<HelixID, { forward: Map<Offset, NtID>, backward: Map<Offset, NtID> }>
        const checkMap = new Map();
        let conflicts = 0;
        for (const [ntId, pos] of grid.entries()) {
            // 1. Initialize Helix Bucket if missing
            if (!checkMap.has(pos.helixId)) {
                checkMap.set(pos.helixId, {
                    forward: new Map(),
                    backward: new Map()
                });
            }
            const helixBuckets = checkMap.get(pos.helixId);
            const strandMap = helixBuckets[pos.direction];
            // 2. Check for collision
            if (strandMap.has(pos.offset)) {
                const existingNt = strandMap.get(pos.offset);
                console.error(`❌ CONFLICT DETECTED:\n` +
                    `   Helix: ${pos.helixId}\n` +
                    `   Strand: ${pos.direction}\n` +
                    `   Offset: ${pos.offset}\n` +
                    `   Fighting Nucleotides: IDs ${existingNt} vs ${ntId}`);
                conflicts++;
            }
            else {
                // 3. Register valid position
                strandMap.set(pos.offset, ntId);
            }
        }
        if (conflicts === 0) {
            console.log(`✅ Grid Validated: ${grid.size} nucleotides assigned with 0 overlapping offsets.`);
        }
        else {
            console.warn(`⚠️ Grid Validation Failed: Found ${conflicts} offset collisions.`);
        }
        if (grid.size !== elements.size) {
            console.log("⚠️ INVALID. Grid does not include all nucleotides from the original element set.");
        }
    }
    toscad.validateGrid = validateGrid;
})(toscad || (toscad = {}));
