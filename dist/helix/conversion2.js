/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
var toscad;
(function (toscad) {
    function helixEndpoints(helices) {
        // Find the two most distant endpoints in the helix using BFS.
        // first we remove duplicates
        // best hope is that there never should be. All of the duplicates must necessarily be removed by findhelix2.ts.
        const nodes = Array.from(new Map(helices.map((n) => [n.id, n])).values());
        if (!nodes.length)
            return null;
        if (nodes.length !== helices.length) {
            console.log("Holy shit the world is doomed");
            console.log("Just kidding, there are duplicates in the helices");
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
    toscad.getScaffoldStrand = getScaffoldStrand;
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
    // TODO: needs more testing.
    function directionAlign2(grid) {
        // ── Initial scan to discover all helices ────────────────────────
        const { helixIds } = collectCrossovers(grid);
        const anchored = new Set();
        const flippedHelices = [];
        // Anchor helix 0
        anchored.add(0);
        const queue = [0];
        let qIdx = 0;
        // Unoptimized, brute force code.
        // Every flip, it rescans all crossovers/
        // Need to make it more efficient.
        while (qIdx < queue.length) {
            const currHelix = queue[qIdx++];
            // Re-scan crossovers from the CURRENT grid state (post-flips)
            const { crossovers } = collectCrossovers(grid);
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
                const { crossovers } = collectCrossovers(grid);
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
                    }
                    anchored.add(neighborHelix);
                    subQueue.push(neighborHelix);
                }
            }
        }
        console.log(`[directionAlign2] Flipped ${flippedHelices.length} helices: [${flippedHelices.sort((a, b) => a - b).join(', ')}]`);
        return {
            flippedHelices: flippedHelices.sort((a, b) => a - b),
            edgeCount: 0
        };
    }
    toscad.directionAlign2 = directionAlign2;
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
        console.log(`[alignGridPrim] Aligned ${helixList.length} helices. Shifts:`, Object.fromEntries(Array.from(cumulativeShift.entries()).sort((a, b) => a[0] - b[0])));
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
            const { crossovers, helixIds } = collectCrossovers(grid);
            const adjacency = buildAdjacency(crossovers);
            const offsetSets = buildOffsetSets();
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
                const compatibleWithHub = neighbors.filter(n => offsetsDisjoint(offsetSets.get(hub), offsetSets.get(n))
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
                            if (!offsetsDisjoint(offsetSets.get(a), offsetSets.get(b)))
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
            const { crossovers, helixIds } = collectCrossovers(grid);
            const adjacency = buildAdjacency(crossovers);
            const offsetSets = buildOffsetSets();
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
                    if (!offsetsDisjoint(offsetSets.get(a), offsetSets.get(b)))
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
        const averageHelixCenter = (helix) => {
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
        const estimateHexBasis3D = (allHelices) => {
            const centers = [];
            const axisSamples = [];
            for (const helix of allHelices) {
                if (!helix || helix.length === 0)
                    continue;
                const center = averageHelixCenter(helix);
                if (center)
                    centers.push(center);
                const ep = helixEndpoints(helix);
                if (!ep)
                    continue;
                const axis = ep.end2.getPos().clone().sub(ep.end1.getPos());
                if (axis.lengthSq() > 1e-8)
                    axisSamples.push(axis.normalize());
            }
            if (centers.length === 0)
                return null;
            const origin = new THREE.Vector3();
            for (const c of centers)
                origin.add(c);
            origin.divideScalar(centers.length);
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
            const projectedCenters = centers.map((c) => {
                const rel = c.clone().sub(origin);
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
                const c = averageHelixCenter(helix);
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
    // ── Relative Position Calculation Methods ─────────────────────────
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
        // 1. Z-axis (Normal): The physical direction of Helix 0 itself.
        const endPts = helixEndpoints(helices[0]);
        let longAxis = new THREE.Vector3(0, 0, 1);
        if (endPts && endPts.end1 && endPts.end2) {
            longAxis.subVectors(endPts.end1.getPos(), endPts.end2.getPos()).normalize();
        }
        // We need to find the first valid neighbor connection to establish the X-axis (u0)
        let u0 = null;
        let u1 = null;
        // BFS
        while (queue.length > 0) {
            const curr = queue.shift();
            const currPos = positions.get(curr);
            for (let i = 0; i < helices.length; i++) {
                if (i === curr)
                    continue;
                // getCrossoverVector gives us the true 3D spatial step between core axes
                const vec = getCrossoverVector(curr, i, helices);
                if (vec) { // connection exists
                    // If we haven't established our flat 2D plane yet, do it on the very first connection!
                    if (!u0 || !u1) {
                        // Project the crossover vector so it's perfectly orthogonal to Helix 0's Z-axis
                        const proj = vec.clone().projectOnPlane(longAxis);
                        u0 = proj.clone().normalize();
                        u1 = new THREE.Vector3().crossVectors(longAxis, u0).normalize();
                    }
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
