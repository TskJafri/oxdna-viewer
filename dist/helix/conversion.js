"use strict";
/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
/*
Very important to know the following before you read the code further:
- scadnano files use strands, not individual nucleotides, and each strand requires a direction. Thus, we use "fwd" and "bwd" for directions.
- 'forward'-labeled nt has offsets that increase along its own 5'->3'.
- buildScadnano3 builds strands from 5'->3', and is fully and only influenced by topology and setGrid.

The full canonical pipeline lives in ScadnanoExportManager.prepareScadnanoLayout
(ts/file_handling/scadnano_export.ts). If you need to drive the stages by hand
from the console, the shape is roughly:

    const { helices } = helix.findHelices(elements, 3);
    const { grid, binderHelices } = toscad.setGrid(helices);
    toscad.directionAlign2(grid);
    toscad.alignGridPrim(grid, binderHelices);
    const angles = toscad.getAngles(grid, helices, 'honeycomb');
    // ... anglecomb3 (needs partials + hashMergePairs) → anglecorr2 → ...
    // see prepareScadnanoLayout for the current wiring.
*/
var toscad;
(function (toscad) {
    // Find the two most distant endpoints in the helix using BFS.
    function helixEndpoints(helix) {
        // first we remove duplicates
        // best hope is that there never should be. All of the duplicates must necessarily be removed by findhelix2.ts.
        const nodes = Array.from(new Map(helix.map((n) => [n.id, n])).values());
        if (!nodes.length)
            return null;
        if (nodes.length !== helix.length) {
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
        // Actual iteration
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
    function setGrid(helices, preserveGrid, preservedNtIds, mergedGroups) {
        // Initialize the map
        const grid = new Map();
        // Check if we want to "preserve" any helices.
        if (preserveGrid && preserveGrid.size > 0) {
            // Map each ntId to a helix in helices[][]
            const ntToCurrentHelixId = new Map();
            for (let slotIdx = 0; slotIdx < helices.length; slotIdx++) {
                const slot = helices[slotIdx];
                if (!slot)
                    continue;
                for (const nt of slot) {
                    ntToCurrentHelixId.set(nt.id, slotIdx);
                }
            }
            // For the existing preservedNtIds, remap them to the current helixId and copy them into the new grid.
            for (const [ntId, markData] of preserveGrid.entries()) {
                if (preservedNtIds && !preservedNtIds.has(ntId))
                    continue;
                const currentHelixId = ntToCurrentHelixId.get(ntId);
                if (currentHelixId === undefined)
                    continue; // nt not in current helices
                grid.set(ntId, { ...markData, helixId: currentHelixId });
            }
        }
        // Helper to mark the grid.
        const mark = (nt, helixId, offset, dir) => {
            if (!grid.has(nt.id)) {
                grid.set(nt.id, { helixId, offset, direction: dir });
            }
        };
        // Is the nucleotide inside the helix we ask for?
        const isInHelix = (set, nt) => !!nt && set.has(nt.id);
        // get the pair WITHIN the set.
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
        // Few other subtle differences.
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
            const sides = [
                { label: 'fwd', curr: start.fwd, dir: dirs.fwd, visited: new Set(), otherStartId: start.bwd?.id ?? null },
                { label: 'bwd', curr: start.bwd, dir: dirs.bwd, visited: new Set(), otherStartId: start.fwd?.id ?? null },
            ];
            for (const s of sides)
                if (s.curr)
                    s.visited.add(s.curr.id);
            // Main loop
            for (let steps = 0; steps < 200; steps++) {
                // Advance both sides
                let anyAdvanced = false;
                for (const s of sides) {
                    if (!s.curr)
                        continue;
                    const next = s.curr[s.dir];
                    if (!isInHelix(set, next) || s.visited.has(next.id)) {
                        s.curr = null;
                        continue;
                    }
                    s.curr = next;
                    s.visited.add(s.curr.id);
                    anyAdvanced = true;
                }
                if (!anyAdvanced)
                    return null;
                // Check the added sides for a pair within the set.
                for (const s of sides) {
                    if (!s.curr)
                        continue;
                    const p = getPair(set, s.curr);
                    if (!p)
                        continue;
                    if (s.otherStartId !== null && p.id === s.otherStartId)
                        continue;
                    return {
                        anchor: s.label === 'fwd' ? s.curr : p,
                        steps: steps + 1,
                        source: s.label === 'fwd' ? 'fwd' : 'bwd_pair',
                    };
                }
            }
            return null;
        };
        // Track which helices are binder-only (no internal base-pairing)
        const binderHelices = [];
        // Main loop for setting the grid.
        helices.forEach((helix, helixId) => {
            if (!helix.length)
                return;
            // nucleotide -> helixId
            const helixSet = new Set(helix.map(n => n.id));
            const endpoints = helixEndpoints(helix);
            // Detect binder helix: no nucleotide has a pair *within* the helix
            const hasPairInHelix = helix.some(n => n.pair && n.pair instanceof Nucleotide && helixSet.has(n.pair.id));
            if (!hasPairInHelix) {
                binderHelices.push(helixId);
            }
            // If the helix is a merged helix, find which ones were the originals. Disconnected helices CANNOT be walked by this walker.
            const originGroups = preserveGrid ? mergedGroups?.get(helixId) : undefined;
            const isMergedHelix = !!(originGroups && originGroups.length > 0);
            if (isMergedHelix) {
                for (const group of originGroups) {
                    for (const ntId of group) {
                        if (grid.has(ntId))
                            continue;
                        const prevMark = preserveGrid.get(ntId);
                        if (prevMark)
                            grid.set(ntId, { ...prevMark, helixId });
                    }
                }
            }
            let offset = 0; // Local offset for the main backbone
            // main body of setting the grid. Only runs across non-merged helices.
            if (!isMergedHelix && endpoints) {
                // start "forward" from any endpoint. They will be oriented later. Our main priority is to generate a grid without overlap and sufficient details.
                const helixFwd = endpoints.end1;
                const helixFwdDir = (isInHelix(helixSet, helixFwd.n3) ? 'n3' : 'n5');
                const helixBwdDir = (helixFwdDir === 'n3' ? 'n5' : 'n3');
                const revFwdDir = helixFwdDir === 'n3' ? 'n5' : 'n3';
                const revBwdDir = helixBwdDir === 'n3' ? 'n5' : 'n3';
                const walkLabel = helixFwdDir === 'n3' ? 'forward' : 'backward';
                const pairLabel = helixFwdDir === 'n3' ? 'backward' : 'forward';
                // Find Head
                let firstAnchor = null;
                let firstAnchorPair = null;
                // if it has a pair, set it as a head otherwise find a new anchorpoint.
                // findNextPaired finds an anchorpoint which does have a valid pair within the helix (to anchor the other direction at some offset)
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
                            mark(n, helixId, (startOffset - fwdHeadLen) + i, walkLabel);
                    });
                    [...headBwd].reverse().forEach((n, i) => {
                        if (i < headBwd.length - 1)
                            mark(n, helixId, (startOffset - bwdHeadLen) + i, pairLabel);
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
                    mark(currFwd, helixId, offset, walkLabel);
                if (currBwd)
                    mark(currBwd, helixId, offset, pairLabel);
                // Track visited anchors to break out of circular helices.
                const visitedAnchors = new Set();
                if (firstAnchor)
                    visitedAnchors.add(firstAnchor.id);
                while (currFwd) {
                    const nextStep = findNextPaired({ fwd: currFwd, bwd: currBwd }, { fwd: helixFwdDir, bwd: helixBwdDir }, helixSet);
                    if (!nextStep)
                        break;
                    const nextAnchor = nextStep.anchor;
                    const nextPair = getPair(helixSet, nextAnchor);
                    // probably doesn't need this check but can happen due to cross/double pairing?
                    if (nextAnchor.id === currFwd.id)
                        break;
                    // Circular helix guard: stop if we've already processed this anchor.
                    if (visitedAnchors.has(nextAnchor.id)) {
                        console.log("CIRCULAR STRAND DETECTED");
                        break;
                    }
                    visitedAnchors.add(nextAnchor.id);
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
                    // this backward trace is actually very significant. Without it, cross-pairing becomes a real issue.
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
                    fwdTail.forEach((n, i) => mark(n, helixId, offset + i + 1, walkLabel));
                    const fwdHeadStart = offset + 1 + (gapLength - fwdHead.length);
                    fwdHead.forEach((n, i) => mark(n, helixId, fwdHeadStart + i, walkLabel));
                    bwdTail.forEach((n, i) => mark(n, helixId, offset + i + 1, pairLabel));
                    const bwdHeadStart = offset + 1 + (gapLength - bwdHead.length);
                    bwdHead.forEach((n, i) => mark(n, helixId, bwdHeadStart + i, pairLabel));
                    offset += gapLength + 1;
                    mark(nextAnchor, helixId, offset, walkLabel);
                    if (nextPair)
                        mark(nextPair, helixId, offset, pairLabel);
                    currFwd = nextAnchor;
                    currBwd = nextPair;
                }
                // The Tail
                if (currFwd) {
                    const fwdTail = tracePath(currFwd, helixFwdDir, helixSet).slice(1);
                    fwdTail.forEach((n, i) => mark(n, helixId, offset + i + 1, walkLabel));
                }
                if (currBwd) {
                    const bwdTail = tracePath(currBwd, helixBwdDir, helixSet).slice(1);
                    bwdTail.forEach((n, i) => mark(n, helixId, offset + i + 1, pairLabel));
                }
            }
            // set containing unvisited nucleotides.
            const unvisitedSet = new Set();
            for (const n of helix) {
                if (!grid.has(n.id))
                    unvisitedSet.add(n.id);
            }
            /* The main walker will miss:
                - pure binder helices (binders always forced to face the forward direction)
                - disconnected helices
                - disconnected nucleotides
            Thus, the following part aims to mark those remaining nucleotides into the grid. It iterates over all the helices.
            */
            if (unvisitedSet.size > 0) {
                console.log(`[setGrid] For binders, processing ${unvisitedSet.size} disconnected items on Helix ${helixId}`);
                // Determine "True" Start Offset from Grid State
                let currentBinderOffset = 0;
                let maxFoundOffset = -1;
                for (const n of helix) {
                    const m = grid.get(n.id);
                    if (m && m.helixId === helixId && m.offset > maxFoundOffset) {
                        maxFoundOffset = m.offset;
                    }
                }
                // If the grid has content, start after it. If empty, start at 0.
                if (maxFoundOffset > -1) {
                    currentBinderOffset = maxFoundOffset + 4; // Add visual buffer
                }
                // track the nucleotides that have been dealt with
                const processedExtras = new Set();
                for (const node of helix) {
                    if (!unvisitedSet.has(node.id))
                        continue;
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
        // Note: The alignment of the merged helices is left upto alignMergedGroups()
        return { grid, binderHelices };
    }
    toscad.setGrid = setGrid;
    ;
    // grid flipper. Great helper function for the final output.
    function gridFlip(grid, helixId) {
        // Collect the offset range of the helix
        let min = Infinity;
        let max = -Infinity;
        for (const [, mark] of grid.entries()) {
            if (mark.helixId !== helixId)
                continue;
            if (mark.offset < min)
                min = mark.offset;
            if (mark.offset > max)
                max = mark.offset;
        }
        // Helix has no marks in the grid — nothing to flip.
        if (min === Infinity)
            return;
        for (const [, mark] of grid.entries()) {
            if (mark.helixId !== helixId)
                continue;
            mark.offset = max + min - mark.offset;
            mark.direction = mark.direction === 'forward' ? 'backward' : 'forward';
        }
    }
    toscad.gridFlip = gridFlip;
    // Collect crossovers between different helices.
    // returns aggregate crossover info, as opposed to crossoverNts(), which returns detailed info per crossover.
    // TODO: Allow this to find direction trends between merged helices to flip one helix-part of the merged helix to align with the other merged helix.
    function collectCrossovers(grid) {
        // collect all nucleotides tht have a grid mark. By design, no nucleotide should be skipped from this, so it should be safe to use all nucleotides as a set instead of doing this...
        const allNtIds = new Set();
        for (const [ntId] of grid.entries())
            allNtIds.add(ntId);
        const visited = new Set();
        // crossovers[fromHelix][toHelix] = { sameWalk: n, diffWalk: n }
        //   sameWalk  = both runs have same offset trend -> need flip
        //   diffWalk  = runs have opposite offset trend -> already correct
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
            // This directly helps directionAlign2 flip the helices to align them properly.
            for (let i = 0; i < runs.length - 1; i++) {
                const runA = runs[i];
                const runB = runs[i + 1];
                if (runA.helixId === runB.helixId)
                    continue;
                // Skip the ones with <2 nts.
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
    // Uses collectCrossovers to determine which helices should be flipped to align the directions of all helices in the grid.
    function directionAlign2(grid) {
        // collect crossover and trend info
        const { crossovers, helixIds } = collectCrossovers(grid);
        // track the "anchored" helices as they get aligned.
        const anchored = new Set();
        const flippedHelices = [];
        let edgeCount = 0;
        // see which trends don't match
        for (const [fromHelix, neighbors] of crossovers.entries()) {
            for (const [toHelix, stats] of neighbors.entries()) {
                if (fromHelix >= toHelix)
                    continue;
                if (stats.sameWalk + stats.diffWalk <= 0)
                    continue;
                edgeCount++;
            }
        }
        // helper to apply the flip to the crossover stats.
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
                // sameWalk = both sides increase (or both decrease) -> flip
                // diffWalk = they alternate -> already correct`
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
        // Handle disconnected helices (by doing a separate BFS for each unanchored helix)
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
    // Now that the grid is in place, we can convert to the scadnano format.
    function buildScadnano3(grid, helices, gridType, helixPositions) {
        // scaffold stuff. Also need to assign colors, so scaffold is blue and the staples are red/green/black.
        const scaffoldStrand = helix.getScaffoldStrand();
        const SCAFFOLD_COLOR = '#0066cc';
        const STAPLE_COLORS = ['#f74308', '#57bb00', '#000000'];
        // helix info.
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
        const scadStrands = [];
        const allSystems = [];
        if (Array.isArray(systems))
            for (const s of systems)
                if (s)
                    allSystems.push(s);
        if (typeof tmpSystems !== 'undefined' && Array.isArray(tmpSystems)) {
            for (const s of tmpSystems)
                if (s)
                    allSystems.push(s);
        }
        // walker to build the scadnano strands. This is a 5' -> 3' walker, so it starts at end5 and walks via n3.
        for (const sys of allSystems) {
            const sysStrands = (sys && Array.isArray(sys.strands)) ? sys.strands : [];
            for (const strand of sysStrands) {
                if (!strand)
                    continue;
                // Start at the 5' end of the strand
                const start = strand.end5;
                if (!(start instanceof Nucleotide))
                    continue;
                // Walk this strand 5' -> 3' via n3
                let sequence = '';
                const domains = [];
                // in this case, isCircular is indicating circularity within strand, not the helix. 
                let isCircular = false;
                let openDomain = null;
                const visited = new Set();
                const closeDomain = () => {
                    if (!openDomain)
                        return;
                    domains.push({
                        helix: openDomain.helixId,
                        forward: openDomain.forward,
                        start: openDomain.minOffset,
                        end: openDomain.maxOffset + 1
                    });
                    openDomain = null;
                };
                const openAt = (nt, mark) => {
                    openDomain = {
                        helixId: mark.helixId,
                        direction: mark.direction,
                        step: mark.direction === 'forward' ? 1 : -1,
                        forward: mark.direction === 'forward',
                        minOffset: mark.offset,
                        maxOffset: mark.offset,
                        lastOffset: mark.offset
                    };
                    sequence += nt.type || 'N';
                };
                let curr = start;
                // Walker loop.
                while (curr instanceof Nucleotide) {
                    if (visited.has(curr.id)) {
                        // Walked back to a node we already emitted => circular.
                        isCircular = true;
                        break;
                    }
                    visited.add(curr.id);
                    const mark = grid.get(curr.id);
                    if (!mark) {
                        // Unplaced nt — close any open domain, skip until we land on a placed nt again.
                        closeDomain();
                    }
                    else if (!openDomain) {
                        openAt(curr, mark);
                    }
                    else {
                        const continues = openDomain.helixId === mark.helixId &&
                            openDomain.direction === mark.direction &&
                            openDomain.lastOffset + openDomain.step === mark.offset;
                        if (continues) {
                            openDomain.lastOffset = mark.offset;
                            if (mark.offset < openDomain.minOffset)
                                openDomain.minOffset = mark.offset;
                            if (mark.offset > openDomain.maxOffset)
                                openDomain.maxOffset = mark.offset;
                            sequence += curr.type || 'N';
                        }
                        else {
                            closeDomain();
                            openAt(curr, mark);
                        }
                    }
                    // Advance via n3. Closed-loop strands have n3 of the 3' end
                    // pointing back to end5, so detect that before stepping.
                    const nextRef = curr.n3;
                    if (nextRef instanceof Nucleotide && nextRef.id === start.id) {
                        isCircular = true;
                        break;
                    }
                    curr = (nextRef instanceof Nucleotide) ? nextRef : null;
                }
                closeDomain();
                if (domains.length === 0)
                    continue;
                const isScaffold = scaffoldStrand !== null && strand === scaffoldStrand;
                const color = isScaffold
                    ? SCAFFOLD_COLOR
                    : STAPLE_COLORS[Math.floor(Math.random() * STAPLE_COLORS.length)];
                const out = { color, sequence, domains };
                if (isScaffold)
                    out.is_scaffold = true;
                if (isCircular)
                    out.circular = true;
                scadStrands.push(out);
            }
        }
        return {
            version: '0.20.1',
            grid: gridType,
            helices: scadHelices,
            strands: scadStrands
        };
    }
    toscad.buildScadnano3 = buildScadnano3;
    ;
    // Confirms whether every offset and its direction is unique. 
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
