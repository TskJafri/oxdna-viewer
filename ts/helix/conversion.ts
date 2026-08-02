/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />

/*
Here's an easy way to use this code:

    const {helices} = helix.findHelices(elements,3)

    const { grid, binderHelices } = toscad.setGrid(helices);
    toscad.directionAlign2(grid);
    toscad.alignGridPrim(grid, binderHelices);

    const angles = toscad.getAngles(grid, helices, 'honeycomb');
    const corrected = toscad.anglecomb2(grid, helices, 'honeycomb', angles);
    const correct = toscad.anglecorr2(grid, helices, 'honeycomb', corrected.networkMap);
*/

namespace toscad {
    export function helixEndpoints(helix: Nucleotide[]) {
        // Find the two most distant endpoints in the helix using BFS.

        // first we remove duplicates
        // best hope is that there never should be. All of the duplicates must necessarily be removed by findhelix2.ts.
        const nodes: Nucleotide[] = Array.from(new Map<number, Nucleotide>(helix.map((n: Nucleotide) => [n.id, n])).values());
        if (!nodes.length) return null;
        if (nodes.length !== helix.length) {
            console.log("Did you give 1 helix as input or all of them? This function only takes 1.")
            console.warn("Pay attention something went wrong")
        }

        const nodeById = new Map<number, Nucleotide>(nodes.map((n: Nucleotide) => [n.id, n]));

        const neighbors = (nt: Nucleotide) => {
            const list: Nucleotide[] = [];
            const n5 = nt.n5; if (n5 instanceof Nucleotide && nodeById.has(n5.id)) list.push(n5);
            const n3 = nt.n3; if (n3 instanceof Nucleotide && nodeById.has(n3.id)) list.push(n3);
            const pair = nt.pair; if (pair instanceof Nucleotide && nodeById.has(pair.id)) list.push(pair);
            return list;
        };

        // BFS algo mentioned earlier
        const bfs = (start: Nucleotide) => {
            // initialize
            const q: Nucleotide[] = [start];
            const dist = new Map<number, number>();
            dist.set(start.id, 0);
            let far = start;

            for (let i = 0; i < q.length; i++) {
                const cur = q[i];
                const d = dist.get(cur.id);
                if (d === undefined) continue;

                const farthestDist = dist.get(far.id);
                if (farthestDist === undefined || d > farthestDist) far = cur;

                for (const nb of neighbors(cur)) {
                    if (!dist.has(nb.id)) {
                        dist.set(nb.id, d + 1);
                        q.push(nb);
                    }
                }
            }
            return { far, dist };
        };

        const visited = new Set<number>();
        let best = { end1: nodes[0], end2: nodes[0], diameter: 0 };

        for (const n of nodes) {
            if (visited.has(n.id)) continue;
            const first = bfs(n);
            first.dist.forEach((_, id) => visited.add(id));

            const second = bfs(first.far);
            let farId = first.far.id;
            let diameter = 0;
            second.dist.forEach((d, id) => {
                if (d > diameter) { diameter = d; farId = id; }
            });

            if (diameter > best.diameter) {
                const farNode = nodeById.get(farId);
                if (farNode) {
                    best = { end1: first.far, end2: farNode, diameter };
                }
            }
        }
        return best;
    };

    export function showHelixEndpoints(helices: Nucleotide[][]) {
        // const helices = await helix.findHelices(elements, 2);
        const endpoints = helices.map((helix, i) => {
            const res = helixEndpoints(helix);
            return {
                helixIndex: i,
                endpointA: res?.end1?.id,
                endpointB: res?.end2?.id,
                diameter: res?.diameter
            };
        })
        console.log(endpoints);
        return endpoints;
    };

    export type GridMark = { helixId: number; offset: number; direction: 'forward' | 'backward' };
    export type GridMap = Map<number, GridMark>;
    export type Direction = 'n3' | 'n5';

    export function crossoverEndpointsHelix(
        grid: GridMap,
        helix: Nucleotide[],
        helixId: number
    ) {
        const nodes: Nucleotide[] = Array.from(
            new Map<number, Nucleotide>(helix.map((nt: Nucleotide) => [nt.id, nt])).values()
        );
        if (!nodes.length) return null;

        const nodeById = new Map<number, Nucleotide>(nodes.map((nt: Nucleotide) => [nt.id, nt]));
        const neighbors = (nt: Nucleotide): Nucleotide[] => {
            const list: Nucleotide[] = [];
            const n5 = nt.n5; if (n5 instanceof Nucleotide && nodeById.has(n5.id)) list.push(n5);
            const n3 = nt.n3; if (n3 instanceof Nucleotide && nodeById.has(n3.id)) list.push(n3);
            const pair = nt.pair; if (pair instanceof Nucleotide && nodeById.has(pair.id)) list.push(pair);
            return list;
        };

        const bfsDistances = (start: Nucleotide) => {
            const q: Nucleotide[] = [start];
            const dist = new Map<number, number>();
            dist.set(start.id, 0);

            for (let i = 0; i < q.length; i++) {
                const cur = q[i];
                const d = dist.get(cur.id);
                if (d === undefined) continue;

                for (const nb of neighbors(cur)) {
                    if (!dist.has(nb.id)) {
                        dist.set(nb.id, d + 1);
                        q.push(nb);
                    }
                }
            }

            return dist;
        };

        const crossoverNtIdsByHelix = new Map<number, Set<number>>();
        const ensureSet = (hId: number) => {
            if (!crossoverNtIdsByHelix.has(hId)) crossoverNtIdsByHelix.set(hId, new Set<number>());
            return crossoverNtIdsByHelix.get(hId)!;
        };

        for (const crossover of crossoverNts(grid)) {
            ensureSet(crossover.fromHelix).add(crossover.fromNt.id);
            ensureSet(crossover.toHelix).add(crossover.toNt.id);
        }

        const crossoverIds = crossoverNtIdsByHelix.get(helixId);
        if (!crossoverIds || crossoverIds.size === 0) return null;

        const helixEnds = helixEndpoints(nodes);
        if (!helixEnds) return null;

        const closestCrossoverFrom = (start: Nucleotide): Nucleotide | null => {
            const dist = bfsDistances(start);
            let bestNt: Nucleotide | null = null;
            let bestDist = Infinity;

            for (const ntId of crossoverIds) {
                const d = dist.get(ntId);
                if (d === undefined) continue;
                if (d < bestDist) {
                    const nt = nodeById.get(ntId);
                    if (!nt) continue;
                    bestDist = d;
                    bestNt = nt;
                }
            }

            return bestNt;
        };

        const end1 = closestCrossoverFrom(helixEnds.end1);
        const end2 = closestCrossoverFrom(helixEnds.end2);
        if (!end1 || !end2) return null;

        const diameter = bfsDistances(end1).get(end2.id) ?? 0;
        return { end1, end2, diameter };
    }


    // Helper interface to keep track of where we are on both strands
    interface DualCursor {
        fwd: Nucleotide | null;
        bwd: Nucleotide | null;
    }

    export function setGrid(
        helices: Nucleotide[][],
        preserveGrid?: GridMap,
        preservedNtIds?: Set<number>,
        mergedGroups?: Map<number, number[][]>
    ): { grid: GridMap; binderHelices: number[] } {
        // mergedGroups: merge provenance from the previous pipeline pass.
        // Key = current helix slot; value = one nt-id array per pre-merge
        // helix (an "origin-group"). Helices listed here are NOT re-walked —
        // the backbone walk can't handle a helix made of disconnected
        // components. Instead each origin-group keeps its internal grid marks
        // from preserveGrid, and after the main loop the groups are re-aligned
        // relative to each other (and the whole lattice) via alignGridPrim.
        const grid: GridMap = new Map();

        // If a previous grid is provided for preservation, copy its marks into
        // the new grid — but ONLY for nucleotides in preservedNtIds (if given).
        // This lets helices that haven't been merged keep their grid positions,
        // while nucleotides in merged/altered helices get fresh assignments.
        // 
        // IMPORTANT: remap the helixId of each preserved mark to match the
        // current helices array slot. The preserved mark's helixId is from the
        // previous iteration's slot ordering, which may differ from the current
        // ordering after renumbering.
        if (preserveGrid && preserveGrid.size > 0) {
            // Build a lookup: ntId → current helixId (slot index)
            const ntToCurrentHelixId = new Map<number, number>();
            for (let slotIdx = 0; slotIdx < helices.length; slotIdx++) {
                const slot = helices[slotIdx];
                if (!slot) continue;
                for (const nt of slot) {
                    ntToCurrentHelixId.set(nt.id, slotIdx);
                }
            }

            for (const [ntId, markData] of preserveGrid.entries()) {
                if (preservedNtIds && !preservedNtIds.has(ntId)) continue;
                const currentHelixId = ntToCurrentHelixId.get(ntId);
                if (currentHelixId === undefined) continue; // nt not in current helices
                grid.set(ntId, { ...markData, helixId: currentHelixId });
            }
        }

        // --- Helpers ---
        const mark = (nt: Nucleotide, helixId: number, offset: number, dir: 'forward' | 'backward') => {
            if (!grid.has(nt.id)) {
                grid.set(nt.id, { helixId, offset, direction: dir });
            }
        };

        const isInHelix = (set: Set<number>, nt: Nucleotide | null): nt is Nucleotide => !!nt && set.has(nt.id);

        const getPair = (set: Set<number>, nt: Nucleotide | null): Nucleotide | null =>
            (nt && nt.pair && set.has(nt.pair.id)) ? (nt.pair as Nucleotide) : null;

        // note: tracePath does NOT include the stopAt nucleotide... 
        // returns a single-segment path in the given direction.
        const tracePath = (start: Nucleotide, dir: Direction, set: Set<number>, maxSteps: number = -1, stopAt?: Nucleotide): Nucleotide[] => {
            const path: Nucleotide[] = [];
            const visited = new Set<number>();
            let curr: Nucleotide | null = start;
            while (curr && isInHelix(set, curr)) {
                if (visited.has(curr.id)) break;
                if (stopAt && curr.id === stopAt.id) break;
                visited.add(curr.id);
                path.push(curr);
                if (maxSteps !== -1 && path.length >= maxSteps) break;
                curr = curr[dir] as Nucleotide | null;
            }
            return path;
        };

        // note: this one DOES include the stopAt nucleotide...
        const tracePathWithStop = (start: Nucleotide, dir: Direction, set: Set<number>, stopAt: Nucleotide): { path: Nucleotide[]; reachedStop: boolean } => {
            const path: Nucleotide[] = [];
            const visited = new Set<number>();
            let curr: Nucleotide | null = start;
            let safety = 0;
            while (curr && isInHelix(set, curr) && safety++ < 500) {
                if (visited.has(curr.id)) break;
                visited.add(curr.id);
                if (curr.id === stopAt.id) {
                    path.push(curr);
                    return { path, reachedStop: true };
                }
                path.push(curr);
                curr = curr[dir] as Nucleotide | null;
            }
            return { path, reachedStop: false };
        };

        // same as tracePath but at the end, the path[] is reversed.
        const traceBack = (start: Nucleotide, revDir: Direction, set: Set<number>, stopAtId: number, maxSteps: number): Nucleotide[] => {
            const path: Nucleotide[] = [];
            let curr: Nucleotide | null = start[revDir] as Nucleotide | null;
            let steps = 0;
            while (curr && isInHelix(set, curr) && steps++ < maxSteps) {
                if (curr.id === stopAtId) break;
                path.push(curr);
                curr = curr[revDir] as Nucleotide | null;
            }
            return path.reverse();
        };

        const findNextPaired = (start: DualCursor, dirs: { fwd: Direction, bwd: Direction }, set: Set<number>): { anchor: Nucleotide, steps: number, source: string } | null => {
            let fwdCurr = start.fwd;
            let bwdCurr = start.bwd;
            let steps = 0;
            const visitedFwd = new Set<number>();
            const visitedBwd = new Set<number>();
            if (fwdCurr) visitedFwd.add(fwdCurr.id);
            if (bwdCurr) visitedBwd.add(bwdCurr.id);

            while (steps < 200) {
                const nextFwd = fwdCurr ? (fwdCurr[dirs.fwd] as Nucleotide | null) : null;
                const nextBwd = bwdCurr ? (bwdCurr[dirs.bwd] as Nucleotide | null) : null;
                const validFwd = isInHelix(set, nextFwd) && nextFwd && !visitedFwd.has(nextFwd.id);
                const validBwd = isInHelix(set, nextBwd) && nextBwd && !visitedBwd.has(nextBwd.id);

                if (!validFwd && !validBwd) return null;
                steps++;
                if (validFwd && nextFwd) { fwdCurr = nextFwd; visitedFwd.add(fwdCurr.id); } else { fwdCurr = null; }
                if (validBwd && nextBwd) { bwdCurr = nextBwd; visitedBwd.add(bwdCurr.id); } else { bwdCurr = null; }

                if (fwdCurr) {
                    const p = getPair(set, fwdCurr);
                    if (p && !(start.bwd && p.id === start.bwd.id)) return { anchor: fwdCurr, steps, source: 'fwd' };
                }
                if (bwdCurr) {
                    const p = getPair(set, bwdCurr);
                    if (p && !(start.fwd && p.id === start.fwd.id)) return { anchor: p, steps, source: 'bwd_pair' };
                }
            }
            return null;
        };

        // Track which helices are binder-only (no internal base-pairing)
        const binderHelices: number[] = [];

        // --- Main Loop ---
        helices.forEach((helix, helixId) => {
            if (!helix.length) return;
            const helixSet = new Set(helix.map(n => n.id));
            const endpoints = helixEndpoints(helix);

            // Detect binder helix: no nucleotide has a pair within the helix
            const hasPairInHelix = helix.some(n =>
                n.pair && n.pair instanceof Nucleotide && helixSet.has(n.pair.id)
            );
            if (!hasPairInHelix) {
                binderHelices.push(helixId);
            }

            // --- MERGED-HELIX FAST PATH ---
            // This helix is the product of merges from the previous pipeline
            // pass. Do NOT re-walk it — the walk can't handle disconnected
            // components. Carry forward each origin-group's internal marks
            // from preserveGrid (remapped to this slot). The groups are
            // re-aligned relative to each other in the post-pass after this
            // loop. Any nt missing from preserveGrid falls through to the
            // binder sweeper below.
            const originGroups = preserveGrid ? mergedGroups?.get(helixId) : undefined;
            const isMergedHelix = !!(originGroups && originGroups.length > 0);
            if (isMergedHelix) {
                for (const group of originGroups!) {
                    for (const ntId of group) {
                        if (grid.has(ntId)) continue;
                        const prevMark = preserveGrid!.get(ntId);
                        if (prevMark) grid.set(ntId, { ...prevMark, helixId });
                    }
                }
            }

            let offset = 0; // Local offset for the main backbone

            // --- A-D. MAIN BACKBONE LOGIC ---
            if (!isMergedHelix && endpoints) {
                // start "forward" from any endpoint. They will be oriented later. Our main priority is to generate a grid without overlap and sufficient details.
                const helixFwd = endpoints.end1;
                const helixFwdDir = (isInHelix(helixSet, helixFwd.n3 as Nucleotide) ? 'n3' : 'n5');
                const helixBwdDir = (helixFwdDir === 'n3' ? 'n5' : 'n3');
                const revFwdDir = helixFwdDir === 'n3' ? 'n5' : 'n3';
                const revBwdDir = helixBwdDir === 'n3' ? 'n5' : 'n3';

                // Convention: a 'forward'-labeled nt has offsets that increase
                // along its own 5'→3'. This is what buildScadnano3 assumes
                // (step = +1 for forward) and what directionAlign2's trend
                // criterion needs in order to be equivalent to label
                // alternation at crossovers.
                //
                // We stamp offsets in the helixFwdDir direction. When
                // helixFwdDir === 'n3' that direction is the walked strand's
                // own 5'→3', so the walked strand naturally has forward=inc
                // and gets the 'forward' label. When helixFwdDir === 'n5'
                // we're stamping along the walked strand's 3'→5', meaning
                // its offsets decrease along its own 5'→3'; the *pair*
                // strand is the one whose offsets increase along its 5'→3'.
                // So we swap the labels in that case.
                const walkLabel: 'forward' | 'backward' = helixFwdDir === 'n3' ? 'forward' : 'backward';
                const pairLabel: 'forward' | 'backward' = helixFwdDir === 'n3' ? 'backward' : 'forward';

                // Find Head
                let firstAnchor: Nucleotide | null = null;
                let firstAnchorPair: Nucleotide | null = null;
                // if it has a pair, set it as a head otherwise find a new anchorpoint.
                // findNextPaired finds an anchorpoint which does have a valid pair within the helix.
                if (getPair(helixSet, helixFwd)) {
                    firstAnchor = helixFwd;
                } else {
                    const result = findNextPaired({ fwd: helixFwd, bwd: null }, { fwd: helixFwdDir, bwd: helixBwdDir }, helixSet);
                    if (result) firstAnchor = result.anchor;
                }

                if (firstAnchor) {
                    // by definition anchor's pair exists.
                    firstAnchorPair = getPair(helixSet, firstAnchor);
                    const headFwd = tracePath(firstAnchor, revFwdDir, helixSet);
                    const headBwd = tracePath(firstAnchorPair!, revBwdDir, helixSet);
                    const fwdHeadLen = headFwd.length - 1;
                    const bwdHeadLen = headBwd.length - 1;
                    const startOffset = Math.max(fwdHeadLen, bwdHeadLen);

                    // remember we won't mark the anchor and it's pair, those will be marked to the grid later.
                    [...headFwd].reverse().forEach((n, i) => {
                        if (i < headFwd.length - 1) mark(n, helixId, (startOffset - fwdHeadLen) + i, walkLabel);
                    });
                    [...headBwd].reverse().forEach((n, i) => {
                        if (i < headBwd.length - 1) mark(n, helixId, (startOffset - bwdHeadLen) + i, pairLabel);
                    });
                    offset = startOffset;
                } else {
                    firstAnchor = helixFwd;
                }
                // By here we have the correct offset for the first anchorpoint.

                // The Body
                let currFwd = firstAnchor;
                let currBwd = firstAnchorPair;
                if (currFwd) mark(currFwd, helixId, offset, walkLabel);
                if (currBwd) mark(currBwd, helixId, offset, pairLabel);

                // Track visited anchors to break out of circular helices.
                const visitedAnchors = new Set<number>();
                if (firstAnchor) visitedAnchors.add(firstAnchor.id);

                while (currFwd) {
                    const nextStep = findNextPaired({ fwd: currFwd, bwd: currBwd }, { fwd: helixFwdDir, bwd: helixBwdDir }, helixSet);
                    if (!nextStep) break;
                    const nextAnchor = nextStep.anchor;
                    const nextPair = getPair(helixSet, nextAnchor);
                    // probably doesn't need this check but can happen due to cross/double pairing?
                    if (nextAnchor.id === currFwd.id) break;
                    // Circular helix guard: stop if we've already processed this anchor.
                    // Only mark circular when the source is 'fwd' — that means the forward
                    // strand's n3 chain physically looped back. A 'bwd_pair' revisit just
                    // means the bwd cursor walked past the fwd strand's end, which is normal.
                    if (visitedAnchors.has(nextAnchor.id)) {
                        break;
                    }
                    visitedAnchors.add(nextAnchor.id);

                    // Gap Detect
                    const fwdTrace = tracePathWithStop(currFwd, helixFwdDir, helixSet, nextAnchor);
                    // exclude the end points (the pairs themselves)
                    let fwdTail = fwdTrace.path.slice(1);
                    if (fwdTrace.reachedStop) fwdTail.pop();
                    let fwdHead: Nucleotide[] = [];
                    if (!fwdTrace.reachedStop || fwdTail.length === 0) {
                        const potentialHead = traceBack(nextAnchor, revFwdDir, helixSet, currFwd.id, 10);
                        if (potentialHead.length > fwdTail.length) { fwdHead = potentialHead; fwdTail = []; }
                    }

                    // this backward trace is actually very significant. 
                    // without it, cross-pairing becomes a real issue.
                    let bwdTail: Nucleotide[] = [];
                    let bwdHead: Nucleotide[] = [];
                    if (currBwd && nextPair) {
                        const bwdTrace = tracePathWithStop(currBwd, helixBwdDir, helixSet, nextPair);
                        bwdTail = bwdTrace.path.slice(1);
                        if (bwdTrace.reachedStop) bwdTail.pop();
                        if (!bwdTrace.reachedStop || bwdTail.length === 0) {
                            const potentialHead = traceBack(nextPair, revBwdDir, helixSet, currBwd.id, 10);
                            if (potentialHead.length > bwdTail.length) { bwdHead = potentialHead; bwdTail = []; }
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
                    if (nextPair) mark(nextPair, helixId, offset, pairLabel);

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
                        if (val.offset > maxFoundOffset) maxFoundOffset = val.offset;
                    }
                }

                // If the grid has content, start after it. If empty, start at 0.
                if (maxFoundOffset > -1) {
                    currentBinderOffset = maxFoundOffset + 4; // Add visual buffer
                }

                // 3. Group and Grid
                const unvisitedSet = new Set(unvisited.map(n => n.id));
                const processedExtras = new Set<number>();

                for (const node of unvisited) {
                    if (processedExtras.has(node.id)) continue;

                    // Trace Back to find segment start
                    let startOfSegment = node;
                    // Walk 5' until we hit something that is NOT in the unvisited set (either in grid or null)
                    let bwdSearch = node.n5 as Nucleotide | null;
                    while (bwdSearch && unvisitedSet.has(bwdSearch.id)) {
                        startOfSegment = bwdSearch;
                        bwdSearch = bwdSearch.n5 as Nucleotide | null;
                    }

                    // Trace Forward (3') to grab the full segment
                    let currSeg: Nucleotide | null = startOfSegment;

                    while (currSeg && unvisitedSet.has(currSeg.id)) {
                        if (processedExtras.has(currSeg.id)) break;

                        processedExtras.add(currSeg.id);
                        mark(currSeg, helixId, currentBinderOffset++, 'forward');

                        currSeg = currSeg.n3 as Nucleotide | null;
                    }

                    // Add small gap between distinct binder segments
                    currentBinderOffset += 4;
                }
            }
        });

        // --- MERGED-HELIX GROUP ALIGNMENT POST-PASS ---
        // Each merged helix's origin-groups were preserved with their
        // pre-merge (separate-helix) offsets. Temporarily split every merged
        // helix into pseudo-helices — one per origin-group — and run
        // alignGridPrim over the whole grid, so each group is anchored by its
        // own crossovers to the full (already-placed) lattice. The groups were
        // mutually disjoint on the offset axis at merge time (disjointness is
        // a merge precondition), so once aligned they can be folded back under
        // the real helixId without overlap. Afterwards, crossovers into the
        // merged helix observe ~0 shift, so the pipeline's own alignGridPrim
        // pass leaves the merged helix in place instead of shifting it as a
        // unit and misaligning one of the groups.
        if (preserveGrid && mergedGroups && mergedGroups.size > 0) {
            let maxHelixId = -1;
            for (const [, m] of grid) {
                if (m.helixId > maxHelixId) maxHelixId = m.helixId;
            }

            // pseudo helixId -> real (merged) helixId
            const pseudoToReal = new Map<number, number>();
            let nextPseudo = maxHelixId + 1;
            mergedGroups.forEach((groups, helixId) => {
                if (groups.length < 2) return;
                // Group 0 keeps the real helixId; the rest become pseudo-helices.
                for (let gi = 1; gi < groups.length; gi++) {
                    const pseudo = nextPseudo++;
                    pseudoToReal.set(pseudo, helixId);
                    for (const ntId of groups[gi]) {
                        const m = grid.get(ntId);
                        if (m) m.helixId = pseudo;
                    }
                }
            });

            if (pseudoToReal.size > 0) {
                alignGridPrim(grid, binderHelices);
                for (const [, m] of grid) {
                    const real = pseudoToReal.get(m.helixId);
                    if (real !== undefined) m.helixId = real;
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
                const shiftGroup = (helixId: number, gi: number, delta: number) => {
                    if (delta === 0) return;
                    const group = mergedGroups.get(helixId)?.[gi];
                    if (!group) return;
                    for (const ntId of group) {
                        const m = grid.get(ntId);
                        if (m && m.helixId === helixId) m.offset += delta;
                    }
                };

                // For (helixId, gi), find the min-|delta| shift such that
                // no mark of gi collides (same direction + same offset)
                // with any mark of ANY OTHER group on the same helix.
                // Returns 0 if already disjoint, or null if no delta
                // within the search bound resolves the conflict.
                const findCleanShift = (helixId: number, gi: number): number | null => {
                    const groups = mergedGroups.get(helixId);
                    if (!groups) return 0;
                    const myGroup = groups[gi];
                    if (!myGroup || myGroup.length === 0) return 0;
                    const mySet = new Set<number>(myGroup);

                    // My occupied cells, indexed by direction.
                    const myCells = new Map<string, Set<number>>();
                    let myMin = Infinity, myMax = -Infinity;
                    for (const ntId of myGroup) {
                        const m = grid.get(ntId);
                        if (!m || m.helixId !== helixId) continue;
                        let s = myCells.get(m.direction);
                        if (!s) { s = new Set(); myCells.set(m.direction, s); }
                        s.add(m.offset);
                        if (m.offset < myMin) myMin = m.offset;
                        if (m.offset > myMax) myMax = m.offset;
                    }
                    if (myCells.size === 0) return 0;

                    // Cells occupied by OTHER groups on the same helix.
                    const otherCells = new Map<string, Set<number>>();
                    let otherMin = Infinity, otherMax = -Infinity;
                    for (const [ntId, m] of grid) {
                        if (m.helixId !== helixId) continue;
                        if (mySet.has(ntId)) continue;
                        let s = otherCells.get(m.direction);
                        if (!s) { s = new Set(); otherCells.set(m.direction, s); }
                        s.add(m.offset);
                        if (m.offset < otherMin) otherMin = m.offset;
                        if (m.offset > otherMax) otherMax = m.offset;
                    }
                    if (otherCells.size === 0) return 0;

                    const conflictsAt = (delta: number): boolean => {
                        for (const [dir, myOffs] of myCells) {
                            const others = otherCells.get(dir);
                            if (!others) continue;
                            for (const off of myOffs) {
                                if (others.has(off + delta)) return true;
                            }
                        }
                        return false;
                    };

                    if (!conflictsAt(0)) return 0;

                    // Bound: any collision requires (my_off + delta) to
                    // land on an other-off, so |delta| ≤ (other-span +
                    // my-span). Add a small pad so we can step JUST past
                    // the far edge in either direction.
                    const mySpan = (myMax - myMin) || 0;
                    const otherSpan = (otherMax - otherMin) || 0;
                    const maxSearch = Math.max(mySpan + otherSpan + 8, 32);

                    for (let mag = 1; mag <= maxSearch; mag++) {
                        if (!conflictsAt(mag)) return mag;
                        if (!conflictsAt(-mag)) return -mag;
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
                            if (delta === null) { unresolvable++; continue; }
                            if (delta === 0) continue;
                            shiftGroup(helixId, gi, delta);
                            didMove = true;
                            movedAny = true;
                        }
                    });
                    if (!didMove) break;
                }

                if (movedAny) {
                    console.warn(
                        `[setGrid] Merged-helix inter-group overlap resolved by min-|delta| shifts` +
                        ` (sweeps=${sweepGuard - 1}${unresolvable > 0 ? `, unresolved=${unresolvable}` : ''}).`
                    );
                }
                if (unresolvable > 0) {
                    console.error(
                        `[setGrid] ${unresolvable} merged-helix group(s) had no collision-free placement` +
                        ` within the search bound. Grid will fail validation downstream.`
                    );
                }
            }
        }

        return { grid, binderHelices };
    };

    // grid flipper. Great helper function for the final output.
    export function gridFlip(grid: GridMap, helixId: number) {
        const ranges = new Map<number, { min: number; max: number }>();

        for (const [, mark] of grid.entries()) {
            if (mark.helixId !== helixId) continue;
            const range = ranges.get(mark.helixId);
            if (!range) {
                ranges.set(mark.helixId, { min: mark.offset, max: mark.offset });
                continue;
            }
            if (mark.offset < range.min) range.min = mark.offset;
            if (mark.offset > range.max) range.max = mark.offset;
        }

        for (const [, mark] of grid.entries()) {
            if (mark.helixId !== helixId) continue;
            const range = ranges.get(mark.helixId);
            if (!range) continue;
            mark.offset = range.max + range.min - mark.offset;
            mark.direction = mark.direction === 'forward' ? 'backward' : 'forward';
        }
    }

    function getScaffoldStrand() {
        let maxLen = 0;
        let scaffold: Strand | null = null;
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


    export function collectCrossovers(grid: GridMap) {
        const allNtIds = new Set<number>();
        for (const [ntId] of grid.entries()) allNtIds.add(ntId);

        const visited = new Set<number>();

        // crossovers[fromHelix][toHelix] = { sameWalk: n, diffWalk: n }
        //   sameWalk  = both runs have same offset trend → need flip
        //   diffWalk  = runs have opposite offset trend → already correct
        const crossovers = new Map<number, Map<number, { sameWalk: number; diffWalk: number }>>();
        const helixIds = new Set<number>();

        const ensureEntry = (from: number, to: number) => {
            if (!crossovers.has(from)) crossovers.set(from, new Map());
            const inner = crossovers.get(from)!;
            if (!inner.has(to)) inner.set(to, { sameWalk: 0, diffWalk: 0 });
            return inner.get(to)!;
        };

        for (const [ntId] of grid.entries()) {
            if (visited.has(ntId)) continue;

            const startNt = elements.get(ntId) as Nucleotide | undefined;
            if (!startNt || !(startNt instanceof Nucleotide)) continue;

            // Find 5' end
            let fivePrime: Nucleotide = startNt;
            const walkBack = new Set<number>();
            walkBack.add(fivePrime.id);
            while (true) {
                const prev = fivePrime.n5;
                if (!prev || !(prev instanceof Nucleotide)) break;
                if (!allNtIds.has(prev.id)) break;
                if (walkBack.has(prev.id)) break;
                walkBack.add(prev.id);
                fivePrime = prev;
            }

            // Walk 5'→3' and split into runs by helixId
            type Run = { helixId: number; offsets: number[] };
            const runs: Run[] = [];
            let currentRun: Run | null = null;

            let curr: Nucleotide | null = fivePrime;
            const walkForward = new Set<number>();

            while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                if (walkForward.has(curr.id)) break;
                walkForward.add(curr.id);
                visited.add(curr.id);

                const mark = grid.get(curr.id);
                if (mark) {
                    helixIds.add(mark.helixId);

                    if (currentRun && currentRun.helixId === mark.helixId) {
                        currentRun.offsets.push(mark.offset);
                    } else {
                        currentRun = { helixId: mark.helixId, offsets: [mark.offset] };
                        runs.push(currentRun);
                    }
                } else {
                    currentRun = null;
                }

                const n3ref: any = curr.n3;
                curr = (n3ref && n3ref instanceof Nucleotide) ? (n3ref as Nucleotide) : null;
            }

            // Now examine consecutive runs for crossovers
            for (let i = 0; i < runs.length - 1; i++) {
                const runA = runs[i];
                const runB = runs[i + 1];
                if (runA.helixId === runB.helixId) continue;

                // Determine offset trend for each run.
                // For runs with ≥2 nts, compare first and last offset.
                // For single-nt runs, skip (can't determine trend).
                if (runA.offsets.length < 2 && runB.offsets.length < 2) continue;

                // Use the trend near the crossover point:
                // runA trend: compare second-to-last offset to last offset
                // runB trend: compare first offset to second offset
                let trendA: 'inc' | 'dec' | null = null;
                let trendB: 'inc' | 'dec' | null = null;

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
                if (!trendA && !trendB) continue;

                // If only one side is known, we still can't compare — skip
                if (!trendA || !trendB) continue;

                const entry = ensureEntry(runA.helixId, runB.helixId);
                const entryRev = ensureEntry(runB.helixId, runA.helixId);

                if (trendA === trendB) {
                    // Same walk direction on both helices → need to flip one
                    entry.sameWalk++;
                    entryRev.sameWalk++;
                } else {
                    // Opposite walk direction → already correct
                    entry.diffWalk++;
                    entryRev.diffWalk++;
                }
            }
        }

        return { crossovers, helixIds };
    };
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

    export function directionAlign2(grid: GridMap) {

        // ── Initial scan to discover all helices and crossover stats ────
        const { crossovers, helixIds } = collectCrossovers(grid);

        const anchored = new Set<number>();
        const flippedHelices: number[] = [];
        let edgeCount = 0;

        for (const [fromHelix, neighbors] of crossovers.entries()) {
            for (const [toHelix, stats] of neighbors.entries()) {
                if (fromHelix >= toHelix) continue;
                if (stats.sameWalk + stats.diffWalk <= 0) continue;
                edgeCount++;
            }
        }

        const applyFlipToCrossoverStats = (helixId: number) => {
            const neighbors = crossovers.get(helixId);
            if (!neighbors) return;

            for (const [neighborHelix, stats] of neighbors.entries()) {
                const reverseStats = crossovers.get(neighborHelix)?.get(helixId);
                if (!reverseStats) continue;

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
        const queue: number[] = [0];
        let qIdx = 0;

        while (qIdx < queue.length) {
            const currHelix = queue[qIdx++];

            const neighbors = crossovers.get(currHelix);
            if (!neighbors) continue;

            for (const [neighborHelix] of neighbors.entries()) {
                if (anchored.has(neighborHelix)) continue;

                // Collect votes from ALL anchored helices to this neighbor
                let totalSameWalk = 0;
                let totalDiffWalk = 0;

                for (const anchoredHelix of anchored) {
                    const anchoredNeighbors = crossovers.get(anchoredHelix);
                    if (!anchoredNeighbors) continue;
                    const s = anchoredNeighbors.get(neighborHelix);
                    if (!s) continue;
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
            if (anchored.has(hId)) continue;

            anchored.add(hId);
            const subQueue: number[] = [hId];
            let subIdx = 0;

            while (subIdx < subQueue.length) {
                const currHelix = subQueue[subIdx++];
                const neighbors = crossovers.get(currHelix);
                if (!neighbors) continue;

                for (const [neighborHelix] of neighbors.entries()) {
                    if (anchored.has(neighborHelix)) continue;

                    let totalSameWalk = 0;
                    let totalDiffWalk = 0;

                    for (const anchoredHelix of anchored) {
                        const anchoredNeighbors = crossovers.get(anchoredHelix);
                        if (!anchoredNeighbors) continue;
                        const s = anchoredNeighbors.get(neighborHelix);
                        if (!s) continue;
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

    // Helper function to collect all backbone crossovers with their helix and offset info.
    // Slightly lengthy but quite useful.
    export function crossoverNts(
        grid: GridMap
    ): Array<{
        fromHelix: number;
        toHelix: number;
        fromOffset: number;
        toOffset: number;
        fromNt: Nucleotide;
        toNt: Nucleotide;
    }> {
        const allNtIds = new Set<number>();
        for (const [ntId] of grid.entries()) allNtIds.add(ntId);

        const visited = new Set<number>();
        const crossovers: Array<{
            fromHelix: number;
            toHelix: number;
            fromOffset: number;
            toOffset: number;
            fromNt: Nucleotide;
            toNt: Nucleotide;
        }> = [];

        for (const [ntId] of grid.entries()) {
            if (visited.has(ntId)) continue;

            const startNt = elements.get(ntId) as Nucleotide | undefined;
            if (!startNt || !(startNt instanceof Nucleotide)) continue;

            // Find 5' end
            let fivePrime: Nucleotide = startNt;
            const walkBack = new Set<number>();
            walkBack.add(fivePrime.id);

            while (true) {
                const prev = fivePrime.n5;
                if (!prev || !(prev instanceof Nucleotide)) break;
                if (!allNtIds.has(prev.id)) break;
                if (walkBack.has(prev.id)) break;
                walkBack.add(prev.id);
                fivePrime = prev;
            }

            // Walk 5' -> 3' and record backbone helix transitions
            let curr: Nucleotide | null = fivePrime;
            const walkForward = new Set<number>();
            let prevNt: Nucleotide | null = null;
            let prevMark: GridMark | null = null;

            while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                if (walkForward.has(curr.id)) break;
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
                } else {
                    prevNt = null;
                    prevMark = null;
                }

                const n3ref: any = curr.n3;
                curr = (n3ref && n3ref instanceof Nucleotide) ? (n3ref as Nucleotide) : null;
            }
        }

        return crossovers;
    }
    
    /**
     * The intent is for the export to track grid edits (combine, move,
     * flip) faithfully — every coordinate read goes through grid.get(),
     * so post-construction edits propagate to the export for free.
     * Topology pointers (n3) are used only to decide "where does this
     * strand go next"; everything emitted is read off the grid.
     *
     *  Domain split rules (close current domain, open a new one) on each
     *  step from currNt to currNt.n3 = nextNt:
     *    - nextNt has no grid mark             → close, skip until placed
     *    - nextNt is on a different helixId    → close + open (crossover)
     *    - nextNt has a different direction    → close + open (reversal)
     *    - nextNt offset != prevOffset + step  → close + open (gap)
     *  Otherwise the open domain extends by one slot.
     *
     *  Same gap example (offsets 1..5, _, _, 8..10 on one helix forward,
     *  then crossover) emits three domains: [1,6), [8,11), and one on
     *  the destination helix.
     */
    export function buildScadnano3(
        grid: GridMap,
        helices: Nucleotide[][],
        gridType?: string,
        helixPositions?: Map<number, [number, number]>
    ) {
        // ── Scaffold detection ──────────────────────────────────────────
        const scaffoldStrand: Strand | null = getScaffoldStrand();
        const SCAFFOLD_COLOR = '#0066cc';
        const STAPLE_COLORS = ['#f74308', '#57bb00', '#000000'];

        // ── Helix metadata ──────────────────────────────────────────────
        const helixCount = helices.length || Math.max(0, ...Array.from(grid.values()).map(m => m.helixId + 1));
        const helixMaxOffsets = new Map<number, number>();
        for (const [, mark] of grid.entries()) {
            const current = helixMaxOffsets.get(mark.helixId) ?? -1;
            if (mark.offset > current) helixMaxOffsets.set(mark.helixId, mark.offset);
        }
        const scadHelices = Array.from({ length: helixCount }, (_, i) => ({
            max_offset: (helixMaxOffsets.get(i) ?? 0) + 1,
            grid_position: helixPositions?.get(i) ?? [0, i]
        }));

        // ── Strand iteration ────────────────────────────────────────────
        type Domain = { helix: number; forward: boolean; start: number; end: number };
        type ScadStrand = {
            color: string;
            sequence: string;
            domains: Domain[];
            is_scaffold?: boolean;
            circular?: boolean;
        };
        type OpenDomain = {
            helixId: number;
            direction: 'forward' | 'backward';
            step: 1 | -1;
            forward: boolean;     // scadnano "forward" flag
            minOffset: number;
            maxOffset: number;
            lastOffset: number;
        };

        const scadStrands: ScadStrand[] = [];

        const allSystems: System[] = [];
        if (Array.isArray(systems)) for (const s of systems) if (s) allSystems.push(s);
        if (typeof tmpSystems !== 'undefined' && Array.isArray(tmpSystems)) {
            for (const s of tmpSystems) if (s) allSystems.push(s);
        }

        for (const sys of allSystems) {
            const sysStrands = (sys && Array.isArray((sys as any).strands)) ? (sys as any).strands as Strand[] : [];

            for (const strand of sysStrands) {
                if (!strand) continue;

                const start: any = (strand as any).end5;
                if (!(start instanceof Nucleotide)) continue;

                // ── Walk this strand 5' → 3' via n3 ─────────────────────
                let sequence = '';
                const domains: Domain[] = [];
                let isCircular = false;
                let openDomain: OpenDomain | null = null;
                const visited = new Set<number>();

                const closeDomain = () => {
                    if (!openDomain) return;
                    domains.push({
                        helix: openDomain.helixId,
                        forward: openDomain.forward,
                        start: openDomain.minOffset,
                        end: openDomain.maxOffset + 1
                    });
                    openDomain = null;
                };

                const openAt = (nt: Nucleotide, mark: GridMark) => {
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

                let curr: Nucleotide | null = start;

                while (curr instanceof Nucleotide) {
                    if (visited.has(curr.id)) {
                        // Walked back to a node we already emitted — circular.
                        isCircular = true;
                        break;
                    }
                    visited.add(curr.id);

                    const mark = grid.get(curr.id);

                    if (!mark) {
                        // Unplaced nt — close any open domain, skip until we
                        // land on a placed nt again.
                        closeDomain();
                    } else if (!openDomain) {
                        openAt(curr, mark);
                    } else {
                        const continues =
                            openDomain.helixId === mark.helixId &&
                            openDomain.direction === mark.direction &&
                            openDomain.lastOffset + openDomain.step === mark.offset;

                        if (continues) {
                            openDomain.lastOffset = mark.offset;
                            if (mark.offset < openDomain.minOffset) openDomain.minOffset = mark.offset;
                            if (mark.offset > openDomain.maxOffset) openDomain.maxOffset = mark.offset;
                            sequence += curr.type || 'N';
                        } else {
                            closeDomain();
                            openAt(curr, mark);
                        }
                    }

                    // Advance via n3. Closed-loop strands have n3 of the 3' end
                    // pointing back to end5, so detect that before stepping.
                    const nextRef: any = (curr as any).n3;
                    if (nextRef instanceof Nucleotide && nextRef.id === start.id) {
                        isCircular = true;
                        break;
                    }
                    curr = (nextRef instanceof Nucleotide) ? nextRef : null;
                }

                closeDomain();

                if (domains.length === 0) continue;

                const isScaffold = scaffoldStrand !== null && strand === scaffoldStrand;
                const color = isScaffold
                    ? SCAFFOLD_COLOR
                    : STAPLE_COLORS[Math.floor(Math.random() * STAPLE_COLORS.length)];

                const out: ScadStrand = { color, sequence, domains };
                if (isScaffold) out.is_scaffold = true;
                if (isCircular) out.circular = true;
                scadStrands.push(out);
            }
        }

        return {
            version: '0.20.1',
            grid: gridType,
            helices: scadHelices,
            strands: scadStrands
        };
    };

    // Confirms whether every offset -> direction is unique. 
    export function validateGrid(grid: Map<number, GridMark>) {
        // Structure: Map<HelixID, { forward: Map<Offset, NtID>, backward: Map<Offset, NtID> }>
        const checkMap = new Map<number, {
            forward: Map<number, number>;
            backward: Map<number, number>;
        }>();

        let conflicts = 0;

        for (const [ntId, pos] of grid.entries()) {
            // 1. Initialize Helix Bucket if missing
            if (!checkMap.has(pos.helixId)) {
                checkMap.set(pos.helixId, {
                    forward: new Map(),
                    backward: new Map()
                });
            }

            const helixBuckets = checkMap.get(pos.helixId)!;
            const strandMap = helixBuckets[pos.direction];

            // 2. Check for collision
            if (strandMap.has(pos.offset)) {
                const existingNt = strandMap.get(pos.offset);
                console.error(
                    `❌ CONFLICT DETECTED:\n` +
                    `   Helix: ${pos.helixId}\n` +
                    `   Strand: ${pos.direction}\n` +
                    `   Offset: ${pos.offset}\n` +
                    `   Fighting Nucleotides: IDs ${existingNt} vs ${ntId}`
                );
                conflicts++;
            } else {
                // 3. Register valid position
                strandMap.set(pos.offset, ntId);
            }
        }

        if (conflicts === 0) {
            console.log(`✅ Grid Validated: ${grid.size} nucleotides assigned with 0 overlapping offsets.`);
        } else {
            console.warn(`⚠️ Grid Validation Failed: Found ${conflicts} offset collisions.`);
        }

        if (grid.size !== elements.size) {
            console.log("⚠️ INVALID. Grid does not include all nucleotides from the original element set.");
        }
    }
}
