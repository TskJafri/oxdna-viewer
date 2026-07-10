"use strict";
/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
/*
// for ease of use, and to prevent dumb mistakes as I code and test things out, here is exactly the commands to use this in console:
findBasepairs(); // use 3 because there are 2 versions of findBasepairs, and 3 is the fastest. Dont ask why i named it that.
helix.dropIntraStrandPairs();
let {partials, unpaired} = helix.findHelixPartials2(elements, 2);
let {ssdna, stubs, longssScaffold} = helix.ssdnaPartials(unpaired);
let ssScaffold = helix.longssScaffoldfunc(longssScaffold, stubs);
let {helices, lastScraps, binders, binder2, disconnected, unhandled} = helix.generateHelix(partials, ssdna, ssScaffold, stubs);
// and helices are what you want!
// This code has been completed (polishing required but sure).
// After running this, check for helix.flat().length == elements.size
// If false, then something went wrong! RIP.
*/
/*
// Check for double-pairing or cross-pairing.
let pairTally = new Map();
let overloadedNucleotides = [];
elements.forEach(nt=>{
    let targetPairId = nt.pair;
    if (targetPairId !== undefined && targetPairId !== null) {
        let currentCount = pairTally.get(targetPairId) || 0;
        pairTally.set(targetPairId, currentCount + 1);
        if (currentCount + 1 === 2) {
            overloadedNucleotides.push(targetPairId);
        }
    }
})
*/
// For even easier use, just run:
// let {helices, partials, usedSides} = await helix.findHelices(elements, 2);
var helix;
(function (helix_1) {
    // helper function cuz didnt want to type this every time
    function checkAngle(n1 = null, n2 = null) {
        if (!n1 || !n2)
            return 0;
        return Math.acos(n1.getA3().dot(n2.getA3())) * (180 / Math.PI);
    }
    helix_1.checkAngle = checkAngle;
    // Returns the longest strand in the system as scaffold...
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
    helix_1.getScaffoldStrand = getScaffoldStrand;
    // Removes intra-strand pairings across the whole structure
    // ALWAYS run this function after findBasepairs3(). Very important.
    function dropIntraStrandPairs() {
        elements.forEach((nt) => {
            if (!(nt instanceof Nucleotide))
                return;
            if (!nt.pair)
                return;
            if (nt.pair.strand === nt.strand) {
                const mate = nt.pair;
                nt.pair = null;
                if (mate.pair === nt) {
                    mate.pair = null;
                }
            }
        });
    }
    helix_1.dropIntraStrandPairs = dropIntraStrandPairs;
    // Finds helix parts using destructive consumption of a working copy of elements (called elmts).
    // tolerance 2 is good enough for most cases. Higher tolerances seem to have no negative consequences, however.
    // go to terminatingConditions for what tolerance is.
    function findHelixPartials2(inputMap, tolerance = 2) {
        const elmts = new Map(inputMap); // copy of the elements map
        const elmts2 = new Map(inputMap); // backup copy for duplication, used later
        const unpaired = new Map(); // unpaired / skipped nts
        let partials = [];
        const record = (list, set, nt = null) => {
            if (nt && !set.has(nt.id)) {
                set.add(nt.id);
                list.push(nt);
            }
        };
        const nextStart = () => elmts.values().next().value;
        while (true) {
            const start = nextStart();
            if (!start)
                break;
            // Collect unpaired/binder nts for downstream ssDNA processing instead of discarding silently.
            // Additionally, only walk if the pair is also present in the current pool (elmts).
            const pairInPool = start.pair ? elmts.get(start.pair.id) : undefined;
            if (!start.pair || !pairInPool) {
                unpaired.set(start.id, start);
                elmts.delete(start.id);
                continue;
            }
            // initialize the walker.
            let curr = start;
            let currPair = pairInPool;
            const strandA = curr.strand; // required to ensure we don't cross strands. This is how we know partials is actually a partial helix.
            const strandB = currPair.strand;
            // this is the partial being built.
            const partial = [];
            const seen = new Set();
            const terminatingConditions = (nucA, nucB, dir) => {
                const nucAc = elmts.get(nucA.id + dir);
                const nucBc = elmts.get(nucB.id - dir);
                if (!nucAc || !nucBc)
                    return null;
                if (nucAc.strand !== strandA || nucBc.strand !== strandB)
                    return null;
                // Tolerance window: OR logic. Continue if ANY offset in [1, tolerance] has forward paired to backward.
                // (allows rescue by walking farther along topological neighbors on both strands).
                let forwardCursor = nucA;
                let backwardCursor = nucB;
                for (let offset = 1; offset <= tolerance; offset++) {
                    const forward = elmts.get(nucA.id + dir * offset);
                    const backward = elmts.get(nucB.id - dir * offset);
                    if (!forward || !backward)
                        continue; // check if either of them even exist.
                    if (forward.strand !== strandA || backward.strand !== strandB)
                        continue;
                    if (forward.pair === backward) {
                        return { nextA: nucAc, nextB: nucBc };
                    }
                }
                return null;
            };
            // actual traversal loop.
            while (curr && currPair) {
                record(partial, seen, curr);
                record(partial, seen, currPair);
                // destructive consumption. This is why we make a copy of elements, and not use the original map directly.
                elmts.delete(curr.id);
                elmts.delete(currPair.id);
                const step = terminatingConditions(curr, currPair, 1) || terminatingConditions(curr, currPair, -1);
                if (!step)
                    break;
                console.log("Circular strand detected");
                // Circular helix guard: if the next step loops back into this partial, stop here.
                if (seen.has(step.nextA.id) || seen.has(step.nextB.id))
                    break;
                console.log("Circular strand stopped");
                // Always consume the immediate neighbors (curr+1 and a-1) even if mismatched.
                record(partial, seen, step.nextA);
                record(partial, seen, step.nextB);
                elmts.delete(step.nextA.id);
                elmts.delete(step.nextB.id);
                // Advance walker by one along each strand.
                curr = step.nextA;
                currPair = step.nextB;
            }
            if (partial.length) {
                partials.push(partial);
            }
        }
        // Deduplicate across all partials: any nucleotide that appears more than once
        // is moved to unpaired along with its pair and any nucleotide paired to it.
        // they will be handled as either ssDNA or just directly added to helix later.
        const seenfordups = new Map();
        const duplicates = new Set();
        partials.forEach(helix => {
            helix.forEach(nt => {
                if (seenfordups.has(nt.id)) {
                    duplicates.add(nt.id);
                }
                else {
                    seenfordups.set(nt.id, nt);
                }
            });
        });
        // had to add the multi-pairing logic (down below) to fix the problem. But still logging for curiosity.
        console.log('Duplicates set: ', duplicates);
        if (duplicates.size) {
            // Collect pair ids for duplicates
            console.log('duplicates found: ', duplicates);
            console.log('total duplicates: ', duplicates.size);
            const pairIds = new Set();
            duplicates.forEach(id => {
                const a = elmts2.get(id);
                if (a) {
                    unpaired.set(a.id, a);
                    if (a.pair) {
                        pairIds.add(a.pair.id);
                        unpaired.set(a.pair.id, a.pair);
                    }
                }
            });
            // Remove duplicates, their pairs, and any nucleotide whose pair is a duplicate
            partials.forEach((helix, i) => {
                partials[i] = helix.filter(nt => {
                    if (duplicates.has(nt.id) || pairIds.has(nt.id)) {
                        unpaired.set(nt.id, nt);
                        return false;
                    }
                    if (nt.pair && duplicates.has(nt.pair.id)) {
                        unpaired.set(nt.id, nt);
                        return false;
                    }
                    return true;
                });
            });
            // Drop empty partial lists
            partials = partials.filter(helix => helix.length > 0);
        }
        return { partials, unpaired: Array.from(unpaired.values()) };
    }
    helix_1.findHelixPartials2 = findHelixPartials2;
    // Groups unpaired/binder nucleotides (unpaired) into ssDNA partials by strand.
    // Only contiguous runs (>2) along a strand are kept; shorter runs go to stubs.
    function ssdnaPartials(unpaired) {
        const unpairStrand = new Map();
        const ssdna = [];
        const stubs = [];
        const longssScaffold = [];
        const scaffold = getScaffoldStrand();
        unpaired.forEach(nt => {
            const arr = unpairStrand.get(nt.strand) || [];
            arr.push(nt);
            unpairStrand.set(nt.strand, arr);
        });
        // For each strand, sort it, find it's 5' ends and walk down n3 to build runs.
        unpairStrand.forEach((list, strand) => {
            const inSet = new Set(list.map(n => n.id));
            const visited = new Set();
            const isScaffoldStrand = scaffold && strand === scaffold;
            list.sort((a, b) => a.id - b.id);
            for (const nt of list) {
                if (visited.has(nt.id))
                    continue;
                const isStart = !nt.n5 || !inSet.has(nt.n5.id);
                if (isStart) {
                    const run = [];
                    let curr = nt;
                    while (curr && inSet.has(curr.id)) {
                        run.push(curr);
                        visited.add(curr.id);
                        curr = curr.n3;
                    }
                    if (run.length > 2) {
                        if (!isScaffoldStrand) {
                            ssdna.push(run);
                        }
                        else {
                            longssScaffold.push(...run);
                        }
                    }
                    else {
                        run.forEach(r => stubs.push(r));
                    }
                }
                // // expand run both directions along n5/n3 within unpaired set
                // const run: Nucleotide[] = [];
                // const pushRun = (node: Nucleotide | null, dir: 'n5' | 'n3') => {
                // 	let curr = node;
                // 	while (curr && inSet.has(curr.id) && !visited.has(curr.id)) {
                // 		run.push(curr);
                // 		visited.add(curr.id);
                // 		curr = curr[dir] as Nucleotide | null;
                // 	}
                // };
                // // walk n5 then n3 from seed to capture contiguous block
                // pushRun(nt, 'n5');
                // // pushRun added seed and upstream; now extend downstream from the last added toward n3
                // // ensure we start from the n3 of the seed to avoid duplicate seed
                // const seedN3 = nt.n3 as Nucleotide | null;
                // pushRun(seedN3, 'n3');
                // if (run.length > 2) {
                // 	ssdna.push(run);
                // } else {
                // 	run.forEach(r => stubs.push(r));
                // }
            }
        });
        return { ssdna, stubs, longssScaffold };
    }
    helix_1.ssdnaPartials = ssdnaPartials;
    // helper function for adding the average a3 vector in canvas. Really should not be here.
    function averageA3a(list) {
        if (!list.length)
            return new THREE.Vector3(0, 0, 0);
        // Align all A3 vectors so they point in a consistent direction before averaging.
        const ref = list[0].getA3().clone().normalize();
        const acc = ref.clone();
        for (let i = 1; i < list.length; i++) {
            const v = list[i].getA3().clone().normalize();
            acc.add(v.dot(ref) < 0 ? v.multiplyScalar(-1) : v);
        }
        acc.divideScalar(list.length);
        const avg = acc.normalize();
        // Visualize the averaged orientation from the first nucleotide origin when possible.
        const origin = list[0]?.getPos();
        if (origin && typeof THREE !== 'undefined' && typeof scene !== 'undefined' && scene?.add) {
            const helper = new THREE.ArrowHelper(avg.clone(), origin, 5);
            scene.add(helper);
        }
        return avg;
    }
    helix_1.averageA3a = averageA3a;
    ;
    // helper to consolidate the nucleotides into contiguous segments, and enforce equal halves for scaffold segments.
    function longssScaffoldfunc(longssScaffold, stubs = []) {
        const ssScaffold = [];
        if (!longssScaffold.length)
            return ssScaffold;
        // Work on a sorted copy so numeric contiguity is easy to detect.
        const sorted = [...longssScaffold].sort((a, b) => a.id - b.id);
        let runAway = [];
        const flushRun = () => {
            if (!runAway.length)
                return;
            if (runAway.length < 3) {
                runAway.forEach(nt => stubs.push(nt));
                runAway = [];
                return;
            }
            // Enforce equal halves; trim one nucleotide if odd-length to satisfy the requirement.
            const evenLen = runAway.length - (runAway.length % 2);
            if (evenLen !== runAway.length) {
                const dropped = runAway[evenLen];
                if (dropped)
                    stubs.push(dropped);
            }
            if (evenLen === 0) {
                runAway = [];
                return;
            }
            const half = evenLen / 2;
            ssScaffold.push(runAway.slice(0, half));
            ssScaffold.push(runAway.slice(half, evenLen));
            runAway = [];
        };
        sorted.forEach(nt => {
            const last = runAway[runAway.length - 1];
            if (!last || nt.id === last.id + 1) {
                runAway.push(nt);
                return;
            }
            flushRun();
            runAway.push(nt);
        });
        flushRun();
        return ssScaffold;
    }
    helix_1.longssScaffoldfunc = longssScaffoldfunc;
    // let partialStrandMap = new Map<number, Map<number, Nucleotide[]>>();
    function mapPartialEnds(partials) {
        const partialEndsMap = new Map();
        partials.forEach((partial, index) => {
            const inSet = new Set(partial.map(n => n.id));
            // Find 3' ends (n3 is missing or outside the partial) and sort them by ID
            const ends3 = partial
                .filter(n => !n.n3 || !inSet.has(n.n3.id))
                .sort((a, b) => a.id - b.id);
            // Find 5' ends (n5 is missing or outside the partial)
            const ends5 = partial.filter(n => !n.n5 || !inSet.has(n.n5.id));
            if (ends3.length >= 2 && ends5.length >= 2) {
                // strand1 is assigned to the one with the lowest ID at the 3' end (ends3[0])
                const start1 = ends5.find(n => n.strand === ends3[0].strand);
                const start2 = ends5.find(n => n.strand === ends3[1].strand);
                // Skip 1-bp partials: when a strand has length 1 inside the partial,
                // start1 == end1 (or start2 == end2). There's no linear extent, so the
                // "two sides" model doesn't apply. These partials are merged by attach-count
                // cap in generateHelix instead of by side uniqueness.
                if (start1.id === ends3[0].id || start2.id === ends3[1].id)
                    return;
                partialEndsMap.set(index, {
                    start1: start1, end1: ends3[0],
                    start2: start2, end2: ends3[1]
                });
            }
        });
        return partialEndsMap;
    }
    helix_1.mapPartialEnds = mapPartialEnds;
    //  Returns the indices of partials that have exactly 1 free side.
    //  Only partials present in `partialEndsMap` (i.e. those with 2 sides) are considered.
    function partialsWithOneFreeSide(partials, partialEndsMap, usedSides) {
        const result = [];
        for (let pIdx = 0; pIdx < partials.length; pIdx++) {
            if (!partialEndsMap.has(pIdx))
                continue;
            const sideUsage = usedSides.get(pIdx);
            const side0Used = (sideUsage?.get(0) ?? 0) > 1e-9;
            const side1Used = (sideUsage?.get(1) ?? 0) > 1e-9;
            if (side0Used !== side1Used) {
                result.push(pIdx);
            }
        }
        return result;
    }
    helix_1.partialsWithOneFreeSide = partialsWithOneFreeSide;
    // Like partialsWithOneFreeSide, but also computes the helical axis for each such partial
    // and flips it so it points TOWARDS the free side. Returns pIdx -> { vector }.
    // planeVector from getPartialAxis points side 0 -> side 1 (guess = midSide1 - midSide0,
    // and planeVector is aligned to guess), so we negate it when the free side is 0.
    function partialAxesTowardFreeSide(partials, partialEndsMap, usedSides) {
        const result = new Map();
        const oneFreeSide = partialsWithOneFreeSide(partials, partialEndsMap, usedSides);
        for (const pIdx of oneFreeSide) {
            const ends = partialEndsMap.get(pIdx);
            if (!ends)
                continue;
            const sideUsage = usedSides.get(pIdx);
            const side0Free = (sideUsage?.get(0) ?? 0) < 1e-9;
            const { planeVector } = getPartialAxis(ends);
            if (side0Free)
                planeVector.negate();
            result.set(pIdx, { vector: planeVector });
        }
        return result;
    }
    helix_1.partialAxesTowardFreeSide = partialAxesTowardFreeSide;
    // Shortest distance between two 3D line segments (P1->P2) and (P3->P4).
    // Eberly's algorithm. Inputs in any consistent unit; output in same unit.
    // Used by hashAxisOverlap for cylinder-vs-cylinder overlap (distance <= 2*radius).
    function segmentDistance3D(P1, P2, P3, P4) {
        const d1x = P2.x - P1.x, d1y = P2.y - P1.y, d1z = P2.z - P1.z;
        const d2x = P4.x - P3.x, d2y = P4.y - P3.y, d2z = P4.z - P3.z;
        const rx = P1.x - P3.x, ry = P1.y - P3.y, rz = P1.z - P3.z;
        const a = d1x * d1x + d1y * d1y + d1z * d1z;
        const e = d2x * d2x + d2y * d2y + d2z * d2z;
        const f = d2x * rx + d2y * ry + d2z * rz;
        const eps = 1e-12;
        let s = 0, t = 0;
        if (a <= eps && e <= eps) {
            return Math.sqrt(rx * rx + ry * ry + rz * rz);
        }
        if (a <= eps) {
            s = 0;
            t = e > eps ? Math.max(0, Math.min(1, f / e)) : 0;
        }
        else {
            const c = d1x * rx + d1y * ry + d1z * rz;
            if (e <= eps) {
                t = 0;
                s = Math.max(0, Math.min(1, -c / a));
            }
            else {
                const b = d1x * d2x + d1y * d2y + d1z * d2z;
                const denom = a * e - b * b;
                if (Math.abs(denom) > eps) {
                    s = Math.max(0, Math.min(1, (b * f - c * e) / denom));
                }
                else {
                    s = 0; // parallel / near-parallel
                }
                t = (b * s + f) / e;
                if (t < 0) {
                    t = 0;
                    s = Math.max(0, Math.min(1, -c / a));
                }
                else if (t > 1) {
                    t = 1;
                    s = Math.max(0, Math.min(1, (b - c) / a));
                }
            }
        }
        const cax = P1.x + s * d1x - (P3.x + t * d2x);
        const cay = P1.y + s * d1y - (P3.y + t * d2y);
        const caz = P1.z + s * d1z - (P3.z + t * d2z);
        return Math.sqrt(cax * cax + cay * cay + caz * caz);
    }
    // Hash-merge using overlapping cylindrical volumes ("hash cells") around each partial's
    // free-side basepair. Each partial contributes one hash-cylinder:
    //   origin = midpoint of the free-side basepair (start1+end2)/2 or (end1+start2)/2
    //   dir    = the partial's helical axis vector pointing TOWARD the free side
    //            (already produced by the caller via partialAxesTowardFreeSide)
    //   length = cylinderLengthAng (default 50 Å)
    //   radius = cylinderRadiusAng (default 12.5 Å)
    // Brute-force pairwise. Two partials are candidate-mergeable iff BOTH:
    //   1. Their hash-cylinders overlap: shortest segment-segment distance between the two
    //      cylinder centerlines <= 2*radius (Eberly's algorithm in segmentDistance3D).
    //   2. axis_a . axis_b < dotThreshold (default -0.9, i.e. strongly anti-parallel).
    //      Rationale: both axes point TOWARD their free sides, so a valid end-to-end merge
    //      (free sides facing each other) makes the axes anti-parallel -> dot close to -1.
    // `dist` is computed but not used as a filter (the overlap test handles spatial rejection).
    //
    // Candidate pairs are greedily resolved in ASCENDING dot order (most negative first = best
    // anti-parallel alignment wins), with each pIdx exclusive (a partial with one free side
    // joins at most one merge). This weighting by `dot` is the tie-break rule.
    //
    // NOTE: `bbOffsets` from Nucleotide.getInstanceParameter3 returns positions in oxView
    // internal units, where 1 unit = 8.518 Å. All distance thresholds given in Å are converted
    // once for internal comparison in oxView units.
    function hashAxisOverlap(partials, partialEndsMap, usedSides, partialAxes, options) {
        const OX_TO_ANG = 8.518;
        const dotThreshold = options?.dotThreshold ?? -0.9;
        const cylRadiusAng = options?.cylinderRadiusAng ?? 12.5;
        const cylLengthAng = options?.cylinderLengthAng ?? 50;
        const cylLengthOx = cylLengthAng / OX_TO_ANG;
        const overlapOx = (2 * cylRadiusAng) / OX_TO_ANG;
        const entries = [];
        partialAxes.forEach(({ vector }, pIdx) => {
            const ends = partialEndsMap.get(pIdx);
            if (!ends)
                return;
            const sideUsage = usedSides.get(pIdx);
            const side0Free = (sideUsage?.get(0) ?? 0) < 1e-9;
            // Side 0 basepair: start1 + end2  ;  Side 1 basepair: end1 + start2
            const bpA = side0Free ? ends.start1 : ends.end1;
            const bpB = side0Free ? ends.end2 : ends.start2;
            const origin = bpA.getInstanceParameter3('bbOffsets')
                .add(bpB.getInstanceParameter3('bbOffsets'))
                .multiplyScalar(0.5);
            const dir = vector.clone().normalize();
            const segEnd = origin.clone().add(dir.clone().multiplyScalar(cylLengthOx));
            entries.push({ pIdx, origin, dir, segEnd });
        });
        // Brute-force: collect all candidate pairs that pass BOTH gates.
        const candidates = [];
        for (let i = 0; i < entries.length; i++) {
            const A = entries[i];
            for (let j = i + 1; j < entries.length; j++) {
                const B = entries[j];
                // Gate 1: cylinder overlap. Each hash-cylinder is a finite segment from its
                // origin to segEnd along the partial's axis. Two cylinders overlap iff the
                // shortest distance between their centerline segments <= 2*radius. This is the
                // essential spatial filter — without it, anti-parallel partials on opposite
                // sides of the scene would falsely match.
                const segDist = segmentDistance3D(A.origin, A.segEnd, B.origin, B.segEnd);
                if (segDist > overlapOx)
                    continue;
                // Gate 2: anti-parallel alignment. Both axes point toward their free sides,
                // so a valid merge (free sides facing each other) requires dot < dotThreshold.
                const dot = A.dir.dot(B.dir);
                if (dot >= dotThreshold)
                    continue;
                // Distance computed for reference only; not used as a filter.
                const distOx = A.origin.distanceTo(B.origin);
                candidates.push({ a: A.pIdx, b: B.pIdx, dot, dist: distOx * OX_TO_ANG });
            }
        }
        // Resolve greedily: lowest dot (most anti-parallel) wins, each pIdx in at most one pair.
        candidates.sort((x, y) => x.dot - y.dot);
        const paired = new Set();
        const result = [];
        for (const c of candidates) {
            if (paired.has(c.a) || paired.has(c.b))
                continue;
            result.push(c);
            paired.add(c.a);
            paired.add(c.b);
        }
        return result;
    }
    helix_1.hashAxisOverlap = hashAxisOverlap;
    // Consume a list of hashAxisOverlap merge pairs and fold any cross-helix pairs
    // into a single helix via combineHelices. For each pair:
    //   1. Look up each partial's current helix via grid.get(partials[x][0].id).helixId.
    //   2. If both partials already belong to the same current helix, skip.
    //   3. If they're in different helices, use the same path as the scadnano
    //      "Combine Helices" button to make the pair disjoint on the grid
    //      (computeCombineShifts + applyCombineShifts) and then combineHelices.
    // Pairs are processed in input order; each iteration re-reads helixIds from
    // the grid, so indices returned by combineHelices' idRemap are tracked
    // correctly for any subsequent pairs.
    function applyAxisOverlapMerge(helices, grid, partials, mergePairs) {
        if (!Array.isArray(helices) || !(grid instanceof Map) ||
            !Array.isArray(partials) || !Array.isArray(mergePairs))
            return;
        if (mergePairs.length === 0)
            return;
        mergePairs.forEach(({ a, b }) => {
            const pa = partials[a];
            const pb = partials[b];
            if (!Array.isArray(pa) || pa.length === 0)
                return;
            if (!Array.isArray(pb) || pb.length === 0)
                return;
            const markA = grid.get(pa[0].id);
            const markB = grid.get(pb[0].id);
            if (!markA || !markB)
                return;
            const helixA = markA.helixId;
            const helixB = markB.helixId;
            if (typeof helixA !== 'number' || typeof helixB !== 'number')
                return;
            if (helixA === helixB)
                return; // same current helix — nothing to do
            if (helixA >= helices.length || helixB >= helices.length)
                return;
            if (!Array.isArray(helices[helixA]) || helices[helixA].length === 0)
                return;
            if (!Array.isArray(helices[helixB]) || helices[helixB].length === 0)
                return;
            // Mirror the Combine Helices button: shift GridMark.offsets so the pair
            // is disjoint on the grid, then perform the logical merge.
            const shifts = toscad.computeCombineShifts(helices, [helixA, helixB], grid);
            if (shifts && shifts.size > 0)
                toscad.applyCombineShifts(grid, shifts);
            combineHelices(helices, [helixA, helixB], grid);
        });
    }
    helix_1.applyAxisOverlapMerge = applyAxisOverlapMerge;
    // Perfected!
    // this one uses average a3 vectors of CONNECTED strands, as opposed to average a3 vectors of the entire partial (which cancels out, due to topology).
    function generateHelix(partials, ssdna, ssScaffold, stubs) {
        // Currently uses partials, ssScaffold and stubs to build perfect (almost) helices.
        // Hence, helices.flat().length and ssdna.flat().length should be the full size of the structure. For any missing piece, check lastScraps[].
        const helices = [];
        // guys for context lastScraps[] basically are the dumb nucleotides that couldnt be placed into helices due to fraying and angle conflicts.
        // Stored as segments so grouped leftovers (e.g. deferred ssScaffold segments) stay together.
        const lastScraps = [];
        if (!partials.length)
            return { helices, lastScraps, usedSides: new Map() }; // surely no helices if no partials.
        const dot = 0.5;
        // quick lookup for id to partial index and stubs index.
        const idToPartial = new Map();
        partials.forEach((list, idx) => {
            list.forEach(nt => idToPartial.set(nt.id, idx));
        });
        const idTostubs = new Map();
        stubs.forEach((nt, idx) => idTostubs.set(nt.id, idx));
        console.log('ID to Partial Map:', idToPartial); // lets check out what the map looks like
        console.log('ID to stubs Map:', idTostubs);
        const averageA3 = (list) => {
            if (!list.length)
                return null;
            const acc = new THREE.Vector3(0, 0, 0);
            list.forEach(nt => {
                acc.add(nt.getA3().clone().normalize());
            });
            const len = acc.length();
            if (len < 1e-6)
                return null;
            return acc.divideScalar(len);
        };
        // within a partial, find the nts that belong to a specific strand within a specific partial. 
        // these will the ones used for finding the average a3 vector, which will later be used for connecting partials.
        const partialStrandMap = new Map();
        const getPartialStrandNts = (partialIdx, strand) => {
            let byStrand = partialStrandMap.get(partialIdx);
            if (!byStrand) {
                byStrand = new Map();
                partialStrandMap.set(partialIdx, byStrand);
            }
            let list = byStrand.get(strand.id);
            if (!list) {
                list = partials[partialIdx].filter(nt => nt.strand === strand);
                byStrand.set(strand.id, list);
            }
            return list;
        };
        // similar to above, but caches average a3 vectors for each partial-strand combo instead of just nucleotide lists.
        const partialStrandA3 = new Map();
        const getPartialStrandA3 = (partialIdx, strand) => {
            let byStrand = partialStrandA3.get(partialIdx);
            if (!byStrand) {
                byStrand = new Map();
                partialStrandA3.set(partialIdx, byStrand);
            }
            let vec = byStrand.get(strand.id);
            if (vec === undefined) {
                const list = getPartialStrandNts(partialIdx, strand);
                vec = list.length ? averageA3(list) : null;
                byStrand.set(strand.id, vec ?? null);
            }
            return vec;
        };
        // stubs are just single nts, so caching their a3 vectors is simpler.
        const stubsA3 = new Map();
        const getstubsA3 = (idx) => {
            let vec = stubsA3.get(idx);
            if (!vec) {
                vec = stubs[idx].getA3().clone().normalize();
                stubsA3.set(idx, vec);
            }
            return vec;
        };
        const totalNodes = partials.length + stubs.length;
        const parent = Array.from({ length: totalNodes }, (_, i) => i);
        const find = (x) => (parent[x] === x ? x : parent[x] = find(parent[x]));
        const unite = (a, b) => {
            const pa = find(a);
            const pb = find(b);
            if (pa !== pb)
                parent[pb] = pa;
        };
        // convert a nucleotide to its corresponding node reference (partial or stubs)...
        const getNodeRef = (nt) => {
            const partialId = idToPartial.get(nt.id);
            if (partialId !== undefined)
                return { node: partialId, kind: 'partial', index: partialId };
            const stubsId = idTostubs.get(nt.id);
            if (stubsId !== undefined)
                return { node: partials.length + stubsId, kind: 'stubs', index: stubsId };
            return null;
        };
        // Each partial has up to 2 sides (from mapPartialEnds). One side gets at most 1 connection to another partial.
        // 1-bp partials are excluded from partialEndsMap upstream — they get no sides here. Their
        // per-partial cap is enforced by noSideAttachCount (max 2 attachments) in the greedy below.
        // Build (partialIdx, ntId) -> sideIdx (0 or 1) so any exit-nt resolves to its side.
        const partialEndsMap = mapPartialEnds(partials);
        const ntToSide = new Map();
        partials.forEach((_, pIdx) => {
            const ends = partialEndsMap.get(pIdx);
            const inner = new Map();
            ntToSide.set(pIdx, inner);
            if (!ends)
                return;
            // Side 0: start1 (5' of strand A) paired with end2 (3' of strand B)
            inner.set(ends.start1.id, 0);
            inner.set(ends.end2.id, 0);
            // Side 1: end1 (3' of strand A) paired with start2 (5' of strand B)
            inner.set(ends.end1.id, 1);
            inner.set(ends.start2.id, 1);
        });
        const getSideForNt = (pIdx, ntId) => {
            return ntToSide.get(pIdx)?.get(ntId);
        };
        // Track direct adjacency between partials (used for stub-bridge safety check).
        const partialAdj = new Map();
        const addPartialAdj = (a, b) => {
            if (a === b)
                return;
            const setA = partialAdj.get(a) || new Set();
            setA.add(b);
            partialAdj.set(a, setA);
            const setB = partialAdj.get(b) || new Set();
            setB.add(a);
            partialAdj.set(b, setB);
        };
        const directEdges = [];
        const addDirectEdge = (a, sideA, b, sideB, dots) => {
            if (a === b)
                return;
            directEdges.push({ a, sideA, b, sideB, dots });
        };
        const stubsLinks = new Map();
        const addstubsLink = (stubNode, partialIdx, dots, strand, partialSide) => {
            const links = stubsLinks.get(stubNode) || new Map();
            const prev = links.get(partialIdx);
            if (!prev || dots > prev.dots) {
                links.set(partialIdx, { partialIdx, dots, strand, partialSide });
            }
            stubsLinks.set(stubNode, links);
        };
        const attachDot = (a, b, strand) => {
            let vecA = null;
            let vecB = null;
            if (a.kind === 'partial')
                vecA = getPartialStrandA3(a.index, strand);
            else
                vecA = getstubsA3(a.index);
            if (b.kind === 'partial')
                vecB = getPartialStrandA3(b.index, strand);
            else
                vecB = getstubsA3(b.index);
            if (!vecA || !vecB)
                return -1;
            return vecA.dot(vecB);
        };
        // First pass. Collects all data for connections (partial-partial or partial-stubs)
        systems.forEach(system => {
            system.strands.forEach(strand => {
                let prev = null;
                strand.forEach(elem => {
                    const nt = elem;
                    if (prev) {
                        const nodeA = getNodeRef(prev);
                        const nodeB = getNodeRef(nt);
                        if (nodeA && nodeB && nodeA.node !== nodeB.node) {
                            if (nodeA.kind === 'partial' && nodeB.kind === 'partial') {
                                addPartialAdj(nodeA.index, nodeB.index);
                            }
                            const d = attachDot(nodeA, nodeB, strand);
                            if (d > dot) {
                                if (nodeA.kind === 'partial' && nodeB.kind === 'partial') {
                                    // prev is the exit-nt of nodeA's partial; nt is the entry-nt of nodeB's partial.
                                    // code does not account for any partials that don't go into mapPartialEnds().
                                    const sideA = getSideForNt(nodeA.index, prev.id);
                                    const sideB = getSideForNt(nodeB.index, nt.id);
                                    addDirectEdge(nodeA.index, sideA, nodeB.index, sideB, d);
                                }
                                else if (nodeA.kind === 'stubs' || nodeB.kind === 'stubs') {
                                    const stubNode = nodeA.kind === 'stubs' ? nodeA : nodeB;
                                    const otherNode = nodeA.kind === 'stubs' ? nodeB : nodeA;
                                    const otherNt = nodeA.kind === 'stubs' ? nt : prev;
                                    if (otherNode.kind === 'partial') {
                                        const partialSide = getSideForNt(otherNode.index, otherNt.id);
                                        addstubsLink(stubNode.node, otherNode.index, d, strand, partialSide);
                                    }
                                }
                            }
                        }
                    }
                    prev = nt;
                });
            });
        });
        // Per-partial used-side amounts. For multi-bp partials (side is always defined when
        // used) this tracks how much of each side has been consumed. Partial-partial edges and
        // stub bridges consume a full side (1.0); ssDNA overhangs consume half a side (0.5),
        // so two ssDNA overhangs may share the same side. A side is free iff its total is < 1.
        const usedSides = new Map();
        // Per-partial attachment count for 0-side partials only. 1-bp partials have no
        // `partialEndsMap` entry, so `side` is undefined and nothing ever lands in
        // `usedSides` for them — this counter is the only way to cap them at 2 neighbors.
        const noSideAttachCount = new Map();
        const PER_PARTIAL_CAP = 2;
        const getNoSideCount = (pIdx) => noSideAttachCount.get(pIdx) ?? 0;
        const slotAvailable = (pIdx, side, amount = 1) => {
            if (side === undefined)
                return getNoSideCount(pIdx) + amount <= PER_PARTIAL_CAP;
            const used = usedSides.get(pIdx)?.get(side) ?? 0;
            return used + amount <= 1 + 1e-9;
        };
        const reserveSlot = (pIdx, side, amount = 1) => {
            if (side === undefined) {
                noSideAttachCount.set(pIdx, getNoSideCount(pIdx) + amount);
                return;
            }
            let bySide = usedSides.get(pIdx);
            if (!bySide) {
                bySide = new Map();
                usedSides.set(pIdx, bySide);
            }
            bySide.set(side, (bySide.get(side) ?? 0) + amount);
        };
        // Union-find helpers needed before the greedy pass for partial-group operations.
        const partialParent = Array.from({ length: partials.length }, (_, i) => i);
        const findPartial = (x) => (partialParent[x] === x ? x : partialParent[x] = findPartial(partialParent[x]));
        const partialMembers = new Map();
        for (let i = 0; i < partials.length; i++) {
            const set = partialMembers.get(i) || new Set();
            set.add(i);
            partialMembers.set(i, set);
        }
        const mergePartialGroups = (a, b) => {
            let ra = findPartial(a);
            let rb = findPartial(b);
            if (ra === rb)
                return ra;
            const setA = partialMembers.get(ra);
            const setB = partialMembers.get(rb);
            if (setA.size < setB.size) {
                const tmp = ra;
                ra = rb;
                rb = tmp;
            }
            const keep = partialMembers.get(ra);
            const drop = partialMembers.get(rb);
            drop.forEach(idx => keep.add(idx));
            partialMembers.set(ra, keep);
            partialMembers.delete(rb);
            partialParent[rb] = ra;
            return ra;
        };
        const hasDirectConnection = (rootA, rootB) => {
            if (rootA === rootB)
                return true;
            const setA = partialMembers.get(rootA);
            const setB = partialMembers.get(rootB);
            if (!setA || !setB)
                return false;
            const [small, large] = setA.size <= setB.size ? [setA, setB] : [setB, setA];
            for (const idx of small) {
                const neighbors = partialAdj.get(idx);
                if (!neighbors)
                    continue;
                for (const n of neighbors) {
                    if (large.has(n))
                        return true;
                }
            }
            return false;
        };
        const stubEdges = [];
        const partialsBridgeDot = (a, b) => {
            const vecA = getPartialStrandA3(a.partialIdx, a.strand);
            const vecB = getPartialStrandA3(b.partialIdx, b.strand);
            if (!vecA || !vecB)
                return -1;
            return vecA.dot(vecB);
        };
        stubsLinks.forEach((linksByPartial, stubNode) => {
            const candidates = Array.from(linksByPartial.values());
            if (candidates.length < 2)
                return;
            for (let i = 0; i < candidates.length; i++) {
                for (let j = i + 1; j < candidates.length; j++) {
                    const a = candidates[i];
                    const b = candidates[j];
                    const d = partialsBridgeDot(a, b);
                    if (d > dot) {
                        stubEdges.push({
                            a: a.partialIdx, sideA: a.partialSide,
                            b: b.partialIdx, sideB: b.partialSide,
                            dots: d,
                            stubNode
                        });
                    }
                }
            }
        });
        // Two-stage greedy merge:
        //
        // Stage A — direct partial<->partial edges only, sorted by dots desc.
        //   Backbone strand continuity (a direct strand-edge between two partials) is the
        //   strongest signal that those partials belong to the same helix. Run this first
        //   so a direct edge can claim a partial's side before any stub bridge can.
        //
        // Stage B — stub bridges only, on whatever sides remain unconsumed.
        //   A stub bridge legitimately merges two partials only when their respective sides
        //   are not already claimed by direct backbone neighbors. If a partial's side was
        //   consumed in Stage A, any stub bridge targeting that side is rejected here — that
        //   prevents a stub from yanking a partial out of its real (direct-edge) helix.
        //
        // Why the same-group check matters: as edges are accepted, partials get merged via
        // union-find. A later edge between two partials that are already in the same group
        // would be redundant — connecting them again does nothing structurally, but it would
        // still consume two sides, blocking those sides from a real cross-group merge.
        const sortedDirect = directEdges.slice().sort((x, y) => y.dots - x.dots);
        const sortedStub = stubEdges.slice().sort((x, y) => y.dots - x.dots);
        // Stage A: direct partial<->partial.
        for (const c of sortedDirect) {
            if (!slotAvailable(c.a, c.sideA))
                continue;
            if (!slotAvailable(c.b, c.sideB))
                continue;
            const rootA = findPartial(c.a);
            const rootB = findPartial(c.b);
            if (rootA === rootB)
                continue;
            unite(c.a, c.b);
            mergePartialGroups(c.a, c.b);
            reserveSlot(c.a, c.sideA, 1);
            reserveSlot(c.b, c.sideB, 1);
        }
        // Stage B: stub bridges. Same checks plus the existing safeguard against bridging
        // two groups that already have a direct partial-partial connection.
        for (const c of sortedStub) {
            if (!slotAvailable(c.a, c.sideA))
                continue;
            if (!slotAvailable(c.b, c.sideB))
                continue;
            const rootA = findPartial(c.a);
            const rootB = findPartial(c.b);
            if (rootA === rootB)
                continue;
            if (hasDirectConnection(rootA, rootB))
                continue;
            unite(c.stubNode, c.a);
            unite(c.stubNode, c.b);
            mergePartialGroups(c.a, c.b);
            reserveSlot(c.a, c.sideA, 1);
            reserveSlot(c.b, c.sideB, 1);
        }
        // Stubs join the best-aligned partial through A3 dots.
        stubsLinks.forEach((linksByPartial, stubNode) => {
            const ranked = Array.from(linksByPartial.values()).sort((x, y) => y.dots - x.dots);
            if (!ranked.length)
                return;
            const primary = ranked[0];
            unite(stubNode, primary.partialIdx);
        });
        const groups = new Map();
        partials.forEach((list, idx) => {
            const root = find(idx);
            const arr = groups.get(root) || [];
            arr.push(...list);
            groups.set(root, arr);
        });
        stubs.forEach((nt, idx) => {
            const node = partials.length + idx;
            const root = find(node);
            if (!groups.has(root)) {
                lastScraps.push([nt]);
                return;
            }
            const arr = groups.get(root) || [];
            arr.push(nt);
            groups.set(root, arr);
        });
        console.log('Merged Groups:', groups);
        // if everything goes right, you should NEVER have duplicates, helices or nucleotides. But this is good for safety.
        groups.forEach(group => {
            const seen = new Set();
            const unique = [];
            group.forEach(nt => {
                if (seen.has(nt.id)) {
                    console.log('Duplicate nucleotide found in helix grouping:', nt);
                    return;
                }
                seen.add(nt.id);
                unique.push(nt);
            });
            if (unique.length)
                helices.push(unique);
        });
        // After helices are built, attach ssScaffold segments to the helix they connect to.
        if (ssScaffold && ssScaffold.length) {
            const idToHelix = new Map();
            helices.forEach((list, idx) => {
                list.forEach(nt => idToHelix.set(nt.id, idx));
            });
            const ssScaffoldIds = new Set();
            ssScaffold.forEach(segment => segment.forEach(nt => ssScaffoldIds.add(nt.id)));
            // fairly obvious. Adds the segment of nucleotides (from ssScaffold) to the target helix index.
            const addToHelix = (targetIdx, segment) => {
                const helix = helices[targetIdx];
                const seen = new Set(helix.map(nt => nt.id));
                segment.forEach(nt => {
                    if (seen.has(nt.id))
                        return;
                    helix.push(nt);
                    seen.add(nt.id);
                    idToHelix.set(nt.id, targetIdx);
                });
            };
            const findSsScaffoldTargets = (segment) => {
                const segmentIds = new Set(segment.map(nt => nt.id));
                const helixIndices = new Set();
                segment.forEach(nt => {
                    const n5 = nt.n5;
                    const n3 = nt.n3;
                    if (n5 && !segmentIds.has(n5.id) && !ssScaffoldIds.has(n5.id)) {
                        const hIdx = idToHelix.get(n5.id);
                        if (hIdx !== undefined)
                            helixIndices.add(hIdx);
                    }
                    if (n3 && !segmentIds.has(n3.id) && !ssScaffoldIds.has(n3.id)) {
                        const hIdx = idToHelix.get(n3.id);
                        if (hIdx !== undefined)
                            helixIndices.add(hIdx);
                    }
                });
                return Array.from(helixIndices.values());
            };
            // For an ssScaffold segment, find the partial side through which it connects to a
            // specific helix. If the connection is not through a partial side (e.g. through a
            // stub or internal partial nucleotide), no side needs to be reserved.
            const findSsScaffoldConnectionSide = (segment, helixIdx) => {
                const segmentIds = new Set(segment.map(nt => nt.id));
                for (const nt of segment) {
                    for (const dir of ['n5', 'n3']) {
                        const neighbor = nt[dir];
                        if (!neighbor || segmentIds.has(neighbor.id) || ssScaffoldIds.has(neighbor.id))
                            continue;
                        if (idToHelix.get(neighbor.id) !== helixIdx)
                            continue;
                        const pIdx = idToPartial.get(neighbor.id);
                        if (pIdx === undefined)
                            continue;
                        return { pIdx, side: getSideForNt(pIdx, neighbor.id) };
                    }
                }
                return null;
            };
            let pending = ssScaffold.filter(segment => segment.length > 0);
            const maxRounds = Math.max(1, pending.length * 2);
            let round = 0;
            while (pending.length) {
                round += 1;
                let attachedThisRound = 0;
                const nextPending = [];
                pending.forEach(segment => {
                    const targets = findSsScaffoldTargets(segment);
                    if (!targets.length) {
                        nextPending.push(segment);
                        return;
                    }
                    // IN CASE that the ssScaffold segment connects to multiple helices, warn the user.
                    if (targets.length > 1) {
                        console.warn('ssScaffold segment connects to multiple helices; attaching to first available.', {
                            helices: targets,
                            segmentLength: segment.length,
                            round
                        });
                    }
                    // Try targets in order and attach to the first one whose partial side is available.
                    let attached = false;
                    for (const hIdx of targets) {
                        const conn = findSsScaffoldConnectionSide(segment, hIdx);
                        if (conn && !slotAvailable(conn.pIdx, conn.side, 0.5))
                            continue;
                        if (conn)
                            reserveSlot(conn.pIdx, conn.side, 0.5);
                        addToHelix(hIdx, segment);
                        attached = true;
                        break;
                    }
                    if (attached) {
                        attachedThisRound += 1;
                    }
                    else {
                        nextPending.push(segment);
                    }
                });
                if (!nextPending.length)
                    break;
                if (!attachedThisRound) {
                    console.warn('[ssScaffold] No attach progress in retry round; moving unresolved segments to lastScraps as grouped segments.', {
                        round,
                        unresolvedSegments: nextPending.length
                    });
                    nextPending.forEach(segment => lastScraps.push(segment.slice()));
                    break;
                }
                if (round >= maxRounds) {
                    console.warn('[ssScaffold] Retry limit reached; moving unresolved segments to lastScraps as grouped segments.', {
                        round,
                        unresolvedSegments: nextPending.length,
                        maxRounds
                    });
                    nextPending.forEach(segment => lastScraps.push(segment.slice()));
                    break;
                }
                pending = nextPending;
            }
        }
        const partialToHelix = new Map();
        helices.forEach((list, hIdx) => {
            list.forEach(nt => {
                const pIdx = idToPartial.get(nt.id);
                if (pIdx !== undefined && !partialToHelix.has(pIdx)) {
                    partialToHelix.set(pIdx, hIdx);
                }
            });
        });
        const addSegmentToHelix = (targetIdx, segment) => {
            const helix = helices[targetIdx];
            if (!helix)
                return;
            const seen = new Set(helix.map(nt => nt.id));
            segment.forEach(nt => {
                if (seen.has(nt.id))
                    return;
                helix.push(nt);
                seen.add(nt.id);
            });
        };
        // const oppositeDir = (dir: 'n5' | 'n3') => (dir === 'n5' ? 'n3' : 'n5');
        const findEndOnSide = (segment, segmentSet, dir) => {
            for (const nt of segment) {
                const neighbor = nt[dir];
                if (!neighbor || !segmentSet.has(neighbor.id))
                    return nt;
            }
            return null;
        };
        const walkForPartial = (start, dir, segmentSet) => {
            let curr = start;
            while (curr) {
                if (segmentSet.has(curr.id))
                    return undefined;
                const pIdx = idToPartial.get(curr.id);
                if (pIdx !== undefined)
                    return { pIdx, node: curr };
                curr = curr[dir];
            }
            return undefined;
        };
        const stepWithinSamePartial = (start, dir, steps, pIdx) => {
            let curr = start;
            let last = start;
            for (let i = 0; i < steps; i++) {
                curr = curr?.[dir];
                if (!curr)
                    return last;
                const idx = idToPartial.get(curr.id);
                if (idx !== pIdx)
                    return last;
                last = curr;
            }
            return last;
        };
        const stepN = (start, dir, steps) => {
            let curr = start;
            for (let i = 0; i < steps; i++) {
                curr = curr?.[dir];
                if (!curr)
                    return undefined;
            }
            return curr;
        };
        const classifySegment = (segment) => {
            const segmentSet = new Set(segment.map(nt => nt.id));
            const analyzeSide = (dir) => {
                const end = findEndOnSide(segment, segmentSet, dir);
                if (!end)
                    return { side: dir, result: 'overhang' };
                const anchor = end[dir];
                if (!anchor)
                    return { side: dir, result: 'overhang' };
                const first = walkForPartial(anchor, dir, segmentSet);
                if (!first)
                    return { side: dir, result: 'overhang' };
                const firstHelixId = partialToHelix.get(first.pIdx);
                const firstPartialSide = getSideForNt(first.pIdx, first.node.id);
                // note: using more than 1 step might look fine, but it can cause issues in edge cases.
                // specifically, structure 51 from nanobase (Dumbbell structure) has issues with this.
                // Sticking to 1 step has NOT shown ANY problems so far.
                const lastInPartial = stepWithinSamePartial(first.node, dir, 1, first.pIdx);
                const pair = lastInPartial.pair;
                if (!pair) {
                    return {
                        side: dir,
                        result: 'overhang',
                        firstPartialId: first.pIdx,
                        firstPartialSide,
                        firstHelixId
                    };
                }
                const oppositeNode = stepN(pair, dir, 3);
                if (!oppositeNode) {
                    return {
                        side: dir,
                        result: 'overhang',
                        firstPartialId: first.pIdx,
                        firstPartialSide,
                        firstHelixId
                    };
                }
                const oppositePartialId = idToPartial.get(oppositeNode.id);
                const oppositeHelixId = oppositePartialId !== undefined ? partialToHelix.get(oppositePartialId) : undefined;
                const binder = oppositePartialId !== undefined &&
                    firstHelixId !== undefined &&
                    oppositeHelixId !== undefined &&
                    firstHelixId === oppositeHelixId;
                return {
                    side: dir,
                    result: binder ? 'binder' : 'overhang',
                    firstPartialId: first.pIdx,
                    firstPartialSide,
                    oppositePartialId,
                    firstHelixId,
                    oppositeHelixId
                };
            };
            const res5 = analyzeSide('n5');
            const res3 = analyzeSide('n3');
            return { res5, res3 };
        };
        const binders = [];
        const binder2 = [];
        const disconnected = [];
        const unhandled = [];
        const isBinder = (res) => res?.result === 'binder';
        const isOverhang = (res) => res?.result === 'overhang';
        const hasPartial = (res) => res?.firstPartialId !== undefined;
        // Helper to claim a partial side for an ssDNA overhang. If the side is already
        // reserved (by a partial-partial edge or a stub bridge), the overhang is reclassified
        // as a binder so it follows binder routing rules instead.
        const tryReserveOverhangSide = (res) => {
            if (!res || res.result !== 'overhang')
                return res;
            if (res.firstPartialId === undefined)
                return res;
            if (!slotAvailable(res.firstPartialId, res.firstPartialSide, 0.5)) {
                return { ...res, result: 'binder' };
            }
            reserveSlot(res.firstPartialId, res.firstPartialSide, 0.5);
            return res;
        };
        // The lot of if statements are required (unless you can figure out a better way).
        // You can read through these, but they mostly comprise of cases where the segment is connected to helices on both ends, and has different types of such connections.
        // example, if overhang on one end and binder on the other, then it will connect to the helix on overhang side.
        ssdna.forEach(segment => {
            if (!segment.length)
                return;
            const raw = classifySegment(segment);
            const res5 = tryReserveOverhangSide(raw.res5);
            const res3 = tryReserveOverhangSide(raw.res3);
            const res5HasPartial = hasPartial(res5);
            const res3HasPartial = hasPartial(res3);
            if (!res5HasPartial && !res3HasPartial) {
                disconnected.push(segment);
                return;
            }
            if (isOverhang(res5) && hasPartial(res5) && isOverhang(res3) && !hasPartial(res3)) {
                if (res5.firstHelixId !== undefined)
                    addSegmentToHelix(res5.firstHelixId, segment);
                return;
            }
            if (isOverhang(res3) && hasPartial(res3) && isOverhang(res5) && !hasPartial(res5)) {
                if (res3.firstHelixId !== undefined)
                    addSegmentToHelix(res3.firstHelixId, segment);
                return;
            }
            if (isOverhang(res5) && res5HasPartial && isOverhang(res3) && res3HasPartial) {
                if (res5.firstHelixId !== undefined && res3.firstHelixId !== undefined) {
                    const half = Math.floor(segment.length / 2);
                    const left = segment.slice(0, half);
                    const right = segment.slice(half);
                    addSegmentToHelix(res5.firstHelixId, left);
                    addSegmentToHelix(res3.firstHelixId, right);
                    return;
                }
            }
            if (isOverhang(res5) && hasPartial(res5) && isBinder(res3)) {
                if (res5.firstHelixId !== undefined)
                    addSegmentToHelix(res5.firstHelixId, segment);
                return;
            }
            if (isOverhang(res3) && hasPartial(res3) && isBinder(res5)) {
                if (res3.firstHelixId !== undefined)
                    addSegmentToHelix(res3.firstHelixId, segment);
                return;
            }
            if (isBinder(res5) && !res5HasPartial && isOverhang(res3) && res3HasPartial) {
                if (res3.firstHelixId !== undefined)
                    addSegmentToHelix(res3.firstHelixId, segment);
                return;
            }
            if (isBinder(res3) && !res3HasPartial && isOverhang(res5) && res5HasPartial) {
                if (res5.firstHelixId !== undefined)
                    addSegmentToHelix(res5.firstHelixId, segment);
                return;
            }
            if (isBinder(res5) && isOverhang(res3) && !hasPartial(res3)) {
                binders.push({ segment, res5, res3 });
                return;
            }
            if (isBinder(res3) && isOverhang(res5) && !hasPartial(res5)) {
                binders.push({ segment, res5, res3 });
                return;
            }
            if (isBinder(res5) && isBinder(res3)) {
                binder2.push({ segment, res5, res3 });
                return;
            }
            unhandled.push(segment);
        });
        // For any binder/binder2 segments, group them by which helix they connect to.
        // If multiple binder segments connect to the same helix, they form a new helix.
        const resolveBinderHelix = (entry) => {
            const { res5, res3 } = entry;
            const helixIds = new Set();
            const collect = (res) => {
                if (!res || !isBinder(res))
                    return;
                if (res.firstHelixId !== undefined)
                    helixIds.add(res.firstHelixId);
                if (res.oppositeHelixId !== undefined)
                    helixIds.add(res.oppositeHelixId);
            };
            collect(res5);
            collect(res3);
            if (helixIds.size === 1)
                return Array.from(helixIds.values())[0];
            return undefined;
        };
        const binderGroups = new Map();
        const addBinderToGroup = (helixId, segment) => {
            const list = binderGroups.get(helixId) || [];
            list.push(segment);
            binderGroups.set(helixId, list);
        };
        // For binder2 segments that span two distinct helices, group by the unordered helix pair.
        // All binder2 segments connecting the SAME two helices get merged into a single new helix.
        // Binder2 segments connecting a DIFFERENT pair get their own new helix.
        const pairKey = (a, b) => (a < b ? `${a}:${b}` : `${b}:${a}`);
        const binderPairGroups = new Map();
        const addBinderToPairGroup = (a, b, segment) => {
            const key = pairKey(a, b);
            const list = binderPairGroups.get(key) || [];
            list.push(segment);
            binderPairGroups.set(key, list);
        };
        // Collect helix ids touched by binder-classified sides of a segment.
        // For binder2 the result has size 1 (both sides agree) or 2 (sides resolve to different helices).
        const getBinderHelixIds = (entry) => {
            const { res5, res3 } = entry;
            const ids = new Set();
            const collect = (res) => {
                if (!res || !isBinder(res))
                    return;
                if (res.firstHelixId !== undefined)
                    ids.add(res.firstHelixId);
                if (res.oppositeHelixId !== undefined)
                    ids.add(res.oppositeHelixId);
            };
            collect(res5);
            collect(res3);
            return Array.from(ids.values());
        };
        binders.forEach(entry => {
            const helixId = resolveBinderHelix(entry);
            if (helixId === undefined)
                return;
            addBinderToGroup(helixId, entry.segment);
        });
        binder2.forEach(entry => {
            const ids = getBinderHelixIds(entry);
            if (ids.length === 1) {
                // Single host helix: same path as a normal binder.
                addBinderToGroup(ids[0], entry.segment);
            }
            else if (ids.length === 2) {
                // Two host helices: group by the unordered pair so all binder2 segments
                // spanning the same {A, B} pair fuse into one new helix together.
                addBinderToPairGroup(ids[0], ids[1], entry.segment);
            }
            // ids.length === 0 or > 2: should be unreachable for binder2; silently dropped.
        });
        const materializeBinderHelix = (segments) => {
            if (!segments.length)
                return;
            const seen = new Set();
            const newHelix = [];
            segments.forEach(segment => {
                segment.forEach(nt => {
                    if (seen.has(nt.id))
                        return;
                    seen.add(nt.id);
                    newHelix.push(nt);
                });
            });
            if (newHelix.length)
                helices.push(newHelix);
        };
        binderGroups.forEach(segments => materializeBinderHelix(segments));
        binderPairGroups.forEach(segments => materializeBinderHelix(segments));
        // This is the code to re-attach lastScraps[] segments to the closest helix by partial connection (n5/n3).
        // For grouped segments (e.g. deferred ssScaffold), we attach the full segment to one chosen helix.
        if (lastScraps.length && helices.length) {
            const idToHelix = new Map();
            helices.forEach((list, idx) => {
                list.forEach(nt => idToHelix.set(nt.id, idx));
            });
            const walkToHelix = (start, dir, owner) => {
                if (!start) {
                    console.warn('[walkToHelix] Side', dir, 'is null for nucleotide', owner.id, '; using opposite side if available.');
                    return null;
                }
                let curr = start;
                while (curr) {
                    // which helix does this current nt belong to?
                    const hIdx = idToHelix.get(curr.id);
                    if (hIdx !== undefined)
                        return { helixIdx: hIdx, anchor: curr };
                    curr = curr[dir];
                }
                console.warn('[walkToHelix] Side', dir, 'for nucleotide', owner.id, 'started at', start.id, 'but did not reach any existing helix.');
                return null;
            };
            const addlastScrapsSegmentToHelix = (targetIdx, segment) => {
                const helix = helices[targetIdx];
                if (!helix)
                    return false;
                const seen = new Set(helix.map(nt => nt.id));
                segment.forEach(nt => {
                    if (seen.has(nt.id))
                        return;
                    helix.push(nt);
                    seen.add(nt.id);
                    idToHelix.set(nt.id, targetIdx);
                });
                return true;
            };
            const pickSegmentTarget = (segment) => {
                let best = null;
                segment.forEach(nt => {
                    const via5 = walkToHelix((nt.n5 ?? null), 'n5', nt);
                    const via3 = walkToHelix((nt.n3 ?? null), 'n3', nt);
                    if (!via5 && via3) {
                        console.log('[walkToHelix] Nucleotide', nt.id, ': n5 lookup failed; using n3 fallback candidate to helix', via3.helixIdx, 'via anchor', via3.anchor.id);
                    }
                    if (!via3 && via5) {
                        console.log('[walkToHelix] Nucleotide', nt.id, ': n3 lookup failed; using n5 fallback candidate to helix', via5.helixIdx, 'via anchor', via5.anchor.id);
                    }
                    if (!via5 && !via3) {
                        console.warn('[walkToHelix] Nucleotide', nt.id, ': both n5 and n3 lookups failed while resolving segment target.');
                        return;
                    }
                    const pos = nt.getPos();
                    const consider = (hit) => {
                        const distance = pos.distanceTo(hit.anchor.getPos());
                        if (!best || distance < best.distance) {
                            best = { hit, distance };
                        }
                    };
                    if (via5)
                        consider(via5);
                    if (via3)
                        consider(via3);
                });
                return best ? best.hit : null;
            };
            const remaining = [];
            lastScraps.forEach(segment => {
                if (!segment.length)
                    return;
                const target = pickSegmentTarget(segment);
                if (!target) {
                    console.warn('[walkToHelix] Could not resolve target helix for lastScraps segment; keeping grouped segment in lastScraps.', {
                        segmentLength: segment.length,
                        segmentIds: segment.map(nt => nt.id)
                    });
                    remaining.push(segment);
                    return;
                }
                if (!addlastScrapsSegmentToHelix(target.helixIdx, segment)) {
                    remaining.push(segment);
                    return;
                }
                console.log('[walkToHelix] Attached lastScraps segment to helix', target.helixIdx, 'segmentLength', segment.length);
            });
            // push the remaining grouped segments back to lastScraps[].
            lastScraps.length = 0;
            lastScraps.push(...remaining);
        }
        // const finalHelices = helices.filter(h => h.length > 0);
        return { helices, lastScraps, binders, binder2, disconnected, unhandled, usedSides };
    }
    helix_1.generateHelix = generateHelix;
    function findHelices(inputMap, tolerance = 2) {
        findBasepairsOptim2();
        dropIntraStrandPairs();
        // ok now we can do the rest of the stuff.
        let { partials, unpaired } = findHelixPartials2(inputMap, tolerance);
        let { ssdna, stubs, longssScaffold } = ssdnaPartials(unpaired);
        let ssScaffold = longssScaffoldfunc(longssScaffold, stubs);
        let { helices, lastScraps, binders, binder2, disconnected, unhandled, usedSides } = generateHelix(partials, ssdna, ssScaffold, stubs);
        console.log("Helices size:", helices.flat().length);
        console.log("Total elements:", inputMap.size);
        return { helices, partials, usedSides };
    }
    helix_1.findHelices = findHelices;
    // Merge two or more helices into the one with the lowest index.
    // Mutates `helices` in place: pushes nucleotides from higher-indexed entries into the kept helix
    // and SPLICES those entries out, so the array shrinks. Also mutates `grid` in place: every
    // GridMark whose helixId points at a merged-away helix is rewritten to the kept id, and every
    // surviving helixId is shifted down through the same remap that's returned. Returns an
    // `idRemap` (oldIdx -> newIdx) that callers still need for any external references that key
    // off helix index — gridview node ids, crossover connection ids, etc.
    // No checks for grid layout yet
    function combineHelices(helices, indices, grid) {
        if (!Array.isArray(helices) || !Array.isArray(indices) || !(grid instanceof Map))
            return null;
        const valid = [];
        const seenIdx = new Set();
        indices.forEach(i => {
            if (i < 0 || i >= helices.length)
                return;
            if (!Array.isArray(helices[i]) || helices[i].length === 0)
                return;
            if (seenIdx.has(i))
                return;
            seenIdx.add(i);
            valid.push(i);
        });
        if (valid.length < 2)
            return null;
        // merge by keeping the lowest id only
        valid.sort((a, b) => a - b);
        const keptIdxOld = valid[0];
        const mergedIdxOld = valid.slice(1);
        // Move nucleotides into the kept helix, deduping by id
        const seenNts = new Set(helices[keptIdxOld].map(nt => nt.id));
        mergedIdxOld.forEach(idx => {
            helices[idx].forEach(nt => {
                if (seenNts.has(nt.id))
                    return;
                seenNts.add(nt.id);
                helices[keptIdxOld].push(nt);
            });
        });
        // Build an oldIdx -> newIdx remap for every helix that survives the splice.
        const removed = new Set(mergedIdxOld);
        const idRemap = new Map();
        let shift = 0;
        for (let i = 0; i < helices.length; i++) {
            if (removed.has(i)) {
                shift += 1;
                continue;
            }
            idRemap.set(i, i - shift);
        }
        // Apply the same remap to the grid so per-nucleotide helixIds stay consistent with the helices array. 
        // Marks pointing at a merged-away helix collapse onto the kept id
        // marks on survivors shift down through idRemap. keptIdxOld is the lowest valid index, so its new id equals its old id
        grid.forEach(mark => {
            if (removed.has(mark.helixId)) {
                mark.helixId = keptIdxOld;
                return;
            }
            const next = idRemap.get(mark.helixId);
            if (next !== undefined)
                mark.helixId = next;
        });
        // Splice in reverse so earlier indices stay valid during removal.
        for (let i = helices.length - 1; i >= 0; i--) {
            if (removed.has(i))
                helices.splice(i, 1);
        }
        // keptIdx is the lowest valid index, so nothing in front of it was removed: its new index
        // is the same as its old one. Look it up via idRemap to stay correct if this invariant ever changes.
        const keptIdx = idRemap.get(keptIdxOld) ?? keptIdxOld;
        return { keptIdx, mergedIdx: mergedIdxOld, idRemap };
    }
    helix_1.combineHelices = combineHelices;
    // Reverse a combineHelices call using a snapshot recorded at merge time.
    // Mutates `helices` in place: rebuilds at full pre-merge length, places each survivor back at
    // its old index, and repopulates each removed slot from its captured `ntIds`. Mutates `grid`
    // in place: marks for nucleotides in any removed slot get their original `helixId` back;
    // every other mark shifts up through the inverse remap.
    // Returns an `inverseRemap(currentId) -> oldId` helper so callers can fix up external state
    // (gridview node ids, connection endpoints) that's still keyed in post-merge numbering.
    function splitHelices(helices, grid, snapshot) {
        if (!Array.isArray(helices) || !(grid instanceof Map) || !snapshot)
            return null;
        if (!Array.isArray(snapshot.idRemap) || !Array.isArray(snapshot.removed))
            return null;
        // Inverse remap: post-merge (current) index -> pre-merge (old) index. Only survivors are
        // in this map; ids outside it are passed through unchanged.
        const survivorInverse = new Map();
        snapshot.idRemap.forEach(([oldIdx, newIdx]) => survivorInverse.set(newIdx, oldIdx));
        const inverseRemap = (currentId) => {
            const oldIdx = survivorInverse.get(currentId);
            return oldIdx !== undefined ? oldIdx : currentId;
        };
        // Locate the kept helix in its current (post-merge) position.
        const keptOldIdx = snapshot.kept;
        let keptCurrentIdx = keptOldIdx;
        for (const [newIdx, oldIdx] of survivorInverse.entries()) {
            if (oldIdx === keptOldIdx) {
                keptCurrentIdx = newIdx;
                break;
            }
        }
        // Pull merged-away nucleotides out of the kept helix.
        const ntsToExtract = new Set();
        snapshot.removed.forEach(slot => slot.ntIds.forEach(id => ntsToExtract.add(id)));
        const keptArr = helices[keptCurrentIdx] || [];
        const keptKept = [];
        keptArr.forEach(nt => { if (!ntsToExtract.has(nt.id))
            keptKept.push(nt); });
        // Rebuild helices at full pre-merge length: survivors back at their old indices, removed
        // slots repopulated by looking up nucleotides in the global elements map.
        const total = helices.length + snapshot.removed.length;
        const rebuilt = new Array(total);
        for (let curIdx = 0; curIdx < helices.length; curIdx++) {
            const oldIdx = survivorInverse.get(curIdx);
            if (oldIdx === undefined)
                continue;
            rebuilt[oldIdx] = curIdx === keptCurrentIdx ? keptKept : helices[curIdx];
        }
        snapshot.removed.forEach(slot => {
            const restored = [];
            slot.ntIds.forEach(id => {
                const nt = elements.get(id);
                if (nt instanceof Nucleotide)
                    restored.push(nt);
            });
            rebuilt[slot.oldIdx] = restored;
        });
        // Mutate the original array in place so existing references stay valid.
        helices.length = 0;
        rebuilt.forEach((slot, i) => { helices[i] = slot || []; });
        // Restore GridMap helixIds. A nucleotide that was in a removed slot snaps back to that
        // slot's oldIdx; every other mark goes through the inverse remap.
        const ntToOldHelix = new Map();
        snapshot.removed.forEach(slot => slot.ntIds.forEach(id => ntToOldHelix.set(id, slot.oldIdx)));
        grid.forEach((mark, ntId) => {
            const overriden = ntToOldHelix.get(ntId);
            if (overriden !== undefined) {
                mark.helixId = overriden;
            }
            else {
                mark.helixId = inverseRemap(mark.helixId);
            }
        });
        return { inverseRemap };
    }
    helix_1.splitHelices = splitHelices;
    // Fits a plane through the given points and returns the plane normal
    // (the eigenvector of the covariance-like matrix with the smallest eigenvalue).
    function fitPlane(points) {
        // centroid
        const rc = new THREE.Vector3(0, 0, 0);
        points.forEach(p => rc.add(p));
        rc.divideScalar(points.length);
        // 3x3 symmetric accumulator A[i][j] += (p-rc)[i] * (p-rc)[j]
        const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        points.forEach(p => {
            const q = [p.x - rc.x, p.y - rc.y, p.z - rc.z];
            for (let i = 0; i < 3; i++) {
                for (let j = 0; j < 3; j++) {
                    A[i][j] += q[i] * q[j];
                }
            }
        });
        // Jacobi eigen-decomposition on 3x3 symmetric A (analog of numpy.linalg.eigh).
        const a = A.map(r => r.slice());
        const v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        const maxIter = 100;
        const eps = 1e-12;
        for (let iter = 0; iter < maxIter; iter++) {
            // find largest off-diagonal entry
            let p = 0, q = 1;
            let maxOff = Math.abs(a[0][1]);
            if (Math.abs(a[0][2]) > maxOff) {
                p = 0;
                q = 2;
                maxOff = Math.abs(a[0][2]);
            }
            if (Math.abs(a[1][2]) > maxOff) {
                p = 1;
                q = 2;
                maxOff = Math.abs(a[1][2]);
            }
            if (maxOff < eps)
                break;
            const app = a[p][p], aqq = a[q][q], apq = a[p][q];
            const theta = (aqq - app) / (2 * apq);
            const t = theta >= 0
                ? 1 / (theta + Math.sqrt(theta * theta + 1))
                : 1 / (theta - Math.sqrt(theta * theta + 1));
            const c = 1 / Math.sqrt(t * t + 1);
            const s = t * c;
            a[p][p] = app - t * apq;
            a[q][q] = aqq + t * apq;
            a[p][q] = 0;
            a[q][p] = 0;
            for (let i = 0; i < 3; i++) {
                if (i !== p && i !== q) {
                    const aip = a[i][p], aiq = a[i][q];
                    a[i][p] = c * aip - s * aiq;
                    a[p][i] = a[i][p];
                    a[i][q] = c * aiq + s * aip;
                    a[q][i] = a[i][q];
                }
                const vip = v[i][p], viq = v[i][q];
                v[i][p] = c * vip - s * viq;
                v[i][q] = c * viq + s * vip;
            }
        }
        const vals = [a[0][0], a[1][1], a[2][2]];
        // numpy.linalg.eigh returns eigenvalues in ascending order; the plane normal
        // is the eigenvector for the smallest eigenvalue (vecs[:, 0]).
        let minIdx = 0;
        if (vals[1] < vals[minIdx])
            minIdx = 1;
        if (vals[2] < vals[minIdx])
            minIdx = 2;
        return new THREE.Vector3(v[0][minIdx], v[1][minIdx], v[2][minIdx]);
    }
    helix_1.fitPlane = fitPlane;
    // Returns the axis of an DNA duplex given the four end nucleotides of the two strands.
    // - start1/end1 are the 5'/3' ends of strand A
    // - start2/end2 are the 5'/3' ends of strand B (start1 pairs with end2, end1 pairs with start2)
    function getPartialAxis(d) {
        const backboneSite = (nt) => nt.getInstanceParameter3('bbOffsets');
        // initial guess vector from the midpoint of the start1-end2 pair to end1-start2 pair
        const midA0 = backboneSite(d.start1).add(backboneSite(d.end2)).multiplyScalar(0.5);
        const midAc0 = backboneSite(d.end1).add(backboneSite(d.start2)).multiplyScalar(0.5);
        const guess = midAc0.clone().sub(midA0);
        if (guess.length() > 0)
            guess.normalize();
        // Walk pairs (nucA on strand A via n3, nucB on strand B via n5) in lockstep.
        // nucAc/nucBc are the next pair along the walk.
        const posAs = [];
        const posBs = [];
        const backPoses = [];
        let nucA = d.start1;
        let nucB = d.end2;
        while (nucA && nucB && nucA !== d.end1) {
            const nucAc = nucA.n3;
            const nucBc = nucB.n5;
            if (!nucAc || !nucBc)
                break;
            posAs.push(backboneSite(nucA));
            posBs.push(backboneSite(nucB));
            // on the last iteration, also push the trailing pair (mirrors Python's `if i == end1-1`)
            if (nucAc === d.end1) {
                posAs.push(backboneSite(nucAc));
                posBs.push(backboneSite(nucBc));
            }
            backPoses.push(backboneSite(nucAc).sub(backboneSite(nucA)));
            backPoses.push(backboneSite(nucBc).sub(backboneSite(nucB)));
            nucA = nucAc;
            nucB = nucBc;
        }
        const planeVector = fitPlane(backPoses);
        if (guess.dot(planeVector) < 0)
            planeVector.multiplyScalar(-1);
        // Find where the helical axis originates by intersecting per-base-pair perpendiculars
        // projected onto the plane.
        const helPos = [];
        for (let i = 0; i < posAs.length - 1; i++) {
            // project current base pair to plane
            let apos = posAs[i].clone();
            let bpos = posBs[i].clone();
            apos.sub(planeVector.clone().multiplyScalar(apos.dot(planeVector)));
            bpos.sub(planeVector.clone().multiplyScalar(bpos.dot(planeVector)));
            const bpVecA = bpos.clone().sub(apos);
            if (bpVecA.length() === 0)
                continue;
            const midpointA = apos.clone().add(bpos).multiplyScalar(0.5);
            const perpA = bpVecA.clone().cross(planeVector).normalize();
            // project next base pair to plane
            let apos2 = posAs[i + 1].clone();
            let bpos2 = posBs[i + 1].clone();
            apos2.sub(planeVector.clone().multiplyScalar(apos2.dot(planeVector)));
            bpos2.sub(planeVector.clone().multiplyScalar(bpos2.dot(planeVector)));
            const bpVecB = bpos2.clone().sub(apos2);
            if (bpVecB.length() === 0)
                continue;
            const midpointB = apos2.clone().add(bpos2).multiplyScalar(0.5);
            const perpB = bpVecB.clone().cross(planeVector).normalize();
            // Solve the 3x2 least-squares system:
            //   [perpA, -perpB] [t, c]^T = midpointB - midpointA
            // via 2x2 normal equations. perpA, perpB are unit vectors so their self-dots are 1.
            const y = midpointB.clone().sub(midpointA);
            const m01 = -perpA.dot(perpB); // = m10
            const b0 = perpA.dot(y);
            const b1 = -perpB.dot(y);
            const det = 1 - m01 * m01;
            if (Math.abs(det) < 1e-12)
                continue;
            const t = (b0 - m01 * b1) / det;
            const c = (b1 - m01 * b0) / det;
            const pointA = midpointA.clone().add(perpA.clone().multiplyScalar(t));
            const pointB = midpointB.clone().add(perpB.clone().multiplyScalar(c));
            if (pointA.distanceTo(pointB) > 1e-6) {
                console.log('Error in finding common intersection point', pointA, pointB);
            }
            helPos.push(pointA);
        }
        const finalHelPos = new THREE.Vector3(0, 0, 0);
        helPos.forEach(p => finalHelPos.add(p));
        if (helPos.length)
            finalHelPos.divideScalar(helPos.length);
        return { planeVector, finalHelPos };
    }
    helix_1.getPartialAxis = getPartialAxis;
    // Draws the duplex axis vector in the scene, anchored at the start1 nucleotide.
    function addPartialAxisToScene(d) {
        const { planeVector } = getPartialAxis(d);
        const origin = d.start1.getInstanceParameter3('bbOffsets')
            .add(d.end2.getInstanceParameter3('bbOffsets'))
            .multiplyScalar(0.5);
        if (typeof THREE !== 'undefined' && typeof scene !== 'undefined' && scene?.add) {
            const arrow = new THREE.ArrowHelper(planeVector.clone().normalize(), origin, 10);
            scene.add(arrow);
        }
    }
    helix_1.addPartialAxisToScene = addPartialAxisToScene;
})(helix || (helix = {}));
