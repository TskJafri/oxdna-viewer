/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
/*
// for ease of use, and to prevent dumb mistakes as I code and test things out, here is exactly the commands to use this in console:
findBasepairs3(); // use 3 because there are 2 versions of findBasepairs, and 3 is the fastest. Dont ask why i named it that.
honda.dropIntraStrandPairs();
let {partials, fishies} = honda.findHelixPartials(elements, 2);
let {ssdna, deadfishies, longssScaffold} = honda.ssdnaPartials(fishies);
let ssScaffold = honda.longssScaffoldfunc(longssScaffold, deadfishies);
let {helices, murdered, binders, binder2, disconnected, unhandled} = honda.generateHelix(partials, ssdna, ssScaffold, deadfishies);
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
// let helices = await honda.findHelices(elements, 2);
var honda;
(function (honda) {
    // helper function cuz didnt want to type this every time
    function checkAngle(n1 = null, n2 = null) {
        if (!n1 || !n2)
            return 0;
        return Math.acos(n1.getA3().dot(n2.getA3())) * (180 / Math.PI);
    }
    honda.checkAngle = checkAngle;
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
    honda.getScaffoldStrand = getScaffoldStrand;
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
    honda.dropIntraStrandPairs = dropIntraStrandPairs;
    // Finds helix parts using destructive consumption of a working copy of elements (called mermaid).
    // tolerance 2 is good enough for most cases. Higher tolerances seem to have no negative consequences, however.
    function findHelixPartials(inputMap, tolerance = 2) {
        const mermaid = new Map(inputMap);
        const mermaid2 = new Map(inputMap); // backup copy for duplication, used later
        const fishies = new Map(); // unpaired / skipped nts
        let partials = [];
        const record = (list, set, nt = null) => {
            if (nt && !set.has(nt.id)) {
                set.add(nt.id);
                list.push(nt);
            }
        };
        const nextStart = () => mermaid.values().next().value;
        while (true) {
            const start = nextStart();
            // console.log('Starting of partials:', start);
            if (!start)
                break;
            // Collect unpaired/binder nts for downstream ssDNA processing instead of discarding silently.
            // Additionally, only walk if the pair is also present in the current pool (mermaid).
            const pairInPool = start.pair ? mermaid.get(start.pair.id) : undefined;
            if (!start.pair || !pairInPool) {
                fishies.set(start.id, start);
                mermaid.delete(start.id);
                continue;
            }
            // initialize the walker.
            let curr = start;
            let ally = pairInPool;
            const strandA = curr.strand; // required to ensure we don't cross strands. This is how we know partials is actually a partial helix.
            const strandB = ally.strand;
            // this is the partial being built.
            const partial = [];
            const seen = new Set();
            const tryDirection = (dir) => {
                const nextCurr = mermaid.get(curr.id + dir);
                const nextAlly = mermaid.get(ally.id - dir);
                if (!nextCurr || !nextAlly)
                    return null;
                if (nextCurr.strand !== strandA || nextAlly.strand !== strandB)
                    return null;
                // Tolerance window: OR logic. Continue if ANY offset in [1, tolerance]
                // has forward paired to backward (allows rescue by curr+2/a-2, curr+3/a-3, etc.).
                let onwards = false;
                for (let offset = 1; offset <= tolerance; offset++) {
                    const forward = mermaid.get(curr.id + dir * offset);
                    const backward = mermaid.get(ally.id - dir * offset);
                    if (!forward || !backward)
                        continue; // check if either of them even exist.
                    if (forward.strand !== strandA || backward.strand !== strandB)
                        continue;
                    if (forward.pair === backward) {
                        onwards = true;
                        break;
                    }
                }
                if (!onwards)
                    return null;
                return { nextCurr, nextAlly };
            };
            // actual traversal loop.
            while (curr && ally) {
                record(partial, seen, curr);
                record(partial, seen, ally);
                // destructive consumption. This is why we make a copy of elements, and not use the original map directly.
                mermaid.delete(curr.id);
                mermaid.delete(ally.id);
                const step = tryDirection(1) || tryDirection(-1);
                if (!step)
                    break;
                // Always consume the immediate neighbors (curr+1 and a-1) even if mismatched.
                record(partial, seen, step.nextCurr);
                record(partial, seen, step.nextAlly);
                mermaid.delete(step.nextCurr.id);
                mermaid.delete(step.nextAlly.id);
                // Advance walker by one along each strand.
                curr = step.nextCurr;
                ally = step.nextAlly;
            }
            if (partial.length) {
                partials.push(partial);
            }
        }
        // Deduplicate across all partials: any nucleotide that appears more than once
        // is moved to fishies along with its pair and any nucleotide paired to it.
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
        // wait lately its being doing some weird crap, lets log it out.
        // had to add the multi-pairing logic (down below) to fix the problem. But still logging for curiosity.
        console.log('Duplicates set: ', duplicates);
        // const pairIds = new Set<number>();
        // partials.forEach((helix,i) => {
        // 	partials[i] = helix.filter(nt => {
        // 		if (duplicates.has(nt.id) || pairIds.has(nt.id)) {
        // 			fishies.set(nt.id, nt);
        // 			return false;
        // 		}
        // 		if (nt.pair && duplicates.has(nt.pair.id)) {
        // 			fishies.set(nt.id, nt);
        // 			return false;
        // 		}
        // 		return true;
        // 	});
        // });
        if (duplicates.size) {
            // Collect pair ids for duplicates
            console.log('duplicates found: ', duplicates);
            console.log('total duplicates: ', duplicates.size);
            const pairIds = new Set();
            duplicates.forEach(id => {
                const a = mermaid2.get(id);
                if (a) {
                    fishies.set(a.id, a);
                    if (a.pair) {
                        pairIds.add(a.pair.id);
                        fishies.set(a.pair.id, a.pair);
                    }
                }
            });
            // Remove duplicates, their pairs, and any nucleotide whose pair is a duplicate
            partials.forEach((helix, i) => {
                partials[i] = helix.filter(nt => {
                    if (duplicates.has(nt.id) || pairIds.has(nt.id)) {
                        fishies.set(nt.id, nt);
                        return false;
                    }
                    if (nt.pair && duplicates.has(nt.pair.id)) {
                        fishies.set(nt.id, nt);
                        return false;
                    }
                    return true;
                });
            });
            // // Drop empty partial lists
            partials = partials.filter(helix => helix.length > 0);
        }
        // Remove any nucleotide that is part of a multi-pairing (2+ nts paired to the same nt)
        // dang it. This doesnt work either, specificially at the ends of the helices, where the fraying causes multi-pairing.
        // const pairedTo = new Map<number, number[]>();
        // mermaid2.forEach(nt => {
        // 	if (!nt.pair) return;
        // 	const list = pairedTo.get(nt.pair.id) || [];
        // 	list.push(nt.id);
        // 	pairedTo.set(nt.pair.id, list);
        // });
        // const multiPairIds = new Set<number>();
        // pairedTo.forEach((list, targetId) => {
        // 	if (list.length >= 2) {
        // 		multiPairIds.add(targetId);
        // 		list.forEach(id => multiPairIds.add(id));
        // 	}
        // });
        // if (multiPairIds.size) {
        // 	partials.forEach((helix, i) => {
        // 		partials[i] = helix.filter(nt => {
        // 			if (multiPairIds.has(nt.id)) {
        // 				fishies.set(nt.id, nt);
        // 				return false;
        // 			}
        // 			return true;
        // 		});
        // 	});
        // 	partials = partials.filter(helix => helix.length > 0);
        // }
        return { partials, fishies: Array.from(fishies.values()) };
    }
    honda.findHelixPartials = findHelixPartials;
    // Groups unpaired/binder nucleotides (fishies) into ssDNA partials by strand.
    // Only contiguous runs (>2) along a strand are kept; shorter runs go to deadfishies.
    function ssdnaPartials(fishies) {
        const byStrand = new Map();
        const ssdna = [];
        const deadfishies = [];
        const longssScaffold = [];
        const scaffold = getScaffoldStrand();
        fishies.forEach(nt => {
            const arr = byStrand.get(nt.strand) || [];
            arr.push(nt);
            byStrand.set(nt.strand, arr);
        });
        // For each strand, sort it, find it's 5' ends and walk down n3 to build runs.
        byStrand.forEach((list, strand) => {
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
                        run.forEach(r => deadfishies.push(r));
                    }
                }
                // // expand run both directions along n5/n3 within fishies set
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
                // 	run.forEach(r => deadfishies.push(r));
                // }
            }
        });
        return { ssdna, deadfishies, longssScaffold };
    }
    honda.ssdnaPartials = ssdnaPartials;
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
    honda.averageA3a = averageA3a;
    ;
    function longssScaffoldfunc(longssScaffold, deadfishies = []) {
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
                runAway.forEach(nt => deadfishies.push(nt));
                runAway = [];
                return;
            }
            // Enforce equal halves; trim one nucleotide if odd-length to satisfy the requirement.
            const evenLen = runAway.length - (runAway.length % 2);
            if (evenLen !== runAway.length) {
                const dropped = runAway[evenLen];
                if (dropped)
                    deadfishies.push(dropped);
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
    honda.longssScaffoldfunc = longssScaffoldfunc;
    // let partialStrandMap = new Map<number, Map<number, Nucleotide[]>>();
    // Perfected!
    // this one uses average a3 vectors of CONNECTED strands, as opposed to average a3 vectors of the entire partial (which cancels out, due to topology).
    function generateHelix(partials, ssdna, ssScaffold, deadfishies) {
        // Currently uses partials, ssScaffold and deadfishies to build perfect (almost) helices.
        // Hence, helices.flat().length and ssdna.flat().length should be the full size of the structure. For any missing piece, check murdered[].
        const helices = [];
        // guys for context murdered[] basically are the dumb nucleotides that couldnt be placed into helices due to fraying and angle conflicts.
        // #theyDeserveIt
        const murdered = [];
        if (!partials.length)
            return { helices, murdered }; // surely no helices if no partials.
        // quick lookup for id to partial index and deadfishies index.
        const idToPartial = new Map();
        partials.forEach((list, idx) => {
            list.forEach(nt => idToPartial.set(nt.id, idx));
        });
        const idToDeadfishy = new Map();
        deadfishies.forEach((nt, idx) => idToDeadfishy.set(nt.id, idx));
        console.log('ID to Partial Map:', idToPartial); // lets check out what the map looks like
        console.log('ID to Deadfishy Map:', idToDeadfishy);
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
        // deadfishies are just single nts, so caching their a3 vectors is simpler.
        const deadfishyA3 = new Map();
        const getDeadfishyA3 = (idx) => {
            let vec = deadfishyA3.get(idx);
            if (!vec) {
                vec = deadfishies[idx].getA3().clone().normalize();
                deadfishyA3.set(idx, vec);
            }
            return vec;
        };
        // https://www.youtube.com/watch?v=ayW5B2W9hfo
        // hopefully i remember to delete the youtube link before i commit or publish this lmao
        // is this O(logN) or O(N)??
        const totalNodes = partials.length + deadfishies.length;
        const parent = Array.from({ length: totalNodes }, (_, i) => i);
        const find = (x) => (parent[x] === x ? x : parent[x] = find(parent[x]));
        const unite = (a, b) => {
            const pa = find(a);
            const pb = find(b);
            if (pa !== pb)
                parent[pb] = pa;
        };
        // convert a nucleotide to its corresponding node reference (partial or deadfishy)...
        const getNodeRef = (nt) => {
            const partialId = idToPartial.get(nt.id);
            if (partialId !== undefined)
                return { node: partialId, kind: 'partial', index: partialId };
            const deadfishyId = idToDeadfishy.get(nt.id);
            if (deadfishyId !== undefined)
                return { node: partials.length + deadfishyId, kind: 'deadfishy', index: deadfishyId };
            return null;
        };
        // Track direct adjacency between partials and collect deadfishy links for a second pass.
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
        const deadfishyLinks = new Map();
        const addDeadfishyLink = (deadNode, partialIdx) => {
            const set = deadfishyLinks.get(deadNode) || new Set();
            set.add(partialIdx);
            deadfishyLinks.set(deadNode, set);
        };
        const canAttach = (a, b, strand) => {
            let vecA = null;
            let vecB = null;
            if (a.kind === 'partial')
                vecA = getPartialStrandA3(a.index, strand);
            else
                vecA = getDeadfishyA3(a.index);
            if (b.kind === 'partial')
                vecB = getPartialStrandA3(b.index, strand);
            else
                vecB = getDeadfishyA3(b.index);
            if (!vecA || !vecB)
                return false;
            return vecA.dot(vecB) > 0;
        };
        // Only connect nodes that are adjacent on the same strand and whose A3 vectors align.
        systems.forEach(system => {
            system.strands.forEach(strand => {
                let prev = null;
                strand.forEach(elem => {
                    const nt = elem;
                    if (prev) {
                        const nodeA = getNodeRef(prev);
                        const nodeB = getNodeRef(nt);
                        if (nodeA && nodeB && nodeA.node !== nodeB.node) {
                            // find adjacency between partials. Reason being, this will then be used for deadfishy connection checks later.
                            // without this, funny unintended behavior CAN happen.
                            // Example scenario: comment this part out and try running this code on Dumbbell Structure (nanobase 51). Helix 34/35 will be merged, unfortunately.
                            if (nodeA.kind === 'partial' && nodeB.kind === 'partial') {
                                addPartialAdj(nodeA.index, nodeB.index);
                            }
                            if (canAttach(nodeA, nodeB, strand)) {
                                if (nodeA.kind === 'partial' && nodeB.kind === 'partial') {
                                    // unionize the partials without question
                                    unite(nodeA.node, nodeB.node);
                                }
                                else if (nodeA.kind === 'deadfishy' || nodeB.kind === 'deadfishy') {
                                    // if either of the nodes are deadfishies, then checks are necessary.
                                    // add this to a "link". They will be processed in the 2nd pass. Slows down but much more accurate.
                                    const deadNode = nodeA.kind === 'deadfishy' ? nodeA : nodeB;
                                    const otherNode = nodeA.kind === 'deadfishy' ? nodeB : nodeA;
                                    // does not do anything if both nodes are deadfishies.
                                    if (otherNode.kind === 'partial') {
                                        addDeadfishyLink(deadNode.node, otherNode.index);
                                    }
                                }
                            }
                        }
                    }
                    prev = nt;
                });
            });
        });
        // Second pass: attach deadfishies after partial unions are finalized.
        const partialParent = Array.from({ length: partials.length }, (_, i) => find(i));
        const findPartial = (x) => (partialParent[x] === x ? x : partialParent[x] = findPartial(partialParent[x]));
        const partialMembers = new Map();
        for (let i = 0; i < partials.length; i++) {
            const root = findPartial(i);
            const set = partialMembers.get(root) || new Set();
            set.add(i);
            partialMembers.set(root, set);
        }
        // merging partials connected via a deadfishy.
        const mergePartialGroups = (a, b) => {
            let ra = findPartial(a);
            let rb = findPartial(b);
            if (ra === rb)
                return ra; // they are already united
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
        // if the 2 "helices" (partial groups) are to be merged, then they can NOT have connections amongst each other. 
        // if they do, then the one of the helix "turned around" to connect somewhere, and therefore should not be merged because it is now a different helix.
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
        deadfishyLinks.forEach((partialsSet, deadNode) => {
            const roots = new Set();
            partialsSet.forEach(pIdx => roots.add(findPartial(pIdx)));
            if (!roots.size)
                return;
            const rootsArr = Array.from(roots.values());
            let primaryRoot = rootsArr[0];
            unite(deadNode, primaryRoot);
            for (let i = 1; i < rootsArr.length; i++) {
                const otherRoot = rootsArr[i];
                if (findPartial(primaryRoot) === findPartial(otherRoot))
                    continue;
                if (hasDirectConnection(findPartial(primaryRoot), findPartial(otherRoot))) {
                    continue; // block deadfishy merge across already-connected helices
                }
                unite(deadNode, otherRoot);
                primaryRoot = mergePartialGroups(primaryRoot, otherRoot);
            }
        });
        const groups = new Map();
        partials.forEach((list, idx) => {
            const root = find(idx);
            const arr = groups.get(root) || [];
            arr.push(...list);
            groups.set(root, arr);
        });
        deadfishies.forEach((nt, idx) => {
            const node = partials.length + idx;
            const root = find(node);
            if (!groups.has(root)) {
                murdered.push(nt);
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
            // // Unnecessary. Might actually just remove this part later.
            // const mergeHelices = (targetIdx: number, sourceIdx: number) => {
            // 	if (targetIdx === sourceIdx) return;
            // 	const source = helices[sourceIdx];
            // 	if (!source.length) return;
            // 	addToHelix(targetIdx, source);
            // 	helices[sourceIdx] = [];
            // };
            // very unnecessarily, copilot has added logic so that if a ssScaffold segment connects to multiple helices, they get merged together.
            // is this even good? Does it cause problems instead of solving them? 
            // it is currently not causing issues AND it is not solving anything either, so whatever.
            // This part does the good stuff too, so dont just delete this lol.
            ssScaffold.forEach(segment => {
                if (!segment.length)
                    return;
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
                const targets = Array.from(helixIndices.values());
                if (!targets.length)
                    return;
                // IN CASE that the ssScaffold segment connects to multiple helices, warn the user.
                // This should NEVER happen, since the longssScaffold segments get divided into equal halves before being processed here...
                if (targets.length > 1) {
                    console.warn('ssScaffold segment connects to multiple helices; attaching to first.', {
                        helices: targets,
                        segmentLength: segment.length
                    });
                }
                const primary = targets[0];
                // for (let i = 1; i < targets.length; i++) {
                // 	mergeHelices(primary, targets[i]);
                // }
                addToHelix(primary, segment);
            });
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
                        firstHelixId
                    };
                }
                const oppositeNode = stepN(pair, dir, 3);
                if (!oppositeNode) {
                    return {
                        side: dir,
                        result: 'overhang',
                        firstPartialId: first.pIdx,
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
        // The lot of if statements are required (unless you can figure out a better way).
        // You can read through these, but they mostly comprise of cases where the segment is connected to helices on both ends, and has different types of such connections.
        // example, if overhang on one end and binder on the other, then it will connect to the helix on overhang side.
        ssdna.forEach(segment => {
            if (!segment.length)
                return;
            const { res5, res3 } = classifySegment(segment);
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
                binders.push(segment);
                return;
            }
            if (isBinder(res3) && isOverhang(res5) && !hasPartial(res5)) {
                binders.push(segment);
                return;
            }
            if (isBinder(res5) && isBinder(res3)) {
                binder2.push(segment);
                return;
            }
            unhandled.push(segment);
        });
        // For any binder/binder2 segments, group them by which helix they connect to.
        // If multiple binder segments connect to the same helix, they form a new helix.
        const resolveBinderHelix = (segment) => {
            const { res5, res3 } = classifySegment(segment);
            const helixIds = new Set();
            const collect = (res) => {
                if (!isBinder(res))
                    return;
                if (res.firstHelixId !== undefined)
                    helixIds.add(res.firstHelixId);
                if (res.oppositeHelixId !== undefined)
                    helixIds.add(res.oppositeHelixId);
            };
            if (res5)
                collect(res5);
            if (res3)
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
        binders.forEach(segment => {
            const helixId = resolveBinderHelix(segment);
            if (helixId === undefined)
                return;
            addBinderToGroup(helixId, segment);
        });
        binder2.forEach(segment => {
            const helixId = resolveBinderHelix(segment);
            if (helixId === undefined)
                return;
            addBinderToGroup(helixId, segment);
        });
        binderGroups.forEach(segments => {
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
        });
        // This is the code to re-attach murdered[] nts to the closest helix by partial connection (n5/n3).
        // murdered[] nts are the ones that got rejected from everything so far, but they still belong to the structure somehow...
        // "obviously", they belong to the closest helix!!! What could go wrong?
        if (murdered.length && helices.length) {
            const idToHelix = new Map();
            helices.forEach((list, idx) => {
                list.forEach(nt => idToHelix.set(nt.id, idx));
            });
            const walkToHelix = (start, dir) => {
                if (!start) {
                    console.warn('Tried to find the nearest connected helix, but start was', start);
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
                console.warn('Tried to find the nearest connected helix, but could not find any');
                console.log('Problematic nucleotide:', start);
                return null;
            };
            const remaining = [];
            murdered.forEach(nt => {
                const via5 = walkToHelix((nt.n5 ?? null), 'n5');
                const via3 = walkToHelix((nt.n3 ?? null), 'n3');
                // hopefully they are connected to the structure, otherwise its a relaxation issue (probably).
                if (!via5 && !via3) {
                    remaining.push(nt);
                    return;
                }
                let target = via5 || via3;
                // if the murdered nucleotide is connected to helices on both sides to 2 different helices, then 
                // it will pick the closest one.
                if (via5 && via3) {
                    const pos = nt.getPos();
                    const d5 = pos.distanceTo(via5.anchor.getPos());
                    const d3 = pos.distanceTo(via3.anchor.getPos());
                    target = d5 <= d3 ? via5 : via3;
                }
                // oops has happened
                if (!target) {
                    remaining.push(nt);
                    return;
                }
                const helix = helices[target.helixIdx];
                if (!helix) {
                    remaining.push(nt);
                    return;
                }
                helix.push(nt);
                idToHelix.set(nt.id, target.helixIdx);
            });
            // push the remaining nts back to murdered[].
            murdered.length = 0;
            murdered.push(...remaining);
        }
        // const finalHelices = helices.filter(h => h.length > 0);
        return { helices, murdered, binders, binder2, disconnected, unhandled };
    }
    honda.generateHelix = generateHelix;
    // export function generateTotal(inputMap: Map<number, Nucleotide>, tolerance = 2) {
    // 	let {partials, fishies} = findHelixPartials(inputMap, tolerance);
    // 	let {ssdna, deadfishies, longssScaffold} = ssdnaPartials(fishies);
    // 	let ssScaffold = longssScaffoldfunc(longssScaffold);
    // }
    async function findHelices(inputMap, tolerance = 2) {
        findBasepairsOptim2();
        dropIntraStrandPairs();
        // ok now we can do the rest of the stuff.
        let { partials, fishies } = findHelixPartials(inputMap, tolerance);
        let { ssdna, deadfishies, longssScaffold } = ssdnaPartials(fishies);
        let ssScaffold = longssScaffoldfunc(longssScaffold, deadfishies);
        let { helices, murdered, binders, binder2, disconnected, unhandled } = generateHelix(partials, ssdna, ssScaffold, deadfishies);
        console.log("Helices size:", helices.flat().length);
        console.log("Total elements:", inputMap.size);
        if (helices.flat().length !== inputMap.size) {
            console.log("Oops has occurred!");
            const helixIds = new Set();
            helices.forEach(list => list.forEach(nt => helixIds.add(nt.id)));
            const missing = [];
            inputMap.forEach((nt) => {
                if (!helixIds.has(nt.id))
                    missing.push(nt.id);
            });
            const listHits = new Map();
            const addHits = (label, list) => {
                const add = (nt) => {
                    if (!missingSet.has(nt.id))
                        return;
                    const set = listHits.get(nt.id) || new Set();
                    set.add(label);
                    listHits.set(nt.id, set);
                };
                if (Array.isArray(list[0])) {
                    list.forEach(segment => segment.forEach(add));
                }
                else {
                    list.forEach(add);
                }
            };
            const missingSet = new Set(missing);
            addHits('partials', partials);
            addHits('ssdna', ssdna);
            addHits('deadfishies', deadfishies);
            addHits('ssScaffold', ssScaffold);
            if (binders)
                addHits('binders', binders);
            if (binder2)
                addHits('binder2', binder2);
            if (disconnected)
                addHits('disconnected', disconnected);
            if (unhandled)
                addHits('unhandled', unhandled);
            addHits('murdered', murdered);
            const report = missing.map(id => ({
                id,
                lists: Array.from(listHits.get(id) || [])
            }));
            console.log('Missing nucleotide IDs:', missing);
            console.log('Missing membership report:', report);
        }
        return helices;
    }
    honda.findHelices = findHelices;
})(honda || (honda = {}));
