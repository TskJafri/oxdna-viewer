"use strict";
/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
/*
// for ease of use, and to prevent dumb mistakes as I code and test things out, here is exactly the commands to use this in console:
findBasepairs3(); // use 3 because there are 2 versions of findBasepairs, and 3 is the fastest. Dont ask why i named it that.
let {partials, unpaired} = airport.findHelixPartials(elements, 2);

let {ssdna, stubs, longssScaffold} = airport.ssdnaPartials(unpaired);
let ssScaffold = airport.longssScaffoldfunc(longssScaffold, stubs);
let {helices, lastScraps, binders, binder2, disconnected, unhandled} = airport.generateHelix(partials, ssdna, ssScaffold, stubs);
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
// let helices = await airport.findHelices(elements, 2);
var helix;
(function (helix_1) {
    // // helper function cuz didnt want to type this every time
    // export function checkAngle(n1: Nucleotide | null = null, n2: Nucleotide | null = null) {
    // 	if (!n1 || !n2) return 0;
    // 	return Math.acos(n1.getA3().dot(n2.getA3())) * (180 / Math.PI);
    // }
    // // Returns the longest strand in the system as scaffold...
    // export function getScaffoldStrand() {
    // 	let maxLen = 0;
    // 	let scaffold: Strand | null = null;
    // 	systems.forEach(s => {
    // 		s.strands.forEach(strand => {
    // 			if (strand.getLength() > maxLen) {
    // 				maxLen = strand.getLength();
    // 				scaffold = strand;
    // 			}
    // 		});
    // 	});
    // 	return scaffold;
    // }
    // Finds helix parts using destructive consumption of a working copy of elements (called elmts).
    // tolerance 2 is good enough for most cases. Higher tolerances seem to have no negative consequences, however.
    // TODO: Remove the strand parameter, maybe later? Check how this works on multi-strand structures. But this is specifically for RNA, so check strand is not needed.
    function findrnapartials(inputMap, tolerance = 2) {
        const elmts = new Map(inputMap);
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
            // console.log('Starting of partials:', start);
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
            let ally = pairInPool;
            const strandA = curr.strand; // required to ensure we don't cross strands. This is how we know partials is actually a partial helix.
            const strandB = ally.strand;
            // this is the partial being built.
            const partial = [];
            const seen = new Set();
            const tryDirection = (dir) => {
                const nextCurr = elmts.get(curr.id + dir);
                const nextAlly = elmts.get(ally.id - dir);
                if (!nextCurr || !nextAlly)
                    return null;
                if (nextCurr.strand !== strandA || nextAlly.strand !== strandB)
                    return null;
                // Tolerance window: OR logic. Continue if ANY offset in [1, tolerance]
                // has forward paired to backward (allows rescue by curr+2/a-2, curr+3/a-3, etc.).
                let onwards = false;
                for (let offset = 1; offset <= tolerance; offset++) {
                    const forward = elmts.get(curr.id + dir * offset);
                    const backward = elmts.get(ally.id - dir * offset);
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
                elmts.delete(curr.id);
                elmts.delete(ally.id);
                const step = tryDirection(1) || tryDirection(-1);
                if (!step)
                    break;
                // Always consume the immediate neighbors (curr+1 and a-1) even if mismatched.
                record(partial, seen, step.nextCurr);
                record(partial, seen, step.nextAlly);
                elmts.delete(step.nextCurr.id);
                elmts.delete(step.nextAlly.id);
                // Advance walker by one along each strand.
                curr = step.nextCurr;
                ally = step.nextAlly;
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
        // wait lately its being doing some weird crap, lets log it out.
        // had to add the multi-pairing logic (down below) to fix the problem. But still logging for curiosity.
        console.log('Duplicates set: ', duplicates);
        // const pairIds = new Set<number>();
        // partials.forEach((helix,i) => {
        // 	partials[i] = helix.filter(nt => {
        // 		if (duplicates.has(nt.id) || pairIds.has(nt.id)) {
        // 			unpaired.set(nt.id, nt);
        // 			return false;
        // 		}
        // 		if (nt.pair && duplicates.has(nt.pair.id)) {
        // 			unpaired.set(nt.id, nt);
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
            // // Drop empty partial lists
            partials = partials.filter(helix => helix.length > 0);
        }
        // Remove any nucleotide that is part of a multi-pairing (2+ nts paired to the same nt)
        // dang it. This doesnt work either, specificially at the ends of the helices, where the fraying causes multi-pairing.
        // const pairedTo = new Map<number, number[]>();
        // elmts2.forEach(nt => {
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
        // 				unpaired.set(nt.id, nt);
        // 				return false;
        // 			}
        // 			return true;
        // 		});
        // 	});
        // 	partials = partials.filter(helix => helix.length > 0);
        // }
        return { partials, unpaired: Array.from(unpaired.values()) };
    }
    helix_1.findrnapartials = findrnapartials;
    function averageA3(list) {
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
    helix_1.averageA3 = averageA3;
    ;
    function averageA1(list) {
        if (!list.length)
            return new THREE.Vector3(0, 0, 0);
        // Align all A1 vectors so they point in a consistent direction before averaging.
        const ref = list[0].getA1().clone().normalize();
        const acc = ref.clone();
        for (let i = 1; i < list.length; i++) {
            const v = list[i].getA1().clone().normalize();
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
    helix_1.averageA1 = averageA1;
    ;
    function rnaseparation(unpaired) {
        const not_sodead_unpaired = [];
        const stubs = [];
        if (!unpaired.length)
            return { not_sodead_unpaired, stubs };
        const fishSet = new Set(unpaired.map(nt => nt.id));
        const idToNt = new Map();
        unpaired.forEach(nt => idToNt.set(nt.id, nt));
        const visited = new Set();
        for (const seed of unpaired) {
            if (visited.has(seed.id))
                continue;
            const component = [];
            const stack = [seed];
            while (stack.length) {
                const curr = stack.pop();
                if (visited.has(curr.id))
                    continue;
                visited.add(curr.id);
                component.push(curr);
                const n5 = curr.n5;
                const n3 = curr.n3;
                if (n5 && fishSet.has(n5.id) && !visited.has(n5.id)) {
                    const next = idToNt.get(n5.id) || n5;
                    stack.push(next);
                }
                if (n3 && fishSet.has(n3.id) && !visited.has(n3.id)) {
                    const next = idToNt.get(n3.id) || n3;
                    stack.push(next);
                }
            }
            if (component.length >= 4) {
                not_sodead_unpaired.push(component);
            }
            else {
                component.forEach(nt => stubs.push(nt));
            }
        }
        return { not_sodead_unpaired, stubs };
    }
    helix_1.rnaseparation = rnaseparation;
    function rnaGenerateHelix(partials, stubs = []) {
        const helices = [];
        if (!partials.length)
            return { helices, eligibleEndpoints: [] };
        const idToPartial = new Map();
        partials.forEach((list, idx) => {
            list.forEach(nt => idToPartial.set(nt.id, idx));
        });
        const idToStub = new Map();
        stubs.forEach((nt, idx) => idToStub.set(nt.id, idx));
        const getNodeRef = (nt) => {
            const pIdx = idToPartial.get(nt.id);
            if (pIdx !== undefined)
                return { kind: 'partial', index: pIdx };
            const dIdx = idToStub.get(nt.id);
            if (dIdx !== undefined)
                return { kind: 'stub', index: dIdx };
            return null;
        };
        // endpointId -> neighboring partial indices found by prev/nt crossing across different partial nodes
        const endpointLinks = new Map();
        const addEndpointLink = (endpoint, neighborPartialIdx) => {
            const set = endpointLinks.get(endpoint.id) || new Set();
            set.add(neighborPartialIdx);
            endpointLinks.set(endpoint.id, set);
        };
        // stubIdx -> connected partial indices (A1-filtered adjacency only)
        const stubLinks = new Map();
        const addStubLink = (stubIdx, partialIdx) => {
            const set = stubLinks.get(stubIdx) || new Set();
            set.add(partialIdx);
            stubLinks.set(stubIdx, set);
        };
        // Endpoint detection method mirrored from dsdnahelix: walk each strand and inspect prev/nt crossings.
        systems.forEach(system => {
            system.strands.forEach(strand => {
                let prev = null;
                strand.forEach(elem => {
                    const nt = elem;
                    if (prev) {
                        const nodeA = getNodeRef(prev);
                        const nodeB = getNodeRef(nt);
                        if (nodeA && nodeB && (nodeA.kind !== nodeB.kind || nodeA.index !== nodeB.index)) {
                            const a1Prev = prev.getA1().clone().normalize();
                            const a1Nt = nt.getA1().clone().normalize();
                            if (a1Prev.dot(a1Nt) > 0.3) {
                                if (nodeA.kind === 'partial' && nodeB.kind === 'partial') {
                                    addEndpointLink(prev, nodeB.index);
                                    addEndpointLink(nt, nodeA.index);
                                }
                                else if (nodeA.kind === 'partial' && nodeB.kind === 'stub') {
                                    addStubLink(nodeB.index, nodeA.index);
                                }
                                else if (nodeA.kind === 'stub' && nodeB.kind === 'partial') {
                                    addStubLink(nodeA.index, nodeB.index);
                                }
                            }
                        }
                    }
                    prev = nt;
                });
            });
        });
        // Union-find over partial indices to actually join partials into helices.
        const parent = Array.from({ length: partials.length }, (_, i) => i);
        const find = (x) => (parent[x] === x ? x : (parent[x] = find(parent[x])));
        const unite = (a, b) => {
            const pa = find(a);
            const pb = find(b);
            if (pa !== pb)
                parent[pb] = pa;
        };
        const eligibleEndpoints = [];
        const eligibleEndpointIds = new Set();
        partials.forEach((partial, pIdx) => {
            if (!partial.length)
                return;
            const partialSet = new Set(partial.map(nt => nt.id));
            for (const endpoint of partial) {
                const linkedPartials = endpointLinks.get(endpoint.id);
                if (!linkedPartials || !linkedPartials.size)
                    continue;
                // Endpoint is eligible only if it has immediate n5/n3 neighbor(s) outside this partial.
                const outsideNeighbors = [];
                const n5 = endpoint.n5;
                const n3 = endpoint.n3;
                if (n5 && !partialSet.has(n5.id))
                    outsideNeighbors.push(n5);
                if (n3 && !partialSet.has(n3.id))
                    outsideNeighbors.push(n3);
                if (!outsideNeighbors.length)
                    continue;
                const endpointA1 = endpoint.getA1().clone().normalize();
                let acceptedAny = false;
                for (const neighbor of outsideNeighbors) {
                    const neighborPartialIdx = idToPartial.get(neighbor.id);
                    if (neighborPartialIdx === undefined || neighborPartialIdx === pIdx)
                        continue;
                    if (!linkedPartials.has(neighborPartialIdx))
                        continue;
                    const neighborA1 = neighbor.getA1().clone().normalize();
                    if (endpointA1.dot(neighborA1) <= 0.5)
                        continue;
                    unite(pIdx, neighborPartialIdx);
                    acceptedAny = true;
                }
                if (acceptedAny && !eligibleEndpointIds.has(endpoint.id)) {
                    eligibleEndpointIds.add(endpoint.id);
                    eligibleEndpoints.push(endpoint);
                }
            }
        });
        // Optional stub-mediated merges across linked partials (A1-filtered during link collection).
        stubLinks.forEach(linkedPartials => {
            const indices = Array.from(linkedPartials.values());
            if (indices.length < 2)
                return;
            const first = indices[0];
            for (let i = 1; i < indices.length; i++) {
                unite(first, indices[i]);
            }
        });
        // Build merged helices from connected partial components.
        const groups = new Map();
        partials.forEach((_, pIdx) => {
            const root = find(pIdx);
            const arr = groups.get(root) || [];
            arr.push(...partials[pIdx]);
            groups.set(root, arr);
        });
        // Attach stubs to their linked helix groups (A1-filtered during adjacency collection).
        stubs.forEach((nt, dIdx) => {
            const linked = stubLinks.get(dIdx);
            if (!linked || !linked.size)
                return;
            const firstPartial = Array.from(linked.values())[0];
            const root = find(firstPartial);
            const arr = groups.get(root) || [];
            arr.push(nt);
            groups.set(root, arr);
        });
        groups.forEach(group => {
            const seen = new Set();
            const unique = [];
            group.forEach(nt => {
                if (seen.has(nt.id))
                    return;
                seen.add(nt.id);
                unique.push(nt);
            });
            if (unique.length)
                helices.push(unique);
        });
        return { helices, eligibleEndpoints };
    }
    helix_1.rnaGenerateHelix = rnaGenerateHelix;
})(helix || (helix = {}));
