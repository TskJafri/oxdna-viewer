/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />

/* 
// for ease of use, and to prevent dumb mistakes as I code and test things out, here is exactly the commands to use this in console:
findBasepairs3(); // use 3 because there are 2 versions of findBasepairs, and 3 is the fastest. Dont ask why i named it that.
helix.dropIntraStrandPairs();
let {partials, unpaired} = helix.findHelixPartials(elements, 2);
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
// let helices = await helix.findHelices(elements, 2);

namespace helix {
	// helper function cuz didnt want to type this every time
	export function checkAngle(n1: Nucleotide | null = null, n2: Nucleotide | null = null) {
		if (!n1 || !n2) return 0;
		return Math.acos(n1.getA3().dot(n2.getA3())) * (180 / Math.PI);
	}

	// Returns the longest strand in the system as scaffold...
	export function getScaffoldStrand() {
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

	// Removes intra-strand pairings across the whole structure
	// ALWAYS run this function after findBasepairs3(). Very important.
	export function dropIntraStrandPairs() {
		elements.forEach((nt: any) => {
			if (!(nt instanceof Nucleotide)) return;
			if (!nt.pair) return;
			if (nt.pair.strand === nt.strand) {
				const mate = nt.pair as Nucleotide;
				nt.pair = null;
				if (mate.pair === nt) {
					mate.pair = null;
				}
			}
		});
	}

	// Finds helix parts using destructive consumption of a working copy of elements (called elmts).
	// tolerance 2 is good enough for most cases. Higher tolerances seem to have no negative consequences, however.
	// go to terminatingConditions for what tolerance is.
	export function findHelixPartials2(inputMap: Map<number, Nucleotide>, tolerance = 2) {
		const elmts = new Map<number, Nucleotide>(inputMap); // copy of the elements map
		const elmts2 = new Map<number, Nucleotide>(inputMap); // backup copy for duplication, used later
		const unpaired = new Map<number, Nucleotide>(); // unpaired / skipped nts

		let partials: Nucleotide[][] = [];
		const record = (list: Nucleotide[], set: Set<number>, nt: Nucleotide | null = null) => {
			if (nt && !set.has(nt.id)) {
				set.add(nt.id);
				list.push(nt);
			}
		};

		const nextStart = () => elmts.values().next().value;

		while (true) {
			const start = nextStart();
			if (!start) break;

			// Collect unpaired/binder nts for downstream ssDNA processing instead of discarding silently.
			// Additionally, only walk if the pair is also present in the current pool (elmts).
			const pairInPool = start.pair ? elmts.get(start.pair.id) as Nucleotide | undefined : undefined;
			if (!start.pair || !pairInPool) {
				unpaired.set(start.id, start as Nucleotide);
				elmts.delete(start.id);
				continue;
			}

			// initialize the walker.
			let curr: Nucleotide | null = start;
			let currPair: Nucleotide | null = pairInPool;
			const strandA = curr.strand; // required to ensure we don't cross strands. This is how we know partials is actually a partial helix.
			const strandB = currPair.strand;

			// this is the partial being built.
			const partial: Nucleotide[] = [];
			const seen = new Set<number>();

			type DirectionStep = { nextA: Nucleotide; nextB: Nucleotide } | null;

			const terminatingConditions = (nucA: Nucleotide, nucB: Nucleotide, dir: 1 | -1): DirectionStep => {
				const nucAc = elmts.get(nucA.id + dir) as Nucleotide | undefined;
				const nucBc = elmts.get(nucB.id - dir) as Nucleotide | undefined;
				if (!nucAc || !nucBc) return null;
				if (nucAc.strand !== strandA || nucBc.strand !== strandB) return null;

				// Tolerance window: OR logic. Continue if ANY offset in [1, tolerance] has forward paired to backward.
				// (allows rescue by walking farther along topological neighbors on both strands).
				let forwardCursor: Nucleotide | null = nucA;
				let backwardCursor: Nucleotide | null = nucB;
				for (let offset = 1; offset <= tolerance; offset++) {
					const forward = elmts.get(nucA.id + dir * offset) as Nucleotide | undefined;
					const backward = elmts.get(nucB.id - dir * offset) as Nucleotide | undefined;
					if (!forward || !backward) continue; // check if either of them even exist.
					if (forward.strand !== strandA || backward.strand !== strandB) continue;
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
				if (!step) break;

				console.log("Circular strand detected")
				// Circular helix guard: if the next step loops back into this partial, stop here.
				if (seen.has(step.nextA.id) || seen.has(step.nextB.id)) break;
				console.log("Circular strand stopped")

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
		const seenfordups = new Map<number, Nucleotide>();
		const duplicates = new Set<number>();
		partials.forEach(helix => {
			helix.forEach(nt => {
				if (seenfordups.has(nt.id)) {
					duplicates.add(nt.id);
				} else {
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
			const pairIds = new Set<number>();
			duplicates.forEach(id => {
				const a = elmts2.get(id) as Nucleotide | undefined;
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

	// Groups unpaired/binder nucleotides (unpaired) into ssDNA partials by strand.
	// Only contiguous runs (>2) along a strand are kept; shorter runs go to stubs.
	export function ssdnaPartials(unpaired: Nucleotide[]) {
		const unpairStrand = new Map<Strand, Nucleotide[]>();
		const ssdna: Nucleotide[][] = [];
		const stubs: Nucleotide[] = [];
		const longssScaffold: Nucleotide[] = [];
		const scaffold = getScaffoldStrand();

		unpaired.forEach(nt => {
			const arr = unpairStrand.get(nt.strand) || [];
			arr.push(nt);
			unpairStrand.set(nt.strand, arr);
		});

		// For each strand, sort it, find it's 5' ends and walk down n3 to build runs.
		unpairStrand.forEach((list, strand) => {
			const inSet = new Set(list.map(n => n.id));
			const visited = new Set<number>();
			const isScaffoldStrand = scaffold && strand === scaffold;

			list.sort((a, b) => a.id - b.id);

			for (const nt of list) {
				if (visited.has(nt.id)) continue;

				const isStart = !nt.n5 || !inSet.has((nt.n5 as Nucleotide).id);

				if (isStart) {
					const run: Nucleotide[] = [];
					let curr: Nucleotide | null = nt;

					while (curr && inSet.has(curr.id)) {
						run.push(curr);
						visited.add(curr.id);
						curr = curr.n3 as Nucleotide | null;
					}

					if (run.length > 2) {
						if (!isScaffoldStrand) {
							ssdna.push(run);
						} else {
							longssScaffold.push(...run);
						}
					} else {
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

	// helper function for adding the average a3 vector in canvas. Really should not be here.
	export function averageA3a(list: Nucleotide[]) {
		if (!list.length) return new THREE.Vector3(0, 0, 0);

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
		if (origin && typeof THREE !== 'undefined' && typeof scene !== 'undefined' && (scene as any)?.add) {
			const helper = new THREE.ArrowHelper(avg.clone(), origin, 5);
			(scene as any).add(helper);
		}

		return avg;
	};

	// helper to consolidate the nucleotides into contiguous segments, and enforce equal halves for scaffold segments.
	export function longssScaffoldfunc(longssScaffold: Nucleotide[], stubs: Nucleotide[] = []) {
		const ssScaffold: Nucleotide[][] = [];
		if (!longssScaffold.length) return ssScaffold;

		// Work on a sorted copy so numeric contiguity is easy to detect.
		const sorted = [...longssScaffold].sort((a, b) => a.id - b.id);
		let runAway: Nucleotide[] = [];

		const flushRun = () => {
			if (!runAway.length) return;
			if (runAway.length < 3) {
				runAway.forEach(nt => stubs.push(nt));
				runAway = [];
				return;
			}
			// Enforce equal halves; trim one nucleotide if odd-length to satisfy the requirement.
			const evenLen = runAway.length - (runAway.length % 2);
			if (evenLen !== runAway.length) {
				const dropped = runAway[evenLen];
				if (dropped) stubs.push(dropped);
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

	// let partialStrandMap = new Map<number, Map<number, Nucleotide[]>>();
	export function mapPartialEnds(partials: Nucleotide[][]) {
		const partialEndsMap = new Map<number, { start1: Nucleotide, end1: Nucleotide, start2: Nucleotide, end2: Nucleotide }>();

		partials.forEach((partial, index) => {
			const inSet = new Set(partial.map(n => n.id));
			
			// Find 3' ends (n3 is missing or outside the partial) and sort them by ID
			const ends3 = partial
				.filter(n => !n.n3 || !inSet.has((n.n3 as Nucleotide).id))
				.sort((a, b) => a.id - b.id);
				
			// Find 5' ends (n5 is missing or outside the partial)
			const ends5 = partial.filter(n => !n.n5 || !inSet.has((n.n5 as Nucleotide).id));

			if (ends3.length >= 2 && ends5.length >= 2) {
				// strand1 is assigned to the one with the lowest ID at the 3' end (ends3[0])
				const start1 = ends5.find(n => n.strand === ends3[0].strand)!;
				const start2 = ends5.find(n => n.strand === ends3[1].strand)!;

				// Skip 1-bp partials: when a strand has length 1 inside the partial,
				// start1 == end1 (or start2 == end2). There's no linear extent, so the
				// "two sides" model doesn't apply. These partials are merged by attach-count
				// cap in generateHelix instead of by side uniqueness.
				if (start1.id === ends3[0].id || start2.id === ends3[1].id) return;

				partialEndsMap.set(index, {
					start1: start1, end1: ends3[0],
					start2: start2, end2: ends3[1]
				});
			}
		});

		return partialEndsMap;
	}


	// Perfected!
	// this one uses average a3 vectors of CONNECTED strands, as opposed to average a3 vectors of the entire partial (which cancels out, due to topology).
	export function generateHelix(partials: Nucleotide[][], ssdna: Nucleotide[][], ssScaffold: Nucleotide[][], stubs: Nucleotide[]) {
		// Currently uses partials, ssScaffold and stubs to build perfect (almost) helices.
		// Hence, helices.flat().length and ssdna.flat().length should be the full size of the structure. For any missing piece, check lastScraps[].
		const helices: Nucleotide[][] = [];
		// guys for context lastScraps[] basically are the dumb nucleotides that couldnt be placed into helices due to fraying and angle conflicts.
		// Stored as segments so grouped leftovers (e.g. deferred ssScaffold segments) stay together.
		const lastScraps: Nucleotide[][] = [];
		if (!partials.length) return { helices, lastScraps }; // surely no helices if no partials.
		const dot = 0.5;

		// quick lookup for id to partial index and stubs index.
		const idToPartial = new Map<number, number>();
		partials.forEach((list, idx) => {
			list.forEach(nt => idToPartial.set(nt.id, idx));
		});
		const idTostubs = new Map<number, number>();
		stubs.forEach((nt, idx) => idTostubs.set(nt.id, idx));
		console.log('ID to Partial Map:', idToPartial); // lets check out what the map looks like
		console.log('ID to stubs Map:', idTostubs);

		const averageA3 = (list: Nucleotide[]) => {
			if (!list.length) return null;
			const acc = new THREE.Vector3(0, 0, 0);
			list.forEach(nt => {
				acc.add(nt.getA3().clone().normalize());
			});
			const len = acc.length();
			if (len < 1e-6) return null;
			return acc.divideScalar(len);
		};

		// within a partial, find the nts that belong to a specific strand within a specific partial. 
		// these will the ones used for finding the average a3 vector, which will later be used for connecting partials.
		const partialStrandMap = new Map<number, Map<number, Nucleotide[]>>();
		const getPartialStrandNts = (partialIdx: number, strand: Strand) => {
			let byStrand = partialStrandMap.get(partialIdx);
			if (!byStrand) {
				byStrand = new Map<number, Nucleotide[]>();
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
		const partialStrandA3 = new Map<number, Map<number, THREE.Vector3 | null>>();
		const getPartialStrandA3 = (partialIdx: number, strand: Strand) => {
			let byStrand = partialStrandA3.get(partialIdx);
			if (!byStrand) {
				byStrand = new Map<number, THREE.Vector3 | null>();
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
		const stubsA3 = new Map<number, THREE.Vector3>();
		const getstubsA3 = (idx: number) => {
			let vec = stubsA3.get(idx);
			if (!vec) {
				vec = stubs[idx].getA3().clone().normalize();
				stubsA3.set(idx, vec);
			}
			return vec;
		};

		const totalNodes = partials.length + stubs.length;
		const parent = Array.from({ length: totalNodes }, (_, i) => i);
		const find = (x: number): number => (parent[x] === x ? x : parent[x] = find(parent[x]));
		const unite = (a: number, b: number) => {
			const pa = find(a);
			const pb = find(b);
			if (pa !== pb) parent[pb] = pa;
		};

		type NodeRef = { node: number; kind: 'partial' | 'stubs'; index: number };
		// convert a nucleotide to its corresponding node reference (partial or stubs)...
		const getNodeRef = (nt: Nucleotide): NodeRef | null => {
			const partialId = idToPartial.get(nt.id);
			if (partialId !== undefined) return { node: partialId, kind: 'partial', index: partialId };
			const stubsId = idTostubs.get(nt.id);
			if (stubsId !== undefined) return { node: partials.length + stubsId, kind: 'stubs', index: stubsId };
			return null;
		};

		// Each partial has up to 2 sides (from mapPartialEnds). One side gets at most 1 connection to another partial.
		// 1-bp partials are excluded from partialEndsMap upstream — they get no sides here. Their
		// per-partial cap is enforced by attachCount (max 2 attachments) in the greedy below.
		// Build (partialIdx, ntId) -> sideIdx (0 or 1) so any exit-nt resolves to its side.
		const partialEndsMap = mapPartialEnds(partials);
		const ntToSide = new Map<number, Map<number, number>>();
		const sideCount = new Map<number, number>(); // partialIdx -> number of usable sides (0,1,2)
		partials.forEach((_, pIdx) => {
			const ends = partialEndsMap.get(pIdx);
			const inner = new Map<number, number>();
			ntToSide.set(pIdx, inner);
			if (!ends) {
				sideCount.set(pIdx, 0);
				return;
			}
			// Side 0: start1 (5' of strand A) paired with end2 (3' of strand B)
			inner.set(ends.start1.id, 0);
			inner.set(ends.end2.id, 0);
			// Side 1: end1 (3' of strand A) paired with start2 (5' of strand B)
			inner.set(ends.end1.id, 1);
			inner.set(ends.start2.id, 1);
			sideCount.set(pIdx, 2);
		});

		const getSideForNt = (pIdx: number, ntId: number): number | undefined => {
			return ntToSide.get(pIdx)?.get(ntId);
		};

		// Track direct adjacency between partials (used for stub-bridge safety check).
		const partialAdj = new Map<number, Set<number>>();
		const addPartialAdj = (a: number, b: number) => {
			if (a === b) return;
			const setA = partialAdj.get(a) || new Set<number>();
			setA.add(b);
			partialAdj.set(a, setA);
			const setB = partialAdj.get(b) || new Set<number>();
			setB.add(a);
			partialAdj.set(b, setB);
		};

		// Direct partial<->partial edges with side info.
		// Each strand crossing between two partials produces one edge. By construction, two paired
		// nucleotides on the same side cannot both cross into the same partial (they would have
		// been part of that partial), so per-strand crossings won't collide on the same side-pair.
		type DirectEdge = { a: number; sideA: number; b: number; sideB: number; dots: number };
		const directEdges: DirectEdge[] = [];
		const addDirectEdge = (a: number, sideA: number, b: number, sideB: number, dots: number) => {
			if (a === b) return;
			directEdges.push({ a, sideA, b, sideB, dots });
		};

		// Stub -> partial links. Track which side of the partial the stub connects to,
		// so stub-bridge edges between two partials know which sides they would consume.
		type stubsLink = {
			partialIdx: number;
			dots: number;
			strand: Strand;
			partialSide: number;
		};
		const stubsLinks = new Map<number, Map<number, stubsLink>>();
		const addstubsLink = (stubNode: number, partialIdx: number, dots: number, strand: Strand, partialSide: number) => {
			const links = stubsLinks.get(stubNode) || new Map<number, stubsLink>();
			const prev = links.get(partialIdx);
			if (!prev || dots > prev.dots) {
				links.set(partialIdx, { partialIdx, dots, strand, partialSide });
			}
			stubsLinks.set(stubNode, links);
		};

		const attachDot = (a: NodeRef, b: NodeRef, strand: Strand) => {
			let vecA: THREE.Vector3 | null = null;
			let vecB: THREE.Vector3 | null = null;

			if (a.kind === 'partial') vecA = getPartialStrandA3(a.index, strand);
			else vecA = getstubsA3(a.index);

			if (b.kind === 'partial') vecB = getPartialStrandA3(b.index, strand);
			else vecB = getstubsA3(b.index);

			if (!vecA || !vecB) return -1;
			return vecA.dot(vecB);
		};

		// First pass. Collects all data for connections (partial-partial or partial-stubs)
		systems.forEach(system => {
			system.strands.forEach(strand => {
				let prev: Nucleotide | null = null;
				strand.forEach(elem => {
					const nt = elem as Nucleotide;
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
									addDirectEdge(nodeA.index, sideA!, nodeB.index, sideB!, d);
								} else if (nodeA.kind === 'stubs' || nodeB.kind === 'stubs') {
									const stubNode = nodeA.kind === 'stubs' ? nodeA : nodeB;
									const otherNode = nodeA.kind === 'stubs' ? nodeB : nodeA;
									const otherNt = nodeA.kind === 'stubs' ? nt : prev;
									if (otherNode.kind === 'partial') {
										const partialSide = getSideForNt(otherNode.index, otherNt.id);
										addstubsLink(stubNode.node, otherNode.index, d, strand, partialSide!);
									}
								}
							}
						}
					}
					prev = nt;
				});
			});
		});

		// Track which (partialIdx, sideIdx) slots are already used by a partial-partial connection.
		const consumedSides = new Set<string>();
		const consumedKey = (pIdx: number, side: number) => `${pIdx}:${side}`;
		const isSideConsumed = (pIdx: number, side: number) => consumedSides.has(consumedKey(pIdx, side));
		const consumeSide = (pIdx: number, side: number) => consumedSides.add(consumedKey(pIdx, side));

		// Per-partial attachment count, capped at 2. For multi-bp partials this is already
		// implicitly enforced by side uniqueness (2 sides, each consumable once = max 2 attachments).
		// For 1-bp partials (no sides), this is the only cap — they can host up to 2 neighbors,
		// chosen greedily by highest dot.
		const PER_PARTIAL_CAP = 2;
		const attachCount = new Map<number, number>();
		const getAttach = (pIdx: number) => attachCount.get(pIdx) ?? 0;
		const slotAvailable = (pIdx: number, side: number | undefined): boolean => {
			if (getAttach(pIdx) >= PER_PARTIAL_CAP) return false;
			if (side === undefined) return true;
			return !isSideConsumed(pIdx, side);
		};
		const reserveSlot = (pIdx: number, side: number | undefined) => {
			attachCount.set(pIdx, getAttach(pIdx) + 1);
			if (side !== undefined) consumeSide(pIdx, side);
		};

		/* 
		Mutual-agreement filter:
		For each (partial P, side σ_P), we look at every direct edge incident to that side and
		keep only edges (P, σ_P) <-> (Q, σ_Q) where Q's side σ_Q would also rank P among its top
		candidates. Since we already cap at "1 per side", a side has effectively top-1 candidates,
		so mutual agreement reduces to: σ_Q's best candidate (by dots) toward this junction is P.
		We approximate this with a sort+greedy pass below: among all surviving direct edges for a
		side, only the highest-dot edge can ever win. If the other side's highest-dot edge also
		names the same counterpart, both will agree naturally during the greedy pass.
		*/

		// Union-find helpers needed before the greedy pass for partial-group operations.
		const partialParent = Array.from({ length: partials.length }, (_, i) => i);
		const findPartial = (x: number): number => (partialParent[x] === x ? x : partialParent[x] = findPartial(partialParent[x]));
		const partialMembers = new Map<number, Set<number>>();
		for (let i = 0; i < partials.length; i++) {
			const set = partialMembers.get(i) || new Set<number>();
			set.add(i);
			partialMembers.set(i, set);
		}
		const mergePartialGroups = (a: number, b: number) => {
			let ra = findPartial(a);
			let rb = findPartial(b);
			if (ra === rb) return ra;
			const setA = partialMembers.get(ra)!;
			const setB = partialMembers.get(rb)!;
			if (setA.size < setB.size) {
				const tmp = ra;
				ra = rb;
				rb = tmp;
			}
			const keep = partialMembers.get(ra)!;
			const drop = partialMembers.get(rb)!;
			drop.forEach(idx => keep.add(idx));
			partialMembers.set(ra, keep);
			partialMembers.delete(rb);
			partialParent[rb] = ra;
			return ra;
		};
		const hasDirectConnection = (rootA: number, rootB: number) => {
			if (rootA === rootB) return true;
			const setA = partialMembers.get(rootA);
			const setB = partialMembers.get(rootB);
			if (!setA || !setB) return false;
			const [small, large] = setA.size <= setB.size ? [setA, setB] : [setB, setA];
			for (const idx of small) {
				const neighbors = partialAdj.get(idx);
				if (!neighbors) continue;
				for (const n of neighbors) {
					if (large.has(n)) return true;
				}
			}
			return false;
		};

		// Build stub-bridge edges. For each stub with >= 2 partial links, generate pairwise edges
		// between its candidate partials. The dot used is the partial<->partial alignment
		type StubEdge = {
			a: number; sideA: number;
			b: number; sideB: number;
			dots: number;
			stubNode: number;
		};
		const stubEdges: StubEdge[] = [];
		const partialsBridgeDot = (a: stubsLink, b: stubsLink) => {
			const vecA = getPartialStrandA3(a.partialIdx, a.strand);
			const vecB = getPartialStrandA3(b.partialIdx, b.strand);
			if (!vecA || !vecB) return -1;
			return vecA.dot(vecB);
		};
		stubsLinks.forEach((linksByPartial, stubNode) => {
			const candidates = Array.from(linksByPartial.values());
			if (candidates.length < 2) return;
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
			if (!slotAvailable(c.a, c.sideA)) continue;
			if (!slotAvailable(c.b, c.sideB)) continue;
			const rootA = findPartial(c.a);
			const rootB = findPartial(c.b);
			if (rootA === rootB) continue;
			unite(c.a, c.b);
			mergePartialGroups(c.a, c.b);
			reserveSlot(c.a, c.sideA);
			reserveSlot(c.b, c.sideB);
		}

		// Stage B: stub bridges. Same checks plus the existing safeguard against bridging
		// two groups that already have a direct partial-partial connection.
		for (const c of sortedStub) {
			if (!slotAvailable(c.a, c.sideA)) continue;
			if (!slotAvailable(c.b, c.sideB)) continue;
			const rootA = findPartial(c.a);
			const rootB = findPartial(c.b);
			if (rootA === rootB) continue;
			if (hasDirectConnection(rootA, rootB)) continue;
			unite(c.stubNode, c.a);
			unite(c.stubNode, c.b);
			mergePartialGroups(c.a, c.b);
			reserveSlot(c.a, c.sideA);
			reserveSlot(c.b, c.sideB);
		}

		// Stubs join the best-aligned partial through A3 dots.
		stubsLinks.forEach((linksByPartial, stubNode) => {
			const ranked = Array.from(linksByPartial.values()).sort((x, y) => y.dots - x.dots);
			if (!ranked.length) return;
			const primary = ranked[0];
			unite(stubNode, primary.partialIdx);
		});

		const groups = new Map<number, Nucleotide[]>();
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
			const seen = new Set<number>();
			const unique: Nucleotide[] = [];
			group.forEach(nt => {
				if (seen.has(nt.id)) {
					console.log('Duplicate nucleotide found in helix grouping:', nt);
					return;
				}
				seen.add(nt.id);
				unique.push(nt);
			});
			if (unique.length) helices.push(unique);
		});

		// After helices are built, attach ssScaffold segments to the helix they connect to.
		if (ssScaffold && ssScaffold.length) {
			const idToHelix = new Map<number, number>();
			helices.forEach((list, idx) => {
				list.forEach(nt => idToHelix.set(nt.id, idx));
			});

			const ssScaffoldIds = new Set<number>();
			ssScaffold.forEach(segment => segment.forEach(nt => ssScaffoldIds.add(nt.id)));

			// fairly obvious. Adds the segment of nucleotides (from ssScaffold) to the target helix index.
			const addToHelix = (targetIdx: number, segment: Nucleotide[]) => {
				const helix = helices[targetIdx];
				const seen = new Set<number>(helix.map(nt => nt.id));
				segment.forEach(nt => {
					if (seen.has(nt.id)) return;
					helix.push(nt);
					seen.add(nt.id);
					idToHelix.set(nt.id, targetIdx);
				});
			};

			const findSsScaffoldTargets = (segment: Nucleotide[]) => {
				const segmentIds = new Set<number>(segment.map(nt => nt.id));
				const helixIndices = new Set<number>();

				segment.forEach(nt => {
					const n5 = nt.n5 as Nucleotide | null;
					const n3 = nt.n3 as Nucleotide | null;

					if (n5 && !segmentIds.has(n5.id) && !ssScaffoldIds.has(n5.id)) {
						const hIdx = idToHelix.get(n5.id);
						if (hIdx !== undefined) helixIndices.add(hIdx);
					}

					if (n3 && !segmentIds.has(n3.id) && !ssScaffoldIds.has(n3.id)) {
						const hIdx = idToHelix.get(n3.id);
						if (hIdx !== undefined) helixIndices.add(hIdx);
					}
				});

				return Array.from(helixIndices.values());
			};

			let pending = ssScaffold.filter(segment => segment.length > 0);
			const maxRounds = Math.max(1, pending.length * 2);
			let round = 0;

			while (pending.length) {
				round += 1;
				let attachedThisRound = 0;
				const nextPending: Nucleotide[][] = [];

				pending.forEach(segment => {
					const targets = findSsScaffoldTargets(segment);
					if (!targets.length) {
						nextPending.push(segment);
						return;
					}

					// IN CASE that the ssScaffold segment connects to multiple helices, warn the user.
					if (targets.length > 1) {
						console.warn('ssScaffold segment connects to multiple helices; attaching to first.', {
							helices: targets,
							segmentLength: segment.length,
							round
						});
					}
					const primary = targets[0];
					addToHelix(primary, segment);
					attachedThisRound += 1;
				});

				if (!nextPending.length) break;

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

		const partialToHelix = new Map<number, number>();
		helices.forEach((list, hIdx) => {
			list.forEach(nt => {
				const pIdx = idToPartial.get(nt.id);
				if (pIdx !== undefined && !partialToHelix.has(pIdx)) {
					partialToHelix.set(pIdx, hIdx);
				}
			});
		});

		const addSegmentToHelix = (targetIdx: number, segment: Nucleotide[]) => {
			const helix = helices[targetIdx];
			if (!helix) return;
			const seen = new Set<number>(helix.map(nt => nt.id));
			segment.forEach(nt => {
				if (seen.has(nt.id)) return;
				helix.push(nt);
				seen.add(nt.id);
			});
		};

		// const oppositeDir = (dir: 'n5' | 'n3') => (dir === 'n5' ? 'n3' : 'n5');

		const findEndOnSide = (segment: Nucleotide[], segmentSet: Set<number>, dir: 'n5' | 'n3') => {
			for (const nt of segment) {
				const neighbor = nt[dir] as Nucleotide | null;
				if (!neighbor || !segmentSet.has(neighbor.id)) return nt;
			}
			return null;
		};

		const walkForPartial = (
			start: Nucleotide | null,
			dir: 'n5' | 'n3',
			segmentSet: Set<number>
		) => {
			let curr = start;
			while (curr) {
				if (segmentSet.has(curr.id)) return undefined;
				const pIdx = idToPartial.get(curr.id);
				if (pIdx !== undefined) return { pIdx, node: curr };
				curr = curr[dir] as Nucleotide | null;
			}
			return undefined;
		};

		const stepWithinSamePartial = (
			start: Nucleotide,
			dir: 'n5' | 'n3',
			steps: number,
			pIdx: number
		) => {
			let curr: Nucleotide | null = start;
			let last: Nucleotide = start;
			for (let i = 0; i < steps; i++) {
				curr = curr?.[dir] as Nucleotide | null;
				if (!curr) return last;
				const idx = idToPartial.get(curr.id);
				if (idx !== pIdx) return last;
				last = curr;
			}
			return last;
		};

		const stepN = (start: Nucleotide, dir: 'n5' | 'n3', steps: number) => {
			let curr: Nucleotide | null = start;
			for (let i = 0; i < steps; i++) {
				curr = curr?.[dir] as Nucleotide | null;
				if (!curr) return undefined;
			}
			return curr;
		};

		type SideResult = {
			side: 'n5' | 'n3';
			result: 'binder' | 'overhang';
			firstPartialId?: number;
			oppositePartialId?: number;
			firstHelixId?: number;
			oppositeHelixId?: number;
		};

		const classifySegment = (segment: Nucleotide[]) => {
			const segmentSet = new Set<number>(segment.map(nt => nt.id));

			const analyzeSide = (dir: 'n5' | 'n3'): SideResult => {
				const end = findEndOnSide(segment, segmentSet, dir);
				if (!end) return { side: dir, result: 'overhang' };

				const anchor = end[dir] as Nucleotide | null;
				if (!anchor) return { side: dir, result: 'overhang' };

				const first = walkForPartial(anchor, dir, segmentSet);
				if (!first) return { side: dir, result: 'overhang' };

				const firstHelixId = partialToHelix.get(first.pIdx);
				// note: using more than 1 step might look fine, but it can cause issues in edge cases.
				// specifically, structure 51 from nanobase (Dumbbell structure) has issues with this.
				// Sticking to 1 step has NOT shown ANY problems so far.
				const lastInPartial = stepWithinSamePartial(first.node, dir, 1, first.pIdx);
				const pair = lastInPartial.pair as Nucleotide | null;
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

				const binder =
					oppositePartialId !== undefined &&
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

		const binders: Nucleotide[][] = [];
		const binder2: Nucleotide[][] = [];
		const disconnected: Nucleotide[][] = [];
		const unhandled: Nucleotide[][] = [];

		const isBinder = (res: SideResult | undefined) => res?.result === 'binder';
		const isOverhang = (res: SideResult | undefined) => res?.result === 'overhang';
		const hasPartial = (res: SideResult | undefined) => res?.firstPartialId !== undefined;

		// The lot of if statements are required (unless you can figure out a better way).
		// You can read through these, but they mostly comprise of cases where the segment is connected to helices on both ends, and has different types of such connections.
		// example, if overhang on one end and binder on the other, then it will connect to the helix on overhang side.
		ssdna.forEach(segment => {
			if (!segment.length) return;
			const { res5, res3 } = classifySegment(segment);
			const res5HasPartial = hasPartial(res5);
			const res3HasPartial = hasPartial(res3);

			if (!res5HasPartial && !res3HasPartial) {
				disconnected.push(segment);
				return;
			}

			if (isOverhang(res5) && hasPartial(res5) && isOverhang(res3) && !hasPartial(res3)) {
				if (res5.firstHelixId !== undefined) addSegmentToHelix(res5.firstHelixId, segment);
				return;
			}
			if (isOverhang(res3) && hasPartial(res3) && isOverhang(res5) && !hasPartial(res5)) {
				if (res3.firstHelixId !== undefined) addSegmentToHelix(res3.firstHelixId, segment);
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
				if (res5.firstHelixId !== undefined) addSegmentToHelix(res5.firstHelixId, segment);
				return;
			}
			if (isOverhang(res3) && hasPartial(res3) && isBinder(res5)) {
				if (res3.firstHelixId !== undefined) addSegmentToHelix(res3.firstHelixId, segment);
				return;
			}

			if (isBinder(res5) && !res5HasPartial && isOverhang(res3) && res3HasPartial) {
				if (res3.firstHelixId !== undefined) addSegmentToHelix(res3.firstHelixId, segment);
				return;
			}
			if (isBinder(res3) && !res3HasPartial && isOverhang(res5) && res5HasPartial) {
				if (res5.firstHelixId !== undefined) addSegmentToHelix(res5.firstHelixId, segment);
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
		const resolveBinderHelix = (segment: Nucleotide[]) => {
			const { res5, res3 } = classifySegment(segment);
			const helixIds = new Set<number>();
			const collect = (res: SideResult) => {
				if (!isBinder(res)) return;
				if (res.firstHelixId !== undefined) helixIds.add(res.firstHelixId);
				if (res.oppositeHelixId !== undefined) helixIds.add(res.oppositeHelixId);
			};
			if (res5) collect(res5);
			if (res3) collect(res3);
			if (helixIds.size === 1) return Array.from(helixIds.values())[0];
			return undefined;
		};

		const binderGroups = new Map<number, Nucleotide[][]>();
		const addBinderToGroup = (helixId: number, segment: Nucleotide[]) => {
			const list = binderGroups.get(helixId) || [];
			list.push(segment);
			binderGroups.set(helixId, list);
		};

		// For binder2 segments that span two distinct helices, group by the unordered helix pair.
		// All binder2 segments connecting the SAME two helices get merged into a single new helix.
		// Binder2 segments connecting a DIFFERENT pair get their own new helix.
		const pairKey = (a: number, b: number) => (a < b ? `${a}:${b}` : `${b}:${a}`);
		const binderPairGroups = new Map<string, Nucleotide[][]>();
		const addBinderToPairGroup = (a: number, b: number, segment: Nucleotide[]) => {
			const key = pairKey(a, b);
			const list = binderPairGroups.get(key) || [];
			list.push(segment);
			binderPairGroups.set(key, list);
		};

		// Collect helix ids touched by binder-classified sides of a segment.
		// For binder2 the result has size 1 (both sides agree) or 2 (sides resolve to different helices).
		const getBinderHelixIds = (segment: Nucleotide[]): number[] => {
			const { res5, res3 } = classifySegment(segment);
			const ids = new Set<number>();
			const collect = (res: SideResult | undefined) => {
				if (!res || !isBinder(res)) return;
				if (res.firstHelixId !== undefined) ids.add(res.firstHelixId);
				if (res.oppositeHelixId !== undefined) ids.add(res.oppositeHelixId);
			};
			collect(res5);
			collect(res3);
			return Array.from(ids.values());
		};

		binders.forEach(segment => {
			const helixId = resolveBinderHelix(segment);
			if (helixId === undefined) return;
			addBinderToGroup(helixId, segment);
		});

		binder2.forEach(segment => {
			const ids = getBinderHelixIds(segment);
			if (ids.length === 1) {
				// Single host helix: same path as a normal binder.
				addBinderToGroup(ids[0], segment);
			} else if (ids.length === 2) {
				// Two host helices: group by the unordered pair so all binder2 segments
				// spanning the same {A, B} pair fuse into one new helix together.
				addBinderToPairGroup(ids[0], ids[1], segment);
			}
			// ids.length === 0 or > 2: should be unreachable for binder2; silently dropped.
		});

		const materializeBinderHelix = (segments: Nucleotide[][]) => {
			if (!segments.length) return;
			const seen = new Set<number>();
			const newHelix: Nucleotide[] = [];
			segments.forEach(segment => {
				segment.forEach(nt => {
					if (seen.has(nt.id)) return;
					seen.add(nt.id);
					newHelix.push(nt);
				});
			});
			if (newHelix.length) helices.push(newHelix);
		};

		binderGroups.forEach(segments => materializeBinderHelix(segments));
		binderPairGroups.forEach(segments => materializeBinderHelix(segments));

		// This is the code to re-attach lastScraps[] segments to the closest helix by partial connection (n5/n3).
		// For grouped segments (e.g. deferred ssScaffold), we attach the full segment to one chosen helix.
		if (lastScraps.length && helices.length) {
			const idToHelix = new Map<number, number>();
			helices.forEach((list, idx) => {
				list.forEach(nt => idToHelix.set(nt.id, idx));
			});

			type WalkHit = { helixIdx: number; anchor: Nucleotide };

			const walkToHelix = (start: Nucleotide | null, dir: 'n5' | 'n3', owner: Nucleotide): WalkHit | null => {
				if (!start) {
					console.warn('[walkToHelix] Side', dir, 'is null for nucleotide', owner.id, '; using opposite side if available.');
					return null;
				}
				let curr: Nucleotide | null = start;
				while (curr) {
					// which helix does this current nt belong to?
					const hIdx = idToHelix.get(curr.id);
					if (hIdx !== undefined) return { helixIdx: hIdx, anchor: curr };
					curr = curr[dir] as Nucleotide | null;
				}
				console.warn('[walkToHelix] Side', dir, 'for nucleotide', owner.id, 'started at', start.id, 'but did not reach any existing helix.');
				return null;
			};

			const addlastScrapsSegmentToHelix = (targetIdx: number, segment: Nucleotide[]) => {
				const helix = helices[targetIdx];
				if (!helix) return false;
				const seen = new Set<number>(helix.map(nt => nt.id));
				segment.forEach(nt => {
					if (seen.has(nt.id)) return;
					helix.push(nt);
					seen.add(nt.id);
					idToHelix.set(nt.id, targetIdx);
				});
				return true;
			};

			const pickSegmentTarget = (segment: Nucleotide[]): WalkHit | null => {
				let best: { hit: WalkHit; distance: number } | null = null;

				segment.forEach(nt => {
					const via5 = walkToHelix((nt.n5 ?? null) as Nucleotide | null, 'n5', nt);
					const via3 = walkToHelix((nt.n3 ?? null) as Nucleotide | null, 'n3', nt);

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
					const consider = (hit: WalkHit) => {
						const distance = pos.distanceTo(hit.anchor.getPos());
						if (!best || distance < best.distance) {
							best = { hit, distance };
						}
					};

					if (via5) consider(via5);
					if (via3) consider(via3);
				});

				return best ? best.hit : null;
			};

			const remaining: Nucleotide[][] = [];
			lastScraps.forEach(segment => {
				if (!segment.length) return;
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
		return { helices, lastScraps, binders, binder2, disconnected, unhandled };
	}

	export function findHelices(inputMap: Map<number, Nucleotide>, tolerance = 2) {
		findBasepairsOptim2();
		dropIntraStrandPairs();
		// ok now we can do the rest of the stuff.
		let { partials, unpaired } = findHelixPartials2(inputMap, tolerance);
		let { ssdna, stubs, longssScaffold } = ssdnaPartials(unpaired);
		let ssScaffold = longssScaffoldfunc(longssScaffold, stubs);
		let { helices, lastScraps, binders, binder2, disconnected, unhandled } = generateHelix(partials, ssdna, ssScaffold, stubs);
		console.log("Helices size:", helices.flat().length);
		console.log("Total elements:", inputMap.size);
		return { helices };
	}

	// Merge two or more helices into the one with the lowest index.
	// Mutates `helices` in place: pushes nucleotides from higher-indexed entries into the kept helix
	// and SPLICES those entries out, so the array shrinks. Also mutates `grid` in place: every
	// GridMark whose helixId points at a merged-away helix is rewritten to the kept id, and every
	// surviving helixId is shifted down through the same remap that's returned. Returns an
	// `idRemap` (oldIdx -> newIdx) that callers still need for any external references that key
	// off helix index — gridview node ids, crossover connection ids, etc.
	// No checks for grid layout yet
	export function combineHelices(
		helices: Nucleotide[][],
		indices: number[],
		grid: toscad.GridMap
	): { keptIdx: number; mergedIdx: number[]; idRemap: Map<number, number> } | null {
		if (!Array.isArray(helices) || !Array.isArray(indices) || !(grid instanceof Map)) return null;

		const valid: number[] = [];
		const seenIdx = new Set<number>();
		indices.forEach(i => {
			if (i < 0 || i >= helices.length) return;
			if (!Array.isArray(helices[i]) || helices[i].length === 0) return;
			if (seenIdx.has(i)) return;
			seenIdx.add(i);
			valid.push(i);
		});

		if (valid.length < 2) return null;

		// merge by keeping the lowest id only
		valid.sort((a, b) => a - b);
		const keptIdxOld = valid[0];
		const mergedIdxOld = valid.slice(1);

		// Move nucleotides into the kept helix, deduping by id
		const seenNts = new Set<number>(helices[keptIdxOld].map(nt => nt.id));
		mergedIdxOld.forEach(idx => {
			helices[idx].forEach(nt => {
				if (seenNts.has(nt.id)) return;
				seenNts.add(nt.id);
				helices[keptIdxOld].push(nt);
			});
		});

		// Build an oldIdx -> newIdx remap for every helix that survives the splice.
		const removed = new Set<number>(mergedIdxOld);
		const idRemap = new Map<number, number>();
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
			if (next !== undefined) mark.helixId = next;
		});

		// Splice in reverse so earlier indices stay valid during removal.
		for (let i = helices.length - 1; i >= 0; i--) {
			if (removed.has(i)) helices.splice(i, 1);
		}

		// keptIdx is the lowest valid index, so nothing in front of it was removed: its new index
		// is the same as its old one. Look it up via idRemap to stay correct if this invariant ever changes.
		const keptIdx = idRemap.get(keptIdxOld) ?? keptIdxOld;

		return { keptIdx, mergedIdx: mergedIdxOld, idRemap };
	}

	// make a snapshot type for ease of undo/redo functionality.
	export type CombineSnapshot = {
		kept: number;                                        // pre-merge index of the kept helix
		idRemap: Array<[number, number]>;                    // [oldIdx, newIdx] for survivors
		removed: Array<{ oldIdx: number; ntIds: number[] }>; // one per merged-away helix
	};

	// Reverse a combineHelices call using a snapshot recorded at merge time.
	// Mutates `helices` in place: rebuilds at full pre-merge length, places each survivor back at
	// its old index, and repopulates each removed slot from its captured `ntIds`. Mutates `grid`
	// in place: marks for nucleotides in any removed slot get their original `helixId` back;
	// every other mark shifts up through the inverse remap.
	// Returns an `inverseRemap(currentId) -> oldId` helper so callers can fix up external state
	// (gridview node ids, connection endpoints) that's still keyed in post-merge numbering.
	export function splitHelices(
		helices: Nucleotide[][],
		grid: toscad.GridMap,
		snapshot: CombineSnapshot
	): { inverseRemap: (currentId: number) => number } | null {
		if (!Array.isArray(helices) || !(grid instanceof Map) || !snapshot) return null;
		if (!Array.isArray(snapshot.idRemap) || !Array.isArray(snapshot.removed)) return null;

		// Inverse remap: post-merge (current) index -> pre-merge (old) index. Only survivors are
		// in this map; ids outside it are passed through unchanged.
		const survivorInverse = new Map<number, number>();
		snapshot.idRemap.forEach(([oldIdx, newIdx]) => survivorInverse.set(newIdx, oldIdx));
		const inverseRemap = (currentId: number): number => {
			const oldIdx = survivorInverse.get(currentId);
			return oldIdx !== undefined ? oldIdx : currentId;
		};

		// Locate the kept helix in its current (post-merge) position.
		const keptOldIdx = snapshot.kept;
		let keptCurrentIdx = keptOldIdx;
		for (const [newIdx, oldIdx] of survivorInverse.entries()) {
			if (oldIdx === keptOldIdx) { keptCurrentIdx = newIdx; break; }
		}

		// Pull merged-away nucleotides out of the kept helix.
		const ntsToExtract = new Set<number>();
		snapshot.removed.forEach(slot => slot.ntIds.forEach(id => ntsToExtract.add(id)));
		const keptArr = helices[keptCurrentIdx] || [];
		const keptKept: Nucleotide[] = [];
		keptArr.forEach(nt => { if (!ntsToExtract.has(nt.id)) keptKept.push(nt); });

		// Rebuild helices at full pre-merge length: survivors back at their old indices, removed
		// slots repopulated by looking up nucleotides in the global elements map.
		const total = helices.length + snapshot.removed.length;
		const rebuilt: Nucleotide[][] = new Array(total);
		for (let curIdx = 0; curIdx < helices.length; curIdx++) {
			const oldIdx = survivorInverse.get(curIdx);
			if (oldIdx === undefined) continue;
			rebuilt[oldIdx] = curIdx === keptCurrentIdx ? keptKept : helices[curIdx];
		}
		snapshot.removed.forEach(slot => {
			const restored: Nucleotide[] = [];
			slot.ntIds.forEach(id => {
				const nt = elements.get(id);
				if (nt instanceof Nucleotide) restored.push(nt);
			});
			rebuilt[slot.oldIdx] = restored;
		});

		// Mutate the original array in place so existing references stay valid.
		helices.length = 0;
		rebuilt.forEach((slot, i) => { helices[i] = slot || []; });

		// Restore GridMap helixIds. A nucleotide that was in a removed slot snaps back to that
		// slot's oldIdx; every other mark goes through the inverse remap.
		const ntToOldHelix = new Map<number, number>();
		snapshot.removed.forEach(slot => slot.ntIds.forEach(id => ntToOldHelix.set(id, slot.oldIdx)));
		grid.forEach((mark, ntId) => {
			const overriden = ntToOldHelix.get(ntId);
			if (overriden !== undefined) {
				mark.helixId = overriden;
			} else {
				mark.helixId = inverseRemap(mark.helixId);
			}
		});

		return { inverseRemap };
	}

	// Fits a plane through the given points and returns the plane normal
	// (the eigenvector of the covariance-like matrix with the smallest eigenvalue).
	export function fitPlane(points: THREE.Vector3[]): THREE.Vector3 {
		// centroid
		const rc = new THREE.Vector3(0, 0, 0);
		points.forEach(p => rc.add(p));
		rc.divideScalar(points.length);

		// 3x3 symmetric accumulator A[i][j] += (p-rc)[i] * (p-rc)[j]
		const A: number[][] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
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
		const v: number[][] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
		const maxIter = 100;
		const eps = 1e-12;
		for (let iter = 0; iter < maxIter; iter++) {
			// find largest off-diagonal entry
			let p = 0, q = 1;
			let maxOff = Math.abs(a[0][1]);
			if (Math.abs(a[0][2]) > maxOff) { p = 0; q = 2; maxOff = Math.abs(a[0][2]); }
			if (Math.abs(a[1][2]) > maxOff) { p = 1; q = 2; maxOff = Math.abs(a[1][2]); }
			if (maxOff < eps) break;

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
		if (vals[1] < vals[minIdx]) minIdx = 1;
		if (vals[2] < vals[minIdx]) minIdx = 2;
		return new THREE.Vector3(v[0][minIdx], v[1][minIdx], v[2][minIdx]);
	}

	// Returns the axis of an RNA/DNA duplex given the four end nucleotides of the two strands.
	// - start1/end1 are the 5'/3' ends of strand A
	// - start2/end2 are the 5'/3' ends of strand B (start1 pairs with end2, end1 pairs with start2)
	export function getRNAAxis(d: {
		start1: Nucleotide;
		end1: Nucleotide;
		start2: Nucleotide;
		end2: Nucleotide;
	}): { planeVector: THREE.Vector3; finalHelPos: THREE.Vector3 } {
		const backboneSite = (nt: Nucleotide): THREE.Vector3 =>
			nt.getInstanceParameter3('bbOffsets');

		// initial guess vector from the midpoint of the start1-end2 pair to end1-start2 pair
		const midA0 = backboneSite(d.start1).add(backboneSite(d.end2)).multiplyScalar(0.5);
		const midAc0 = backboneSite(d.end1).add(backboneSite(d.start2)).multiplyScalar(0.5);
		const guess = midAc0.clone().sub(midA0);
		if (guess.length() > 0) guess.normalize();

		// Walk pairs (nucA on strand A via n3, nucB on strand B via n5) in lockstep.
		// nucAc/nucBc are the next pair along the walk.
		const posAs: THREE.Vector3[] = [];
		const posBs: THREE.Vector3[] = [];
		const backPoses: THREE.Vector3[] = [];

		let nucA: Nucleotide | null = d.start1;
		let nucB: Nucleotide | null = d.end2;
		while (nucA && nucB && nucA !== d.end1) {
			const nucAc = nucA.n3 as Nucleotide | null;
			const nucBc = nucB.n5 as Nucleotide | null;
			if (!nucAc || !nucBc) break;

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
		if (guess.dot(planeVector) < 0) planeVector.multiplyScalar(-1);

		// Find where the helical axis originates by intersecting per-base-pair perpendiculars
		// projected onto the plane.
		const helPos: THREE.Vector3[] = [];
		for (let i = 0; i < posAs.length - 1; i++) {
			// project current base pair to plane
			let apos = posAs[i].clone();
			let bpos = posBs[i].clone();
			apos.sub(planeVector.clone().multiplyScalar(apos.dot(planeVector)));
			bpos.sub(planeVector.clone().multiplyScalar(bpos.dot(planeVector)));
			const bpVecA = bpos.clone().sub(apos);
			if (bpVecA.length() === 0) continue;
			const midpointA = apos.clone().add(bpos).multiplyScalar(0.5);
			const perpA = bpVecA.clone().cross(planeVector).normalize();

			// project next base pair to plane
			let apos2 = posAs[i + 1].clone();
			let bpos2 = posBs[i + 1].clone();
			apos2.sub(planeVector.clone().multiplyScalar(apos2.dot(planeVector)));
			bpos2.sub(planeVector.clone().multiplyScalar(bpos2.dot(planeVector)));
			const bpVecB = bpos2.clone().sub(apos2);
			if (bpVecB.length() === 0) continue;
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
			if (Math.abs(det) < 1e-12) continue;
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
		if (helPos.length) finalHelPos.divideScalar(helPos.length);

		return { planeVector, finalHelPos };
	}

	// Draws the duplex axis vector in the scene, anchored at the start1 nucleotide.
	export function addRNAAxisToScene(d: {
		start1: Nucleotide;
		end1: Nucleotide;
		start2: Nucleotide;
		end2: Nucleotide;
	}) {
		const { planeVector } = getRNAAxis(d);
		const origin = d.start1.getInstanceParameter3('bbOffsets')
			.add(d.end2.getInstanceParameter3('bbOffsets'))
			.multiplyScalar(0.5);
		if (typeof THREE !== 'undefined' && typeof scene !== 'undefined' && (scene as any)?.add) {
			const arrow = new THREE.ArrowHelper(planeVector.clone().normalize(), origin, 10);
			(scene as any).add(arrow);
		}
	}

}