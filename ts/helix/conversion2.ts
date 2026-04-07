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

namespace toscad {
    export function helixEndpoints(helix: Nucleotide[]) {
        // Find the two most distant endpoints in the helix using BFS.

        // first we remove duplicates
        // best hope is that there never should be. All of the duplicates must necessarily be removed by findhelix2.ts.
        const nodes: Nucleotide[] = Array.from(new Map<number, Nucleotide>(helix.map((n: Nucleotide) => [n.id, n])).values());
        if (!nodes.length) return null;
        if (nodes.length !== helix.length) {
            console.log("Holy shit the world is doomed");
            console.log("Just kidding, there are duplicates in the helix");
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
        // const helices = await honda.findHelices(elements, 2);
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

    type GridMark = { helixId: number; offset: number; direction: 'forward' | 'backward' };
    type GridMap = Map<number, GridMark>;
    type Direction = 'n3' | 'n5';
    type LatticeKind = 'honeycomb' | 'square';
    type PhaseParity = 0 | 1;

    interface LatticePhaseConfig {
        basesPerTurn: number;
        phases: Record<PhaseParity, number[]>;
        voidPhase: number | null;
        tieEpsilon: number;
    }

    const LATTICE_CONFIG: Record<LatticeKind, LatticePhaseConfig> = {
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

    const resolveLatticeKind = (lattice: string): LatticeKind =>
        (lattice ?? '').toLowerCase() === 'square' ? 'square' : 'honeycomb';

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

    export function getAngles(
        grid: GridMap,
        helices: Nucleotide[][],
        helixId: number,
        lattice: string
    ): Map<number, {
        helixId: number;
        adj_helix: number;
        angle: number;
    }> {
        const result = new Map<number, {
            helixId: number;
            adj_helix: number;
            angle: number;
        }>();

        const requestedLattice = (lattice ?? '').toLowerCase();
        const latticeType = resolveLatticeKind(lattice);
        const latticeConfig = LATTICE_CONFIG[latticeType];

        const nearestPhase = (value: number, parity: PhaseParity) => {
            const candidates = latticeConfig.phases[parity] ?? [];
            if (!candidates.length) return value;

            return candidates.reduce((prev, curr) => {
                const distPrev = Math.abs(prev - value);
                const distCurr = Math.abs(curr - value);

                // For opposite-direction snaps, avoid void phase ties when configured.
                if (
                    parity === 1
                    && latticeConfig.voidPhase !== null
                    && Math.abs(distPrev - distCurr) < latticeConfig.tieEpsilon
                ) {
                    if (prev === latticeConfig.voidPhase && curr !== latticeConfig.voidPhase) return curr;
                    if (curr === latticeConfig.voidPhase && prev !== latticeConfig.voidPhase) return prev;
                }

                return distCurr < distPrev ? curr : prev;
            }, candidates[0]);
        };

        console.log(
            `[getAngles] Start helix=${helixId}, requestedLattice=${requestedLattice || 'unspecified'}, ` +
            `resolvedLattice=${latticeType}, basesPerTurn=${latticeConfig.basesPerTurn}, snapMode=phaseDictionary`
        );

        const helix = helices[helixId] ?? [];

        const endpointPair = crossoverEndpointsHelix(grid, helix, helixId);
        if (!endpointPair) {
            console.log(`[getAngles] No crossover endpoints found for helix ${helixId}`);
            return result;
        }

        const firstEnd = endpointPair.end1;
        const secondEnd = endpointPair.end2;
        const firstEndMark = grid.get(firstEnd.id);
        const secondEndMark = grid.get(secondEnd.id);

        if (!firstEndMark && !secondEndMark) {
            console.log(
                `[getAngles] Both endpoint marks missing from grid: end1=${firstEnd.id}, end2=${secondEnd.id}`
            );
            return result;
        }

        const useSecondEnd = !firstEndMark
            || (!!firstEndMark && !!secondEndMark && secondEndMark.offset < firstEndMark.offset);

        const end1 = useSecondEnd ? secondEnd : firstEnd;
        const end2 = useSecondEnd ? firstEnd : secondEnd;
        const end1Mark = useSecondEnd ? secondEndMark : firstEndMark;
        const end2Mark = useSecondEnd ? firstEndMark : secondEndMark;

        if (!end1Mark) {
            console.log(`[getAngles] Chosen endpoint nt=${end1.id} is missing from grid`);
            return result;
        }

        console.log(
            `[getAngles] crossoverEndpointsHelix(${helixId}) => end1=${endpointPair.end1.id}, end2=${endpointPair.end2.id}, distance=${endpointPair.diameter}`
        );
        console.log(
            `[getAngles] Reference endpoint selected: end1=${end1.id} (offset=${end1Mark.offset}), ` +
            `end2=${end2.id} (offset=${end2Mark?.offset ?? 'missing'})`
        );

        const nodes: Nucleotide[] = Array.from(
            new Map<number, Nucleotide>(helix.map((nt: Nucleotide) => [nt.id, nt])).values()
        );
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

        // nt -> set of connected helix ids through crossover transitions.
        const ntToConnectedHelices = new Map<number, Set<number>>();
        const addConnection = (ntId: number, otherHelix: number) => {
            if (!ntToConnectedHelices.has(ntId)) ntToConnectedHelices.set(ntId, new Set<number>());
            ntToConnectedHelices.get(ntId)!.add(otherHelix);
        };

        const crossoverList = crossoverNts(grid);
        for (const crossover of crossoverList) {
            addConnection(crossover.fromNt.id, crossover.toHelix);
            addConnection(crossover.toNt.id, crossover.fromHelix);
        }

        const helixCrossNtIds = new Set<number>();
        for (const nt of nodes) {
            if (ntToConnectedHelices.has(nt.id)) helixCrossNtIds.add(nt.id);
        }
        if (!helixCrossNtIds.size) {
            console.log(`[getAngles] Helix ${helixId} has no crossover nucleotides in this grid`);
            return result;
        }

        // Use collectCrossovers to discover which helices this helix connects to.
        const { crossovers } = collectCrossovers(grid);
        const connectedHelices = new Set<number>();
        const forward = crossovers.get(helixId);
        if (forward) {
            for (const [toHelix] of forward.entries()) connectedHelices.add(toHelix);
        }
        for (const [fromHelix, toMap] of crossovers.entries()) {
            if (toMap.has(helixId)) connectedHelices.add(fromHelix);
        }

        // Fallback from per-nt connectivity if collectCrossovers misses anything.
        for (const ntId of helixCrossNtIds) {
            const connected = ntToConnectedHelices.get(ntId);
            if (!connected) continue;
            for (const h of connected) connectedHelices.add(h);
        }
        connectedHelices.delete(helixId);
        if (!connectedHelices.size) {
            console.log(`[getAngles] Helix ${helixId} has no connected helices`);
            return result;
        }

        console.log(
            `[getAngles] Reference nt end1=${end1.id} offset=${end1Mark.offset} direction=${end1Mark.direction}`
        );

        const end1Connected = ntToConnectedHelices.get(end1.id) ?? new Set<number>();
        const seenHelices = new Set<number>();

        console.log(
            `[getAngles] end1 nt=${end1.id} connects to helices: [${Array.from(end1Connected).sort((a, b) => a - b).join(', ')}]`
        );

        // Baseline: helix connected at end1 has zero relative angle.
        for (const h of end1Connected) {
            if (!connectedHelices.has(h)) continue;
            result.set(h, {
                helixId,
                adj_helix: h,
                angle: 0
            });
            seenHelices.add(h);
            console.log(
                `[getAngles] Baseline angle set: nt=${end1.id}, helix ${helixId} -> helix ${h}, angle=0`
            );
        }

        const startNode = nodeById.get(end1.id) ?? end1;
        const distFromEnd1 = bfsDistances(startNode);
        const orderedCrossovers = Array.from(helixCrossNtIds)
            .filter((ntId) => ntId !== end1.id)
            .map((ntId) => ({ ntId, dist: distFromEnd1.get(ntId) }))
            .filter((entry): entry is { ntId: number; dist: number } => entry.dist !== undefined)
            .sort((a, b) => a.dist - b.dist || a.ntId - b.ntId);

        for (const { ntId, dist } of orderedCrossovers) {
            if (seenHelices.size >= connectedHelices.size) break;

            const candidateConnections = ntToConnectedHelices.get(ntId);
            if (!candidateConnections) {
                console.log(`[getAngles] Skip nt=${ntId} (distance=${dist}): no connected helices found`);
                continue;
            }

            const candidateConnectionList = Array.from(candidateConnections).sort((a, b) => a - b);
            console.log(
                `[getAngles] Consider nt=${ntId} (distance=${dist}) connected helices=[${candidateConnectionList.join(', ')}]`
            );

            const freshHelices: number[] = [];
            for (const h of candidateConnections) {
                if (!connectedHelices.has(h)) continue;
                if (seenHelices.has(h)) continue;
                if (end1Connected.has(h)) continue;
                freshHelices.push(h);
            }
            if (freshHelices.length === 0) {
                console.log(
                    `[getAngles] Skip nt=${ntId}: only already-seen or end1-connected helices`
                );
                continue;
            }

            const otherMark = grid.get(ntId);
            if (!otherMark) {
                console.log(`[getAngles] Skip nt=${ntId}: missing grid mark`);
                continue;
            }

            const rawX = otherMark.offset - end1Mark.offset;
            const absX = Math.abs(rawX);
            const sign = Math.sign(rawX) || 1;
            const y: PhaseParity = end1Mark.direction === otherMark.direction ? 0 : 1;
            const phase = absX % latticeConfig.basesPerTurn;
            const phases = latticeConfig.phases[y] ?? [];
            const idealPhase = nearestPhase(phase, y);
            const idealAbsX = (Math.floor(absX / latticeConfig.basesPerTurn) * latticeConfig.basesPerTurn) + idealPhase;
            const idealX = idealAbsX * sign + 1;
            const angleRaw = (360 / latticeConfig.basesPerTurn) * idealX + (y * 215 * sign);
            const angle = Math.round(((angleRaw % 360) + 360) % 360);

            console.log(
                `[getAngles] Use nts end1=${end1.id} and other=${ntId}; ` +
                `offsets=(${end1Mark.offset},${otherMark.offset}) directions=(${end1Mark.direction},${otherMark.direction}) ` +
                `rawX=${rawX} absX=${absX} y=${y} basesPerTurn=${latticeConfig.basesPerTurn} ` +
                `phase=${phase} candidates=[${phases.join(', ')}] idealPhase=${idealPhase} idealX=${idealX} raw=${angleRaw} angle=${angle}`
            );

            for (const h of freshHelices) {
                result.set(h, {
                    helixId,
                    adj_helix: h,
                    angle: angle
                });
                seenHelices.add(h);
                console.log(
                    `[getAngles] Angle set: helix ${helixId} -> helix ${h} via nt ${ntId} = ${angle}`
                );
            }
        }

        console.log(
            `[getAngles] Final angle map for helix ${helixId}: ` +
            `${JSON.stringify(Array.from(result.entries()).sort((a, b) => a[0] - b[0]))}`
        );

        return result;
    }

    // Helper interface to keep track of where we are on both strands
    interface DualCursor {
        fwd: Nucleotide | null;
        bwd: Nucleotide | null;
    }

    export function setGrid(helices: Nucleotide[][]): { grid: GridMap; binderHelices: number[] } {
        const grid: GridMap = new Map();

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

            let offset = 0; // Local offset for the main backbone

            // --- A-D. MAIN BACKBONE LOGIC ---
            if (endpoints) {
                // start "forward" from any endpoint. They will be oriented later. Our main priority is to generate a grid without overlap and sufficient details.
                const helixFwd = endpoints.end1;
                const helixFwdDir = (isInHelix(helixSet, helixFwd.n3 as Nucleotide) ? 'n3' : 'n5');
                const helixBwdDir = (helixFwdDir === 'n3' ? 'n5' : 'n3');
                const revFwdDir = helixFwdDir === 'n3' ? 'n5' : 'n3';
                const revBwdDir = helixBwdDir === 'n3' ? 'n5' : 'n3';

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
                        if (i < headFwd.length - 1) mark(n, helixId, (startOffset - fwdHeadLen) + i, 'forward');
                    });
                    [...headBwd].reverse().forEach((n, i) => {
                        if (i < headBwd.length - 1) mark(n, helixId, (startOffset - bwdHeadLen) + i, 'backward');
                    });
                    offset = startOffset;
                } else {
                    firstAnchor = helixFwd;
                }
                // By here we have the correct offset for the first anchorpoint.

                // The Body
                let currFwd = firstAnchor;
                let currBwd = firstAnchorPair;
                if (currFwd) mark(currFwd, helixId, offset, 'forward');
                if (currBwd) mark(currBwd, helixId, offset, 'backward');

                while (currFwd) {
                    const nextStep = findNextPaired({ fwd: currFwd, bwd: currBwd }, { fwd: helixFwdDir, bwd: helixBwdDir }, helixSet);
                    if (!nextStep) break;
                    const nextAnchor = nextStep.anchor;
                    const nextPair = getPair(helixSet, nextAnchor);
                    // probably doesn't need this check but can happen due to cross/double pairing?
                    if (nextAnchor.id === currFwd.id) break;

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
                    fwdTail.forEach((n, i) => mark(n, helixId, offset + i + 1, 'forward'));
                    const fwdHeadStart = offset + 1 + (gapLength - fwdHead.length);
                    fwdHead.forEach((n, i) => mark(n, helixId, fwdHeadStart + i, 'forward'));
                    bwdTail.forEach((n, i) => mark(n, helixId, offset + i + 1, 'backward'));
                    const bwdHeadStart = offset + 1 + (gapLength - bwdHead.length);
                    bwdHead.forEach((n, i) => mark(n, helixId, bwdHeadStart + i, 'backward'));

                    offset += gapLength + 1;
                    mark(nextAnchor, helixId, offset, 'forward');
                    if (nextPair) mark(nextPair, helixId, offset, 'backward');

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

    // ── Shared helper: collect all crossover shift observations ──────
    // For each pair of helices connected by backbone crossovers, returns
    // the list of observed shifts (offsetA - offsetB for each crossover
    // from A→B). Used by both alignGridPrim and alignGridDP.
    export function collectShiftObservations(grid: GridMap) {
        const allNtIds = new Set<number>();
        for (const [ntId] of grid.entries()) allNtIds.add(ntId);

        const visited = new Set<number>();
        const helixIds = new Set<number>();

        // shifts[a][b] = array of (offsetA - offsetB) values
        const shifts = new Map<number, Map<number, number[]>>();

        const ensurePair = (a: number, b: number) => {
            if (!shifts.has(a)) shifts.set(a, new Map());
            if (!shifts.get(a)!.has(b)) shifts.get(a)!.set(b, []);
            return shifts.get(a)!.get(b)!;
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

            // Walk 5'→3'
            let curr: Nucleotide | null = fivePrime;
            const walkForward = new Set<number>();
            let prevMark: GridMark | null = null;

            while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                if (walkForward.has(curr.id)) break;
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
                } else {
                    prevMark = null;
                }

                const n3ref: any = curr.n3;
                curr = (n3ref && n3ref instanceof Nucleotide) ? (n3ref as Nucleotide) : null;
            }
        }

        return { shifts, helixIds };
    }

    // ── Helper: compute median of a sorted-or-unsorted number array ─
    function median(arr: number[]): number {
        if (arr.length === 0) return 0;
        const sorted = arr.slice().sort((a, b) => a - b);
        const mid = Math.floor(sorted.length / 2);
        return sorted.length % 2 === 0
            ? Math.round((sorted[mid - 1] + sorted[mid]) / 2)
            : sorted[mid];
    }

    // ── Helper: apply shift to every nt on a helix, then normalize ──
    function applyHelixShifts(grid: GridMap, shiftMap: Map<number, number>) {
        for (const [, mark] of grid.entries()) {
            const s = shiftMap.get(mark.helixId);
            if (s !== undefined && s !== 0) {
                mark.offset += s;
            }
        }

        // Normalize: find global minimum, shift everything so min = 0
        let globalMin = Infinity;
        for (const [, mark] of grid.entries()) {
            if (mark.offset < globalMin) globalMin = mark.offset;
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
    export function alignGridPrim(grid: GridMap, binderHelices?: number[]) {
        const { shifts, helixIds } = collectShiftObservations(grid);

        const countGridConflicts = (currentGrid: GridMap): number => {
            const checkMap = new Map<number, {
                forward: Map<number, number>;
                backward: Map<number, number>;
            }>();

            let conflicts = 0;
            for (const [ntId, pos] of currentGrid.entries()) {
                if (!checkMap.has(pos.helixId)) {
                    checkMap.set(pos.helixId, {
                        forward: new Map(),
                        backward: new Map()
                    });
                }

                const strandMap = checkMap.get(pos.helixId)![pos.direction];
                if (strandMap.has(pos.offset) && strandMap.get(pos.offset) !== ntId) {
                    conflicts++;
                } else {
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
        const inMST = new Set<number>();
        // mstEdges: parent → child with median shift
        const mstParent = new Map<number, { parent: number; shift: number }>();
        // Priority: pick the edge with the highest weight (most observations)
        inMST.add(0);

        while (inMST.size < helixList.length) {
            let bestNeighbor = -1;
            let bestFrom = -1;
            let bestWeight = 0;

            for (const inNode of inMST) {
                const neighbors = shifts.get(inNode);
                if (!neighbors) continue;

                for (const [neighbor, observations] of neighbors.entries()) {
                    if (inMST.has(neighbor)) continue;
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
            const observations = shifts.get(bestFrom)!.get(bestNeighbor)!;
            const medianShift = median(observations);

            mstParent.set(bestNeighbor, { parent: bestFrom, shift: medianShift });
            inMST.add(bestNeighbor);
        }

        // ── Walk MST from helix 0 to compute cumulative shifts ──────────
        const cumulativeShift = new Map<number, number>();
        cumulativeShift.set(0, 0); // anchor

        // BFS order: process nodes so parent's cumulative shift is known
        const bfsQueue: number[] = [0];
        let qi = 0;

        // Build children adjacency from mstParent
        const children = new Map<number, number[]>();
        for (const [child, { parent }] of mstParent.entries()) {
            if (!children.has(parent)) children.set(parent, []);
            children.get(parent)!.push(child);
        }

        while (qi < bfsQueue.length) {
            const node = bfsQueue[qi++];
            const nodeShift = cumulativeShift.get(node) ?? 0;
            const kids = children.get(node) ?? [];

            for (const child of kids) {
                const edge = mstParent.get(child)!;
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
            const allNtIds = new Set<number>();
            for (const [ntId] of grid.entries()) allNtIds.add(ntId);
            const binderVisited = new Set<number>();

            for (const [ntId] of grid.entries()) {
                if (binderVisited.has(ntId)) continue;

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

                // Walk 5'→3', split into runs by helixId
                type BinderRun = { helixId: number; ntIds: number[] };
                const runs: BinderRun[] = [];
                let currentRun: BinderRun | null = null;
                let curr: Nucleotide | null = fivePrime;
                const walkFwd = new Set<number>();

                while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                    if (walkFwd.has(curr.id)) break;
                    walkFwd.add(curr.id);
                    binderVisited.add(curr.id);

                    const mark = grid.get(curr.id);
                    if (mark) {
                        if (currentRun && currentRun.helixId === mark.helixId) {
                            currentRun.ntIds.push(curr.id);
                        } else {
                            currentRun = { helixId: mark.helixId, ntIds: [curr.id] };
                            runs.push(currentRun);
                        }
                    } else {
                        currentRun = null;
                    }

                    const n3ref: any = curr.n3;
                    curr = (n3ref && n3ref instanceof Nucleotide) ? (n3ref as Nucleotide) : null;
                }

                // For each run on a binder helix, find the crossover offset
                // from its adjacent non-binder run and align.
                for (let i = 0; i < runs.length; i++) {
                    const run = runs[i];
                    if (!binderSet.has(run.helixId)) continue;

                    // Look for the adjacent non-binder run to get the
                    // crossover offset. Check the run before and after.
                    let parentOffset: number | null = null;
                    let binderCrossoverOffset: number | null = null;

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
                                if (mark) mark.offset += segmentShift;
                            }
                        }
                    }
                }
            }
        }

        // ── Final normalization ──────────────────────────────────────────
        let globalMin2 = Infinity;
        for (const [, mark] of grid.entries()) {
            if (mark.offset < globalMin2) globalMin2 = mark.offset;
        }
        if (globalMin2 !== 0 && globalMin2 !== Infinity) {
            for (const [, mark] of grid.entries()) {
                mark.offset -= globalMin2;
            }
        }

        // ── Binder post-pass (after alignGridPrim) ─────────────────────
        const binderPostSet = new Set<number>(binderHelices ?? []);
        if (binderPostSet.size > 0) {
            const binderList = Array.from(binderPostSet).sort((a, b) => a - b);
            console.log(`[alignGridPrim] Binder helices noted: [${binderList.join(', ')}]`);

            type BinderRun = { helixId: number; ntIds: number[] };
            const allNtIds = new Set<number>();
            for (const [ntId] of grid.entries()) allNtIds.add(ntId);
            const visited = new Set<number>();
            const binderRuns: BinderRun[] = [];
            const ntToRun = new Map<number, number>();

            // Rebuild strand runs and keep only runs on binder helices.
            for (const [ntId] of grid.entries()) {
                if (visited.has(ntId)) continue;

                const startNt = elements.get(ntId) as Nucleotide | undefined;
                if (!startNt || !(startNt instanceof Nucleotide)) continue;

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

                let curr: Nucleotide | null = fivePrime;
                const walkFwd = new Set<number>();
                let currentRunNtIds: number[] = [];
                let currentHelix: number | null = null;

                while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                    if (walkFwd.has(curr.id)) break;
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
                    } else if (currentHelix === mark.helixId || currentHelix === null) {
                        currentHelix = mark.helixId;
                        currentRunNtIds.push(curr.id);
                    } else {
                        if (binderPostSet.has(currentHelix) && currentRunNtIds.length > 0) {
                            const runIdx = binderRuns.length;
                            binderRuns.push({ helixId: currentHelix, ntIds: currentRunNtIds });
                            currentRunNtIds.forEach((id) => ntToRun.set(id, runIdx));
                        }
                        currentHelix = mark.helixId;
                        currentRunNtIds = [curr.id];
                    }

                    const n3ref: any = curr.n3;
                    curr = (n3ref && n3ref instanceof Nucleotide) ? (n3ref as Nucleotide) : null;
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
                    if (mark) mark.direction = 'forward';
                }
            }

            const runCenter = (run: BinderRun): number => {
                let sum = 0;
                let count = 0;
                for (const ntId of run.ntIds) {
                    const mark = grid.get(ntId);
                    if (!mark) continue;
                    sum += mark.offset;
                    count++;
                }
                return count > 0 ? sum / count : 0;
            };

            const shiftRunBy = (run: BinderRun, delta: number) => {
                if (delta === 0) return;
                for (const ntId of run.ntIds) {
                    const mark = grid.get(ntId);
                    if (mark) mark.offset += delta;
                }
            };

            // Keep pushing the front run by +5 until no binder run overlaps remain.
            let guard = 0;
            while (guard++ < 2000) {
                const overlapBuckets = new Map<string, Set<number>>();

                for (const run of binderRuns) {
                    for (const ntId of run.ntIds) {
                        const mark = grid.get(ntId);
                        if (!mark || !binderPostSet.has(mark.helixId)) continue;
                        const key = `${mark.helixId}|${mark.offset}`;
                        if (!overlapBuckets.has(key)) overlapBuckets.set(key, new Set<number>());
                        const runIdx = ntToRun.get(ntId);
                        if (runIdx !== undefined) overlapBuckets.get(key)!.add(runIdx);
                    }
                }

                let moved = false;
                for (const [, runSet] of overlapBuckets.entries()) {
                    if (runSet.size <= 1) continue;

                    let frontRun: BinderRun | null = null;
                    let frontCenter = -Infinity;
                    for (const runIdx of runSet.values()) {
                        const run = binderRuns[runIdx];
                        if (!run) continue;
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

                if (!moved) break;
            }
        }

        console.log(`[alignGridPrim] Aligned ${helixList.length} helices. Shifts:`,
            Object.fromEntries(Array.from(cumulativeShift.entries()).sort((a, b) => a[0] - b[0]))
        );

        validateGrid(grid);
        const conflictsAfter = countGridConflicts(grid);
        if (conflictsAfter > 0) {
            throw new Error(`[alignGridPrim] Overlaps remain after post-pass: ${conflictsAfter}`);
        }

        return { shifts: cumulativeShift };
    }

    function computeHelixPcaAxis(helix: Nucleotide[] | undefined): THREE.Vector3 | null {
        if (!helix || helix.length < 2) return null;

        const points: THREE.Vector3[] = [];
        for (const nt of helix) {
            if (!nt) continue;
            const pos = nt.getPos();
            if (!pos) continue;
            points.push(pos.clone());
        }
        if (points.length < 2) return null;

        const centroid = new THREE.Vector3();
        for (const p of points) centroid.add(p);
        centroid.multiplyScalar(1 / points.length);

        const centered: THREE.Vector3[] = [];
        for (const p of points) {
            const rel = p.clone().sub(centroid);
            if (rel.lengthSq() > 1e-12) centered.push(rel);
        }
        if (centered.length < 2) return null;

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
        if (axis.lengthSq() <= 1e-12) return null;

        for (let i = 0; i < 16; i++) {
            const next = new THREE.Vector3(
                cov.xx * axis.x + cov.xy * axis.y + cov.xz * axis.z,
                cov.xy * axis.x + cov.yy * axis.y + cov.yz * axis.z,
                cov.xz * axis.x + cov.yz * axis.y + cov.zz * axis.z,
            );
            if (next.lengthSq() <= 1e-12) break;
            next.normalize();
            if (next.dot(axis) < 0) next.multiplyScalar(-1);
            axis = next;
        }

        return axis.lengthSq() > 1e-12 ? axis.normalize() : null;
    }

    function computeHelixCentroid(helix: Nucleotide[] | undefined): THREE.Vector3 | null {
        if (!helix || helix.length === 0) return null;
        const sum = new THREE.Vector3();
        let count = 0;
        for (const nt of helix) {
            if (!nt) continue;
            const pos = nt.getPos();
            if (!pos) continue;
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
    function areHelicesColinear(
        helixAId: number,
        helixBId: number,
        helices: Nucleotide[][] | undefined,
        pcaAxisCache: Map<number, THREE.Vector3 | null>,
        centroidCache: Map<number, THREE.Vector3 | null>,
        maxOffAxisDist: number = 10,
    ): boolean {
        if (!helices) return true;

        const getAxis = (hId: number) => {
            if (!pcaAxisCache.has(hId)) pcaAxisCache.set(hId, computeHelixPcaAxis(helices[hId]));
            return pcaAxisCache.get(hId) ?? null;
        };
        const getCentroid = (hId: number) => {
            if (!centroidCache.has(hId)) centroidCache.set(hId, computeHelixCentroid(helices[hId]));
            return centroidCache.get(hId) ?? null;
        };

        const axisA = getAxis(helixAId);
        const centroidA = getCentroid(helixAId);
        const centroidB = getCentroid(helixBId);

        if (!axisA || !centroidA || !centroidB) return true; // no data — allow

        // Perpendicular (rejection) distance from centroidB to the line
        // through centroidA along axisA:
        //   d = |(centroidB - centroidA) × axisA|
        const diff = centroidB.clone().sub(centroidA);
        const rejection = diff.clone().cross(axisA).length();
        return rejection <= maxOffAxisDist;
    }

    function areHelixPcaAxesCompatible(
        helixAId: number,
        helixBId: number,
        helices: Nucleotide[][] | undefined,
        pcaAxisCache: Map<number, THREE.Vector3 | null>,
        maxAngleDeg: number = 45,
    ): boolean {
        if (!helices) return true;

        const getAxis = (helixId: number) => {
            if (!pcaAxisCache.has(helixId)) {
                pcaAxisCache.set(helixId, computeHelixPcaAxis(helices[helixId]));
            }
            return pcaAxisCache.get(helixId) ?? null;
        };

        const axisA = getAxis(helixAId);
        const axisB = getAxis(helixBId);

        if (!axisA || !axisB) return false;

        const cosThreshold = Math.cos(maxAngleDeg * Math.PI / 180);
        const cosine = Math.min(1, Math.max(-1, Math.abs(axisA.dot(axisB))));
        return cosine >= cosThreshold;
    }

    // 2 helices are candidates for combining when one of their mutual connection has >3 helices connected to it.
    // They get combined only if they have no offset overalp and their PCA axes are compatible (roughly parallel, allowing for some angle).
    export function combinedHelices(maxOffsetDist: number, grid: GridMap, helices?: Nucleotide[][], binderHelices?: number[]) {
        const retiredHelices = new Set<number>();
        const touchedHelices = new Set<number>();
        const combinedPairs: Array<{ keep: number; merged: number }> = [];
        const binderSet = new Set<number>(binderHelices ?? []);
        const pcaAxisCache = new Map<number, THREE.Vector3 | null>();
        const centroidCache = new Map<number, THREE.Vector3 | null>();
        const maxOffs: number = maxOffsetDist;

        const buildAdjacency = (crossovers: Map<number, Map<number, { sameWalk: number; diffWalk: number }>>) => {
            const adjacency = new Map<number, Set<number>>();

            const ensure = (hId: number) => {
                if (!adjacency.has(hId)) adjacency.set(hId, new Set<number>());
                return adjacency.get(hId)!;
            };

            for (const [from, toMap] of crossovers.entries()) {
                const fromSet = ensure(from);
                for (const [to, stats] of toMap.entries()) {
                    const totalConnections = stats.sameWalk + stats.diffWalk;
                    if (totalConnections <= 0) continue;
                    fromSet.add(to);
                    ensure(to).add(from);
                }
            }

            return adjacency;
        };

        const buildOffsetSets = () => {
            const offsetSets = new Map<number, Set<number>>();
            for (const [, mark] of grid.entries()) {
                if (!offsetSets.has(mark.helixId)) offsetSets.set(mark.helixId, new Set<number>());
                offsetSets.get(mark.helixId)!.add(mark.offset);
            }
            return offsetSets;
        };

        // stupid??
        const offsetsDisjoint = (a: Set<number> | undefined, b: Set<number> | undefined) => {
            if (!a || !b || a.size === 0 || b.size === 0) return false;
            const smaller = a.size <= b.size ? a : b;
            const larger = a.size <= b.size ? b : a;
            for (const offs of smaller) {
                if (larger.has(offs)) return false;
            }
            return true;
        };

        const distanceSquared = (a: number, b: number, positions: Map<number, { x: number, y: number }>) => {
            const pa = positions.get(a);
            const pb = positions.get(b);
            if (!pa || !pb) return Number.POSITIVE_INFINITY;
            const dx = pb.x - pa.x;
            const dy = pb.y - pa.y;
            return dx * dx + dy * dy;
        };

        const mergeHelixInto = (keep: number, merged: number) => {
            for (const [, mark] of grid.entries()) {
                if (mark.helixId === merged) {
                    mark.helixId = keep;
                }
            }

            if (helices && helices[merged] && helices[merged].length > 0) {
                if (!helices[keep]) helices[keep] = [];
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

        const positions = helices ? getRelativePositions(helices) : new Map<number, { x: number, y: number }>();

        let iteration = 0;
        while (iteration++ < 100) {
            const { crossovers, helixIds } = collectCrossovers(grid);
            const adjacency = buildAdjacency(crossovers);
            const offsetSets = buildOffsetSets();

            for (const hId of helixIds) {
                if (!adjacency.has(hId)) adjacency.set(hId, new Set<number>());
            }

            const hubs = Array.from(adjacency.entries())
                .filter(([helixId, neighbors]) => neighbors.size > 3 && !retiredHelices.has(helixId) && !binderSet.has(helixId))
                .map(([helixId]) => helixId)
                .sort((a, b) => a - b);

            if (hubs.length === 0) break;

            let mergedInThisIteration = false;

            for (const hub of hubs) {
                if (retiredHelices.has(hub)) continue;
                if (binderSet.has(hub)) continue;

                const neighbors = Array.from(adjacency.get(hub) ?? [])
                    .filter(n => !retiredHelices.has(n) && n !== hub && !binderSet.has(n));

                if (neighbors.length === 0) continue;

                // Step 1: try combining one neighbor directly into the overloaded hub.
                const compatibleWithHub = neighbors.filter(n =>
                    offsetsDisjoint(offsetSets.get(hub), offsetSets.get(n))
                    && areHelixPcaAxesCompatible(hub, n, helices, pcaAxisCache, 45)
                    && areHelicesColinear(hub, n, helices, pcaAxisCache, centroidCache, maxOffs)
                );

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
                        if (!offsetSets.has(hub)) offsetSets.set(hub, new Set<number>());
                        for (const off of neighborOffs) offsetSets.get(hub)!.add(off);
                    }

                    mergeHelixInto(hub, bestNeighbor);
                    mergedInThisIteration = true;
                }

                // Step 2: try combining overloaded hub neighbors with each other.
                const remainingNeighbors = neighbors
                    .filter(n => !retiredHelices.has(n));

                while (remainingNeighbors.length >= 2) {
                    let bestPair: [number, number] | null = null;
                    let bestPairDist = Number.POSITIVE_INFINITY;

                    for (let i = 0; i < remainingNeighbors.length; i++) {
                        for (let j = i + 1; j < remainingNeighbors.length; j++) {
                            const a = remainingNeighbors[i];
                            const b = remainingNeighbors[j];
                            if (!offsetsDisjoint(offsetSets.get(a), offsetSets.get(b))) continue;
                            if (!areHelixPcaAxesCompatible(a, b, helices, pcaAxisCache, 45)) continue;
                            if (!areHelicesColinear(a, b, helices, pcaAxisCache, centroidCache, maxOffs)) continue;
                            const distSq = distanceSquared(a, b, positions);
                            if (distSq < bestPairDist) {
                                bestPairDist = distSq;
                                bestPair = [a, b];
                            }
                        }
                    }

                    if (!bestPair) break;

                    const keep = Math.min(bestPair[0], bestPair[1]);
                    const merged = Math.max(bestPair[0], bestPair[1]);

                    // Refresh keep's offset snapshot so the next pair check in this loop sees the combined range.
                    const mergedOffs = offsetSets.get(merged);
                    if (mergedOffs) {
                        if (!offsetSets.has(keep)) offsetSets.set(keep, new Set<number>());
                        for (const off of mergedOffs) offsetSets.get(keep)!.add(off);
                    }

                    mergeHelixInto(keep, merged);
                    mergedInThisIteration = true;

                    const mergedIdx = remainingNeighbors.indexOf(merged);
                    if (mergedIdx >= 0) remainingNeighbors.splice(mergedIdx, 1);
                }
            }

            if (!mergedInThisIteration) break;
        }

        const activeAfterMerge = new Set<number>();
        for (const [, mark] of grid.entries()) {
            activeAfterMerge.add(mark.helixId);
        }

        const maxHelixIdAfterMerge = activeAfterMerge.size > 0
            ? Math.max(...Array.from(activeAfterMerge))
            : -1;

        const removedHelices: number[] = [];
        for (let hId = 0; hId <= maxHelixIdAfterMerge; hId++) {
            if (!activeAfterMerge.has(hId)) removedHelices.push(hId);
        }

        const sortedActive = Array.from(activeAfterMerge).sort((a, b) => a - b);
        const helixIdRemap = new Map<number, number>();
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
                combinedPairs.map(p => `${p.merged}->${p.keep}`).join(', ')
            );
            if (removedHelices.length > 0) {
                console.log(`[combinedHelices] Removed empty helices: [${removedHelices.join(', ')}]`);
            }
        } else {
            console.log('[combinedHelices] No helix combinations applied.');
        }

        return {
            combinedPairs,
            combinedHelices: Array.from(touchedHelices).sort((a, b) => a - b),
            removedHelices,
            helixIdRemap: Object.fromEntries(helixIdRemap.entries())
        };
    }

    // a much more hand-wavy function, combines helices if they don't overlap and they have the same connections.
    export function trialComb(grid: GridMap, helices?: Nucleotide[][], binderHelices?: number[]) {
        const combinedPairs: Array<{ keep: number; merged: number }> = [];
        const retiredHelices = new Set<number>();
        const touchedHelices = new Set<number>();
        const binderSet = new Set<number>(binderHelices ?? []);
        const pcaAxisCache = new Map<number, THREE.Vector3 | null>();

        const buildAdjacency = (crossovers: Map<number, Map<number, { sameWalk: number; diffWalk: number }>>) => {
            const adjacency = new Map<number, Set<number>>();

            const ensure = (hId: number) => {
                if (!adjacency.has(hId)) adjacency.set(hId, new Set<number>());
                return adjacency.get(hId)!;
            };

            for (const [from, toMap] of crossovers.entries()) {
                const fromSet = ensure(from);
                for (const [to, stats] of toMap.entries()) {
                    const totalConnections = stats.sameWalk + stats.diffWalk;
                    if (totalConnections <= 0) continue;
                    fromSet.add(to);
                    ensure(to).add(from);
                }
            }

            return adjacency;
        };

        const buildOffsetSets = () => {
            const offsetSets = new Map<number, Set<number>>();
            for (const [, mark] of grid.entries()) {
                if (!offsetSets.has(mark.helixId)) offsetSets.set(mark.helixId, new Set<number>());
                offsetSets.get(mark.helixId)!.add(mark.offset);
            }
            return offsetSets;
        };

        const setEqual = (a: Set<number>, b: Set<number>) => {
            if (a.size !== b.size) return false;
            for (const v of a) {
                if (!b.has(v)) return false;
            }
            return true;
        };

        const offsetsDisjoint = (a: Set<number> | undefined, b: Set<number> | undefined) => {
            if (!a || !b || a.size === 0 || b.size === 0) return false;
            const smaller = a.size <= b.size ? a : b;
            const larger = a.size <= b.size ? b : a;
            for (const off of smaller) {
                if (larger.has(off)) return false;
            }
            return true;
        };

        let iteration = 0;
        while (iteration++ < 100) {
            const { crossovers, helixIds } = collectCrossovers(grid);
            const adjacency = buildAdjacency(crossovers);
            const offsetSets = buildOffsetSets();

            for (const hId of helixIds) {
                if (!adjacency.has(hId)) adjacency.set(hId, new Set<number>());
            }

            const candidates = Array.from(adjacency.keys()).sort((a, b) => a - b);
            let mergedInIteration = false;

            for (let i = 0; i < candidates.length; i++) {
                const a = candidates[i];
                if (retiredHelices.has(a)) continue;
                if (binderSet.has(a)) continue;

                for (let j = i + 1; j < candidates.length; j++) {
                    const b = candidates[j];
                    if (retiredHelices.has(b)) continue;
                    if (binderSet.has(b)) continue;

                    const aNeighbors = adjacency.get(a) ?? new Set<number>();
                    const bNeighbors = adjacency.get(b) ?? new Set<number>();

                    if (!setEqual(aNeighbors, bNeighbors)) continue;
                    if (!offsetsDisjoint(offsetSets.get(a), offsetSets.get(b))) continue;
                    if (!areHelixPcaAxesCompatible(a, b, helices, pcaAxisCache, 45)) continue;

                    const keep = Math.min(a, b);
                    const merged = Math.max(a, b);

                    // Refresh keep's offset snapshot before merging so later candidates in this pass
                    // see the full combined range and don't incorrectly pass the disjoint check.
                    const trialMergedOffs = offsetSets.get(merged);
                    if (trialMergedOffs) {
                        if (!offsetSets.has(keep)) offsetSets.set(keep, new Set<number>());
                        for (const off of trialMergedOffs) offsetSets.get(keep)!.add(off);
                    }

                    for (const [, mark] of grid.entries()) {
                        if (mark.helixId === merged) {
                            mark.helixId = keep;
                        }
                    }

                    if (helices && helices[merged] && helices[merged].length > 0) {
                        if (!helices[keep]) helices[keep] = [];
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

            if (!mergedInIteration) break;
        }

        const activeAfterMerge = new Set<number>();
        for (const [, mark] of grid.entries()) {
            activeAfterMerge.add(mark.helixId);
        }

        const maxHelixIdAfterMerge = activeAfterMerge.size > 0
            ? Math.max(...Array.from(activeAfterMerge))
            : -1;

        const removedHelices: number[] = [];
        for (let hId = 0; hId <= maxHelixIdAfterMerge; hId++) {
            if (!activeAfterMerge.has(hId)) removedHelices.push(hId);
        }

        const sortedActive = Array.from(activeAfterMerge).sort((a, b) => a - b);
        const helixIdRemap = new Map<number, number>();
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
                combinedPairs.map(p => `${p.merged}->${p.keep}`).join(', ')
            );
            if (removedHelices.length > 0) {
                console.log(`[trialComb] Removed empty helices: [${removedHelices.join(', ')}]`);
            }
        } else {
            console.log('[trialComb] No helix combinations applied.');
        }

        return {
            combinedPairs,
            combinedHelices: Array.from(touchedHelices).sort((a, b) => a - b),
            removedHelices,
            helixIdRemap: Object.fromEntries(helixIdRemap.entries())
        };
    }

    function subtractPositions(map: Map<number, { x: number, y: number }>, keyA: number, keyB: number) {
        const posA = map.get(keyA);
        const posB = map.get(keyB);
        if (!posA || !posB) return null;

        return {
            x: posA.x - posB.x,
            y: posA.y - posB.y
        };
    }

    export function HelixPosByRelativeBfs(grid: GridMap, helices: Nucleotide[][]): Map<number, [number, number]> {
        type GridPosition = [number, number];
        type CartesianCoord = { x: number; y: number };

        const placed = new Map<number, CartesianCoord>();
        const occupied = new Set<string>();
        const placedByCoord = new Map<string, number>();
        const orphaned: Array<{ helixId: number; parentId: number }> = [];
        const orphanedIds = new Set<number>();
        const relativePositions = getRelativePositions(helices);
        const { crossovers, helixIds } = collectCrossovers(grid);
        const helixCount = Math.max(
            helices.length,
            ...Array.from(grid.values()).map((mark) => mark.helixId + 1),
            ...Array.from(helixIds).map((helixId) => helixId + 1),
            0
        );

        const positionKey = (coord: CartesianCoord) => `${coord.x},${coord.y}`;
        const isOpen = (coord: CartesianCoord) => !occupied.has(positionKey(coord));
        const getOccupant = (coord: CartesianCoord) => placedByCoord.get(positionKey(coord));
        const tryPlace = (helixId: number, coord: CartesianCoord) => {
            if (placed.has(helixId) || !isOpen(coord)) return false;
            placed.set(helixId, coord);
            occupied.add(positionKey(coord));
            placedByCoord.set(positionKey(coord), helixId);
            return true;
        };

        const enqueueOrphan = (helixId: number, parentId: number) => {
            if (placed.has(helixId) || orphanedIds.has(helixId)) return;
            orphaned.push({ helixId, parentId });
            orphanedIds.add(helixId);
        };

        const buildAdjacency = () => {
            const adjacency = new Map<number, Set<number>>();

            const ensure = (helixId: number) => {
                if (!adjacency.has(helixId)) adjacency.set(helixId, new Set<number>());
                return adjacency.get(helixId)!;
            };

            for (let helixId = 0; helixId < helixCount; helixId++) ensure(helixId);

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
                if (!nt || !(nt instanceof Nucleotide)) continue;
                const pairNt = nt.pair;
                if (!pairNt || !(pairNt instanceof Nucleotide)) continue;
                const pairMark = grid.get(pairNt.id);
                if (!pairMark) continue;
                if (pairMark.helixId === mark.helixId) continue;
                ensure(mark.helixId).add(pairMark.helixId);
                ensure(pairMark.helixId).add(mark.helixId);
            }

            return adjacency;
        };

        const adjacency = buildAdjacency();

        // Spiral outward from `origin` in true Manhattan order (distance 1, 2, 3, …)
        // so orphans are always placed as close to their parent as possible.
        const findNearestOpenAround = (origin: CartesianCoord): CartesianCoord => {
            if (isOpen(origin)) return origin;

            const maxDist = helixCount + 32;
            for (let dist = 1; dist <= maxDist; dist++) {
                // Enumerate all integer coords at exact Manhattan distance = dist
                // and pick the first open one, preferring low |y| then low |x|
                // (i.e. stay close to the same row as origin first).
                let best: CartesianCoord | null = null;
                let bestKey = Infinity;

                for (let dx = -dist; dx <= dist; dx++) {
                    const dyAbs = dist - Math.abs(dx);
                    for (const dy of dyAbs === 0 ? [0] : [-dyAbs, dyAbs]) {
                        const cand = { x: origin.x + dx, y: origin.y + dy };
                        if (!isOpen(cand)) continue;
                        // Sort key: |dy| first (favour same row), then |dx|
                        const key = Math.abs(dy) * (maxDist * 2 + 1) + Math.abs(dx);
                        if (key < bestKey) { bestKey = key; best = cand; }
                    }
                }

                if (best) return best;
            }

            return { x: origin.x + maxDist + 1, y: origin.y };
        };

        const chooseRoot = () => {
            for (let helixId = 0; helixId < helixCount; helixId++) {
                if ((adjacency.get(helixId)?.size ?? 0) === 3) return helixId;
            }
            for (let helixId = 0; helixId < helixCount; helixId++) {
                if ((adjacency.get(helixId)?.size ?? 0) > 0) return helixId;
            }
            return 0;
        };

        const getSeedCoordinate = (helixId: number): CartesianCoord => {
            const degree = adjacency.get(helixId)?.size ?? 0;
            if (degree === 3) return { x: 0, y: 0 };
            // Do NOT use raw relativePositions as grid coords — those are in
            // simulation units (2-3 oxDNA units per helix spacing) and would
            // place the seed 10-30+ cells away from the origin.
            // Just find the nearest open cell to the origin.
            return findNearestOpenAround({ x: 0, y: 0 });
        };

        // Returns the nearest already-placed neighbor's grid coord, or falls
        // back to the given default.  Used to anchor orphan/straggler placement.
        const nearestPlacedNeighbor = (helixId: number, fallback: CartesianCoord): CartesianCoord => {
            for (const nbId of adjacency.get(helixId) ?? []) {
                const nb = placed.get(nbId);
                if (nb) return nb;
            }
            return fallback;
        };

        const placeNeighborGroup = (parentId: number, queue: number[]) => {
            const parentCoord = placed.get(parentId);
            if (!parentCoord) return;

            const neighbors = Array.from(adjacency.get(parentId) ?? [])
                .map((helixId) => ({
                    helixId,
                    delta: subtractPositions(relativePositions, helixId, parentId)
                }))
                .filter((entry): entry is { helixId: number; delta: { x: number; y: number } } => !!entry.delta)
                .sort((a, b) => {
                    const dyDelta = Math.abs(b.delta.y) - Math.abs(a.delta.y);
                    if (dyDelta !== 0) return dyDelta;
                    const dxDelta = Math.abs(b.delta.x) - Math.abs(a.delta.x);
                    if (dxDelta !== 0) return dxDelta;
                    return a.helixId - b.helixId;
                });

            if (neighbors.length === 0) return;

            const planned = new Set<number>();
            const assignedTargets = new Map<number, CartesianCoord>();

            // Only unplaced neighbors are candidates for new slots.
            const unplacedNeighbors = neighbors.filter((entry) => !placed.has(entry.helixId));

            // Assign Y only when at least one candidate is truly vertical,
            // i.e. |dy| > |dx|. Otherwise skip Y and keep placements on ±X.
            const verticalCandidates = unplacedNeighbors
                .filter((entry) => Math.abs(entry.delta.y) > Math.abs(entry.delta.x))
                .sort((a, b) => {
                    const dyDelta = Math.abs(b.delta.y) - Math.abs(a.delta.y);
                    if (dyDelta !== 0) return dyDelta;
                    const ratioA = Math.abs(a.delta.y) / (Math.abs(a.delta.x) + 1e-9);
                    const ratioB = Math.abs(b.delta.y) / (Math.abs(b.delta.x) + 1e-9);
                    if (ratioB !== ratioA) return ratioB - ratioA;
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
                const yStep: 1 | -1 = totalConnections <= 3
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

            const assignXSlot = (entry: { helixId: number; delta: { x: number; y: number } }, xStep: -1 | 1) => {
                const target = { x: parentCoord.x + xStep, y: parentCoord.y };
                assignedTargets.set(entry.helixId, target);
                planned.add(entry.helixId);
            };

            if (positive.length > 0 && negative.length > 0) {
                assignXSlot(positive[0], 1);
                assignXSlot(negative[0], -1);
            } else if (positive.length > 0) {
                assignXSlot(positive[0], 1);
            } else if (negative.length > 0) {
                assignXSlot(negative[0], -1);
            }

            for (const entry of remaining) {
                if (!planned.has(entry.helixId)) enqueueOrphan(entry.helixId, parentId);
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
                } else {
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

        const queue: number[] = [root];
        let qIdx = 0;

        while (qIdx < queue.length) {
            placeNeighborGroup(queue[qIdx++], queue);
        }

        for (const orphan of orphaned) {
            if (placed.has(orphan.helixId)) continue;
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
            if (placed.has(helixId)) continue;
            // Anchor stragglers to a placed neighbor, not raw relativePositions
            // (which are in simulation units and would scatter them far away).
            const anchor = nearestPlacedNeighbor(helixId, { x: 0, y: 0 });
            tryPlace(helixId, findNearestOpenAround(anchor));
        }

        const result = new Map<number, GridPosition>();
        for (let helixId = 0; helixId < helixCount; helixId++) {
            const coord = placed.get(helixId) ?? { x: helixId, y: 0 };
            result.set(helixId, [coord.x, coord.y]);
        }

        return result;
    }

    export function HelixPos(grid: GridMap, helices: Nucleotide[][]): Map<number, [number, number]> {
        type GridPosition = [number, number];
        type AxialCoord = { q: number; r: number };

        interface HexBasis3D {
            origin: THREE.Vector3;
            qVec: THREE.Vector3;
            rVec: THREE.Vector3;
            axisVec: THREE.Vector3;
        }

        const HEX_AXIAL_DIRS: AxialCoord[] = [
            { q: 1, r: 0 },
            { q: 1, r: -1 },
            { q: 0, r: -1 },
            { q: -1, r: 0 },
            { q: -1, r: 1 },
            { q: 0, r: 1 }
        ];

        const axialAdd = (a: AxialCoord, b: AxialCoord): AxialCoord => ({ q: a.q + b.q, r: a.r + b.r });

        const axialDistance = (a: AxialCoord, b: AxialCoord): number => {
            const dq = a.q - b.q;
            const dr = a.r - b.r;
            const ds = (-a.q - a.r) - (-b.q - b.r);
            return Math.max(Math.abs(dq), Math.abs(dr), Math.abs(ds));
        };

        const axialKey = (a: AxialCoord): string => `${a.q},${a.r}`;
        const edgeKey = (a: number, b: number): string => `${a}|${b}`;
        const axialRToOddQRow = (q: number, r: number): number => r + ((q - (q & 1)) / 2);

        // yes
        const helixCOM = (helix: Nucleotide[]): THREE.Vector3 | null => {
            if (!helix || helix.length === 0) return null;
            const sum = new THREE.Vector3();
            let count = 0;
            for (const nt of helix) {
                if (!(nt instanceof Nucleotide)) continue;
                sum.add(nt.getPos());
                count++;
            }
            if (count === 0) return null;
            return sum.divideScalar(count);
        };

        // Creates a coordinate system based on 
        const estimateHexBasis3D = (allHelices: Nucleotide[][]): HexBasis3D | null => {
            const centers: THREE.Vector3[] = [];
            const axisSamples: THREE.Vector3[] = [];

            // "PCA" is basically using endpoints 1 and 2 as a rough axis direction, and then averaging the centroids to find the origin.
            for (const helix of allHelices) {
                if (!helix || helix.length === 0) continue;
                const center = helixCOM(helix);
                if (center) centers.push(center);

                const ep = helixEndpoints(helix);
                if (!ep) continue;
                const axis = ep.end2.getPos().clone().sub(ep.end1.getPos());
                if (axis.lengthSq() > 1e-8) axisSamples.push(axis.normalize());
            }

            if (centers.length === 0) return null;

            // find the average center to use as the origin
            const origin = new THREE.Vector3();
            for (const c of centers) origin.add(c);
            origin.divideScalar(centers.length);

            // Average helix direction to get a rough axis. Not really for any rigorous reason, it works fine for this and thus, we will use it.
            let axisVec = new THREE.Vector3(0, 0, 1);
            if (axisSamples.length > 0) {
                axisVec.set(0, 0, 0);
                for (const s of axisSamples) axisVec.add(s);
                if (axisVec.lengthSq() < 1e-8) axisVec.set(0, 0, 1);
                else axisVec.normalize();
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

            const nearestDistances: number[] = [];
            for (let i = 0; i < projectedCenters.length; i++) {
                let nearest = Infinity;
                for (let j = 0; j < projectedCenters.length; j++) {
                    if (i === j) continue;
                    const d = projectedCenters[j].clone().sub(projectedCenters[i]).length();
                    if (d > 1e-6 && d < nearest) nearest = d;
                }
                if (nearest < Infinity) nearestDistances.push(nearest);
            }

            let spacing = 1;
            if (nearestDistances.length > 0) {
                nearestDistances.sort((a, b) => a - b);
                spacing = nearestDistances[Math.floor(nearestDistances.length / 2)] || 1;
                if (spacing <= 1e-6) spacing = 1;
            }

            const qVec = e1.clone().multiplyScalar(spacing);
            const rVec = e1.clone().multiplyScalar(0.5 * spacing)
                .add(e2.clone().multiplyScalar((Math.sqrt(3) / 2) * spacing));

            return { origin, qVec, rVec, axisVec };
        };

        const vectorToAxialContinuous = (v: THREE.Vector3, basis: HexBasis3D): AxialCoord | null => {
            const inPlane = v.clone().sub(basis.axisVec.clone().multiplyScalar(v.dot(basis.axisVec)));
            const aa = basis.qVec.dot(basis.qVec);
            const ab = basis.qVec.dot(basis.rVec);
            const bb = basis.rVec.dot(basis.rVec);
            const ap = basis.qVec.dot(inPlane);
            const bp = basis.rVec.dot(inPlane);
            const det = aa * bb - ab * ab;
            if (Math.abs(det) < 1e-10) return null;
            return {
                q: (ap * bb - bp * ab) / det,
                r: (bp * aa - ap * ab) / det
            };
        };

        const quantizeToHexDirection = (v: THREE.Vector3, basis: HexBasis3D): AxialCoord => {
            let best = HEX_AXIAL_DIRS[0];
            let bestScore = -Infinity;
            const dirVec = v.clone().sub(basis.axisVec.clone().multiplyScalar(v.dot(basis.axisVec)));
            if (dirVec.lengthSq() < 1e-12) return best;
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

        const collectHelixAdjacency = (currentGrid: GridMap): Map<number, Map<number, number>> => {
            const adjacency = new Map<number, Map<number, number>>();
            const ensure = (a: number, b: number) => {
                if (!adjacency.has(a)) adjacency.set(a, new Map<number, number>());
                const row = adjacency.get(a)!;
                row.set(b, (row.get(b) ?? 0) + 1);
            };

            const allNtIds = new Set<number>();
            for (const [ntId] of currentGrid.entries()) allNtIds.add(ntId);
            const visited = new Set<number>();

            for (const [ntId] of currentGrid.entries()) {
                if (visited.has(ntId)) continue;

                const startNt = elements.get(ntId) as Nucleotide | undefined;
                if (!startNt || !(startNt instanceof Nucleotide)) continue;

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

                let curr: Nucleotide | null = fivePrime;
                const walkForward = new Set<number>();
                let prevMark: GridMark | null = null;

                while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                    if (walkForward.has(curr.id)) break;
                    walkForward.add(curr.id);
                    visited.add(curr.id);

                    const mark = currentGrid.get(curr.id);
                    if (mark) {
                        if (prevMark && prevMark.helixId !== mark.helixId) {
                            ensure(prevMark.helixId, mark.helixId);
                            ensure(mark.helixId, prevMark.helixId);
                        }
                        prevMark = mark;
                    } else {
                        prevMark = null;
                    }

                    const n3ref: any = curr.n3;
                    curr = (n3ref && n3ref instanceof Nucleotide) ? (n3ref as Nucleotide) : null;
                }
            }

            return adjacency;
        };

        const nearestOpenAround = (target: AxialCoord, occupied: Set<string>, maxRadius: number = 16): AxialCoord => {
            if (!occupied.has(axialKey(target))) return target;

            for (let radius = 1; radius <= maxRadius; radius++) {
                let best: AxialCoord | null = null;
                let bestDist = Infinity;
                for (let dq = -radius; dq <= radius; dq++) {
                    const rMin = Math.max(-radius, -dq - radius);
                    const rMax = Math.min(radius, -dq + radius);
                    for (let dr = rMin; dr <= rMax; dr++) {
                        const cand = { q: target.q + dq, r: target.r + dr };
                        const k = axialKey(cand);
                        if (occupied.has(k)) continue;
                        const d = axialDistance(cand, target);
                        if (d < bestDist) {
                            bestDist = d;
                            best = cand;
                        }
                    }
                }
                if (best) return best;
            }

            return { q: target.q + maxRadius + 1, r: target.r };
        };

        const computeOddQGridPositions = (currentGrid: GridMap, allHelices: Nucleotide[][]): Map<number, GridPosition> => {
            const positions = new Map<number, GridPosition>();
            const helixCount = Math.max(0, ...Array.from(currentGrid.values()).map((m) => m.helixId + 1));
            if (helixCount === 0) return positions;

            const basis = estimateHexBasis3D(allHelices);
            const adjacency = collectHelixAdjacency(currentGrid);
            const centers = new Map<number, THREE.Vector3>();
            const projectedAxial = new Map<number, AxialCoord>();

            for (let h = 0; h < helixCount; h++) {
                const helix = allHelices[h] ?? [];
                const c = helixCOM(helix);
                if (!c) continue;
                centers.set(h, c);
                if (basis) {
                    const rel = c.clone().sub(basis.origin);
                    const uv = vectorToAxialContinuous(rel, basis);
                    if (uv) projectedAxial.set(h, uv);
                }
            }

            const preferredDelta = new Map<string, AxialCoord>();
            if (basis) {
                for (let a = 0; a < helixCount; a++) {
                    const row = adjacency.get(a);
                    if (!row) continue;
                    const ca = centers.get(a);
                    if (!ca) continue;
                    for (const [b] of row.entries()) {
                        const cb = centers.get(b);
                        if (!cb) continue;
                        const d = quantizeToHexDirection(cb.clone().sub(ca), basis);
                        preferredDelta.set(edgeKey(a, b), d);
                        preferredDelta.set(edgeKey(b, a), { q: -d.q, r: -d.r });
                    }
                }
            }

            const placed = new Map<number, AxialCoord>();
            const occupied = new Set<string>();
            const place = (helixId: number, coord: AxialCoord): boolean => {
                const k = axialKey(coord);
                if (occupied.has(k)) return false;
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
                if (placed.has(root)) continue;

                const preferredRoot = projectedAxial.get(root)
                    ? { q: Math.round(projectedAxial.get(root)!.q), r: Math.round(projectedAxial.get(root)!.r) }
                    : { q: 0, r: 0 };
                const rootCoord = (root === 0 && !occupied.has(axialKey({ q: 0, r: 0 })))
                    ? { q: 0, r: 0 }
                    : nearestOpenAround(preferredRoot, occupied);

                place(root, rootCoord);

                const queue: number[] = [root];
                let qi = 0;

                while (qi < queue.length) {
                    const current = queue[qi++];
                    const currentPos = placed.get(current);
                    if (!currentPos) continue;

                    const neighbors = Array.from((adjacency.get(current) ?? new Map<number, number>()).entries())
                        .sort((a, b) => b[1] - a[1])
                        .map(([id]) => id);

                    for (const nb of neighbors) {
                        if (placed.has(nb)) continue;

                        const pref = preferredDelta.get(edgeKey(current, nb)) ?? HEX_AXIAL_DIRS[0];
                        const base = axialAdd(currentPos, pref);

                        const candidates = HEX_AXIAL_DIRS
                            .map((d) => axialAdd(currentPos, d))
                            .sort((a, b) => {
                                const score = (coord: AxialCoord) => {
                                    let s = 0;
                                    if (coord.q === base.q && coord.r === base.r) s -= 5;

                                    const nbRow = adjacency.get(nb);
                                    if (nbRow) {
                                        for (const [p] of nbRow.entries()) {
                                            const placedP = placed.get(p);
                                            if (!placedP) continue;
                                            const pd = preferredDelta.get(edgeKey(p, nb));
                                            if (!pd) continue;
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
                if (placed.has(h)) continue;
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
    export function buildScadnano2(
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

        // ── Step 1: Discover all strands via backbone topology ──────────
        // Build a set of all nucleotide ids that exist in the grid so we
        // only emit nucleotides that were actually placed.
        const allNtIds = new Set<number>();
        for (const [ntId] of grid.entries()) {
            allNtIds.add(ntId);
        }

        // Track which nucleotides have been assigned to a strand already.
        const visited = new Set<number>();

        // We'll collect strand data here.
        const scadStrands: Array<{
            color: string;
            sequence: string;
            domains: Array<{ helix: number; forward: boolean; start: number; end: number }>;
            is_scaffold?: boolean;
        }> = [];

        // Iterate over every nucleotide in the grid and discover strands.
        for (const [ntId] of grid.entries()) {
            if (visited.has(ntId)) continue;

            const startNt = elements.get(ntId) as Nucleotide | undefined;
            if (!startNt || !(startNt instanceof Nucleotide)) continue;

            // ── 1a. Find the 5' end of this strand ──────────────────────
            // Walk n5 until we can't anymore (the node with no n5, or
            // whose n5 is not in the grid, is the 5' end).
            let fivePrime: Nucleotide = startNt;
            const walkBack = new Set<number>();
            walkBack.add(fivePrime.id);
            while (true) {
                const prev = fivePrime.n5;
                if (!prev || !(prev instanceof Nucleotide)) break;
                if (!allNtIds.has(prev.id)) break;   // not in grid
                if (walkBack.has(prev.id)) break;     // circular — stop
                walkBack.add(prev.id);
                fivePrime = prev;
            }

            // Detect circular: if fivePrime still has a valid n5 that
            // we stopped on because of the visited guard, it's circular.
            const n5OfFive = fivePrime.n5;
            const isCircular =
                n5OfFive instanceof Nucleotide &&
                allNtIds.has(n5OfFive.id) &&
                walkBack.has(n5OfFive.id);

            // ── 1b. Walk n3 from 5' end to build ordered nt list ────────
            const orderedNts: Nucleotide[] = [];
            let curr: Nucleotide | null = fivePrime;
            const walkForward = new Set<number>();

            while (curr && curr instanceof Nucleotide && allNtIds.has(curr.id)) {
                if (walkForward.has(curr.id)) break;  // full circle
                walkForward.add(curr.id);
                visited.add(curr.id);
                orderedNts.push(curr);
                const n3ref: any = curr.n3;
                curr = (n3ref && n3ref instanceof Nucleotide) ? (n3ref as Nucleotide) : null;
            }

            if (orderedNts.length === 0) continue;

            // ── Step 2: Split into domains by helixId AND direction ────
            // Each "run" is a maximal contiguous subsequence on the same
            // helix with the same grid direction. Break whenever either changes.
            type Run = { helixId: number; direction: 'forward' | 'backward'; nts: Nucleotide[] };
            const runs: Run[] = [];
            let currentRun: Run | null = null;

            for (const nt of orderedNts) {
                const mark = grid.get(nt.id);
                if (!mark) {
                    currentRun = null;
                    continue;
                }

                if (
                    currentRun &&
                    currentRun.helixId === mark.helixId &&
                    currentRun.direction === mark.direction
                ) {
                    currentRun.nts.push(nt);
                } else {
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
            const domains: Array<{ helix: number; forward: boolean; start: number; end: number }> = [];

            for (const run of runs) {
                // Collect (offset, base, nt) tuples in walk order
                const entries: { offset: number; base: string; nt: Nucleotide }[] = [];
                for (const nt of run.nts) {
                    const mark = grid.get(nt.id)!;
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
                type SubRun = typeof entries;
                const subRuns: SubRun[] = [];
                let currentSub: SubRun = [entries[0]];

                for (let i = 1; i < entries.length; i++) {
                    const prev = entries[i - 1].offset;
                    const curr = entries[i].offset;
                    if (curr === prev + step) {
                        currentSub.push(entries[i]);
                    } else {
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
                        end: maxOff + 1   // exclusive end
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

                const strandObj: {
                    color: string;
                    sequence: string;
                    domains: Array<{ helix: number; forward: boolean; start: number; end: number }>;
                    is_scaffold?: boolean;
                } = { color, sequence, domains };

                if (isScaffold) {
                    strandObj.is_scaffold = true;
                }
                if (isCircular) {
                    (strandObj as any).circular = true;
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

    // ── Relative Position Calculation Methods ─────────────────────────

    /**
     * Finds all crossovers from helix1 to helix2, computes the COM of the 4 nucleotides
     * involved in each crossover, and averages them to return a single 3D vector
     * representing the relative connection from helix1 to helix2.
     */
    export function getCrossoverVector(helix1: number, helix2: number, helices: Nucleotide[][]): THREE.Vector3 | null {
        // Collect all nucleotides in helix1 into a Set for fast lookup
        // const h1Set = new Set(helices[helix1].map(n => n.id));
        const h2Set = new Set(helices[helix2].map(n => n.id));

        const crossoverVectors: THREE.Vector3[] = [];

        // Scan all nucleotides in helix 1 to find connections to helix 2
        for (const n1 of helices[helix1]) {
            // Check 5' backbone connection
            if (n1.n5 && n1.n5 instanceof Nucleotide && h2Set.has(n1.n5.id)) {
                // We found a backbone step from helix1 to helix2!
                const n2 = n1.n5 as Nucleotide;
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
                const n2 = n1.n3 as Nucleotide;
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

        if (crossoverVectors.length === 0) return null;

        // Average all crossover displacement vectors into a final definitive step vector
        const avgVector = new THREE.Vector3(0, 0, 0);
        for (const v of crossoverVectors) {
            avgVector.add(v);
        }
        avgVector.divideScalar(crossoverVectors.length);

        return avgVector;
    }

    export function getCrossoverVectorWithBinders(helix1: number, helix2: number, helices: Nucleotide[][]): THREE.Vector3 | null {
        // Collect all nucleotides in helix1 into a Set for fast lookup
        // const h1Set = new Set(helices[helix1].map(n => n.id));
        const h2Set = new Set(helices[helix2].map(n => n.id));

        const crossoverVectors: THREE.Vector3[] = [];

        // For binder-like regions, pair may be missing on one/both sides.
        // Use the nucleotide position itself as fallback so these links still
        // contribute to relative placement.
        const localCenter = (nt: Nucleotide): THREE.Vector3 => {
            const p = nt.pair;
            if (p && p instanceof Nucleotide) {
                return nt.getPos().clone().add(p.getPos()).multiplyScalar(0.5);
            }
            return nt.getPos().clone();
        };

        const pushVector = (n1: Nucleotide, n2: Nucleotide) => {
            const c1 = localCenter(n1);
            const c2 = localCenter(n2);
            crossoverVectors.push(new THREE.Vector3().subVectors(c2, c1));
        };

        // Scan all nucleotides in helix 1 to find connections to helix 2
        for (const n1 of helices[helix1]) {
            // Check 5' backbone connection
            if (n1.n5 && n1.n5 instanceof Nucleotide && h2Set.has(n1.n5.id)) {
                const n2 = n1.n5 as Nucleotide;
                pushVector(n1, n2);
            }

            // Check 3' backbone connection
            if (n1.n3 && n1.n3 instanceof Nucleotide && h2Set.has(n1.n3.id)) {
                const n2 = n1.n3 as Nucleotide;
                pushVector(n1, n2);
            }

            // Check direct pair bridge across helices. This captures binder-like
            // attachments where backbone crossover signatures are sparse/absent.
            if (n1.pair && n1.pair instanceof Nucleotide && h2Set.has(n1.pair.id)) {
                const n2 = n1.pair as Nucleotide;
                pushVector(n1, n2);
            }
            // same as running an average over all 4 nucleotides and then running an average over THOSE vectors
        }

        if (crossoverVectors.length === 0) return null;

        // Average all crossover displacement vectors into a final definitive step vector
        const avgVector = new THREE.Vector3(0, 0, 0);
        for (const v of crossoverVectors) {
            avgVector.add(v);
        }
        avgVector.divideScalar(crossoverVectors.length);

        return avgVector;
    }

    /**
     * Traverses the connections starting from Helix 0, and plots every connected helix
     * onto a 2D coordinate plane locally aligned relative to Helix 0's axis.
     * Returns a Map of HelixId -> { x, y }
     */
    export function getRelativePositions(helices: Nucleotide[][]): Map<number, { x: number, y: number }> {
        const positions = new Map<number, { x: number, y: number }>();
        const visited = new Set<number>();
        const queue: number[] = [];

        // Find anchor (Helix 0)
        positions.set(0, { x: 0, y: 0 });
        visited.add(0);
        queue.push(0);

        // 1. Helix long-axis: the direction the helices run along (their "Z").
        const endPts = helixEndpoints(helices[0]);
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
        const deriveU0 = (axis: THREE.Vector3): THREE.Vector3 => {
            const candidates = [
                new THREE.Vector3(1, 0, 0),
                new THREE.Vector3(0, 1, 0),
                new THREE.Vector3(0, 0, 1),
            ];
            for (const c of candidates) {
                const proj = c.clone().projectOnPlane(axis);
                if (proj.lengthSq() > 0.01) return proj.normalize();
            }
            return new THREE.Vector3(1, 0, 0); // degenerate fallback
        };
        const u0: THREE.Vector3 = deriveU0(longAxis);
        const u1: THREE.Vector3 = new THREE.Vector3().crossVectors(longAxis, u0).normalize();

        // BFS
        while (queue.length > 0) {
            const curr = queue.shift()!;
            const currPos = positions.get(curr)!;

            for (let i = 0; i < helices.length; i++) {
                if (i === curr) continue;

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
}
