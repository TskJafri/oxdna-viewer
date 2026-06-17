/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />

namespace toscad {

    type LatticeKind = 'honeycomb' | 'square';
    type PhaseParity = 0 | 1;

    interface LatticePhaseConfig {
        basesPerTurn: number;
        phases: Record<PhaseParity, number[]>;
        voidPhase: number | null;
        tieEpsilon: number;
    }


    // easy to use constant for later reference.
    const resolveLatticeKind = (lattice: string): LatticeKind =>
        (lattice ?? '').toLowerCase() === 'square' ? 'square' : 'honeycomb';

    const LATTICE_CONFIG: Record<LatticeKind, LatticePhaseConfig> = {
        honeycomb: {
            basesPerTurn: 10.5,
            phases: {
                0: [0.0, 3.5, 7.0, 10.5],
                // 1: [1.75, 5.25, 8.75]
                // 1: [2.77, 6.27, 9.77]
                // numbers found by using 215/145 asymmetry in major/minor groove and allowing void phase to break ties.
                1: [0.7292, 4.229, 7.7292]
            },
            voidPhase: 5.25,
            tieEpsilon: 0.01
        },
        square: {
            basesPerTurn: 32 / 3,
            phases: {
                // Symmetric and therefore simple.
                0: [0.0, 8 / 3, 16 / 3, 8.0, 32 / 3],
                1: [4 / 3, 4.0, 20 / 3, 28 / 3]
            },
            voidPhase: null,
            tieEpsilon: 0.01
        }
    };

    export function findHelixID(targetId: number, helices: Nucleotide[][]): number | null {
        for (let i = 0; i < helices.length; i++) {
            const helix = helices[i];
            for (const nt of helix) {
                if (nt.id === targetId) return i;
            }
        }
        return null;
    }

    // helper function
    // Finds the local angle distribution of a given helix relative to its neighbors.
    // Uses crossover offsets defined by LATTICE_CONFIG to determine ideal angle buckets. 
    function getAngleHelix(
        grid: GridMap,
        helices: Nucleotide[][],
        helixId: number,
        lattice: string
    ): Map<number, {
        helixId: number;
        adj_helix: number;
        angle: number;
    }> {
        void helices;

        const result = new Map<number, {
            helixId: number;
            adj_helix: number;
            angle: number;
        }>();

        const latticeType = resolveLatticeKind(lattice);
        const latticeConfig = LATTICE_CONFIG[latticeType];

        const nearestPhase = (value: number, parity: PhaseParity) => {
            const candidates = latticeConfig.phases[parity] ?? [];
            if (!candidates.length) return value;

            return candidates.reduce((prev, curr) => {
                const distPrev = Math.abs(prev - value);
                const distCurr = Math.abs(curr - value);

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

        type HubCrossover = {
            adj_helix: number;
            offset: number;
            direction: 'forward' | 'backward';
        };

        const hubCrossovers: HubCrossover[] = [];
        for (const crossover of crossoverNts(grid)) {
            if (crossover.fromHelix === helixId) {
                const hubMark = grid.get(crossover.fromNt.id);
                if (!hubMark) continue;
                hubCrossovers.push({
                    adj_helix: crossover.toHelix,
                    offset: crossover.fromOffset,
                    direction: hubMark.direction
                });
            } else if (crossover.toHelix === helixId) {
                const hubMark = grid.get(crossover.toNt.id);
                if (!hubMark) continue;
                hubCrossovers.push({
                    adj_helix: crossover.fromHelix,
                    offset: crossover.toOffset,
                    direction: hubMark.direction
                });
            }
        }

        if (!hubCrossovers.length) return result;

        const groupedByNeighbor = new Map<number, HubCrossover[]>();
        for (const crossover of hubCrossovers) {
            if (!groupedByNeighbor.has(crossover.adj_helix)) {
                groupedByNeighbor.set(crossover.adj_helix, []);
            }
            groupedByNeighbor.get(crossover.adj_helix)!.push(crossover);
        }

        const neighbors = Array.from(groupedByNeighbor.keys()).sort((a, b) => a - b);
        if (!neighbors.length) return result;

        const pairTallies = new Map<string, Map<number, number>>();
        for (let i = 0; i < neighbors.length; i++) {
            const neighborA = neighbors[i];
            const groupA = groupedByNeighbor.get(neighborA) ?? [];

            for (let j = i + 1; j < neighbors.length; j++) {
                const neighborB = neighbors[j];
                const groupB = groupedByNeighbor.get(neighborB) ?? [];
                const pairKey = `${neighborA}|${neighborB}`;
                const bucket = new Map<number, number>();

                for (const crossoverA of groupA) {
                    for (const crossoverB of groupB) {
                        const offsetA = crossoverA.offset;
                        const dirA = crossoverA.direction;
                        const offsetB = crossoverB.offset;
                        const dirB = crossoverB.direction;

                        // Here's the math i spent so long figuring out:
                        const rawX = offsetB - offsetA;
                        const y: PhaseParity = dirA === dirB ? 0 : 1;

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

        const pairConsensus = new Map<string, number>();
        for (const [pairKey, bucket] of pairTallies.entries()) {
            if (!bucket.size) continue;

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

    // Now run getAngleHelix for every helix to get a full network map.
    // Note that every helix is has RELATIVE angles to its neighbors, they are not globally aligned to anything yet.
    // The global alignment is done in calculateGridPositions. 
    export function getAngles(grid: GridMap, helices: Nucleotide[][], lattice: string = 'honeycomb'): Map<number, Map<number, number>> {
        const networkMap = new Map<number, Map<number, number>>();
        const helixIds = new Set<number>();

        for (const [, mark] of grid.entries()) {
            helixIds.add(mark.helixId);
        }

        const sortedhids = Array.from(helixIds).sort((a, b) => a - b);
        for (const currentHID of sortedhids) {
            const helixAngles = getAngleHelix(grid, helices, currentHID, lattice);
            const angleMap = new Map<number, number>();

            for (const [adjHelixId, angleInfo] of helixAngles.entries()) {
                angleMap.set(adjHelixId, angleInfo.angle);
            }

            networkMap.set(currentHID, angleMap);
        }

        return networkMap;
    }

    // For every helix, count how many backbone crossovers it shares with each neighbor. Useful to remove bad combinations (such as end-only, which happens when a helix is broken in 2 pieces)
    export function getConnectionCounts(grid: GridMap): Map<number, Map<number, number>> {
        const counts = new Map<number, Map<number, number>>();

        const bump = (a: number, b: number) => {
            if (!counts.has(a)) counts.set(a, new Map<number, number>());
            const inner = counts.get(a)!;
            inner.set(b, (inner.get(b) ?? 0) + 1);
        };

        for (const crossover of crossoverNts(grid)) {
            if (crossover.fromHelix === crossover.toHelix) continue;
            bump(crossover.fromHelix, crossover.toHelix);
            bump(crossover.toHelix, crossover.fromHelix);
        }

        return counts;
    }

    // helper function to check for angle collisions in the map.
    export function angleCollisions(networkMap: Map<number, Map<number, number>>) {
        const overlappingHelices = [];

        for (const [helixId, angleMap] of networkMap.entries()) {
            if (!angleMap || angleMap.size === 0) continue;

            const anglesToNeighbors = new Map<number, number[]>();
            for (const [adjHelix, angle] of angleMap.entries()) {
                if (!anglesToNeighbors.has(angle)) anglesToNeighbors.set(angle, []);
                anglesToNeighbors.get(angle)!.push(adjHelix);
            }

            const localConflicts = [];
            for (const [angle, collidedHelices] of anglesToNeighbors.entries()) {
                const uniqueCollided = Array.from(new Set(collidedHelices)).sort((a, b) => a - b);
                if (uniqueCollided.length < 2) continue;
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

    // helper function. Checks if 2 helices are mutually disjoint (through offsets in the grid). 
    // typically used AFTER aligning the grid otherwise this is nonsense.
    export function disjoint(helices: Nucleotide[][], h1: number, h2: number, grid: GridMap): boolean {
        const buildSignedOffsets = (helixId: number): Set<number> => {
            const signedOffsets = new Set<number>();
            const nts = helices[helixId] ?? [];

            for (const nt of nts) {
                const mark = grid.get(nt.id);
                if (!mark || mark.helixId !== helixId) continue;

                const x = mark.offset;
                const y = mark.direction === 'backward' ? 1 : 0;
                signedOffsets.add((2 * x) + y);
            }

            return signedOffsets;
        };

        const set1 = buildSignedOffsets(h1);
        const set2 = buildSignedOffsets(h2);

        if (set1.size === 0 || set2.size === 0) return false;

        const smaller = set1.size <= set2.size ? set1 : set2;
        const larger = set1.size <= set2.size ? set2 : set1;

        for (const offs of smaller) {
            if (larger.has(offs)) return false;
        }

        return true;
    }

    /* If 2 helices have the same angle as given by getAngleHelix, then check if they are disjoint (big)? If yes, then merge them.
    If they have the same angle and are NOT disjoint, then check anglecorr() function.
    */
    export function anglecomb(grid: GridMap, helices: Nucleotide[][], lattice: string = 'honeycomb',
        angleMap: Map<number, Map<number, number>> = getAngles(grid, helices, lattice)
    ): {
        networkMap: Map<number, Map<number, number>>;
        mergedPairs: Array<{ sourceHelix: number; keepHelix: number; mergedHelix: number; angle: number }>;
    } {
        // angle requirement for colinearity; this is so they don't merge some extremely non-linear axes.
        const AXIS_DOT_THRESHOLD = 0.935;

        let networkMap = angleMap;
        const mergedPairs: Array<{ sourceHelix: number; keepHelix: number; mergedHelix: number; angle: number }> = [];

        const mergeHelixInto = (keepHelix: number, mergedHelix: number) => {
            const remapHelixId = (helixId: number): number => {
                if (helixId === mergedHelix) return keepHelix;
                if (helixId > mergedHelix) return helixId - 1;
                return helixId;
            };

            for (const [, mark] of grid.entries()) {
                mark.helixId = remapHelixId(mark.helixId);
            }

            const keepNts = helices[keepHelix] ?? [];
            const mergedNts = helices[mergedHelix] ?? [];
            helices[keepHelix] = keepNts.concat(mergedNts);

            // Remove the merged helix node and keep helix indexing compact.
            helices.splice(mergedHelix, 1);
        };

        // Diagnostic helper: summarise a helix's grid footprint. Useful in
        // the merge / rejection logs so we can tell at a glance whether two
        // helices are plausibly pieces of the same physical strand or are
        // far apart in offset space (= probably not the same helix).
        const summarizeHelix = (helixId: number): string => {
            const nts = helices[helixId] ?? [];
            let minOffset = Number.POSITIVE_INFINITY;
            let maxOffset = Number.NEGATIVE_INFINITY;
            let fwdCount = 0;
            let backCount = 0;
            let counted = 0;

            for (const nt of nts) {
                if (!(nt instanceof Nucleotide)) continue;
                const mark = grid.get(nt.id);
                if (!mark || mark.helixId !== helixId) continue;

                if (mark.offset < minOffset) minOffset = mark.offset;
                if (mark.offset > maxOffset) maxOffset = mark.offset;
                if (mark.direction === 'forward') fwdCount++;
                else backCount++;
                counted++;
            }

            if (counted === 0) return `[empty]`;
            return `[offsets=${minOffset}..${maxOffset}, fwd=${fwdCount}, back=${backCount}, total=${counted}]`;
        };

        // Helix axis (normalised). Picks the nucleotide at the helix's min
        // and max grid offsets and returns (pos(min) - pos(max)).normalize().
        // Two helices that are physically pieces of the same line have
        // |axisA · axisB| ≈ 1; two unrelated helices oriented differently
        // give a much smaller dot. Returns null when the helix is a single
        // point or otherwise has no definable axis.
        const helixAxis = (helixId: number): THREE.Vector3 | null => {
            const nts = helices[helixId] ?? [];
            let minOffset = Number.POSITIVE_INFINITY;
            let maxOffset = Number.NEGATIVE_INFINITY;
            let minNt: Nucleotide | null = null;
            let maxNt: Nucleotide | null = null;

            for (const nt of nts) {
                if (!(nt instanceof Nucleotide)) continue;
                const mark = grid.get(nt.id);
                if (!mark || mark.helixId !== helixId) continue;

                if (mark.offset < minOffset) {
                    minOffset = mark.offset;
                    minNt = nt;
                }
                if (mark.offset > maxOffset) {
                    maxOffset = mark.offset;
                    maxNt = nt;
                }
            }

            if (!minNt || !maxNt || minNt === maxNt) return null;

            const a = minNt.getPos();
            const b = maxNt.getPos();
            const dir = a.clone().sub(b);
            const len = dir.length();
            if (!isFinite(len) || len === 0) return null;
            return dir.divideScalar(len);
        };

        let pass = 0;
        while (pass++ < 200) {
            const collisionReports = angleCollisions(networkMap);

            if (collisionReports.length === 0) {
                console.log(`[anglecomb] pass=${pass} no collisions remain`);
                break;
            }

            // Rebuild crossover-count map each pass — merges shift helix IDs.
            const connectionCounts = getConnectionCounts(grid);

            let mergedInThisPass = false;

            outer:
            for (const report of collisionReports) {
                const sourceHelix = report.helixId;
                const sourceCounts = connectionCounts.get(sourceHelix);

                for (const conflict of report.conflicts) {
                    const collidedHelices = conflict.colliding_adj_helices
                        .filter((helixId) => helixId !== sourceHelix)
                        .sort((a, b) => a - b);

                    if (collidedHelices.length < 2) continue;

                    for (let i = 0; i < collidedHelices.length; i++) {
                        const helixA = collidedHelices[i];

                        for (let j = i + 1; j < collidedHelices.length; j++) {
                            const helixB = collidedHelices[j];

                            // Gate: angle inference between A and B rests on
                            // crossovers source→A × source→B. If either side has
                            // <2 crossovers, the consensus is built off a single
                            // (possibly accidental) link — discard the merge.
                            const countA = sourceCounts?.get(helixA) ?? 0;
                            const countB = sourceCounts?.get(helixB) ?? 0;
                            const directAB = connectionCounts.get(helixA)?.get(helixB) ?? 0;
                            if (countA < 2 || countB < 2) {
                                console.log(
                                    `[anglecomb] rejected potential merge: source=${sourceHelix}, ` +
                                    `helixA=${helixA} ${summarizeHelix(helixA)} (crossovers=${countA}), ` +
                                    `helixB=${helixB} ${summarizeHelix(helixB)} (crossovers=${countB}), ` +
                                    `directAB=${directAB}, angle=${conflict.angle} — weakConnection`
                                );
                                continue;
                            }

                            // Gate: helices being merged should be roughly
                            // collinear in 3D. Two pieces of the same broken
                            // helix share an axis (parallel or antiparallel,
                            // hence absolute dot). Unrelated helices typically
                            // have very different axes.
                            const axisA = helixAxis(helixA);
                            const axisB = helixAxis(helixB);
                            const axisDot = (axisA && axisB) ? Math.abs(axisA.dot(axisB)) : NaN;
                            if (!axisA || !axisB || axisDot < AXIS_DOT_THRESHOLD) {
                                console.log(
                                    `[anglecomb] rejected potential merge: source=${sourceHelix}, ` +
                                    `helixA=${helixA} ${summarizeHelix(helixA)}, ` +
                                    `helixB=${helixB} ${summarizeHelix(helixB)}, ` +
                                    `|axisA·axisB|=${isFinite(axisDot) ? axisDot.toFixed(3) : 'undefined'}, ` +
                                    `angle=${conflict.angle} — axesNotCollinear`
                                );
                                continue;
                            }

                            if (!disjoint(helices, helixA, helixB, grid)) continue;

                            const keepHelix = Math.min(helixA, helixB);
                            const mergedHelix = Math.max(helixA, helixB);

                            // Capture summaries BEFORE the merge mutates helices[].
                            const summaryKeep = summarizeHelix(keepHelix);
                            const summaryMerged = summarizeHelix(mergedHelix);

                            mergeHelixInto(keepHelix, mergedHelix);

                            mergedPairs.push({
                                sourceHelix,
                                keepHelix,
                                mergedHelix,
                                angle: conflict.angle
                            });

                            console.log(
                                `[anglecomb] combined helix ${mergedHelix} ${summaryMerged} into ${keepHelix} ${summaryKeep} ` +
                                `(source=${sourceHelix}, angle=${conflict.angle}, ` +
                                `count(source→keep)=${sourceCounts?.get(keepHelix) ?? 0}, ` +
                                `count(source→merged)=${sourceCounts?.get(mergedHelix) ?? 0}, ` +
                                `directAB=${directAB}, |axisDot|=${axisDot.toFixed(3)})`
                            );

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

    // If 2 helices have the same angle as given by getAngleHelix, and they are NOT disjoint, then MOVE one of the helices to an empty spot. For them to physically exist, one of these 2 checks must be true.
    export function anglecorr(grid: GridMap, helices: Nucleotide[][], lattice: string = 'honeycomb',
        angleMap: Map<number, Map<number, number>> = getAngles(grid, helices, lattice)
    ): {
        networkMap: Map<number, Map<number, number>>;
        correctedPairs: Array<{
            sourceHelix: number;
            adjustedHelix: number;
            oldAngle: number;
            newAngle: number;
            conflictAngle: number;
        }>;
    } {
        let networkMap = angleMap;
        const correctedPairs: Array<{
            sourceHelix: number;
            adjustedHelix: number;
            oldAngle: number;
            newAngle: number;
            conflictAngle: number;
        }> = [];

        const normalizeAngle = (angle: number) => ((angle % 360) + 360) % 360;
        const latticeType = resolveLatticeKind(lattice);
        const correctionStep = latticeType === 'square' ? 90 : 120;
        const maxLatticeAngles = latticeType === 'square' ? 4 : 3;

        let pass = 0;
        while (pass++ < 200) {
            const collisionReports = angleCollisions(networkMap);

            if (collisionReports.length === 0) {
                console.log(`[anglecorr] pass=${pass} no collisions remain`);
                break;
            }

            let correctedInThisPass = false;

            outer:
            for (const report of collisionReports) {
                const sourceHelix = report.helixId;
                const sourceMap = networkMap.get(sourceHelix);
                if (!sourceMap || sourceMap.size === 0) continue;

                const baseHelices = Array.from(sourceMap.entries())
                    .filter(([, angle]) => normalizeAngle(angle) === 0)
                    .map(([adjHelix]) => adjHelix);
                const baseHelix = baseHelices.length > 0 ? Math.min(...baseHelices) : null;

                for (const conflict of report.conflicts) {
                    const collidedHelices = conflict.colliding_adj_helices
                        .filter((helixId) => helixId !== sourceHelix)
                        .sort((a, b) => a - b);

                    if (collidedHelices.length < 2) continue;

                    for (let i = 0; i < collidedHelices.length; i++) {
                        const helixA = collidedHelices[i];

                        for (let j = i + 1; j < collidedHelices.length; j++) {
                            const helixB = collidedHelices[j];

                            if (disjoint(helices, helixA, helixB, grid)) continue;

                            const connectionCount = sourceMap.size;
                            if (connectionCount > maxLatticeAngles) {
                                console.warn(
                                    `[anglecorr] source=${sourceHelix} has ${connectionCount} connections and non-disjoint ` +
                                    `collision at angle=${conflict.angle}; no empty ${latticeType} slots available`
                                );
                                continue;
                            }

                            const pairDescending = [helixA, helixB].sort((a, b) => b - a);
                            let adjustedHelix = pairDescending.find((hId) => hId !== baseHelix) ?? pairDescending[0];

                            const oldAngleValue = sourceMap.get(adjustedHelix);
                            if (oldAngleValue === undefined) continue;

                            const oldAngle = normalizeAngle(oldAngleValue);

                            const usedAngles = new Set<number>();
                            for (const [adjHelix, angle] of sourceMap.entries()) {
                                if (adjHelix === adjustedHelix) continue;
                                usedAngles.add(normalizeAngle(angle));
                            }

                            let newAngle = oldAngle;
                            let foundOpenAngle = false;
                            for (let attempt = 1; attempt < maxLatticeAngles; attempt++) {
                                const candidate = normalizeAngle(oldAngle + (attempt * correctionStep));
                                if (!usedAngles.has(candidate)) {
                                    newAngle = candidate;
                                    foundOpenAngle = true;
                                    break;
                                }
                            }

                            if (!foundOpenAngle) {
                                console.warn(
                                    `[anglecorr] source=${sourceHelix} could not place corrected angle for helix=${adjustedHelix}; ` +
                                    `all ${correctionStep}-degree alternatives occupied`
                                );
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

                            console.log(
                                `[anglecorr] corrected source=${sourceHelix}: helix ${adjustedHelix} angle ${oldAngle} -> ${newAngle} ` +
                                `(conflict angle=${conflict.angle})`
                            );

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

    type CrossoverWeightMap = Map<number, Map<number, number>>;
    type CrossoverWeightEdge =
        | { from: number; to: number; weight: number }
        | { fromHelix: number; toHelix: number; weight: number };

    /**
     * Calculates absolute grid coordinates from local helix-to-helix angles.
     *
     * Phase 1: Build a strict lattice spanning tree from helix 0 using
     * weighted BFS (top-3 children only at each node).
     *
     * Phase 2: Place deferred/artefact helices in nearest open coordinates
     * around their parent once the phase-1 core is locked.
     */
    export function calculateGlobalPositions(
        networkMap: Map<number, Map<number, number>>,
        crossoverWeights?: CrossoverWeightMap | CrossoverWeightEdge[],
        options?: {
            runDeferredSubBfs?: boolean;
            maxSearchRadius?: number;
        },
        lattice: string = 'honeycomb'
    ): Map<number, [number, number]> {
        type GridCoord = { col: number; row: number };
        type QueueNode = { helixId: number; offset: number };
        type DeferredItem = { parentId: number; helixId: number; inheritedOffset: number };
        type Parity = 'even' | 'odd';
        type HoneycombAngle = 0 | 120 | 240;
        type SquareAngle = 0 | 90 | 180 | 270;
        type LatticeAngle = HoneycombAngle | SquareAngle;

        const latticeType = resolveLatticeKind(lattice);

        const HONEYCOMB_ANGLES: HoneycombAngle[] = [0, 120, 240];
        const SQUARE_ANGLES: SquareAngle[] = [0, 90, 180, 270];
        const ANGLES: LatticeAngle[] = latticeType === 'square'
            ? SQUARE_ANGLES.slice()
            : HONEYCOMB_ANGLES.slice();
        const maxChildrenPerNode = ANGLES.length;

        const HONEYCOMB_STEP_BY_PARITY: Record<Parity, Record<HoneycombAngle, { dCol: number; dRow: number }>> = {
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

        const SQUARE_STEPS: Record<SquareAngle, { dCol: number; dRow: number }> = {
            0: { dCol: 1, dRow: 0 },
            90: { dCol: 0, dRow: 1 },
            180: { dCol: -1, dRow: 0 },
            270: { dCol: 0, dRow: -1 }
        };

        const normalizeAngle = (angle: number) => ((angle % 360) + 360) % 360;
        const angleDistance = (a: number, b: number) => {
            const diff = Math.abs(normalizeAngle(a) - normalizeAngle(b));
            return Math.min(diff, 360 - diff);
        };

        const parityAt = (coord: GridCoord): Parity =>
            (((coord.col + coord.row) & 1) === 0 ? 'even' : 'odd');

        const getStep = (coord: GridCoord, angle: LatticeAngle): { dCol: number; dRow: number } => {
            if (latticeType === 'square') {
                return SQUARE_STEPS[angle as SquareAngle];
            }

            return HONEYCOMB_STEP_BY_PARITY[parityAt(coord)][angle as HoneycombAngle];
        };

        const snapToLatticeAngle = (angle: number): LatticeAngle => {
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

        const latticeAngleFromDelta = (coord: GridCoord, dCol: number, dRow: number): LatticeAngle | null => {
            for (const angle of ANGLES) {
                const delta = getStep(coord, angle);
                if (delta.dCol === dCol && delta.dRow === dRow) return angle;
            }
            return null;
        };

        const weightMap: CrossoverWeightMap = new Map();
        const ensureWeightRow = (helixId: number): Map<number, number> => {
            if (!weightMap.has(helixId)) weightMap.set(helixId, new Map<number, number>());
            return weightMap.get(helixId)!;
        };
        const setWeight = (from: number, to: number, weight: number) => {
            ensureWeightRow(from).set(to, weight);
        };

        if (crossoverWeights instanceof Map) {
            for (const [from, row] of crossoverWeights.entries()) {
                for (const [to, weight] of row.entries()) {
                    setWeight(from, to, weight);
                }
            }
        } else if (Array.isArray(crossoverWeights)) {
            for (const edge of crossoverWeights) {
                if ('fromHelix' in edge && 'toHelix' in edge) {
                    setWeight(edge.fromHelix, edge.toHelix, edge.weight);
                } else {
                    setWeight(edge.from, edge.to, edge.weight);
                }
            }
        }

        const getWeight = (a: number, b: number): number => {
            const direct = weightMap.get(a)?.get(b);
            if (direct !== undefined) return direct;

            const reverse = weightMap.get(b)?.get(a);
            if (reverse !== undefined) return reverse;

            return 1;
        };

        const allHelixIds = new Set<number>();
        allHelixIds.add(0);

        for (const [from, row] of networkMap.entries()) {
            allHelixIds.add(from);
            for (const to of row.keys()) allHelixIds.add(to);
        }
        for (const [from, row] of weightMap.entries()) {
            allHelixIds.add(from);
            for (const to of row.keys()) allHelixIds.add(to);
        }

        const positions = new Map<number, GridCoord>();
        const occupied = new Map<string, number>();
        const globalRotationOffsets = new Map<number, number>();
        const deferredQueue: DeferredItem[] = [];
        const deferredSeen = new Set<string>();

        const keyOf = (coord: GridCoord) => `${coord.col},${coord.row}`;

        const enqueueDeferred = (parentId: number, helixId: number, inheritedOffset: number) => {
            const key = `${parentId}|${helixId}`;
            if (deferredSeen.has(key) || positions.has(helixId)) return;
            deferredSeen.add(key);
            deferredQueue.push({ parentId, helixId, inheritedOffset });
        };

        const maxSearchRadius = Math.max(1, options?.maxSearchRadius ?? 256);

        const latticeNeighbors = (coord: GridCoord): GridCoord[] => {
            return ANGLES.map((angle) => {
                const step = getStep(coord, angle);
                return {
                    col: coord.col + step.dCol,
                    row: coord.row + step.dRow
                };
            });
        };

        const findNearestOpen = (anchor: GridCoord): GridCoord => {
            if (!occupied.has(keyOf(anchor))) return { col: anchor.col, row: anchor.row };

            const visited = new Set<string>();
            const queue: Array<{ coord: GridCoord; dist: number }> = [{ coord: anchor, dist: 0 }];
            visited.add(keyOf(anchor));

            let qIdx = 0;
            while (qIdx < queue.length) {
                const { coord, dist } = queue[qIdx++];
                if (dist >= maxSearchRadius) continue;

                const neighbors = latticeNeighbors(coord).sort((a, b) => {
                    const da = Math.abs(a.col - anchor.col) + Math.abs(a.row - anchor.row);
                    const db = Math.abs(b.col - anchor.col) + Math.abs(b.row - anchor.row);
                    if (da !== db) return da - db;
                    if (a.col !== b.col) return a.col - b.col;
                    return a.row - b.row;
                });

                for (const next of neighbors) {
                    const nextKey = keyOf(next);
                    if (visited.has(nextKey)) continue;
                    visited.add(nextKey);

                    if (!occupied.has(nextKey)) return next;

                    queue.push({ coord: next, dist: dist + 1 });
                }
            }

            const fallback = { col: anchor.col + maxSearchRadius + 1, row: anchor.row };
            while (occupied.has(keyOf(fallback))) {
                fallback.col += 1;
            }
            return fallback;
        };

        const computeChildOffset = (
            parentId: number,
            childId: number,
            parentCoord: GridCoord,
            childCoord: GridCoord,
            parentLocalAngle: number,
            snappedForwardAngle: LatticeAngle
        ): number => {
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

        const processNode = (node: QueueNode, targetQueue: QueueNode[]) => {
            const parentCoord = positions.get(node.helixId);
            if (!parentCoord) return;

            const localEdges = networkMap.get(node.helixId);
            if (!localEdges || localEdges.size === 0) return;

            const candidates = Array.from(localEdges.entries())
                .filter(([neighborId]) => neighborId !== node.helixId && !positions.has(neighborId))
                .map(([neighborId, localAngle]) => ({
                    neighborId,
                    localAngle: normalizeAngle(localAngle),
                    weight: getWeight(node.helixId, neighborId)
                }))
                .sort((a, b) => {
                    if (b.weight !== a.weight) return b.weight - a.weight;
                    return a.neighborId - b.neighborId;
                });

            const selected = candidates.slice(0, maxChildrenPerNode);
            const overflow = candidates.slice(maxChildrenPerNode);
            for (const item of overflow) {
                enqueueDeferred(node.helixId, item.neighborId, node.offset);
            }

            for (const item of selected) {
                if (positions.has(item.neighborId)) continue;

                const predictedGlobal = normalizeAngle(item.localAngle + node.offset);
                const snappedAngle = snapToLatticeAngle(predictedGlobal);
                const step = getStep(parentCoord, snappedAngle);

                const childCoord: GridCoord = {
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

                const childOffset = computeChildOffset(
                    node.helixId,
                    item.neighborId,
                    parentCoord,
                    childCoord,
                    item.localAngle,
                    snappedAngle
                );

                globalRotationOffsets.set(item.neighborId, childOffset);
                targetQueue.push({ helixId: item.neighborId, offset: childOffset });
            }
        };

        // Phase 1: strict weighted BFS spanning tree from helix 0.
        positions.set(0, { col: 0, row: 0 });
        occupied.set('0,0', 0);
        globalRotationOffsets.set(0, 0);

        const mainQueue: QueueNode[] = [{ helixId: 0, offset: 0 }];
        let mainIdx = 0;
        while (mainIdx < mainQueue.length) {
            processNode(mainQueue[mainIdx++], mainQueue);
        }

        // Phase 2: place deferred artefacts nearest to their parent.
        const runDeferredSubBfs = options?.runDeferredSubBfs ?? false;
        const deferredSubQueue: QueueNode[] = [];

        let deferredIdx = 0;
        while (deferredIdx < deferredQueue.length) {
            const item = deferredQueue[deferredIdx++];
            if (positions.has(item.helixId)) continue;

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
            if (positions.has(helixId)) continue;
            const coord = findNearestOpen(rootCoord);
            positions.set(helixId, coord);
            occupied.set(keyOf(coord), helixId);
        }

        const result = new Map<number, [number, number]>();
        for (const [helixId, coord] of positions.entries()) {
            result.set(helixId, [coord.col, coord.row]);
        }

        return result;
    }


    // ─────────────────────────────────────────────────────────────────────
    //  Helix renumbering: MST-DFS + 2-Opt
    //
    //  Goal: choose a permutation of helix IDs so consecutive numbers
    //  share a crossover whenever possible, and otherwise are physically
    //  close on the grid.
    //
    //  Cost(a, b) = (jump ? 1000 : 0) + euclideanGridDistance(a, b)
    //
    //  Scadnano lattice distances stay well under 1000 in the units used
    //  by scadnano.honeycombToWorld / squareToWorld, so this single scalar
    //  is automatically lex-correct: jumps dominate, distance breaks ties.
    //
    //  Pipeline:
    //    1. Build adjacency (crossovers) + distance helpers from grid.
    //    2. Connected components → lower bound on jumps (just for logging).
    //    3. Prim's MST on the cost graph (same pattern as alignGridPrim,
    //       cost-minimizing).
    //    4. Pre-order DFS down the MST → starting Hamiltonian path.
    //    5. 2-Opt cleanup, two phases:
    //         (a) accept only swaps that strictly reduce jumps,
    //         (b) freeze jumps, accept swaps that reduce distance.
    //
    //  Returns order[] (newIdx → oldHelixId) and remap (oldId → newId).
    //  This function does NOT mutate the grid or helices array; the caller
    //  applies the remap wherever helix IDs live.
    // ─────────────────────────────────────────────────────────────────────

    const RENUMBER_JUMP_PENALTY = 1000;

    export interface RenumberStats {
        n: number;
        components: number;
        lowerBoundJumps: number;
        initialJumps: number;
        initialDistance: number;
        afterPhase1Jumps: number;
        finalJumps: number;
        finalDistance: number;
    }

    export interface RenumberResult {
        order: number[];                    // order[newIdx] = oldHelixId
        remap: Map<number, number>;         // oldHelixId → newHelixId
        stats: RenumberStats;
    }

    export function renumberHelices(
        grid: GridMap,
        helixPos: Map<number, [number, number]>,
        lattice: string = 'honeycomb'
    ): RenumberResult {

        // ── 1. Adjacency + cost helpers ────────────────────────────────
        const { crossovers, helixIds } = collectCrossovers(grid);
        for (const id of helixPos.keys()) helixIds.add(id);
        const nodes = Array.from(helixIds).sort((a, b) => a - b);
        const n = nodes.length;
        const idx = new Map<number, number>();
        nodes.forEach((id, i) => idx.set(id, i));

        const isSquare = (lattice ?? '').toLowerCase() === 'square';
        const world = nodes.map(id => {
            const pos = helixPos.get(id) ?? [0, 0];
            const v = isSquare
                ? scadnano.squareToWorld(pos[0], pos[1])
                : scadnano.honeycombToWorld(pos[0], pos[1]);
            return [v.x, v.y] as [number, number];
        });

        const dist = (a: number, b: number): number => {
            const [ax, ay] = world[a];
            const [bx, by] = world[b];
            return Math.hypot(ax - bx, ay - by);
        };

        // a–b adjacent iff at least one observed crossover (either walk).
        const adj: boolean[][] = Array.from({ length: n }, () => new Array(n).fill(false));
        for (const [from, inner] of crossovers.entries()) {
            const i = idx.get(from);
            if (i === undefined) continue;
            for (const [to, counts] of inner.entries()) {
                const j = idx.get(to);
                if (j === undefined) continue;
                if ((counts.sameWalk + counts.diffWalk) > 0) {
                    adj[i][j] = true;
                    adj[j][i] = true;
                }
            }
        }

        const isJump = (a: number, b: number) => !adj[a][b];
        const edgeCost = (a: number, b: number) =>
            (isJump(a, b) ? RENUMBER_JUMP_PENALTY : 0) + dist(a, b);

        // Trivial cases: 0 or 1 helices → identity remap.
        if (n <= 1) {
            const order = nodes.slice();
            const remap = new Map<number, number>();
            order.forEach((id, i) => remap.set(id, i));
            console.log(`[renumberHelices] n=${n}, nothing to renumber.`);
            return {
                order, remap,
                stats: {
                    n, components: n, lowerBoundJumps: 0,
                    initialJumps: 0, initialDistance: 0,
                    afterPhase1Jumps: 0, finalJumps: 0, finalDistance: 0
                }
            };
        }

        // ── 2. Connected components (lower bound on jumps) ──────────────
        const compId = new Array<number>(n).fill(-1);
        let nComps = 0;
        for (let s = 0; s < n; s++) {
            if (compId[s] !== -1) continue;
            const stack = [s];
            compId[s] = nComps;
            while (stack.length) {
                const u = stack.pop()!;
                for (let v = 0; v < n; v++) {
                    if (compId[v] !== -1) continue;
                    if (adj[u][v]) { compId[v] = nComps; stack.push(v); }
                }
            }
            nComps++;
        }
        const lowerBoundJumps = Math.max(0, nComps - 1);

        // ── 3. Prim's MST (cost-minimizing) ─────────────────────────────
        const inMST = new Array<boolean>(n).fill(false);
        const mstChildren: number[][] = Array.from({ length: n }, () => []);
        const bestEdge = new Array<number>(n).fill(Infinity);
        const bestParent = new Array<number>(n).fill(-1);

        inMST[0] = true;
        for (let v = 1; v < n; v++) {
            bestEdge[v] = edgeCost(0, v);
            bestParent[v] = 0;
        }
        for (let added = 1; added < n; added++) {
            let u = -1, best = Infinity;
            for (let v = 0; v < n; v++) {
                if (!inMST[v] && bestEdge[v] < best) { best = bestEdge[v]; u = v; }
            }
            if (u === -1) break;
            inMST[u] = true;
            mstChildren[bestParent[u]].push(u);
            for (let v = 0; v < n; v++) {
                if (!inMST[v]) {
                    const c = edgeCost(u, v);
                    if (c < bestEdge[v]) { bestEdge[v] = c; bestParent[v] = u; }
                }
            }
        }

        // ── 4. Pre-order DFS along the MST → starting path ──────────────
        // Visit the cheapest child first at each branch — slightly tighter
        // initial tour, gives 2-Opt less to clean up.
        const path: number[] = [];
        const dfsStack: number[] = [0];
        const visited = new Array<boolean>(n).fill(false);
        // Iterative pre-order DFS to avoid recursion on very deep trees.
        while (dfsStack.length) {
            const u = dfsStack.pop()!;
            if (visited[u]) continue;
            visited[u] = true;
            path.push(u);
            const kids = mstChildren[u].slice().sort((a, b) => edgeCost(u, b) - edgeCost(u, a));
            // Push in reverse-cost order so cheapest is popped first.
            for (const c of kids) dfsStack.push(c);
        }

        // ── Path cost helpers ──────────────────────────────────────────
        const pathStats = (p: number[]): { jumps: number; distance: number } => {
            let jumps = 0, distance = 0;
            for (let i = 0; i + 1 < p.length; i++) {
                if (isJump(p[i], p[i + 1])) jumps++;
                distance += dist(p[i], p[i + 1]);
            }
            return { jumps, distance };
        };

        const initial = pathStats(path);

        // ── 5. 2-Opt cleanup ────────────────────────────────────────────
        const reverseSegment = (p: number[], i: number, j: number) => {
            while (i < j) { const t = p[i]; p[i] = p[j]; p[j] = t; i++; j--; }
        };

        // 2-Opt reverses path[i..j]. Internal edges within that segment
        // keep the same endpoints (just reversed direction), so only the
        // two boundary edges change. Both deltas reduce to a 4-vertex check.
        const jumpDelta = (p: number[], i: number, j: number): number => {
            const a = p[i - 1], b = p[i], c = p[j], d = p[j + 1];
            const before = (isJump(a, b) ? 1 : 0) + (isJump(c, d) ? 1 : 0);
            const after  = (isJump(a, c) ? 1 : 0) + (isJump(b, d) ? 1 : 0);
            return after - before;
        };
        const distDelta = (p: number[], i: number, j: number): number => {
            const a = p[i - 1], b = p[i], c = p[j], d = p[j + 1];
            return (dist(a, c) + dist(b, d)) - (dist(a, b) + dist(c, d));
        };

        // Phase A: only accept swaps that strictly reduce jumps.
        let improved = true;
        while (improved) {
            improved = false;
            for (let i = 1; i < path.length - 1; i++) {
                for (let j = i + 1; j < path.length - 1; j++) {
                    if (jumpDelta(path, i, j) < 0) {
                        reverseSegment(path, i, j);
                        improved = true;
                    }
                }
            }
        }
        const afterPhase1 = pathStats(path);

        // Phase B: jumps frozen — only accept swaps that don't increase
        // jumps and strictly reduce distance.
        improved = true;
        while (improved) {
            improved = false;
            for (let i = 1; i < path.length - 1; i++) {
                for (let j = i + 1; j < path.length - 1; j++) {
                    if (jumpDelta(path, i, j) !== 0) continue;
                    if (distDelta(path, i, j) < -1e-9) {
                        reverseSegment(path, i, j);
                        improved = true;
                    }
                }
            }
        }
        const final = pathStats(path);

        // ── 6. Build remap and report ──────────────────────────────────
        const order = path.map(i => nodes[i]);
        const remap = new Map<number, number>();
        order.forEach((id, i) => remap.set(id, i));

        const stats: RenumberStats = {
            n,
            components: nComps,
            lowerBoundJumps,
            initialJumps: initial.jumps,
            initialDistance: initial.distance,
            afterPhase1Jumps: afterPhase1.jumps,
            finalJumps: final.jumps,
            finalDistance: final.distance,
        };

        console.log(
            `[renumberHelices] n=${n}, components=${nComps}, ` +
            `jumps lower bound=${lowerBoundJumps} | ` +
            `initial=[${initial.jumps} jumps, dist=${initial.distance.toFixed(2)}] → ` +
            `phase1=[${afterPhase1.jumps} jumps] → ` +
            `final=[${final.jumps} jumps, dist=${final.distance.toFixed(2)}]`
        );

        return { order, remap, stats };
    }

    /**
     * Apply a helix renumber remap to the in-memory layout.
     *
     * Mutates `grid` (rewrites every mark's helixId), reorders `helices`
     * so its index matches the new helix ID, and returns a fresh
     * `helixPos` map keyed by new helix IDs. Pure on the originals
     * otherwise — Nucleotide objects are not touched.
     *
     * The remap must be a bijection over [0, helices.length).
     */
    export function applyHelixRenumber(
        helices: Nucleotide[][],
        grid: GridMap,
        helixPos: Map<number, [number, number]>,
        remap: Map<number, number>
    ): {
        helices: Nucleotide[][];
        helixPos: Map<number, [number, number]>;
    } {
        const n = helices.length;
        const newHelices: Nucleotide[][] = new Array(n);
        for (let oldId = 0; oldId < n; oldId++) {
            const newId = remap.get(oldId);
            const slot = (newId !== undefined && newId >= 0 && newId < n) ? newId : oldId;
            newHelices[slot] = helices[oldId];
        }

        for (const mark of grid.values()) {
            const newId = remap.get(mark.helixId);
            if (newId !== undefined) mark.helixId = newId;
        }

        const newHelixPos = new Map<number, [number, number]>();
        for (const [oldId, pos] of helixPos.entries()) {
            const newId = remap.get(oldId);
            newHelixPos.set(newId !== undefined ? newId : oldId, pos);
        }

        return { helices: newHelices, helixPos: newHelixPos };
    }
}