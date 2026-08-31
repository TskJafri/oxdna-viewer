/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />

namespace toscad {

    // Tuples are [helixA, helixB, crossoverCount], sorted by weight descending.
    export function crossoverPairs(grid: GridMap): [number, number, number][] {
        const pairs: [number, number, number][] = [];
        for (const [a, neighbors] of getConnectionCounts(grid))
            for (const [b, w] of neighbors)
                if (a < b) pairs.push([a, b, w]);
        return pairs.sort((p, q) => q[2] - p[2] || p[0] - q[0] || p[1] - q[1]);
    }

    // Tuples are [helixA, helixB, crossoverCount, neighborhoodSum], sorted by count then sum.
    export function crossoverPairsNbhd(grid: GridMap, pairs = crossoverPairs(grid)): [number, number, number, number][] {
        const counts = getConnectionCounts(grid);
        const deg = (h: number) => [...(counts.get(h)?.values() ?? [])].reduce((s, w) => s + w, 0);
        return pairs.map(([a, b, w]) => [a, b, w, deg(a) + deg(b)] as [number, number, number, number])
                    .sort((p, q) => q[2] - p[2] || q[3] - p[3] || p[0] - q[0] || p[1] - q[1]);
    }

    const norm360 = (a: number) => ((a % 360) + 360) % 360;

    // Valid global direction labels. Array index == rotation slot.
    export function latticeDirs(lat: LatticeKind): number[] {
        return lat === 'square' ? [0, 90, 180, 270] : [0, 120, 240];
    }

    // Nearest valid direction; ties go to the smaller angle.
    export function snapDir(angle: number, lat: LatticeKind): number {
        const dirs = latticeDirs(lat);
        const dist = (x: number, y: number) => {
            const d = Math.abs(norm360(x) - norm360(y));
            return Math.min(d, 360 - d);
        };
        let best = dirs[0];
        for (const c of dirs) {
            const dc = dist(angle, c), db = dist(angle, best);
            if (dc < db || (dc === db && c < best)) best = c;
        }
        return best;
    }

    // Direction label pointing back, as seen from the cell you just stepped into.
    // Square: always +180. Honeycomb: +240 from an even cell, +120 from an odd one.
    export function oppositeDir(dir: number, parity: 0 | 1, lat: LatticeKind): number {
        const d = norm360(dir);
        if (lat === 'square') return (d + 180) % 360;
        return parity === 0 ? (d + 240) % 360 : (d + 120) % 360;
    }

    // (dCol,dRow) for a direction label, given the parity of the cell you stand on.
    // stepFor(oppositeDir(d,p), 1-p) === -stepFor(d,p) for every d,p on both lattices.
    export function stepFor(dir: number, parity: 0 | 1, lat: LatticeKind): { dCol: number; dRow: number } | null {
        const d = norm360(dir);
        if (lat === 'square') {
            if (d === 0) return { dCol: 1, dRow: 0 };
            if (d === 90) return { dCol: 0, dRow: 1 };
            if (d === 180) return { dCol: -1, dRow: 0 };
            if (d === 270) return { dCol: 0, dRow: -1 };
            return null;
        }
        if (parity === 0) {
            if (d === 0) return { dCol: 1, dRow: 0 };
            if (d === 120) return { dCol: -1, dRow: 0 };
            if (d === 240) return { dCol: 0, dRow: -1 };
            return null;
        }
        if (d === 0) return { dCol: 1, dRow: 0 };
        if (d === 120) return { dCol: 0, dRow: 1 };
        if (d === 240) return { dCol: -1, dRow: 0 };
        return null;
    }

    export interface OrientEdge {
        a: number;
        b: number;
        weight: number;
        compat: Uint8Array; // compat[slotA * K + slotB] === 1 when satisfied
    }

    // Run a small sweep to go from local angles in getAngles() to global orientations **of the helices**.
    // First use BFS to seed a spanning tree, and then use votes to repair cycles. The result is then used by Kruskal's.
    // TODO: There's a TODO in here. Go hunt for it.
    export function voteOrientations(
        networkMap: Map<number, Map<number, number>>,
        grid: GridMap,
        lattice: string = 'honeycomb',
        init?: Map<number, number>,
        maxSweeps: number = 50
    ) {
        const lat = resolveLatticeKind(lattice);
        const dirs = latticeDirs(lat);
        const K = dirs.length;
        const snap = (angle: number) => snapDir(angle, lat);
        const opposite = (dir: number, parity: 0 | 1) => oppositeDir(dir, parity, lat);

        // parity based evaluation for honeycomb (a->b is not 180º rotated from b->a; it is rotated either 120º or 240º).
        const parity = new Map<number, 0 | 1>();
        const oddEdges: Array<[number, number]> = [];
        for (const seed of [...networkMap.keys()].sort((a, b) => a - b)) {
            if (parity.has(seed)) continue;
            parity.set(seed, 0);
            const queue = [seed];
            for (let qi = 0; qi < queue.length; qi++) {
                const h = queue[qi];
                const p = parity.get(h)!;
                for (const n of networkMap.get(h)?.keys() ?? []) {
                    if (!parity.has(n)) { parity.set(n, p === 0 ? 1 : 0); queue.push(n); }
                    else if (parity.get(n) === p) oddEdges.push([Math.min(h, n), Math.max(h, n)]);
                }
            }
        }

        // build K×K compatibility table per edge
        const counts = getConnectionCounts(grid);
        const edges: OrientEdge[] = [];
        for (const [a, nbrs] of networkMap) {
            for (const [b, angleAB] of nbrs) {
                if (a >= b) continue;
                const angleBA = networkMap.get(b)?.get(a);
                if (typeof angleBA !== 'number') continue;
                const pa = parity.get(a) ?? 0;
                const compat = new Uint8Array(K * K);
                for (let ra = 0; ra < K; ra++) {
                    const need = opposite(snap(angleAB + dirs[ra]), pa);
                    for (let rb = 0; rb < K; rb++) {
                        if (snap(angleBA + dirs[rb]) === need) compat[ra * K + rb] = 1;
                    }
                }
                edges.push({ a, b, weight: counts.get(a)?.get(b) ?? 1, compat });
            }
        }

        // incident edges per helix
        const incident = new Map<number, Array<{ e: OrientEdge; isA: boolean }>>();
        const touch = (h: number, e: OrientEdge, isA: boolean) => {
            let arr = incident.get(h);
            if (!arr) { arr = []; incident.set(h, arr); }
            arr.push({ e, isA });
        };
        for (const e of edges) { touch(e.a, e, true); touch(e.b, e, false); }

        const order = [...networkMap.keys()].sort((a, b) => a - b);

        // Initialize with realistic global orientations by using BFS (as opposed to starting with 0º). 
        // TODO: Square lattice has 4 neighbors, causing more ambiguity than honeycomb.
        const slot = new Map<number, number>();
        let seededAmbiguous = 0;

        if (init) {
            for (const h of order) {
                const seed = init.get(h);
                slot.set(h, typeof seed === 'number' ? dirs.indexOf(snap(seed)) : 0);
            }
        } else {
            for (const seedHelix of order) {
                if (slot.has(seedHelix)) continue;
                slot.set(seedHelix, 0);              // component gauge: only differences matter
                const q = [seedHelix];
                for (let qi = 0; qi < q.length; qi++) {
                    const h = q[qi];
                    const rh = slot.get(h)!;
                    for (const { e, isA } of incident.get(h) ?? []) {
                        const other = isA ? e.b : e.a;
                        if (slot.has(other)) continue;
                        let pick = -1, hits = 0;
                        for (let r = 0; r < K; r++) {
                            if (isA ? e.compat[rh * K + r] : e.compat[r * K + rh]) {
                                if (pick < 0) pick = r;
                                hits++;
                            }
                        }
                        if (pick < 0) continue;      // unsatisfiable at rh; let another edge seed it
                        if (hits > 1) seededAmbiguous++;
                        slot.set(other, pick);
                        q.push(other);
                    }
                }
            }
            for (const h of order) if (!slot.has(h)) slot.set(h, 0);   // isolated helices
        }
        if (seededAmbiguous > 0) {
            console.log(`[voteOrientations] ${seededAmbiguous} seed step(s) had >1 compatible slot; took the lowest.`);
        }

        // iterative majority vote
        let sweeps = 0;
        let converged = false;

        while (sweeps < maxSweeps) {
            sweeps++;
            let changed = false;
            for (const h of order) {
                const inc = incident.get(h);
                if (!inc || !inc.length) continue;
                const current = slot.get(h) ?? 0;
                let bestSlot = current, bestScore = -1;
                for (let r = 0; r < K; r++) {
                    let score = 0;
                    for (const { e, isA } of inc) {
                        const idx = isA ? r * K + (slot.get(e.b) ?? 0) : (slot.get(e.a) ?? 0) * K + r;
                        if (e.compat[idx]) score += e.weight;
                    }
                    if (score > bestScore || (score === bestScore && r === current)) { bestScore = score; bestSlot = r; }
                }
                if (bestSlot !== current) { slot.set(h, bestSlot); changed = true; }
            }
            if (!changed) { converged = true; break; }
        }

        // ── score & package results ──
        let satisfiedWeight = 0;
        const violated: Array<[number, number]> = [];
        for (const e of edges) {
            if (e.compat[(slot.get(e.a) ?? 0) * K + (slot.get(e.b) ?? 0)]) satisfiedWeight += e.weight;
            else violated.push([e.a, e.b]);
        }

        const orientation = new Map<number, number>();
        for (const [h, r] of slot) orientation.set(h, dirs[r]);

        return {
            slot, orientation, violated, satisfiedWeight,
            totalWeight: edges.reduce((s, e) => s + e.weight, 0),
            sweeps, converged, parity, oddEdges, edges
        };
    }

    export function kruskals(
        networkMap: Map<number, Map<number, number>>,
        grid: GridMap,
        lattice: string = 'honeycomb',
        vote = voteOrientations(networkMap, grid, lattice),
        pinnedHelices: number[] = []
    ) {
        const lat = resolveLatticeKind(lattice);
        const dirs = latticeDirs(lat);
        const K = dirs.length;
        const { slot, parity, edges } = vote;

        const satisfied = (e: OrientEdge) =>
            e.compat[(slot.get(e.a) ?? 0) * K + (slot.get(e.b) ?? 0)] === 1;

        // Step from a to b using a's global slot and a's 2-coloring parity.
        const stepAB = (a: number, b: number): { dCol: number; dRow: number } | null => {
            const raw = networkMap.get(a)?.get(b);
            if (typeof raw !== 'number') return null;
            return stepFor(snapDir(raw + dirs[slot.get(a) ?? 0], lat), parity.get(a) ?? 0, lat);
        };

        // Union-find carrying explicit member lists, because merging has to translate every member of one side.
        const root = new Map<number, number>();
        const members = new Map<number, number[]>();
        const pos = new Map<number, [number, number]>();

        for (const h of [...networkMap.keys()].sort((x, y) => x - y)) {
            root.set(h, h);
            members.set(h, [h]);
            // Seed so (col+row)&1 equals the helix's 2-coloring parity. 
            // Every lattice step flips both, so the invariant then holds across the whole component.
            pos.set(h, (parity.get(h) ?? 0) === 0 ? [0, 0] : [1, 0]);
        }
        const find = (h: number) => root.get(h)!;

        const treeEdges: Array<[number, number]> = [];
        const cycleEdges: Array<[number, number]> = [];
        const usedViolatedEdges: Array<[number, number]> = [];
        const parityConflicts: Array<[number, number]> = [];

        // merging 2 groups (or clumps)... 
        const tryMerge = (e: OrientEdge) => {
            const { a, b } = e;
            const ra = find(a), rb = find(b);
            if (ra === rb) { cycleEdges.push([a, b]); return; }

            const step = stepAB(a, b);
            if (!step) return;

            const pa = pos.get(a)!, pb = pos.get(b)!;
            // Where b must sit relative to a, then the shift that moves b's frame there.
            const T: [number, number] = [pa[0] + step.dCol - pb[0], pa[1] + step.dRow - pb[1]];

            // A parity-flipping shift would invalidate every honeycomb step inside the moved component.
            // With a valid 2-coloring T is always even, so this firing means the colouring is broken and the edge must not define geometry.
            if (lat === 'honeycomb' && (((T[0] + T[1]) & 1) !== 0)) {
                parityConflicts.push([a, b]);
                return;
            }

            const listA = members.get(ra)!, listB = members.get(rb)!;
            const moveB = listB.length <= listA.length;
            const moveList = moveB ? listB : listA;
            const keepRoot = moveB ? ra : rb;
            const dropRoot = moveB ? rb : ra;
            // Moving A instead of B needs -T: pos_a - T = pos_b - step, so pos_a + step = pos_b.
            const d: [number, number] = moveB ? T : [-T[0], -T[1]];

            for (const h of moveList) {
                const p = pos.get(h)!;
                pos.set(h, [p[0] + d[0], p[1] + d[1]]);
                root.set(h, keepRoot);
            }
            members.get(keepRoot)!.push(...moveList);
            members.delete(dropRoot);

            treeEdges.push([a, b]);
            if (!satisfied(e)) usedViolatedEdges.push([a, b]);
        };

        // order: satisfied > heaviest > helixID (lowest ID wins) > violated
        const ordered = [...edges].sort((p, q) => {
            const sp = satisfied(p) ? 1 : 0, sq = satisfied(q) ? 1 : 0;
            return sq - sp || q.weight - p.weight || p.a - q.a || p.b - q.b;
        });

        // Pinned helices jump the queue: every edge touching one is processed first, so pinning
        // h3 pins h3-h2, h3-h5, h3-h8, etc. Pinned helices need not be connected to each other.
        const pinSet = new Set(pinnedHelices);
        const isPinned = (e: OrientEdge) => pinSet.has(e.a) || pinSet.has(e.b);

        // Actual merging...
        for (const e of ordered) if (isPinned(e)) tryMerge(e);
        for (const e of ordered) if (!isPinned(e)) tryMerge(e);

        // Spread components side by side. Overlaps WITHIN a component are kept but separate components all sit near their own origin and would pile up.
        // The shift is nudged to stay parity-even so honeycomb steps remain valid.
        const components = [...members.values()]
            .map(list => [...list].sort((x, y) => x - y))
            .sort((x, y) => y.length - x.length || x[0] - y[0]);

        let cursor = 0;
        for (const list of components) {
            let minCol = Infinity, maxCol = -Infinity, minRow = Infinity;
            for (const h of list) {
                const p = pos.get(h)!;
                minCol = Math.min(minCol, p[0]);
                maxCol = Math.max(maxCol, p[0]);
                minRow = Math.min(minRow, p[1]);
            }
            let dCol = cursor - minCol;
            const dRow = -minRow;
            if (((dCol + dRow) & 1) !== 0) dCol += 1;
            for (const h of list) {
                const p = pos.get(h)!;
                pos.set(h, [p[0] + dCol, p[1] + dRow]);
            }
            cursor = maxCol + dCol + 3;
        }

        // Overlaps: reported, not fixed. Note grid-editor.html keys nodes by "col,row", so
        // it will render only ONE helix per cell -- this list is the only place they show up.
        const byCell = new Map<string, number[]>();
        for (const [h, p] of pos) {
            const k = `${p[0]},${p[1]}`;
            if (!byCell.has(k)) byCell.set(k, []);
            byCell.get(k)!.push(h);
        }
        const overlaps = [...byCell.entries()]
            .filter(([, v]) => v.length > 1)
            .map(([k, v]) => {
                const [c, r] = k.split(',').map(Number);
                return { cell: [c, r] as [number, number], helices: v.sort((x, y) => x - y) };
            });

        // Orientation agreement does NOT imply position closure. This just checks (doesn't fix).
        const translationViolations: Array<[number, number]> = [];
        for (const e of edges) {
            const step = stepAB(e.a, e.b);
            const pa = pos.get(e.a), pb = pos.get(e.b);
            if (!step || !pa || !pb) continue;
            if (pa[0] + step.dCol !== pb[0] || pa[1] + step.dRow !== pb[1]) {
                translationViolations.push([e.a, e.b]);
            }
        }

        console.log(
            `[kruskals] ${pos.size} helices, ${components.length} component(s) | ` +
            `tree=${treeEdges.length} cycle=${cycleEdges.length} | ` +
            `usedViolated=${usedViolatedEdges.length} parityConflicts=${parityConflicts.length} | ` +
            `overlaps=${overlaps.length} translationViolations=${translationViolations.length}`
        );

        return {
            positions: pos,
            components,
            treeEdges,
            cycleEdges,
            usedViolatedEdges,
            parityConflicts,
            overlaps,
            translationViolations,
            slot,
            orientation: vote.orientation
        };
    }

    export function positionsToJSON(positions: Map<number, [number, number]>): string {
        return JSON.stringify([...positions.entries()].sort((a, b) => a[0] - b[0]));
    }
}
