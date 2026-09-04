"use strict";
/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
var toscad;
(function (toscad) {
    // Tuples are [helixA, helixB, crossoverCount], sorted by weight descending.
    function crossoverPairs(grid) {
        const pairs = [];
        for (const [a, neighbors] of toscad.getConnectionCounts(grid))
            for (const [b, w] of neighbors)
                if (a < b)
                    pairs.push([a, b, w]);
        return pairs.sort((p, q) => q[2] - p[2] || p[0] - q[0] || p[1] - q[1]);
    }
    toscad.crossoverPairs = crossoverPairs;
    // Tuples are [helixA, helixB, crossoverCount, neighborhoodSum], sorted by count then sum.
    function crossoverPairsNbhd(grid, pairs = crossoverPairs(grid)) {
        const counts = toscad.getConnectionCounts(grid);
        const deg = (h) => [...(counts.get(h)?.values() ?? [])].reduce((s, w) => s + w, 0);
        return pairs.map(([a, b, w]) => [a, b, w, deg(a) + deg(b)])
            .sort((p, q) => q[2] - p[2] || q[3] - p[3] || p[0] - q[0] || p[1] - q[1]);
    }
    toscad.crossoverPairsNbhd = crossoverPairsNbhd;
    const norm360 = (a) => ((a % 360) + 360) % 360;
    // Valid global direction labels. Array index == rotation slot.
    function latticeDirs(lat) {
        return lat === 'square' ? [0, 90, 180, 270] : [0, 120, 240];
    }
    toscad.latticeDirs = latticeDirs;
    // Nearest valid direction; ties go to the smaller angle.
    function snapDir(angle, lat) {
        const dirs = latticeDirs(lat);
        const dist = (x, y) => {
            const d = Math.abs(norm360(x) - norm360(y));
            return Math.min(d, 360 - d);
        };
        let best = dirs[0];
        for (const c of dirs) {
            const dc = dist(angle, c), db = dist(angle, best);
            if (dc < db || (dc === db && c < best))
                best = c;
        }
        return best;
    }
    toscad.snapDir = snapDir;
    // Direction label pointing back, as seen from the cell you just stepped into.
    // Square: always +180. Honeycomb: +240 from an even cell, +120 from an odd one.
    function oppositeDir(dir, parity, lat) {
        const d = norm360(dir);
        if (lat === 'square')
            return (d + 180) % 360;
        return parity === 0 ? (d + 240) % 360 : (d + 120) % 360;
    }
    toscad.oppositeDir = oppositeDir;
    // (dCol,dRow) for a direction label, given the parity of the cell you stand on.
    // stepFor(oppositeDir(d,p), 1-p) === -stepFor(d,p) for every d,p on both lattices.
    function stepFor(dir, parity, lat) {
        const d = norm360(dir);
        if (lat === 'square') {
            if (d === 0)
                return { dCol: 1, dRow: 0 };
            if (d === 90)
                return { dCol: 0, dRow: 1 };
            if (d === 180)
                return { dCol: -1, dRow: 0 };
            if (d === 270)
                return { dCol: 0, dRow: -1 };
            return null;
        }
        if (parity === 0) {
            if (d === 0)
                return { dCol: 1, dRow: 0 };
            if (d === 120)
                return { dCol: -1, dRow: 0 };
            if (d === 240)
                return { dCol: 0, dRow: -1 };
            return null;
        }
        if (d === 0)
            return { dCol: 1, dRow: 0 };
        if (d === 120)
            return { dCol: 0, dRow: 1 };
        if (d === 240)
            return { dCol: -1, dRow: 0 };
        return null;
    }
    toscad.stepFor = stepFor;
    // Run a small sweep to go from local angles in getAngles() to global orientations **of the helices**.
    // First use BFS to seed a spanning tree, and then use votes to repair cycles. The result is then used by Kruskal's.
    // TODO: There's a TODO in here. Go hunt for it.
    function voteOrientations(networkMap, grid, lattice = 'honeycomb', init, maxSweeps = 50, paritySeed) {
        const lat = toscad.resolveLatticeKind(lattice);
        const dirs = latticeDirs(lat);
        const K = dirs.length;
        const snap = (angle) => snapDir(angle, lat);
        const opposite = (dir, parity) => oppositeDir(dir, parity, lat);
        // parity based evaluation for honeycomb (a->b is not 180º rotated from b->a; it is rotated either 120º or 240º).
        //
        // The 2-coloring is only determined up to a GLOBAL FLIP per component, so seeding every
        // component at 0 picks the absolute parity arbitrarily -- a structure that should come out
        // odd-even-odd can just as easily come out even-odd-even. When the caller has pinned helices
        // at known cells, `paritySeed` fixes that choice, and it has to happen HERE: the compat
        // tables below call opposite(dir, pa), which differs by parity, so flipping labels after the
        // vote would invalidate the entire solve.
        //
        // Seeded helices are visited first so their component is colored outward from them.
        const parity = new Map();
        const oddEdges = [];
        const seedOrder = [
            ...[...(paritySeed?.keys() ?? [])].filter(h => networkMap.has(h)).sort((a, b) => a - b),
            ...[...networkMap.keys()].sort((a, b) => a - b)
        ];
        for (const seed of seedOrder) {
            if (parity.has(seed))
                continue;
            parity.set(seed, paritySeed?.get(seed) ?? 0);
            const queue = [seed];
            for (let qi = 0; qi < queue.length; qi++) {
                const h = queue[qi];
                const p = parity.get(h);
                for (const n of networkMap.get(h)?.keys() ?? []) {
                    if (!parity.has(n)) {
                        parity.set(n, p === 0 ? 1 : 0);
                        queue.push(n);
                    }
                    else if (parity.get(n) === p)
                        oddEdges.push([Math.min(h, n), Math.max(h, n)]);
                }
            }
        }
        // Two seeded helices in the SAME component can demand incompatible parities, since their
        // relative parity is already forced by the path length between them. The first one wins
        // (it seeded the colouring); report the rest rather than silently ignoring them.
        const paritySeedConflicts = [];
        if (paritySeed) {
            for (const [h, want] of paritySeed) {
                const got = parity.get(h);
                if (got !== undefined && got !== want) {
                    paritySeedConflicts.push({ helix: h, requested: want, actual: got });
                }
            }
            if (paritySeedConflicts.length) {
                console.warn(`[voteOrientations] ${paritySeedConflicts.length} parity seed(s) unsatisfiable ` +
                    `(relative parity is fixed by path length):`, paritySeedConflicts.map(c => `h${c.helix} wanted ${c.requested}, got ${c.actual}`).join('; '));
            }
        }
        // build K×K compatibility table per edge
        const counts = toscad.getConnectionCounts(grid);
        const edges = [];
        for (const [a, nbrs] of networkMap) {
            for (const [b, angleAB] of nbrs) {
                if (a >= b)
                    continue;
                const angleBA = networkMap.get(b)?.get(a);
                if (typeof angleBA !== 'number')
                    continue;
                const pa = parity.get(a) ?? 0;
                const compat = new Uint8Array(K * K);
                for (let ra = 0; ra < K; ra++) {
                    const need = opposite(snap(angleAB + dirs[ra]), pa);
                    for (let rb = 0; rb < K; rb++) {
                        if (snap(angleBA + dirs[rb]) === need)
                            compat[ra * K + rb] = 1;
                    }
                }
                edges.push({ a, b, weight: counts.get(a)?.get(b) ?? 1, compat });
            }
        }
        // incident edges per helix
        const incident = new Map();
        const touch = (h, e, isA) => {
            let arr = incident.get(h);
            if (!arr) {
                arr = [];
                incident.set(h, arr);
            }
            arr.push({ e, isA });
        };
        for (const e of edges) {
            touch(e.a, e, true);
            touch(e.b, e, false);
        }
        const order = [...networkMap.keys()].sort((a, b) => a - b);
        // Initialize with realistic global orientations by using BFS (as opposed to starting with 0º). 
        // TODO: Square lattice has 4 neighbors, causing more ambiguity than honeycomb.
        const slot = new Map();
        let seededAmbiguous = 0;
        if (init) {
            for (const h of order) {
                const seed = init.get(h);
                slot.set(h, typeof seed === 'number' ? dirs.indexOf(snap(seed)) : 0);
            }
        }
        else {
            for (const seedHelix of order) {
                if (slot.has(seedHelix))
                    continue;
                slot.set(seedHelix, 0); // component gauge: only differences matter
                const q = [seedHelix];
                for (let qi = 0; qi < q.length; qi++) {
                    const h = q[qi];
                    const rh = slot.get(h);
                    for (const { e, isA } of incident.get(h) ?? []) {
                        const other = isA ? e.b : e.a;
                        if (slot.has(other))
                            continue;
                        let pick = -1, hits = 0;
                        for (let r = 0; r < K; r++) {
                            if (isA ? e.compat[rh * K + r] : e.compat[r * K + rh]) {
                                if (pick < 0)
                                    pick = r;
                                hits++;
                            }
                        }
                        if (pick < 0)
                            continue; // unsatisfiable at rh; let another edge seed it
                        if (hits > 1)
                            seededAmbiguous++;
                        slot.set(other, pick);
                        q.push(other);
                    }
                }
            }
            for (const h of order)
                if (!slot.has(h))
                    slot.set(h, 0); // isolated helices
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
                if (!inc || !inc.length)
                    continue;
                const current = slot.get(h) ?? 0;
                let bestSlot = current, bestScore = -1;
                for (let r = 0; r < K; r++) {
                    let score = 0;
                    for (const { e, isA } of inc) {
                        const idx = isA ? r * K + (slot.get(e.b) ?? 0) : (slot.get(e.a) ?? 0) * K + r;
                        if (e.compat[idx])
                            score += e.weight;
                    }
                    if (score > bestScore || (score === bestScore && r === current)) {
                        bestScore = score;
                        bestSlot = r;
                    }
                }
                if (bestSlot !== current) {
                    slot.set(h, bestSlot);
                    changed = true;
                }
            }
            if (!changed) {
                converged = true;
                break;
            }
        }
        // ── score & package results ──
        let satisfiedWeight = 0;
        const violated = [];
        for (const e of edges) {
            if (e.compat[(slot.get(e.a) ?? 0) * K + (slot.get(e.b) ?? 0)])
                satisfiedWeight += e.weight;
            else
                violated.push([e.a, e.b]);
        }
        const orientation = new Map();
        for (const [h, r] of slot)
            orientation.set(h, dirs[r]);
        return {
            slot, orientation, violated, satisfiedWeight,
            totalWeight: edges.reduce((s, e) => s + e.weight, 0),
            sweeps, converged, parity, oddEdges, edges, paritySeedConflicts
        };
    }
    toscad.voteOrientations = voteOrientations;
    function kruskals(networkMap, grid, lattice = 'honeycomb', vote = voteOrientations(networkMap, grid, lattice), pins = []) {
        const lat = toscad.resolveLatticeKind(lattice);
        const dirs = latticeDirs(lat);
        const K = dirs.length;
        const { slot, parity, edges } = vote;
        const satisfied = (e) => e.compat[(slot.get(e.a) ?? 0) * K + (slot.get(e.b) ?? 0)] === 1;
        // Step from a to b using a's global slot and a's 2-coloring parity.
        const stepAB = (a, b) => {
            const raw = networkMap.get(a)?.get(b);
            if (typeof raw !== 'number')
                return null;
            return stepFor(snapDir(raw + dirs[slot.get(a) ?? 0], lat), parity.get(a) ?? 0, lat);
        };
        // Union-find carrying explicit member lists, because merging has to translate every member of one side.
        const root = new Map();
        const members = new Map();
        const pos = new Map();
        for (const h of [...networkMap.keys()].sort((x, y) => x - y)) {
            root.set(h, h);
            members.set(h, [h]);
            // Seed so (col+row)&1 equals the helix's 2-coloring parity. 
            // Every lattice step flips both, so the invariant then holds across the whole component.
            pos.set(h, (parity.get(h) ?? 0) === 0 ? [0, 0] : [1, 0]);
        }
        const find = (h) => root.get(h);
        const treeEdges = [];
        const cycleEdges = [];
        const usedViolatedEdges = [];
        const parityConflicts = [];
        // merging 2 groups (or clumps)... 
        const tryMerge = (e, forceStep, forcePin = false) => {
            const { a, b } = e;
            const ra = find(a), rb = find(b);
            if (ra === rb) {
                cycleEdges.push([a, b]);
                return;
            }
            // forceStep is how a pin injects a user-supplied offset instead of an angle-derived one.
            const step = forceStep ?? stepAB(a, b);
            if (!step)
                return;
            const pa = pos.get(a), pb = pos.get(b);
            // Where b must sit relative to a, then the shift that moves b's frame there.
            const T = [pa[0] + step.dCol - pb[0], pa[1] + step.dRow - pb[1]];
            // A parity-flipping shift moves the whole component onto the opposite sublattice, so its
            // helices' actual (col+row)&1 stops matching their 2-coloring labels and every honeycomb
            // step derived from them afterwards uses the wrong table. Recorded, NOT blocked: a pin the
            // user got "wrong" is allowed to be wrong, and blocking here would silently fragment the
            // layout instead. Expect these pairs to render non-adjacent.
            if (lat === 'honeycomb' && (((T[0] + T[1]) & 1) !== 0)) {
                parityConflicts.push([a, b]);
                // A pin overrides this and is applied verbatim -- pins are sacred. For an ordinary
                // angle-derived edge we still bail: a parity-violating T there means the 2-coloring
                // is broken (see vote.oddEdges), and corrupting a whole component's parity alignment
                // is worse than leaving the two sides unmerged.
                if (!forcePin)
                    return;
            }
            const listA = members.get(ra), listB = members.get(rb);
            const moveB = listB.length <= listA.length;
            const moveList = moveB ? listB : listA;
            const keepRoot = moveB ? ra : rb;
            const dropRoot = moveB ? rb : ra;
            // Moving A instead of B needs -T: pos_a - T = pos_b - step, so pos_a + step = pos_b.
            const d = moveB ? T : [-T[0], -T[1]];
            for (const h of moveList) {
                const p = pos.get(h);
                pos.set(h, [p[0] + d[0], p[1] + d[1]]);
                root.set(h, keepRoot);
            }
            members.get(keepRoot).push(...moveList);
            members.delete(dropRoot);
            treeEdges.push([a, b]);
            if (!forceStep && !satisfied(e))
                usedViolatedEdges.push([a, b]);
            return true;
        };
        // order: satisfied > heaviest > helixID (lowest ID wins) > violated
        const ordered = [...edges].sort((p, q) => {
            const sp = satisfied(p) ? 1 : 0, sq = satisfied(q) ? 1 : 0;
            return sq - sp || q.weight - p.weight || p.a - q.a || p.b - q.b;
        });
        // Checks which "graph" the helix belongs to. Used for disconnected components, so parity is kept valid.
        const graphComponentOf = new Map();
        {
            let cid = 0;
            for (const seed of [...networkMap.keys()].sort((x, y) => x - y)) {
                if (graphComponentOf.has(seed))
                    continue;
                graphComponentOf.set(seed, cid);
                const q = [seed];
                for (let qi = 0; qi < q.length; qi++) {
                    for (const n of networkMap.get(q[qi])?.keys() ?? []) {
                        if (!graphComponentOf.has(n)) {
                            graphComponentOf.set(n, cid);
                            q.push(n);
                        }
                    }
                }
                cid++;
            }
        }
        const pinFailures = [];
        const pinAdjustments = [];
        for (const pin of pins) {
            const { a, b, dCol, dRow } = pin;
            if (!pos.has(a) || !pos.has(b)) {
                pinFailures.push({ pin, reason: 'helix not present in networkMap' });
                continue;
            }
            const useCol = dCol;
            if (lat === 'honeycomb' && graphComponentOf.get(a) === graphComponentOf.get(b)) {
                const wantOdd = ((dCol + dRow) & 1) !== 0;
                const parityDiffers = (parity.get(a) ?? 0) !== (parity.get(b) ?? 0);
                if (wantOdd !== parityDiffers) {
                    pinAdjustments.push({
                        pin,
                        applied: [dCol, dRow],
                        reason: `(${dCol},${dRow}) has ${wantOdd ? 'odd' : 'even'} parity but h${a}/h${b} ` +
                            `are an ${parityDiffers ? 'odd' : 'even'} number of lattice steps apart. ` +
                            `Applied as requested; h${b} will sit on the opposite sublattice.`
                    });
                }
            }
            const ra = find(a), rb = find(b);
            if (ra === rb) {
                // Already related, either by an earlier pin or a chain of them. Verify rather than move.
                const pa = pos.get(a), pb = pos.get(b);
                if (pb[0] - pa[0] !== useCol || pb[1] - pa[1] !== dRow) {
                    pinFailures.push({
                        pin,
                        reason: `conflicts with an earlier pin: already at ` +
                            `(${pb[0] - pa[0]},${pb[1] - pa[1]}), requested (${useCol},${dRow})`
                    });
                }
                continue;
            }
            const fake = { a, b, weight: Infinity, compat: new Uint8Array(K * K) };
            // forcePin bypasses tryMerge's parity guard so the offset is applied verbatim. fake because it creates a fake edge between pins.
            if (!tryMerge(fake, { dCol: useCol, dRow }, true)) {
                pinFailures.push({ pin, reason: 'merge rejected' });
            }
        }
        // Helices touching a pin get their remaining edges processed first.
        const pinSet = new Set();
        for (const p of pins) {
            pinSet.add(p.a);
            pinSet.add(p.b);
        }
        const touchesPin = (e) => pinSet.has(e.a) || pinSet.has(e.b);
        // Parity check. If the user fixes to the wrong parity (wrong as per voteOrientations), then the resulting placement is "impossible" but the algorithm will still try its best.
        const parityBroken = new Set();
        if (lat === 'honeycomb') {
            for (const [h, p] of pos) {
                if (((p[0] + p[1]) & 1) !== (parity.get(h) ?? 0))
                    parityBroken.add(h);
            }
        }
        const parityBrokenEdges = [];
        const usable = (e) => {
            if (!parityBroken.has(e.a) && !parityBroken.has(e.b))
                return true;
            parityBrokenEdges.push([e.a, e.b]);
            return false;
        };
        // Actual merging...
        for (const e of ordered)
            if (touchesPin(e) && usable(e))
                tryMerge(e);
        for (const e of ordered)
            if (!touchesPin(e) && usable(e))
                tryMerge(e);
        const packedHelices = [];
        for (let progress = true; progress;) {
            progress = false;
            for (const h of [...networkMap.keys()].sort((x, y) => x - y)) {
                const rh = find(h);
                if ((members.get(rh)?.length ?? 0) > 1)
                    continue; // already in a component
                // Pick the neighbour sitting in the largest placed component.
                let anchor = -1, anchorSize = 0;
                for (const n of networkMap.get(h)?.keys() ?? []) {
                    const size = members.get(find(n))?.length ?? 0;
                    if (find(n) !== rh && size > anchorSize) {
                        anchor = n;
                        anchorSize = size;
                    }
                }
                if (anchor < 0 || anchorSize < 2)
                    continue;
                const rootA = find(anchor);
                const occupied = new Set();
                for (const m of members.get(rootA)) {
                    const p = pos.get(m);
                    occupied.add(`${p[0]},${p[1]}`);
                }
                // Ring search, rejecting cells on the wrong parity by marking them occupied and retrying.
                const want = parity.get(h) ?? 0;
                let cell = toscad.findNearestOpenPos(pos.get(anchor), occupied);
                for (let tries = 0; lat === 'honeycomb' && tries < 64; tries++) {
                    if (((cell[0] + cell[1]) & 1) === want)
                        break;
                    occupied.add(`${cell[0]},${cell[1]}`);
                    cell = toscad.findNearestOpenPos(pos.get(anchor), occupied);
                }
                pos.set(h, cell);
                root.set(h, rootA);
                members.get(rootA).push(h);
                members.delete(rh);
                packedHelices.push({ helix: h, near: anchor, at: cell });
                progress = true;
            }
        }
        // Edges that disagree with the pins
        const pinOverruledEdges = [];
        for (const [a, b] of cycleEdges) {
            const step = stepAB(a, b);
            const pa = pos.get(a), pb = pos.get(b);
            if (!step || !pa || !pb)
                continue;
            if ((pa[0] + step.dCol !== pb[0] || pa[1] + step.dRow !== pb[1]) &&
                (pinSet.has(a) || pinSet.has(b))) {
                pinOverruledEdges.push([a, b]);
            }
        }
        if (pinAdjustments.length) {
            console.warn(`[kruskals] ${pinAdjustments.length} pin(s) applied at a geometrically impossible offset:`);
            for (const adj of pinAdjustments) {
                console.warn(`  h${adj.pin.a} -> h${adj.pin.b}: ${adj.reason}`);
            }
        }
        if (pinFailures.length) {
            console.warn(`[kruskals] ${pinFailures.length} pin(s) could not be applied:`);
            for (const f of pinFailures) {
                console.warn(`  h${f.pin.a} -> h${f.pin.b} (${f.pin.dCol},${f.pin.dRow}): ${f.reason}`);
            }
        }
        if (parityBroken.size) {
            console.warn(`[kruskals] ${parityBroken.size} helix/helices sit on a parity-flipped cell ` +
                `[${[...parityBroken].join(',')}]; ${parityBrokenEdges.length} incident edge(s) excluded ` +
                `from the tree, ${packedHelices.length} helix/helices packed by proximity instead.`);
            for (const p of packedHelices) {
                console.warn(`  h${p.helix} packed at [${p.at}] near h${p.near}`);
            }
        }
        // Spread components side by side. Overlaps WITHIN a component are kept but separate components all sit near their own origin and would pile up.
        // The shift is nudged to stay parity-even so honeycomb steps remain valid.
        const components = [...members.values()]
            .map(list => [...list].sort((x, y) => x - y))
            .sort((x, y) => y.length - x.length || x[0] - y[0]);
        let cursor = 0;
        for (const list of components) {
            let minCol = Infinity, maxCol = -Infinity, minRow = Infinity;
            for (const h of list) {
                const p = pos.get(h);
                minCol = Math.min(minCol, p[0]);
                maxCol = Math.max(maxCol, p[0]);
                minRow = Math.min(minRow, p[1]);
            }
            let dCol = cursor - minCol;
            const dRow = -minRow;
            if (((dCol + dRow) & 1) !== 0)
                dCol += 1;
            for (const h of list) {
                const p = pos.get(h);
                pos.set(h, [p[0] + dCol, p[1] + dRow]);
            }
            cursor = maxCol + dCol + 3;
        }
        // Overlaps: reported, not fixed. Note grid-editor.html keys nodes by "col,row", so
        // it will render only ONE helix per cell -- this list is the only place they show up.
        const byCell = new Map();
        for (const [h, p] of pos) {
            const k = `${p[0]},${p[1]}`;
            if (!byCell.has(k))
                byCell.set(k, []);
            byCell.get(k).push(h);
        }
        const overlaps = [...byCell.entries()]
            .filter(([, v]) => v.length > 1)
            .map(([k, v]) => {
            const [c, r] = k.split(',').map(Number);
            return { cell: [c, r], helices: v.sort((x, y) => x - y) };
        });
        // Orientation agreement does NOT imply position closure. This just checks (doesn't fix).
        const translationViolations = [];
        for (const e of edges) {
            const step = stepAB(e.a, e.b);
            const pa = pos.get(e.a), pb = pos.get(e.b);
            if (!step || !pa || !pb)
                continue;
            if (pa[0] + step.dCol !== pb[0] || pa[1] + step.dRow !== pb[1]) {
                translationViolations.push([e.a, e.b]);
            }
        }
        console.log(`[kruskals] ${pos.size} helices, ${components.length} component(s) | ` +
            `tree=${treeEdges.length} cycle=${cycleEdges.length} | ` +
            `usedViolated=${usedViolatedEdges.length} parityConflicts=${parityConflicts.length} | ` +
            `overlaps=${overlaps.length} translationViolations=${translationViolations.length}` +
            (pins.length ? ` | pins=${pins.length} failed=${pinFailures.length} overruledEdges=${pinOverruledEdges.length}` : ''));
        return {
            positions: pos,
            components,
            treeEdges,
            cycleEdges,
            usedViolatedEdges,
            parityConflicts,
            overlaps,
            translationViolations,
            pinFailures,
            pinAdjustments,
            parityBroken: [...parityBroken],
            parityBrokenEdges,
            packedHelices,
            pinOverruledEdges,
            slot,
            orientation: vote.orientation
        };
    }
    toscad.kruskals = kruskals;
    // Helper function to merge 2 helices.
    function kmMergeHelixInto(grid, helices, keepHelix, mergedHelix) {
        const keepNts = helices[keepHelix] ?? [];
        const mergedNts = helices[mergedHelix] ?? [];
        const keepNtIds = keepNts.map(n => n.id);
        const mergedNtIds = mergedNts.map(n => n.id);
        for (const [, mark] of grid.entries()) {
            if (mark.helixId === mergedHelix)
                mark.helixId = keepHelix;
            else if (mark.helixId > mergedHelix)
                mark.helixId = mark.helixId - 1;
        }
        helices[keepHelix] = keepNts.concat(mergedNts);
        helices.splice(mergedHelix, 1);
        return { keepNtIds, mergedNtIds };
    }
    toscad.kmMergeHelixInto = kmMergeHelixInto;
    // Cache the output from kmHelixEnds.
    const kmEndsCache = new WeakMap();
    function kmHelixEnds(hx) {
        if (!Array.isArray(hx) || hx.length < 2)
            return null;
        if (kmEndsCache.has(hx))
            return kmEndsCache.get(hx);
        const ep = toscad.helixEndpoints(hx);
        if (!ep || !ep.end1 || !ep.end2 || ep.end1.id === ep.end2.id) {
            kmEndsCache.set(hx, null);
            return null;
        }
        const ids = new Set();
        for (const nt of hx)
            ids.add(nt.id);
        const backboneSite = (nt) => nt.getInstanceParameter3('bbOffsets');
        const axisPoint = (nt) => {
            const pair = nt.pair;
            if (pair instanceof Nucleotide && ids.has(pair.id)) {
                return backboneSite(nt).add(backboneSite(pair)).multiplyScalar(0.5);
            }
            return backboneSite(nt);
        };
        const out = { p0: axisPoint(ep.end1), p1: axisPoint(ep.end2) };
        kmEndsCache.set(hx, out);
        return out;
    }
    toscad.kmHelixEnds = kmHelixEnds;
    // Approximate axis of a helix. A fitted axis rather than an end-to-end chord, see helix.fitPlane().
    function kmHelixAxis(grid, helices, hid) {
        void grid; // ends are topological now; grid offsets are not consulted
        const ends = kmHelixEnds(helices[hid] ?? []);
        if (!ends)
            return null;
        const dir = ends.p1.clone().sub(ends.p0);
        const len = dir.length();
        if (!isFinite(len) || len === 0)
            return null;
        return dir.divideScalar(len);
    }
    toscad.kmHelixAxis = kmHelixAxis;
    // Checks disjointness.
    function kmDisjoint(grid, helices, a, b) {
        const signed = (hid) => {
            const out = new Set();
            for (const nt of helices[hid] ?? []) {
                const mark = grid.get(nt.id);
                if (!mark || mark.helixId !== hid)
                    continue;
                out.add(2 * mark.offset + (mark.direction === 'backward' ? 1 : 0));
            }
            return out;
        };
        const s1 = signed(a), s2 = signed(b);
        if (s1.size === 0 || s2.size === 0)
            return false;
        const [small, large] = s1.size <= s2.size ? [s1, s2] : [s2, s1];
        for (const o of small)
            if (large.has(o))
                return false;
        return true;
    }
    toscad.kmDisjoint = kmDisjoint;
    // Returns fraction of each helix's axis covered by the other's shadow. 
    function kmAxisShadowOverlap(helices, a, b, threshold = 0.3) {
        const endsA = kmHelixEnds(helices[a] ?? []);
        const endsB = kmHelixEnds(helices[b] ?? []);
        if (!endsA || !endsB)
            return { overlap: false, covA: 0, covB: 0 };
        const epA = [endsA.p0, endsA.p1];
        const epB = [endsB.p0, endsB.p1];
        const coverage = (S0, S1, X0, X1) => {
            const axis = S1.clone().sub(S0);
            const L = axis.length();
            if (L < 1e-9)
                return 0;
            axis.divideScalar(L);
            const t0 = X0.clone().sub(S0).dot(axis);
            const t1 = X1.clone().sub(S0).dot(axis);
            const lo = Math.max(0, Math.min(t0, t1));
            const hi = Math.min(L, Math.max(t0, t1));
            return Math.max(0, hi - lo) / L;
        };
        const covA = coverage(epA[0], epA[1], epB[0], epB[1]);
        const covB = coverage(epB[0], epB[1], epA[0], epA[1]);
        return { overlap: Math.max(covA, covB) >= threshold, covA, covB };
    }
    toscad.kmAxisShadowOverlap = kmAxisShadowOverlap;
    // "Cosine" between A's and B's crossover-count vectors over the union of their neighbors.
    function kmHelixPairCosine(a, b, cc) {
        const nbA = cc.get(a), nbB = cc.get(b);
        if (!nbA || !nbB)
            return 0;
        const universe = new Set();
        for (const n of nbA.keys())
            if (n !== a && n !== b)
                universe.add(n);
        for (const n of nbB.keys())
            if (n !== a && n !== b)
                universe.add(n);
        let dot = 0, sqA = 0, sqB = 0;
        for (const i of universe) {
            const wA = nbA.get(i) ?? 0, wB = nbB.get(i) ?? 0;
            dot += wA * wB;
            sqA += wA * wA;
            sqB += wB * wB;
        }
        if (sqA === 0 || sqB === 0)
            return 0;
        // "cosine"
        return Math.max(0, Math.min(1, dot / (Math.sqrt(sqA) * Math.sqrt(sqB))));
    }
    toscad.kmHelixPairCosine = kmHelixPairCosine;
    // Find the overlapping helices in each cell as given by Kruskal's.
    function overlapPairs(kr) {
        const out = [];
        for (const ov of kr.overlaps) {
            const hs = [...ov.helices].sort((x, y) => x - y);
            for (let i = 0; i < hs.length; i++) {
                for (let j = i + 1; j < hs.length; j++) {
                    out.push({ a: hs[i], b: hs[j], cell: ov.cell, cellOccupancy: hs.length });
                }
            }
        }
        // Tighter cells first: a 2-helix cell is an unambiguous pairing, a 4-helix cell is a pile that needs sorting out.
        return out.sort((p, q) => p.cellOccupancy - q.cellOccupancy || p.a - q.a || p.b - q.b);
    }
    toscad.overlapPairs = overlapPairs;
    // Resolve Kruskal's overlaps by using the findNearestPos()
    // `keep` holds helices that must not move (pinned/locked ones). Without it a pin could be
    // honored by kruskals and then quietly relocated here, which would defeat the whole point.
    function posCorr5(kr, helices, keep) {
        const pos = new Map();
        const occupied = new Set();
        for (const [h, p] of kr.positions) {
            pos.set(h, [p[0], p[1]]);
            occupied.add(`${p[0]},${p[1]}`);
        }
        const size = (h) => (helices[h] ?? []).length;
        for (const ov of kr.overlaps) {
            // Pinned helices win outright; otherwise biggest stays put (ties -> lowest id).
            const [, ...movers] = [...ov.helices].sort((x, y) => {
                const kx = keep?.has(x) ? 1 : 0, ky = keep?.has(y) ? 1 : 0;
                return ky - kx || size(y) - size(x) || x - y;
            });
            for (const h of movers) {
                if (keep?.has(h))
                    continue; // never relocate a pinned helix
                const p = toscad.findNearestOpenPos(ov.cell, occupied);
                pos.set(h, p);
                occupied.add(`${p[0]},${p[1]}`);
            }
        }
        return pos;
    }
    toscad.posCorr5 = posCorr5;
    function axisMerge(grid, helices, partials, usedSides, lattice = 'honeycomb', opts = {}) {
        // Currently, hashAxisOverlap is 100% trusted. enforceGates = true allows the code to reject a merge through the gates, even if hashAxisOverlap allows it.
        const ENFORCE = opts.enforceGates ?? false;
        const AXIS_DOT = opts.axisDotThreshold ?? 0.9;
        const SHADOW = opts.shadowThreshold ?? 0.3;
        const mergedPairs = [];
        if (!Array.isArray(partials) || partials.length === 0) {
            console.log('[axisMerge] no partials — nothing to do');
            return { networkMap: toscad.getAngles(grid, helices, lattice), mergedPairs, unpairedPartials: [] };
        }
        const partialEnds = helix.mapPartialEnds(partials);
        if (partialEnds.size === 0) {
            console.log('[axisMerge] no partial ends resolved — nothing to do');
            return { networkMap: toscad.getAngles(grid, helices, lattice), mergedPairs, unpairedPartials: [] };
        }
        const partialAxes = helix.partialAxesTowardFreeSide(partials, partialEnds, usedSides);
        if (partialAxes.size === 0) {
            console.log('[axisMerge] no free-side axes — nothing to do');
            return { networkMap: toscad.getAngles(grid, helices, lattice), mergedPairs, unpairedPartials: [] };
        }
        const hashPairs = helix.hashAxisOverlap(partials, partialEnds, usedSides, partialAxes, opts.hashOpts ?? {});
        console.log(`[axisMerge] hashAxisOverlap nominated ${hashPairs.length} pair(s) from ${partialAxes.size} axes`);
        // Partial index -> current helixId, resolved through the grid so it stays correct as earlier merges shift helixIds.
        const pIdxToHelixId = (pIdx) => {
            const p = partials[pIdx];
            if (!Array.isArray(p) || p.length === 0)
                return -1;
            const mark = grid.get(p[0].id);
            return mark && typeof mark.helixId === 'number' ? mark.helixId : -1;
        };
        const paired = new Set();
        // merge loop
        for (const mp of hashPairs) {
            const hA = pIdxToHelixId(mp.a);
            const hB = pIdxToHelixId(mp.b);
            if (hA < 0 || hB < 0)
                continue;
            if (hA === hB) {
                console.log(`[axisMerge] pIdx ${mp.a}/${mp.b} already in helix ${hA} — skip`);
                paired.add(mp.a);
                paired.add(mp.b);
                continue;
            }
            if (!helices[hA] || !helices[hB])
                continue;
            const axA = kmHelixAxis(grid, helices, hA);
            const axB = kmHelixAxis(grid, helices, hB);
            const axisDot = (axA && axB) ? Math.abs(axA.dot(axB)) : 0;
            const disj = kmDisjoint(grid, helices, hA, hB);
            const shadow = kmAxisShadowOverlap(helices, hA, hB, SHADOW);
            const gatesPassed = disj && axisDot >= AXIS_DOT && !shadow.overlap;
            // By default, ENFORCE is false.
            if (ENFORCE && !gatesPassed) {
                console.warn(`[axisMerge] REJECT (${hA},${hB}) pIdx ${mp.a}/${mp.b} — ` +
                    `disjoint=${disj} axisDot=${axisDot.toFixed(3)} ` +
                    `shadow=${shadow.overlap} (covA=${shadow.covA.toFixed(3)} covB=${shadow.covB.toFixed(3)})`);
                continue;
            }
            const keep = Math.min(hA, hB);
            const merged = Math.max(hA, hB);
            const ids = kmMergeHelixInto(grid, helices, keep, merged);
            mergedPairs.push({
                keepHelix: keep, mergedHelix: merged,
                keepNtIds: ids.keepNtIds, mergedNtIds: ids.mergedNtIds,
                pIdxA: mp.a, pIdxB: mp.b,
                hashDot: mp.dot, hashDistAng: mp.dist,
                axisDot, disjoint: disj, shadowOverlap: shadow.overlap,
                gatesPassed,
            });
            paired.add(mp.a);
            paired.add(mp.b);
            console.log(`[axisMerge] MERGED ${merged}→${keep} pIdx ${mp.a}/${mp.b} ` +
                `hashDot=${mp.dot.toFixed(3)} dist=${mp.dist.toFixed(1)}A ` +
                `axisDot=${axisDot.toFixed(3)} disjoint=${disj} shadow=${shadow.overlap}`);
            if (!gatesPassed && !ENFORCE) {
                console.warn(`[axisMerge] ⚠️ ${merged}→${keep} merged on 3D trust but would FAIL gates — ` +
                    `disjoint=${disj} axisDot=${axisDot.toFixed(3)} shadow=${shadow.overlap}. ` +
                    (!disj ? `Merged helix will have overlapping offsets. ` : ``) +
                    `Set enforceGates to reject these.`);
            }
        }
        const unpairedPartials = [];
        for (let i = 0; i < partials.length; i++)
            if (!paired.has(i))
                unpairedPartials.push(i);
        const networkMap = toscad.getAngles(grid, helices, lattice);
        console.log(`[axisMerge] done — ${mergedPairs.length} merge(s), ` +
            `${unpairedPartials.length} partial(s) left unpaired`);
        return { networkMap, mergedPairs, unpairedPartials };
    }
    toscad.axisMerge = axisMerge;
    function anglecomb5(grid, helices, lattice = 'honeycomb', angleMap, opts = {}) {
        const AXIS_DOT = opts.axisDotThreshold ?? 0.9;
        const SHADOW = opts.shadowThreshold ?? 0.3;
        const STRONG = opts.strongWeight ?? 2;
        const MIN_COS = opts.minCosine ?? 0;
        const REQUIRE_CONF = opts.requireConfirmation ?? false;
        const MAX_ITERATIONS = opts.maxIterations ?? 200;
        const lat = toscad.resolveLatticeKind(lattice);
        const K = latticeDirs(lat).length;
        let networkMap = angleMap ?? toscad.getAngles(grid, helices, lattice);
        const mergedPairs = [];
        let iteration = 0;
        let mergedThisIteration = true;
        const kmMinNtId = (helices, hid) => {
            let m = Number.POSITIVE_INFINITY;
            for (const nt of helices[hid] ?? [])
                if (nt && nt.id < m)
                    m = nt.id;
            return isFinite(m) ? m : -1;
        };
        while (mergedThisIteration && iteration < MAX_ITERATIONS) {
            iteration++;
            mergedThisIteration = false;
            // Rebuild the graph every iteration.
            const cc = toscad.getConnectionCounts(grid);
            const vote = voteOrientations(networkMap, grid, lattice);
            const kr = kruskals(networkMap, grid, lattice, vote);
            // Per-edge trust straight off the converged vote.
            const pairKey = (a, b) => `${Math.min(a, b)}|${Math.max(a, b)}`;
            const satisfiedEdges = new Set();
            for (const e of vote.edges) {
                const ra = vote.slot.get(e.a) ?? 0, rb = vote.slot.get(e.b) ?? 0;
                if (e.compat[ra * K + rb] === 1)
                    satisfiedEdges.add(pairKey(e.a, e.b));
            }
            // Edges that built pos (so they agree by construction, and prove nothing).
            const treeEdges = new Set();
            for (const [a, b] of kr.treeEdges)
                treeEdges.add(pairKey(a, b));
            // Edges whose endpoints' final positions do not satisfy their step.
            const violatedSteps = new Set();
            for (const [a, b] of kr.translationViolations)
                violatedSteps.add(pairKey(a, b));
            // Gate 4.
            const placementVerdict = (h) => {
                let sawTestable = false;
                for (const [n, w] of cc.get(h) ?? []) {
                    if (n === h || w < STRONG)
                        continue;
                    const key = pairKey(h, n);
                    if (!satisfiedEdges.has(key))
                        continue; // vote says untrustworthy: abstain
                    if (treeEdges.has(key))
                        continue; // built pos: no information
                    sawTestable = true;
                    if (violatedSteps.has(key))
                        return 'contradicted';
                }
                return sawTestable ? 'confirmed' : 'unknown';
            };
            const candidates = [];
            const nominated = overlapPairs(kr);
            // Filter to check which ones are allowed to merge.
            for (const { a, b, cell } of nominated) {
                if (a === b || !helices[a] || !helices[b])
                    continue;
                // Gate 1: Disjointness on the grid offsets.
                if (!kmDisjoint(grid, helices, a, b)) {
                    console.log(`[anglecomb5] iter=${iteration} reject (${a},${b}) gate1 offset-collision`);
                    continue;
                }
                // Gate 2: Check for parallel axes.
                const axA = kmHelixAxis(grid, helices, a);
                const axB = kmHelixAxis(grid, helices, b);
                if (!axA || !axB) {
                    console.log(`[anglecomb5] iter=${iteration} reject (${a},${b}) gate2 missing-axis`);
                    continue;
                }
                const axisDot = Math.abs(axA.dot(axB));
                if (axisDot < AXIS_DOT) {
                    console.log(`[anglecomb5] iter=${iteration} reject (${a},${b}) gate2 axisDot=${axisDot.toFixed(3)} < ${AXIS_DOT}`);
                    continue;
                }
                // Gate 3: Check for "side-by-side" overlap - aka, shadow.
                const shadow = kmAxisShadowOverlap(helices, a, b, SHADOW);
                if (shadow.overlap) {
                    console.log(`[anglecomb5] iter=${iteration} reject (${a},${b}) gate3 side-by-side ` +
                        `covA=${shadow.covA.toFixed(3)} covB=${shadow.covB.toFixed(3)}`);
                    continue;
                }
                // Gate 4: Check for placement contradiction.
                const verdictA = placementVerdict(a);
                const verdictB = placementVerdict(b);
                if (verdictA === 'contradicted' || verdictB === 'contradicted') {
                    console.log(`[anglecomb5] iter=${iteration} reject (${a},${b}) gate4 contradicted ` +
                        `A=${verdictA} B=${verdictB}`);
                    continue;
                }
                if (REQUIRE_CONF && (verdictA !== 'confirmed' || verdictB !== 'confirmed')) {
                    console.log(`[anglecomb5] iter=${iteration} reject (${a},${b}) gate4 unconfirmed ` +
                        `A=${verdictA} B=${verdictB} (requireConfirmation on)`);
                    continue;
                }
                const cos = kmHelixPairCosine(a, b, cc);
                if (MIN_COS > 0 && cos < MIN_COS) {
                    console.log(`[anglecomb5] iter=${iteration} reject (${a},${b}) cos=${cos.toFixed(3)} < ${MIN_COS}`);
                    continue;
                }
                const minNtA = kmMinNtId(helices, a);
                const minNtB = kmMinNtId(helices, b);
                if (minNtA < 0 || minNtB < 0)
                    continue;
                candidates.push({
                    a, b, minNtA, minNtB, cell,
                    cos, axisDot, covA: shadow.covA, covB: shadow.covB,
                    verdictA, verdictB,
                });
            }
            if (candidates.length === 0) {
                console.log(`[anglecomb5] iter=${iteration} no surviving candidates — done`);
                break;
            }
            // Ranking only. Confirmed placements first, then cosine, then axis agreement, then lowest id for determinism.
            const confRank = (v) => v === 'confirmed' ? 0 : 1;
            candidates.sort((x, y) => (confRank(x.verdictA) + confRank(x.verdictB)) - (confRank(y.verdictA) + confRank(y.verdictB)) ||
                y.cos - x.cos ||
                y.axisDot - x.axisDot ||
                x.a - y.a || x.b - y.b);
            console.log(`[anglecomb5] iter=${iteration} ${candidates.length} candidate(s): ` +
                candidates.map(c => `(${c.a},${c.b})cos=${c.cos.toFixed(2)}`).join(' '));
            // Drain. One helix per drain.
            const consumed = new Set();
            for (const c of candidates) {
                if (consumed.has(c.minNtA) || consumed.has(c.minNtB)) {
                    console.log(`[anglecomb5] iter=${iteration} defer (${c.a},${c.b}) — helix consumed this drain`);
                    continue;
                }
                // helixIds shift on every splice; re-resolve through nt ids.
                const currA = grid.get(c.minNtA)?.helixId ?? -1;
                const currB = grid.get(c.minNtB)?.helixId ?? -1;
                if (currA < 0 || currB < 0 || currA === currB)
                    continue;
                if (!helices[currA] || !helices[currB])
                    continue;
                const keep = Math.min(currA, currB);
                const merged = Math.max(currA, currB);
                const ids = kmMergeHelixInto(grid, helices, keep, merged);
                mergedPairs.push({
                    keepHelix: keep, mergedHelix: merged,
                    keepNtIds: ids.keepNtIds, mergedNtIds: ids.mergedNtIds,
                    cell: c.cell, cos: c.cos, axisDot: c.axisDot,
                    covA: c.covA, covB: c.covB,
                    verdictA: c.verdictA, verdictB: c.verdictB,
                    iteration,
                });
                console.log(`[anglecomb5] iter=${iteration} MERGED ${merged}→${keep} ` +
                    `cell=(${c.cell[0]},${c.cell[1]}) cos=${c.cos.toFixed(3)} ` +
                    `axisDot=${c.axisDot.toFixed(3)} cov=(${c.covA.toFixed(2)},${c.covB.toFixed(2)}) ` +
                    `placement=(${c.verdictA},${c.verdictB})`);
                consumed.add(c.minNtA);
                consumed.add(c.minNtB);
                mergedThisIteration = true;
            }
            if (mergedThisIteration)
                networkMap = toscad.getAngles(grid, helices, lattice);
        }
        if (iteration >= MAX_ITERATIONS) {
            console.warn(`[anglecomb5] hit MAX_ITERATIONS=${MAX_ITERATIONS}`);
        }
        const unknowns = mergedPairs.filter(m => m.verdictA === 'unknown' || m.verdictB === 'unknown').length;
        console.log(`[anglecomb5] done — ${mergedPairs.length} merge(s) over ${iteration} iteration(s), ` +
            `${unknowns} involved an unconfirmable placement`);
        return { networkMap, mergedPairs, iterations: iteration };
    }
    toscad.anglecomb5 = anglecomb5;
})(toscad || (toscad = {}));
