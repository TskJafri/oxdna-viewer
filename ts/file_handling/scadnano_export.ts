/// <reference path="../typescript_definitions/index.d.ts" />

type ScadnanoGridType = 'honeycomb' | 'square';
// What the dialog/grid-type dropdown can hand us. 'automatic' is resolved to
// one of the concrete kinds inside prepareScadnanoLayout via toscad.detectLatticeKind.
type ScadnanoRequestedGridType = ScadnanoGridType | 'automatic';
type HelixPosMap = Map<number, [number, number]>;

type ScadnanoPreparedLayout = {
    latticeType: ScadnanoGridType;
    nucleotideCount: number;
    helices: Nucleotide[][];
    grid: any;
    helixPos: HelixPosMap;
    wireframe: boolean;
};

type ScadnanoDialogOptions = {
    name: string;
    gridType: ScadnanoRequestedGridType;
    includeHelixPos: boolean;
    wireframe: boolean;
};

interface Window {
    view?: {
        longCalculation: (calc: () => void, message: string, callback?: () => void) => void;
    };
    api?: {
        selectElements?: (elementsToSelect: any) => void;
    };
    scadnano?: any;

    currentScadnanoHelixPos?: HelixPosMap;
    exportScadnanoFromGridView?: (helixPosInput?: unknown) => void;
    showScadnanoGridFromHelixPos?: (helixPosInput?: unknown, gridTypeInput?: unknown) => void;
    hideScadnanoGridPane?: () => void;
    scadnanoDialogExport?: () => void;
    toggleGridDropdown?: (checkboxElement: HTMLInputElement) => void;
    scadnanoSelectHelixFromNucleotide?: (nucleotideInput?: unknown, additive?: boolean) => void;
    scadnanoGetHelices?: () => Nucleotide[][] | null;
    scadnanoGridUndo?: () => void;
    scadnanoGridRedo?: () => void;
    scadnanoGridGetHistory?: () => any;
}

// ── Edit-history journal types ───────────────────────────────────────────────
// Stored as plain JSON-friendly objects so the journal stays readable / editable.
type ScadnanoMoveEntry = {
    op: 'move';
    id: number;
    from: [number, number];
    to:   [number, number];
};

type ScadnanoCombineRemoved = {
    oldIdx: number;
    col:    number;
    row:    number;
    ntIds:  number[];
};

type ScadnanoCombineEntry = {
    op: 'combine';
    kept: number;
    indices: number[];
    // Map<oldIdx, newIdx> for survivors only, serialized as pairs for JSON.
    idRemap: Array<[number, number]>;
    removed: ScadnanoCombineRemoved[];
    connections: {
        before: Array<[number, number]>;
        after:  Array<[number, number]>;
    };
    // Per-helix offset shifts applied to grid marks BEFORE the merge, keyed by pre-merge helix id.
    // Empty when the selected helices were already disjoint. Undo reverses these AFTER splitting,
    // redo reapplies them BEFORE combining. Serialized as [oldIdx, delta] pairs for JSON safety.
    preCombineShifts: Array<[number, number]>;
};

type ScadnanoHistoryEntry = ScadnanoMoveEntry | ScadnanoCombineEntry;

type ScadnanoHistoryJournal = {
    entries: { [key: string]: ScadnanoHistoryEntry };
    order:   string[];
    cursor:  number;
    nextMove:    number;
    nextCombine: number;
};

class ScadnanoExportManager {
    private currentScadnanoHelices: Nucleotide[][] | null = null;
    private currentScadnanoConnections: Array<[number, number]> = [];
    private currentScadnanoLayout: ScadnanoPreparedLayout | null = null;

    // Cached findHelices() output. Populated lazily by calculateScadnanoHelices and
    // reused by the export pipeline so the partials + usedSides ledger computed by
    // generateHelix is available to applyAxisOverlapMerge without re-walking the strands.
    private currentScadnanoPartials: Nucleotide[][] | null = null;
    private currentScadnanoUsedSides: Map<number, Map<number, number>> | null = null;

    private scadnanoGridEditor: any = null;
    private scadnanoGridEditorType: ScadnanoGridType | null = null;
    private suppressNodeSelectedCallback = false;

    // Edit-history journal. Lives only while the grid pane is open; cleared on export / close /
    // reload. Holds combine and move operations in the readable JSON form discussed.
    private history: ScadnanoHistoryJournal = this.createEmptyHistory();

    private createEmptyHistory(): ScadnanoHistoryJournal {
        return { entries: {}, order: [], cursor: 0, nextMove: 0, nextCombine: 0 };
    }

    // Helix ids that have been locked by the user. Persists for the lifetime of the grid pane
    // session; cleared when the pane is closed or the grid is reloaded.
    private lockedHelices: Set<number> = new Set();

    private clearHistory(): void {
        this.history = this.createEmptyHistory();
        this.refreshHistoryButtons();
    }

    // Drop any redo-tail entries before a new operation extends the journal.
    private truncateRedoTail(): void {
        if (this.history.cursor >= this.history.order.length) return;
        const tail = this.history.order.splice(this.history.cursor);
        tail.forEach(key => { delete this.history.entries[key]; });
    }

    private pushHistoryEntry(entry: ScadnanoHistoryEntry): string {
        this.truncateRedoTail();
        const key = entry.op === 'move'
            ? `m${++this.history.nextMove}`
            : `c${++this.history.nextCombine}`;
        this.history.entries[key] = entry;
        this.history.order.push(key);
        this.history.cursor = this.history.order.length;
        this.refreshHistoryButtons();
        return key;
    }

    private canUndo(): boolean { return this.history.cursor > 0; }
    private canRedo(): boolean { return this.history.cursor < this.history.order.length; }

    private refreshHistoryButtons(): void {
        const undoBtn = document.getElementById('scadnanoGridUndoBtn') as HTMLButtonElement | null;
        const redoBtn = document.getElementById('scadnanoGridRedoBtn') as HTMLButtonElement | null;
        if (undoBtn) undoBtn.disabled = !this.canUndo();
        if (redoBtn) redoBtn.disabled = !this.canRedo();
    }

    constructor() {
        this.initScadnanoGridPaneControls();
    }

    public handleDialogExport(): void {
        const options = this.readDialogOptions();
        if (!options) return;

        this.closeScadnanoDialog();
        requestAnimationFrame(() => {
            requestAnimationFrame(() => {
                this.runDialogExport(options);
            });
        });
    }

    public exportFromGridView(helixPosInput?: unknown): void {
        const { name, gridType, wireframe } = this.readCurrentExportTarget();
        const map = this.normalizeHelixPosMap(helixPosInput ?? window.currentScadnanoHelixPos);

        if (!map || map.size === 0) {
            notify('No edited helix positions available to export.', 'warning');
            return;
        }

        try {
            this.exportToScadnano(name, gridType, map, wireframe);
        } catch (err) {
            notify(`Scadnano export failed: ${err}`, 'alert');
        }
    }

    // Triggered by the "Recalculate Grid" button. Reruns the full layout pipeline on the
    // Triggered by the "Recalculate Grid" button. Reruns the layout pipeline on the
    // existing helices and grid (as the user has combined/edited them) — skipping setGrid.
    // Pipeline: directionAlign2 → alignGridPrim → detectLatticeKind →
    //           getAngles → anglecomb → anglecorr → calculateGlobalPositions →
    //           renumberHelicesGNN → applyHelixRenumber → collectCrossovers.
    // Locked helix positions are tracked by nucleotide membership so they survive renumbering.
    private recalculateGridFromScratch(): void {
        if (!this.currentScadnanoLayout) {
            notify('Open the grid view first before recalculating.', 'warning');
            return;
        }

        // Save one nucleotide id per locked helix as an anchor. Nucleotide ids are stable
        // across any renumbering, so we can find the new helix index after the pipeline runs.
        const lockedNtAnchors: number[] = [];
        if (this.lockedHelices.size > 0) {
            const oldHelices = this.currentScadnanoLayout.helices;
            this.lockedHelices.forEach(helixId => {
                const slot = oldHelices[helixId];
                if (Array.isArray(slot) && slot.length > 0) {
                    lockedNtAnchors.push(slot[0].id);
                }
            });
        }

        const { gridType, wireframe } = this.readCurrentExportTarget();
        let failed = false;

        this.runScadnanoLongCalculation(
            () => {
                try {
                    // Reuse the existing grid — setGrid is skipped entirely.
                    // Shallow-clone helices so anglecomb's in-place splices don't
                    // mutate currentScadnanoLayout.helices, keeping recalculate idempotent.
                    const grid = this.currentScadnanoLayout!.grid;
                    const helices: Nucleotide[][] = this.currentScadnanoLayout!.helices
                        .map(slot => slot.slice());

                    // binderHelices are not stored on the layout; pass empty so
                    // alignGridPrim and detectLatticeKind treat all helices as lattice members.
                    const binderHelices: number[] = [];

                    toscad.directionAlign2(grid);
                    toscad.alignGridPrim(grid, binderHelices);

                    let latticeType: ScadnanoGridType;
                    if (gridType === 'automatic') {
                        latticeType = toscad.detectLatticeKind(grid, binderHelices);
                        notify(`Auto-detected lattice: ${latticeType}`, 'success');
                    } else {
                        latticeType = gridType;
                    }

                    const angles = toscad.getAngles(grid, helices, latticeType);
                    let networkMap = angles;

                    if (!wireframe) {
                        const corrected = toscad.anglecomb(grid, helices, latticeType, angles);
                        const correct = toscad.anglecorr(grid, helices, latticeType, corrected.networkMap);
                        networkMap = correct.networkMap;

                        // Sweep for any axis-overlapping helices the angle chain missed.
                        // Reuses the partials + usedSides ledger built by findHelices → generateHelix.
                        const axisMerges = this.runAxisOverlapMerge(helices, grid, latticeType);
                        if (axisMerges > 0) {
                            networkMap = toscad.getAngles(grid, helices, latticeType);
                        }
                    }

                    let helixPos = toscad.calculateGlobalPositions(networkMap, undefined, undefined, latticeType);

                    const renumber = toscad.renumberHelicesGNN(grid, helixPos, latticeType);
                    const renumbered = toscad.applyHelixRenumber(helices, grid, helixPos, renumber.remap);
                    helixPos = renumbered.helixPos;

                    const { crossovers } = toscad.collectCrossovers(grid);
                    this.currentScadnanoConnections = this.buildScadnanoConnections(crossovers);

                    this.currentScadnanoHelices = renumbered.helices;
                    this.currentScadnanoLayout = {
                        latticeType,
                        nucleotideCount: this.currentScadnanoLayout!.nucleotideCount,
                        helices: renumbered.helices,
                        grid,
                        helixPos,
                        wireframe
                    };
                } catch (err) {
                    failed = true;
                    notify(`Recalculate grid failed: ${err}`, 'alert');
                }
            },
            () => {
                if (failed || !this.currentScadnanoLayout) return;

                const newLayout = this.currentScadnanoLayout;

                // Remap locked helix ids through the new numbering using nucleotide anchors.
                const newLocked = new Set<number>();
                lockedNtAnchors.forEach(ntId => {
                    const newHelixId = toscad.findHelixID(ntId, newLayout.helices);
                    if (newHelixId !== null) newLocked.add(newHelixId);
                });

                const helixPos = this.cloneHelixPosMap(newLayout.helixPos);
                window.currentScadnanoHelixPos = helixPos;

                // Reload the editor. showGridFromHelixPos clears history and lockedHelices,
                // so we restore the remapped set immediately after.
                this.showGridFromHelixPos(helixPos, newLayout.latticeType);
                this.lockedHelices = newLocked;
                this.applyLockedColors();

                notify('Grid recalculated.', 'success');
            }
        );
    }

    // Triggered by the "Lock/Unlock Helices" button in the grid view.
    // Toggles the locked state of every currently-selected helix:
    //   - Unlocked → grey fill, added to lockedHelices.
    //   - Already locked → yellow fill restored, removed from lockedHelices.
    private lockSelectedHelices(): void {
        const editor = this.scadnanoGridEditor;
        if (!editor || typeof editor.setNodeColor !== 'function') return;

        const ids: number[] = typeof editor.getSelectedHelixIds === 'function'
            ? editor.getSelectedHelixIds()
            : [];

        if (ids.length === 0) return;

        const LOCKED_COLOR   = 0x808080;
        const UNLOCKED_COLOR = 0xffd400; // DOT_COLOR from scadnano_gridview.ts

        ids.forEach(id => {
            if (this.lockedHelices.has(id)) {
                this.lockedHelices.delete(id);
                editor.setNodeColor(id, UNLOCKED_COLOR);
            } else {
                this.lockedHelices.add(id);
                editor.setNodeColor(id, LOCKED_COLOR);
            }
        });

        // Keep the editor's own set in sync so its drag/move guards work.
        editor.lockedHelices = new Set(this.lockedHelices);
    }

    // Triggered by the "Combine" button in the grid view.
    // Merges the helices currently selected in the grid editor into the lowest-numbered one,
    public combineSelectedHelicesFromGridView(): void {
        const editor = this.scadnanoGridEditor;
        if (!editor) {
            notify('Open the scadnano grid view before combining helices.', 'warning');
            return;
        }

        const ids = typeof editor.getSelectedHelixIds === 'function'
            ? editor.getSelectedHelixIds() as number[]
            : [];

        if (!Array.isArray(ids) || ids.length < 2) {
            notify('Select two or more helices (cmd/ctrl+click) before combining.', 'warning');
            return;
        }

        // Refuse to combine if any of the selected helices are locked.
        const lockedInSelection = ids.filter(id => this.lockedHelices.has(id));
        if (lockedInSelection.length > 0) {
            notify(`Cannot combine: helix ${lockedInSelection.join(', ')} ${lockedInSelection.length === 1 ? 'is' : 'are'} locked. Unlock ${lockedInSelection.length === 1 ? 'it' : 'them'} first.`, 'warning');
            return;
        }

        const helices = this.ensureScadnanoHelicesCache();
        if (!helices) {
            notify('Helix data is not yet available.', 'alert');
            return;
        }

        if (!this.currentScadnanoLayout) {
            notify('Helix layout is not yet available.', 'alert');
            return;
        }

        // Snapshot the pre-merge state of every helix that's about to be removed. This is what
        // makes combine reversible: we record nucleotide ids per slot plus the editor cell each
        // slot occupied. The kept helix doesn't need a snapshot — only its members from the
        // merged-away helices need to be pulled back out on undo.
        const sortedIds = [...new Set(ids)].sort((a, b) => a - b);

        // If the selection contains overlapping helices, ask the shift-solver to compute an
        // offset delta per helix that (a) makes the group mutually disjoint on the grid and
        // (b) minimizes disruption to crossovers with helices outside the reordered subset.
        // Helices that are already disjoint from everyone selected get no shift and stay put.
        const grid = this.currentScadnanoLayout.grid;
        const preCombineShifts = toscad.computeCombineShifts(helices, sortedIds, grid);
        if (preCombineShifts.size > 0) {
            toscad.applyCombineShifts(grid, preCombineShifts);
        }

        // Verify every pair is now mutually disjoint. If the solver punted (overlap group too
        // large, or fundamentally unresolvable) some pair will still collide; bail with the same
        // error message so the user knows the automatic path didn't fix things.
        for (let i = 0; i < sortedIds.length; i++) {
            for (let j = i + 1; j < sortedIds.length; j++) {
                const h1 = sortedIds[i];
                const h2 = sortedIds[j];
                if (!toscad.disjoint(helices, h1, h2, grid)) {
                    // Roll back any shifts we applied so the grid stays consistent with the
                    // pre-combine state the user was looking at.
                    if (preCombineShifts.size > 0) {
                        const inverse = new Map<number, number>();
                        preCombineShifts.forEach((delta, hid) => inverse.set(hid, -delta));
                        toscad.applyCombineShifts(grid, inverse);
                    }
                    notify(
                        `Cannot combine: helices ${h1} and ${h2} are not disjoint (overlapping nucleotides on the grid).`,
                        'alert'
                    );
                    return;
                }
            }
        }

        const willRemoveOld = sortedIds.slice(1);
        const editorNodesBefore = typeof editor.getNodes === 'function'
            ? editor.getNodes() as Array<{ id: number; col: number; row: number }>
            : [];
        const cellByOldId = new Map<number, [number, number]>();
        editorNodesBefore.forEach(n => cellByOldId.set(Number(n.id), [Number(n.col), Number(n.row)]));

        const removedSnapshots: ScadnanoCombineRemoved[] = willRemoveOld.map(oldIdx => {
            const cell = cellByOldId.get(oldIdx) ?? [0, 0];
            const slot = helices[oldIdx] || [];
            return {
                oldIdx,
                col: cell[0],
                row: cell[1],
                ntIds: slot.map(nt => nt.id)
            };
        });
        const connectionsBefore = this.currentScadnanoConnections.map(([a, b]) => [a, b] as [number, number]);

        // combineHelices mutates `helices` AND `currentScadnanoLayout.grid` in place.
        const result = helix.combineHelices(helices, ids, this.currentScadnanoLayout.grid);
        if (!result) {
            notify('Nothing to combine.', 'warning');
            return;
        }

        const { keptIdx, mergedIdx, idRemap } = result;

        // Cascade the merge through the editor (nodes + connections), helixPos, and selection.
        this.applyCombineCascade(result);

        // Record the combine in the history journal so it's reversible.
        const entry: ScadnanoCombineEntry = {
            op: 'combine',
            kept: keptIdx,
            indices: sortedIds,
            idRemap: Array.from(idRemap.entries()),
            removed: removedSnapshots,
            connections: {
                before: connectionsBefore,
                after:  this.currentScadnanoConnections.map(([a, b]) => [a, b] as [number, number])
            },
            preCombineShifts: Array.from(preCombineShifts.entries())
        };
        this.pushHistoryEntry(entry);

        if (preCombineShifts.size > 0) {
            const shiftSummary = Array.from(preCombineShifts.entries())
                .map(([hid, d]) => `${hid}:${d > 0 ? '+' : ''}${d}`)
                .join(', ');
            notify(`Combined helix ${mergedIdx.join(', ')} into helix ${keptIdx}. Shifted overlapping helices to make them disjoint: ${shiftSummary}.`);
        } else {
            notify(`Combined helix ${mergedIdx.join(', ')} into helix ${keptIdx}.`);
        }
    }

    // ── Undo / Redo ─────────────────────────────────────────────────────────
    public undoFromGridView(): void {
        if (!this.canUndo()) return;
        const key = this.history.order[this.history.cursor - 1];
        const entry = this.history.entries[key];
        if (!entry) return;

        if (entry.op === 'move')    this.undoMove(entry);
        else if (entry.op === 'combine') this.undoCombine(entry);

        this.history.cursor -= 1;
        this.refreshHistoryButtons();
    }

    public redoFromGridView(): void {
        if (!this.canRedo()) return;
        const key = this.history.order[this.history.cursor];
        const entry = this.history.entries[key];
        if (!entry) return;

        if (entry.op === 'move')    this.redoMove(entry);
        else if (entry.op === 'combine') this.redoCombine(entry);

        this.history.cursor += 1;
        this.refreshHistoryButtons();
    }

    public getHistorySnapshot(): ScadnanoHistoryJournal {
        // Returned by reference for inspection / debugging. Don't mutate from outside.
        return this.history;
    }

    private redoMove(entry: ScadnanoMoveEntry): void {
        const editor = this.scadnanoGridEditor;
        if (!editor || typeof editor.moveNodeById !== 'function') return;
        editor.moveNodeById(entry.id, entry.to[0], entry.to[1]);
        this.publishCurrentHelixPosFromEditor();
        this.syncLayoutHelixPosFromEditor();
    }

    private undoMove(entry: ScadnanoMoveEntry): void {
        const editor = this.scadnanoGridEditor;
        if (!editor || typeof editor.moveNodeById !== 'function') return;
        editor.moveNodeById(entry.id, entry.from[0], entry.from[1]);
        this.publishCurrentHelixPosFromEditor();
        this.syncLayoutHelixPosFromEditor();
    }

    // Replay a combine forward by re-running combineHelices on the saved indices, then handing
    // off to the same cascade helper the live combine path uses. The only difference is that
    // redo applies the cached post-merge connection list directly (passed via override), instead
    // of re-deriving it through the remap.
    private redoCombine(entry: ScadnanoCombineEntry): void {
        const editor = this.scadnanoGridEditor;
        if (!editor) return;
        const helices = this.ensureScadnanoHelicesCache();
        if (!helices) return;
        if (!this.currentScadnanoLayout) return;

        // Reapply the pre-combine offset shifts BEFORE combining, mirroring the live path.
        // Keyed by pre-merge helix id, which is what the grid marks still carry at this point.
        if (Array.isArray(entry.preCombineShifts) && entry.preCombineShifts.length > 0) {
            const shifts = new Map<number, number>(entry.preCombineShifts);
            toscad.applyCombineShifts(this.currentScadnanoLayout.grid, shifts);
        }

        const result = helix.combineHelices(helices, entry.indices, this.currentScadnanoLayout.grid);
        if (!result) return;

        this.applyCombineCascade(result, entry.connections.after);
    }

    // Reverse a combine: re-insert the removed helix slots, pull their nucleotides back out of the
    // kept helix, then renumber everything in the GridMap and editor through the inverse remap.
    private undoCombine(entry: ScadnanoCombineEntry): void {
        const editor = this.scadnanoGridEditor;
        if (!editor) return;
        const helices = this.ensureScadnanoHelicesCache();
        if (!helices) return;
        if (!this.currentScadnanoLayout) return;

        // splitHelices mutates `helices` AND `currentScadnanoLayout.grid` in place, mirroring
        // what combineHelices did on the forward path. It hands back the inverse remap the
        // cascade needs to renumber editor nodes that are still in post-merge id space.
        const result = helix.splitHelices(helices, this.currentScadnanoLayout.grid, entry);
        if (!result) return;

        // AFTER splitting, grid marks are back in their pre-combine helix ids, but any offset
        // shifts applied by the shift-solver are still present. Reverse them with negated deltas.
        if (Array.isArray(entry.preCombineShifts) && entry.preCombineShifts.length > 0) {
            const inverse = new Map<number, number>();
            entry.preCombineShifts.forEach(([hid, delta]) => inverse.set(hid, -delta));
            toscad.applyCombineShifts(this.currentScadnanoLayout.grid, inverse);
        }

        this.applyCombineCascadeInverse(entry, result.inverseRemap);
    }

    private remapCachedGridIds(remap: (oldId: number) => number | null): void {
        if (!this.currentScadnanoLayout) return;
        this.currentScadnanoLayout.grid.forEach((mark: { helixId: number; }) => {
            const newId = remap(mark.helixId);
            if (newId === null) return;
            mark.helixId = newId;
        });
    }

    // Run after combineHelices has already mutated the helices array AND the cached grid in place.
    // Cascades the merge through the visual editor (nodes + connections), the cached helixPos,
    // selection state, and the published helix-pos map. Single entry point so a debugger breakpoint
    // here covers every downstream side effect of a combine.
    //
    // `connectionsOverride` is for redo: live combines remap the current connection list through
    // `remapId`, but redo already has the post-merge connection list cached on the journal entry
    // (`entry.connections.after`), so it passes that in directly to skip the remap.
    private applyCombineCascade(
        result: { keptIdx: number; mergedIdx: number[]; idRemap: Map<number, number> },
        connectionsOverride?: Array<[number, number]>
    ): void {
        const editor = this.scadnanoGridEditor;
        if (!editor) return;

        const { keptIdx, mergedIdx, idRemap } = result;
        const removed = new Set<number>(mergedIdx);
        const remapId = (oldId: number): number | null => {
            if (removed.has(oldId)) return keptIdx;
            const next = idRemap.get(oldId);
            return next === undefined ? null : next;
        };

        // Renumber editor nodes through the remap, drop merged-away duplicates. Survivors keep
        // their col/row — the kept helix stays in its original cell.
        if (typeof editor.getNodes === 'function') {
            const nodes = editor.getNodes() as Array<{ id: number; col: number; row: number; label?: string; color?: number }>;
            const remapped: Array<{ id: number; col: number; row: number; label?: string; color?: number }> = [];
            const seenNew = new Set<number>();
            nodes.forEach(node => {
                const oldId = Number(node.id);
                const newId = remapId(oldId);
                if (newId === null) return;
                if (seenNew.has(newId)) return; // skip duplicates (e.g. merged-away nodes mapping to keptIdx)
                seenNew.add(newId);
                remapped.push({ ...node, id: newId, label: String(newId) });
            });
            if (typeof editor.setNodes === 'function') editor.setNodes(remapped);
        }

        // Connections: redo applies the cached post-merge list verbatim; live combines remap the
        // current list (drop self-loops + dedup).
        if (connectionsOverride) {
            this.currentScadnanoConnections = connectionsOverride.map(([a, b]) => [a, b] as [number, number]);
        } else {
            const remappedConnKeys = new Set<string>();
            const updatedConnections: Array<[number, number]> = [];
            this.currentScadnanoConnections.forEach(([from, to]) => {
                const a = remapId(from);
                const b = remapId(to);
                if (a === null || b === null) return;
                if (a === b) return;
                const key = a < b ? `${a}:${b}` : `${b}:${a}`;
                if (remappedConnKeys.has(key)) return;
                remappedConnKeys.add(key);
                updatedConnections.push([a, b]);
            });
            this.currentScadnanoConnections = updatedConnections;
        }
        if (typeof editor.setConnections === 'function') {
            editor.setConnections(this.currentScadnanoConnections);
        }

        // helixPos is rebuilt from the (now post-merge) editor node set.
        this.syncLayoutHelixPosFromEditor();

        if (typeof editor.clearSelection === 'function') editor.clearSelection();

        // Refresh the published helix-pos map from the editor (it just lost a few nodes).
        this.publishCurrentHelixPosFromEditor();

        // Remap lockedHelices through the same id transformation so locks follow their nodes.
        this.remapLockedHelicesForward(remapId);
    }

    // Inverse of applyCombineCascade. Run after splitHelices has already restored helices + grid
    // back to their pre-merge state. Renumbers existing editor nodes through the inverse remap,
    // re-adds the removed nodes at their saved cells, and restores the pre-merge connection list
    // from the journal entry.
    private applyCombineCascadeInverse(
        entry: ScadnanoCombineEntry,
        inverseRemap: (currentId: number) => number
    ): void {
        const editor = this.scadnanoGridEditor;
        if (!editor) return;

        // Renumber existing editor nodes (still in post-merge id space) back to old ids, then
        // re-insert the removed nodes at their captured cells.
        const currentNodes = typeof editor.getNodes === 'function'
            ? editor.getNodes() as Array<{ id: number; col: number; row: number; label?: string; color?: number }>
            : [];
        const restoredNodes: Array<{ id: number; col: number; row: number; label?: string; color?: number }> = [];
        currentNodes.forEach(node => {
            const oldId = inverseRemap(Number(node.id));
            restoredNodes.push({ ...node, id: oldId, label: String(oldId) });
        });
        entry.removed.forEach(slot => {
            restoredNodes.push({ id: slot.oldIdx, col: slot.col, row: slot.row, label: String(slot.oldIdx) });
        });
        if (typeof editor.setNodes === 'function') editor.setNodes(restoredNodes);

        // Restore the pre-merge connection list verbatim.
        this.currentScadnanoConnections = entry.connections.before.map(([a, b]) => [a, b] as [number, number]);
        if (typeof editor.setConnections === 'function') {
            editor.setConnections(this.currentScadnanoConnections);
        }

        this.syncLayoutHelixPosFromEditor();
        this.publishCurrentHelixPosFromEditor();
        if (typeof editor.clearSelection === 'function') editor.clearSelection();

        // Restore the locked set back to pre-merge ids and reapply colors.
        this.remapLockedHelicesInverse(inverseRemap, entry);
    }

    // Remap lockedHelices forward after a combine: each old id is mapped through remapId to its
    // new post-merge id. If a merged-away helix was locked, its lock carries over to the kept
    // helix (keptIdx). Reapplies grey fill to all still-locked nodes.
    private remapLockedHelicesForward(remapId: (oldId: number) => number | null): void {
        if (this.lockedHelices.size === 0) return;

        const next = new Set<number>();
        this.lockedHelices.forEach(id => {
            const newId = remapId(id);
            if (newId !== null) next.add(newId);
        });
        this.lockedHelices = next;
        this.applyLockedColors();
    }

    // Remap lockedHelices back to pre-merge ids after an undo-combine. inverseRemap maps the
    // current post-merge id → the original pre-merge id. Reapplies grey fill to locked nodes
    // and restores yellow on the re-split nodes that are no longer locked.
    private remapLockedHelicesInverse(
        inverseRemap: (currentId: number) => number,
        entry: ScadnanoCombineEntry
    ): void {
        if (this.lockedHelices.size === 0 && entry.removed.length === 0) return;

        const next = new Set<number>();
        this.lockedHelices.forEach(id => {
            next.add(inverseRemap(id));
        });
        this.lockedHelices = next;
        this.applyLockedColors();
    }

    // Reapply the locked (grey) color to every node currently in lockedHelices, and restore
    // the default yellow to every other node. Called after any operation that rebuilds the
    // editor node set (combine, undo, redo).
    private applyLockedColors(): void {
        const editor = this.scadnanoGridEditor;
        if (!editor || typeof editor.setNodeColor !== 'function') return;

        const LOCKED_COLOR   = 0x808080;
        const UNLOCKED_COLOR = 0xffd400;

        const nodes: Array<{ id: number }> = typeof editor.getNodes === 'function'
            ? editor.getNodes()
            : [];

        nodes.forEach(node => {
            const id = Number(node.id);
            editor.setNodeColor(id, this.lockedHelices.has(id) ? LOCKED_COLOR : UNLOCKED_COLOR);
        });

        // Keep the editor's drag/move guard in sync.
        editor.lockedHelices = new Set(this.lockedHelices);
    }

    private syncLayoutHelixPosFromEditor(): void {
        if (!this.currentScadnanoLayout) return;
        const editor = this.scadnanoGridEditor;
        if (!editor || typeof editor.getNodes !== 'function') return;
        const refreshed = new Map<number, [number, number]>();
        (editor.getNodes() as Array<{ id: number; col: number; row: number }>).forEach(node => {
            refreshed.set(Number(node.id), [Number(node.col), Number(node.row)]);
        });
        this.currentScadnanoLayout.helixPos = refreshed;
    }

    public showGridFromHelixPos(helixPosInput?: unknown, gridTypeInput?: unknown): void {
        const pane = this.getScadnanoGridPane();
        if (!pane) {
            notify('Scadnano grid pane is unavailable.', 'alert');
            return;
        }

        const gridType: ScadnanoGridType = gridTypeInput === 'square' ? 'square' : 'honeycomb';

        document.body.classList.add('scadnano-grid-open');
        this.resizeScadnanoGridCanvas();

        const editor = this.ensureScadnanoGridEditor(gridType);
        if (!editor) {
            notify('Unable to open scadnano grid view.', 'alert');
            return;
        }

        if (typeof editor.resize === 'function') {
            editor.resize();
        }

        const map = this.normalizeHelixPosMap(helixPosInput);
        if (!map || map.size === 0) {
            notify('No helix positions available for the grid view.', 'warning');
            return;
        }

        editor.loadFromHelixPos(map);
        if (typeof editor.setConnections === 'function') {
            editor.setConnections(this.currentScadnanoConnections);
        }
        this.publishCurrentHelixPosFromEditor();

        // Fresh grid view = fresh history and no locks.
        this.lockedHelices.clear();
        if (this.scadnanoGridEditor) this.scadnanoGridEditor.lockedHelices = new Set();
        this.clearHistory();
    }

    public hideScadnanoGridPane(): void {
        document.body.classList.remove('scadnano-grid-open');
        this.lockedHelices.clear();
        if (this.scadnanoGridEditor) this.scadnanoGridEditor.lockedHelices = new Set();
        this.clearHistory();
    }

    public toggleGridDropdown(checkboxElement: HTMLInputElement): void {
        const gridDropdown = document.getElementById('scadnanoGrid') as HTMLSelectElement | null;
        if (gridDropdown) {
            gridDropdown.disabled = !checkboxElement.checked;
        }
    }

    public getHelices(): Nucleotide[][] | null {
        return this.ensureScadnanoHelicesCache();
    }

    public selectHelixFromNucleotide(nucleotideInput?: unknown, additive: boolean = false): void {
        if (!document.body.classList.contains('scadnano-grid-open')) return;
        if (!this.scadnanoGridEditor || typeof this.scadnanoGridEditor.selectNodeById !== 'function') return;

        const nucleotide = nucleotideInput instanceof Nucleotide
            ? nucleotideInput
            : (typeof nucleotideInput === 'number' ? elements.get(nucleotideInput) : null);

        if (!nucleotide || !(nucleotide instanceof Nucleotide)) return;

        const helices = this.ensureScadnanoHelicesCache();
        if (!helices) return;

        const helixId = toscad.findHelixID(nucleotide.id, helices);
        if (helixId === null) return;

        this.suppressNodeSelectedCallback = true;
        try {
            this.scadnanoGridEditor.selectNodeById(helixId, additive);
        } finally {
            this.suppressNodeSelectedCallback = false;
        }
    }

    private runDialogExport(options: ScadnanoDialogOptions): void {
        if (!options.includeHelixPos) {
            this.runScadnanoLongCalculation(() => {
                try {
                    this.exportToScadnano(options.name, options.gridType, undefined, options.wireframe);
                } catch (err) {
                    notify(`Scadnano export failed: ${err}`, 'alert');
                }
            });
            return;
        }

        let helixPos: HelixPosMap | null = null;
        let resolvedGridType: ScadnanoGridType = 'honeycomb';
        let failed = false;

        this.runScadnanoLongCalculation(
            () => {
                try {
                    const result = this.calculateScadnanoHelixPos(options.gridType, options.wireframe);
                    helixPos = result.helixPos;
                    resolvedGridType = result.latticeType;
                    if (helixPos) {
                        window.currentScadnanoHelixPos = this.cloneHelixPosMap(helixPos);
                    }
                } catch (err) {
                    failed = true;
                    notify(`Scadnano export failed: ${err}`, 'alert');
                }
            },
            () => {
                if (failed || !helixPos) return;
                this.showGridFromHelixPos(helixPos, resolvedGridType);
            }
        );
    }

    private readDialogOptions(): ScadnanoDialogOptions | null {
        const nameInput = document.getElementById('scadnanoFilename') as HTMLInputElement | null;
        const helixPosCheckbox = document.getElementById('scadnanoIncludeHPos') as HTMLInputElement | null;
        const scadnanoGrid = document.getElementById('scadnanoGrid') as HTMLInputElement | null;
        const wireframeCheckbox = document.getElementById('scadnanoWireframe') as HTMLInputElement | null;

        if (!nameInput || !helixPosCheckbox || !scadnanoGrid || !wireframeCheckbox) {
            console.warn('scadnano export dialog missing inputs');
            return null;
        }

        return {
            name: nameInput.value.trim() || 'output',
            gridType: this.normalizeRequestedGridType(scadnanoGrid.value),
            includeHelixPos: helixPosCheckbox.checked,
            wireframe: wireframeCheckbox.checked,
        };
    }

    private readCurrentExportTarget(): { name: string; gridType: ScadnanoRequestedGridType; wireframe: boolean } {
        const nameInput = document.getElementById('scadnanoFilename') as HTMLInputElement | null;
        const scadnanoGrid = document.getElementById('scadnanoGrid') as HTMLInputElement | null;
        const wireframeCheckbox = document.getElementById('scadnanoWireframe') as HTMLInputElement | null;

        return {
            name: nameInput?.value.trim() || 'output',
            gridType: this.normalizeRequestedGridType(scadnanoGrid?.value),
            wireframe: Boolean(wireframeCheckbox?.checked),
        };
    }

    // Dialog/dropdown value — accepts 'automatic' for downstream resolution.
    private normalizeRequestedGridType(value?: string): ScadnanoRequestedGridType {
        if (value === 'square') return 'square';
        if (value === 'honeycomb') return 'honeycomb';
        return 'automatic';
    }

    private runScadnanoLongCalculation(calc: () => void, callback?: () => void): void {
        const longCalculation = window.view?.longCalculation;
        if (typeof longCalculation === 'function') {
            longCalculation(calc, 'Preparing scadnano export, please be patient...', callback);
            return;
        }

        calc();
        if (callback) callback();
    }

    private closeScadnanoDialog(): void {
        let closedByMetro = false;
        const metroDialog = (window as any)?.Metro?.dialog;

        if (metroDialog && typeof metroDialog.close === 'function') {
            try {
                metroDialog.close('#scadnanoDialog');
                closedByMetro = true;
            } catch (err) {
                console.warn('Failed to close scadnano dialog via Metro API:', err);
            }
        }

        if (!closedByMetro) {
            const closeBtn = document.querySelector('#scadnanoDialog .js-dialog-close') as HTMLElement | null;
            if (closeBtn) closeBtn.click();
        }

        const dialogEl = document.getElementById('scadnanoDialog') as HTMLElement | null;
        if (dialogEl) {
            dialogEl.classList.remove('open');
            dialogEl.setAttribute('aria-hidden', 'true');
        }
    }

    private exportToScadnano(
        name: string,
        gridType: ScadnanoRequestedGridType,
        helixPos?: HelixPosMap,
        wireframe = false
    ): void {
        const layout = this.prepareScadnanoLayout(gridType, false, wireframe);
        const resolvedGridType = layout.latticeType;
        const { helices, grid } = layout;

        // Switch to toscad.buildScadnano2 here to fall back to the old
        // topology-driven export. buildScadnano3 reads boundaries from the
        // grid (helixId / direction / offset step) so post-construction
        // grid edits propagate to the export.
        const scadnano = helixPos
            ? toscad.buildScadnano3(grid, helices, resolvedGridType, helixPos)
            : toscad.buildScadnano3(grid, helices, resolvedGridType);

        const fileName = name ? `${name}.sc` : 'output.sc';
        makeTextFile(fileName, JSON.stringify(scadnano, null, 2));
    }

    private getCurrentNucleotideCount(): number {
        let count = 0;
        elements.forEach((element) => {
            if (element instanceof Nucleotide) count += 1;
        });
        return count;
    }

    private cloneHelixPosMap(input: HelixPosMap): HelixPosMap {
        const out = new Map<number, [number, number]>();
        input.forEach((value, key) => {
            out.set(Number(key), [Number(value[0]), Number(value[1])]);
        });
        return out;
    }

    private calculateScadnanoHelices(): Nucleotide[][] {
        const nucleotideElements = new Map<number, Nucleotide>();
        elements.forEach((element, id) => {
            if (element instanceof Nucleotide) {
                nucleotideElements.set(id, element);
            }
        });

        const result = helix.findHelices(nucleotideElements, 3) as {
            helices: Nucleotide[][];
            partials: Nucleotide[][];
            usedSides: Map<number, Map<number, number>>;
        };
        const helices = result?.helices ?? [];
        this.currentScadnanoPartials = result?.partials ?? null;
        this.currentScadnanoUsedSides = result?.usedSides ?? null;
        this.notifyHelixCoverageMismatch(helices, nucleotideElements);
        return helices;
    }

    // Build the partialEnds + partialAxes inputs for hashAxisOverlap from the
    // cached findHelices output, run hashAxisOverlap, and feed the resulting
    // merge pairs into helix.applyAxisOverlapMerge. Returns the number of
    // merges performed (caller refreshes getAngles only when this is > 0).
    private runAxisOverlapMerge(helices: Nucleotide[][], grid: any, latticeType: ScadnanoGridType): number {
        const partials = this.currentScadnanoPartials;
        const usedSides = this.currentScadnanoUsedSides;
        if (!Array.isArray(partials) || !usedSides || partials.length === 0) return 0;

        const partialEnds = helix.mapPartialEnds(partials);
        if (partialEnds.size === 0) return 0;

        const partialAxes = helix.partialAxesTowardFreeSide(partials, partialEnds, usedSides);
        if (partialAxes.size === 0) return 0;

        const mergePairs = helix.hashAxisOverlap(partials, partialEnds, usedSides, partialAxes);
        if (!Array.isArray(mergePairs) || mergePairs.length === 0) return 0;

        const before = helices.length;
        helix.applyAxisOverlapMerge(helices, grid, partials, mergePairs);
        return before - helices.length;
    }

    private prepareScadnanoLayout(
        requestedLatticeType: ScadnanoRequestedGridType,
        forceRecompute = false,
        wireframe = false
    ): ScadnanoPreparedLayout {
        const nucleotideCount = this.getCurrentNucleotideCount();

        // Cache hit:
        //   - For concrete kinds: cached latticeType must match.
        //   - For 'automatic': any cached latticeType is acceptable (it was
        //     either detected the same way last time or explicitly chosen).
        if (
            !forceRecompute &&
            this.currentScadnanoLayout &&
            this.currentScadnanoLayout.nucleotideCount === nucleotideCount &&
            this.currentScadnanoLayout.wireframe === wireframe &&
            (requestedLatticeType === 'automatic' ||
                this.currentScadnanoLayout.latticeType === requestedLatticeType)
        ) {
            this.currentScadnanoHelices = this.currentScadnanoLayout.helices;
            return this.currentScadnanoLayout;
        }

        const helices = this.calculateScadnanoHelices();
        this.currentScadnanoHelices = helices;

        const { grid, binderHelices } = toscad.setGrid(helices);
        toscad.directionAlign2(grid);
        toscad.alignGridPrim(grid, binderHelices);

        // Resolve 'automatic' once the grid is built — detection reads helixId,
        // offset and direction off the grid via crossoverNts. Binder helices are
        // excluded from the lattice vote because their crossover spacing is
        // non-standard (attach points to scaffold from outside the lattice).
        let latticeType: ScadnanoGridType;
        if (requestedLatticeType === 'automatic') {
            latticeType = toscad.detectLatticeKind(grid, binderHelices);
            notify(`Auto-detected lattice: ${latticeType}`, 'success');
        } else {
            latticeType = requestedLatticeType;
        }

        const angles = toscad.getAngles(grid, helices, latticeType);
        let networkMap = angles;

        if (!wireframe) {
            const corrected = toscad.anglecomb(grid, helices, latticeType, angles);
            const correct = toscad.anglecorr(grid, helices, latticeType, corrected.networkMap);
            networkMap = correct.networkMap;

            // Sweep for any axis-overlapping helices the angle chain missed and merge them.
            // Reuses the partials + usedSides ledger built by findHelices → generateHelix.
            // Refresh getAngles afterwards so networkMap reflects the post-merge state.
            const axisMerges = this.runAxisOverlapMerge(helices, grid, latticeType);
            if (axisMerges > 0) {
                networkMap = toscad.getAngles(grid, helices, latticeType);
            }
        }

        // Initial helix grid positions (still using pre-renumber helix IDs).
        let helixPos = toscad.calculateGlobalPositions(networkMap, undefined, undefined, latticeType);

        // Renumber helices so consecutive IDs share crossovers when possible.
        // Runs after combine (anglecomb), so it operates on the final helix
        // count. Mutates `grid` in place; produces a renumbered helices
        // array and a helixPos map keyed by new IDs.
        const renumber = toscad.renumberHelicesGNN(grid, helixPos, latticeType);
        const renumbered = toscad.applyHelixRenumber(helices, grid, helixPos, renumber.remap);
        const finalHelices = renumbered.helices;
        helixPos = renumbered.helixPos;
        this.currentScadnanoHelices = finalHelices;

        // Build connections AFTER renumber so they reference the new IDs.
        const { crossovers } = toscad.collectCrossovers(grid);
        this.currentScadnanoConnections = this.buildScadnanoConnections(crossovers);

        this.currentScadnanoLayout = {
            latticeType,
            nucleotideCount,
            helices: finalHelices,
            grid,
            helixPos,
            wireframe
        };

        return this.currentScadnanoLayout;
    }

    private calculateScadnanoHelixPos(
        latticeType: ScadnanoRequestedGridType = 'automatic',
        wireframe = false
    ): { helixPos: HelixPosMap; latticeType: ScadnanoGridType } {
        const layout = this.prepareScadnanoLayout(latticeType, false, wireframe);
        return {
            helixPos: this.cloneHelixPosMap(layout.helixPos),
            latticeType: layout.latticeType
        };
    }

    private notifyHelixCoverageMismatch(
        helices: Nucleotide[][],
        inputMap: Map<number, Nucleotide>
    ): void {
        const helixCount = helices.flat().length;
        const totalCount = inputMap.size;
        if (helixCount === totalCount) return;

        const missingCount = totalCount - helixCount;

        notify(
            `Helix mapping error: ${helixCount}/${totalCount} nucleotides were mapped. Missing ${missingCount} nucleotides.`,
            'alert',
            true
        );
    }

    private buildScadnanoConnections(
        crossovers: Map<number, Map<number, { sameWalk: number; diffWalk: number }>>
    ): Array<[number, number]> {
        const uniquePairs = new Set<string>();
        const pairs: Array<[number, number]> = [];

        for (const [fromHelix, toMap] of crossovers.entries()) {
            for (const [toHelix, counts] of toMap.entries()) {
                const totalConnections = Number(counts?.sameWalk ?? 0) + Number(counts?.diffWalk ?? 0);
                if (totalConnections <= 0) continue;

                const a = Math.min(fromHelix, toHelix);
                const b = Math.max(fromHelix, toHelix);
                if (a === b) continue;

                const key = `${a}:${b}`;
                if (uniquePairs.has(key)) continue;
                uniquePairs.add(key);
                pairs.push([a, b]);
            }
        }

        return pairs;
    }

    private normalizeHelixPosMap(input: unknown): HelixPosMap | null {
        if (!input) return null;

        if (input instanceof Map) {
            const out = new Map<number, [number, number]>();
            input.forEach((value, key) => {
                if (!Array.isArray(value) || value.length < 2) return;
                const helixId = Number(key);
                const col = Number(value[0]);
                const row = Number(value[1]);
                if (Number.isFinite(helixId) && Number.isFinite(col) && Number.isFinite(row)) {
                    out.set(helixId, [col, row]);
                }
            });
            return out;
        }

        if (Array.isArray(input)) {
            const out = new Map<number, [number, number]>();
            input.forEach((entry: any) => {
                if (!Array.isArray(entry) || entry.length < 2) return;
                const helixId = Number(entry[0]);
                const value = entry[1];
                if (!Array.isArray(value) || value.length < 2) return;
                const col = Number(value[0]);
                const row = Number(value[1]);
                if (Number.isFinite(helixId) && Number.isFinite(col) && Number.isFinite(row)) {
                    out.set(helixId, [col, row]);
                }
            });
            return out;
        }

        if (typeof input === 'object') {
            const out = new Map<number, [number, number]>();
            Object.keys(input as Record<string, any>).forEach((k) => {
                const value = (input as Record<string, any>)[k];
                if (!Array.isArray(value) || value.length < 2) return;
                const helixId = Number(k);
                const col = Number(value[0]);
                const row = Number(value[1]);
                if (Number.isFinite(helixId) && Number.isFinite(col) && Number.isFinite(row)) {
                    out.set(helixId, [col, row]);
                }
            });
            return out;
        }

        return null;
    }

    private getScadnanoGridPane(): HTMLElement | null {
        return document.getElementById('scadnanoGridPane');
    }

    private getScadnanoGridCanvas(): HTMLCanvasElement | null {
        return document.getElementById('scadnanoGridCanvas') as HTMLCanvasElement | null;
    }

    private setScadnanoPaneWidth(widthPx: number): void {
        const minW = 240;
        const maxW = Math.max(minW, Math.floor(window.innerWidth * 0.75));
        const clamped = Math.max(minW, Math.min(maxW, Math.round(widthPx)));
        document.documentElement.style.setProperty('--scadnano-pane-width', `${clamped}px`);
    }

    private resizeScadnanoGridCanvas(): void {
        const pane = this.getScadnanoGridPane();
        const canvas = this.getScadnanoGridCanvas();
        if (!pane || !canvas) return;
        canvas.width = pane.clientWidth;
        canvas.height = pane.clientHeight;
    }

    private mapFromEditorNodes(editor: any): HelixPosMap {
        const out = new Map<number, [number, number]>();
        const nodes = typeof editor.getNodes === 'function' ? editor.getNodes() : [];
        nodes.forEach((node: any) => {
            out.set(Number(node.id), [Number(node.col), Number(node.row)]);
        });
        return out;
    }

    private publishCurrentHelixPosFromEditor(): void {
        if (!this.scadnanoGridEditor) return;
        window.currentScadnanoHelixPos = this.mapFromEditorNodes(this.scadnanoGridEditor);
    }

    private ensureScadnanoHelicesCache(): Nucleotide[][] | null {
        if (this.currentScadnanoHelices && this.currentScadnanoHelices.length > 0) {
            return this.currentScadnanoHelices;
        }

        if (this.currentScadnanoLayout && this.currentScadnanoLayout.helices.length > 0) {
            this.currentScadnanoHelices = this.currentScadnanoLayout.helices;
            return this.currentScadnanoHelices;
        }

        try {
            this.currentScadnanoHelices = this.calculateScadnanoHelices();
            return this.currentScadnanoHelices;
        } catch (err) {
            notify(`Unable to map grid helix selection: ${err}`, 'warning');
            return null;
        }
    }

    private selectHelixFromGridNode(helixId: number): void {
        const helices = this.ensureScadnanoHelicesCache();
        if (!helices) return;

        const helix = helices[helixId];
        if (!Array.isArray(helix) || helix.length === 0) return;

        const selectElements = window.api?.selectElements;
        if (typeof selectElements !== 'function') return;

        selectElements(helix);
    }

    // Called when the grid editor's multi-select set changes (cmd/ctrl+click).
    // Selects all nucleotides belonging to every currently-selected helix in the 3D scene.
    private selectHelicesFromMultiSelect(ids: number[]): void {
        const helices = this.ensureScadnanoHelicesCache();
        if (!helices) return;

        const selectElements = window.api?.selectElements;
        if (typeof selectElements !== 'function') return;

        const allNucleos: any[] = [];
        ids.forEach(id => {
            const helix = helices[id];
            if (Array.isArray(helix)) allNucleos.push(...helix);
        });

        // Pass the full combined set in one call so the scene clears its old
        // selection and replaces it with exactly the highlighted helices.
        selectElements(allNucleos);
    }

    private ensureScadnanoGridEditor(gridType: ScadnanoGridType): any | null {
        if (this.scadnanoGridEditor && this.scadnanoGridEditorType === gridType) return this.scadnanoGridEditor;

        if (this.scadnanoGridEditor && this.scadnanoGridEditorType !== gridType) {
            if (typeof this.scadnanoGridEditor.dispose === 'function') {
                this.scadnanoGridEditor.dispose();
            }
            this.scadnanoGridEditor = null;
            this.scadnanoGridEditorType = null;
        }

        const canvas = this.getScadnanoGridCanvas();
        if (!canvas) return null;

        const scadnanoNs = window.scadnano;
        if (!scadnanoNs) return null;

        const editorCtorName = gridType === 'square' ? 'SquareEditor' : 'HoneycombEditor';
        const EditorCtor = scadnanoNs[editorCtorName];
        if (typeof EditorCtor !== 'function') return null;

        this.scadnanoGridEditor = new EditorCtor(canvas);
        this.scadnanoGridEditorType = gridType;

        this.scadnanoGridEditor.onNodesChanged = () => {
            this.publishCurrentHelixPosFromEditor();
        };

        this.scadnanoGridEditor.onNodeSelected = (node: any) => {
            if (this.suppressNodeSelectedCallback) return;
            const helixId = Number(node?.id);
            if (!Number.isFinite(helixId)) return;
            this.selectHelixFromGridNode(helixId);
        };

        // Multi-select (cmd/ctrl+click): sync all highlighted grid nodes to the 3D scene.
        this.scadnanoGridEditor.onSelectionChanged = (ids: number[]) => {
            if (this.suppressNodeSelectedCallback) return;
            this.selectHelicesFromMultiSelect(ids);
        };

        // Genuine user-initiated drags push a move entry onto the history journal.
        // Programmatic moves (undo/redo) suppress this callback inside the editor.
        this.scadnanoGridEditor.onNodeMoved = (info: { id: number; from: [number, number]; to: [number, number] }) => {
            const id = Number(info?.id);
            if (!Number.isFinite(id)) return;
            const fromCol = Number(info.from?.[0]);
            const fromRow = Number(info.from?.[1]);
            const toCol = Number(info.to?.[0]);
            const toRow = Number(info.to?.[1]);
            if (!Number.isFinite(fromCol) || !Number.isFinite(fromRow)) return;
            if (!Number.isFinite(toCol) || !Number.isFinite(toRow)) return;
            if (fromCol === toCol && fromRow === toRow) return;
            this.pushHistoryEntry({
                op: 'move',
                id,
                from: [fromCol, fromRow],
                to:   [toCol, toRow]
            });
            this.syncLayoutHelixPosFromEditor();
        };

        return this.scadnanoGridEditor;
    }

    private initScadnanoGridPaneControls(): void {
        const closeBtn = document.getElementById('scdgridClose');
        if (closeBtn) {
            closeBtn.addEventListener('click', () => {
                this.hideScadnanoGridPane();
            });
        }

        const exportBtn = document.getElementById('scadnanoGridExportBtn');
        if (exportBtn) {
            exportBtn.addEventListener('click', () => {
                this.publishCurrentHelixPosFromEditor();
                this.exportFromGridView(window.currentScadnanoHelixPos);
                // Drop the history cache after a successful export, per spec.
                this.clearHistory();
            });
        }

        const combineBtn = document.getElementById('scadnanoGridCombineBtn');
        if (combineBtn) {
            combineBtn.addEventListener('click', () => {
                this.combineSelectedHelicesFromGridView();
            });
        }

        const lockHelicesBtn = document.getElementById('scadnanoGridLockHelicesBtn');
        if (lockHelicesBtn) {
            lockHelicesBtn.addEventListener('click', () => {
                this.lockSelectedHelices();
            });
        }

        const recalcBtn = document.getElementById('scadnanoGridRecalcBtn');
        if (recalcBtn) {
            recalcBtn.addEventListener('click', () => {
                this.recalculateGridFromScratch();
            });
        }

        const undoBtn = document.getElementById('scadnanoGridUndoBtn');
        if (undoBtn) {
            undoBtn.addEventListener('click', () => {
                this.undoFromGridView();
            });
        }

        const redoBtn = document.getElementById('scadnanoGridRedoBtn');
        if (redoBtn) {
            redoBtn.addEventListener('click', () => {
                this.redoFromGridView();
            });
        }

        // Cmd/Ctrl+Z and Cmd/Ctrl+Shift+Z are scoped to the gridview canvas so they don't
        // conflict with the global editHistory bindings on the 3D scene canvas.
        const gridCanvas = this.getScadnanoGridCanvas();
        if (gridCanvas) {
            // Make the canvas focusable so it can receive keydown events.
            if (!gridCanvas.hasAttribute('tabindex')) {
                gridCanvas.setAttribute('tabindex', '0');
            }
            gridCanvas.addEventListener('keydown', (e: KeyboardEvent) => {
                const cmd = e.ctrlKey || e.metaKey;
                if (!cmd) return;
                if (e.key.toLowerCase() === 'z') {
                    e.preventDefault();
                    if (e.shiftKey) this.redoFromGridView();
                    else this.undoFromGridView();
                } else if (e.key.toLowerCase() === 'y') {
                    e.preventDefault();
                    this.redoFromGridView();
                }
            });
        }

        this.refreshHistoryButtons();

        const resizeHandle = document.getElementById('scadnanoGridResizeHandle');
        let resizing = false;

        if (resizeHandle) {
            resizeHandle.addEventListener('mousedown', (e: MouseEvent) => {
                resizing = true;
                document.body.classList.add('scadnano-grid-resizing');
                e.preventDefault();
            });
        }

        window.addEventListener('mousemove', (e: MouseEvent) => {
            if (!resizing) return;
            this.setScadnanoPaneWidth(e.clientX);
            this.resizeScadnanoGridCanvas();
            if (this.scadnanoGridEditor && typeof this.scadnanoGridEditor.resize === 'function') {
                this.scadnanoGridEditor.resize();
            }
        });

        window.addEventListener('mouseup', () => {
            if (!resizing) return;
            resizing = false;
            document.body.classList.remove('scadnano-grid-resizing');
        });

        window.addEventListener('resize', () => {
            this.resizeScadnanoGridCanvas();
            if (this.scadnanoGridEditor && typeof this.scadnanoGridEditor.resize === 'function') {
                this.scadnanoGridEditor.resize();
            }
        });
    }
}

const scadnanoManager = new ScadnanoExportManager();

function registerScadnanoWindowApi(): void {
    window.exportScadnanoFromGridView = (helixPosInput?: unknown) => {
        scadnanoManager.exportFromGridView(helixPosInput);
    };

    window.showScadnanoGridFromHelixPos = (helixPosInput?: unknown, gridTypeInput?: unknown) => {
        scadnanoManager.showGridFromHelixPos(helixPosInput, gridTypeInput);
    };

    window.hideScadnanoGridPane = () => {
        scadnanoManager.hideScadnanoGridPane();
    };

    window.scadnanoDialogExport = () => {
        scadnanoManager.handleDialogExport();
    };

    window.toggleGridDropdown = (checkboxElement: HTMLInputElement) => {
        scadnanoManager.toggleGridDropdown(checkboxElement);
    };

    window.scadnanoSelectHelixFromNucleotide = (nucleotideInput?: unknown, additive?: boolean) => {
        scadnanoManager.selectHelixFromNucleotide(nucleotideInput, additive === true);
    };

    window.scadnanoGetHelices = () => scadnanoManager.getHelices();

    window.scadnanoGridUndo = () => scadnanoManager.undoFromGridView();
    window.scadnanoGridRedo = () => scadnanoManager.redoFromGridView();
    window.scadnanoGridGetHistory = () => scadnanoManager.getHistorySnapshot();
}

registerScadnanoWindowApi();

// Keep these named wrappers for inline HTML handlers and backwards compatibility.
function scadnanoDialogExport(): void {
    scadnanoManager.handleDialogExport();
}

function toggleGridDropdown(checkboxElement: HTMLInputElement): void {
    scadnanoManager.toggleGridDropdown(checkboxElement);
}
