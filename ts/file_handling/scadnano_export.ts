/// <reference path="../typescript_definitions/index.d.ts" />

type ScadnanoGridType = 'honeycomb' | 'square';
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
    gridType: ScadnanoGridType;
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
    scadnanoSelectHelixFromNucleotide?: (nucleotideInput?: unknown) => void;
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

    private scadnanoGridEditor: any = null;
    private scadnanoGridEditorType: ScadnanoGridType | null = null;
    private suppressNodeSelectedCallback = false;

    // Edit-history journal. Lives only while the grid pane is open; cleared on export / close /
    // reload. Holds combine and move operations in the readable JSON form discussed.
    private history: ScadnanoHistoryJournal = this.createEmptyHistory();

    private createEmptyHistory(): ScadnanoHistoryJournal {
        return { entries: {}, order: [], cursor: 0, nextMove: 0, nextCombine: 0 };
    }

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

        // Verify every pair of selected helices is mutually disjoint on the grid. Combining
        // non-disjoint helices would collide nucleotides into the same offset/direction slot,
        // so bail out with an error and let the user resolve the overlap first.
        for (let i = 0; i < sortedIds.length; i++) {
            for (let j = i + 1; j < sortedIds.length; j++) {
                const h1 = sortedIds[i];
                const h2 = sortedIds[j];
                if (!toscad.disjoint(helices, h1, h2, this.currentScadnanoLayout.grid)) {
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
            }
        };
        this.pushHistoryEntry(entry);

        notify(`Combined helix ${mergedIdx.join(', ')} into helix ${keptIdx}.`);
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

        this.applyCombineCascadeInverse(entry, result.inverseRemap);
    }

    private remapCachedGridIds(remap: (oldId: number) => number | null): void {
        if (!this.currentScadnanoLayout) return;
        this.currentScadnanoLayout.grid.forEach(mark => {
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

        // Fresh grid view = fresh history.
        this.clearHistory();
    }

    public hideScadnanoGridPane(): void {
        document.body.classList.remove('scadnano-grid-open');
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

    public selectHelixFromNucleotide(nucleotideInput?: unknown): void {
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
            this.scadnanoGridEditor.selectNodeById(helixId);
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
        let failed = false;

        this.runScadnanoLongCalculation(
            () => {
                try {
                    helixPos = this.calculateScadnanoHelixPos(options.gridType, options.wireframe);
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
                this.showGridFromHelixPos(helixPos, options.gridType);
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
            gridType: this.normalizeGridType(scadnanoGrid.value),
            includeHelixPos: helixPosCheckbox.checked,
            wireframe: wireframeCheckbox.checked,
        };
    }

    private readCurrentExportTarget(): { name: string; gridType: ScadnanoGridType; wireframe: boolean } {
        const nameInput = document.getElementById('scadnanoFilename') as HTMLInputElement | null;
        const scadnanoGrid = document.getElementById('scadnanoGrid') as HTMLInputElement | null;
        const wireframeCheckbox = document.getElementById('scadnanoWireframe') as HTMLInputElement | null;

        return {
            name: nameInput?.value.trim() || 'output',
            gridType: this.normalizeGridType(scadnanoGrid?.value),
            wireframe: Boolean(wireframeCheckbox?.checked),
        };
    }

    private normalizeGridType(value?: string): ScadnanoGridType {
        return value === 'honeycomb' ? 'honeycomb' : 'square';
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
        gridType: ScadnanoGridType,
        helixPos?: HelixPosMap,
        wireframe = false
    ): void {
        const latticeType: ScadnanoGridType = this.normalizeGridType(gridType);
        const { helices, grid } = this.prepareScadnanoLayout(latticeType, false, wireframe);

        // Switch to toscad.buildScadnano2 here to fall back to the old
        // topology-driven export. buildScadnano3 reads boundaries from the
        // grid (helixId / direction / offset step) so post-construction
        // grid edits propagate to the export.
        const scadnano = helixPos
            ? toscad.buildScadnano3(grid, helices, gridType, helixPos)
            : toscad.buildScadnano3(grid, helices, gridType);

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

        const result = helix.findHelices(nucleotideElements, 3) as { helices: Nucleotide[][] };
        const helices = result?.helices ?? [];
        this.notifyHelixCoverageMismatch(helices, nucleotideElements);
        return helices;
    }

    private prepareScadnanoLayout(
        latticeType: ScadnanoGridType,
        forceRecompute = false,
        wireframe = false
    ): ScadnanoPreparedLayout {
        const nucleotideCount = this.getCurrentNucleotideCount();
        if (
            !forceRecompute &&
            this.currentScadnanoLayout &&
            this.currentScadnanoLayout.latticeType === latticeType &&
            this.currentScadnanoLayout.nucleotideCount === nucleotideCount &&
            this.currentScadnanoLayout.wireframe === wireframe
        ) {
            this.currentScadnanoHelices = this.currentScadnanoLayout.helices;
            return this.currentScadnanoLayout;
        }

        const helices = this.calculateScadnanoHelices();
        this.currentScadnanoHelices = helices;

        const { grid, binderHelices } = toscad.setGrid(helices);
        toscad.directionAlign2(grid);
        toscad.alignGridPrim(grid, binderHelices);

        const angles = toscad.getAngles(grid, helices, latticeType);
        let networkMap = angles;

        if (!wireframe) {
            const corrected = toscad.anglecomb(grid, helices, latticeType, angles);
            const correct = toscad.anglecorr(grid, helices, latticeType, corrected.networkMap);
            networkMap = correct.networkMap;
        }

        // Initial helix grid positions (still using pre-renumber helix IDs).
        let helixPos = toscad.calculateGlobalPositions(networkMap, undefined, undefined, latticeType);

        // Renumber helices so consecutive IDs share crossovers when possible.
        // Runs after combine (anglecomb), so it operates on the final helix
        // count. Mutates `grid` in place; produces a renumbered helices
        // array and a helixPos map keyed by new IDs.
        const renumber = toscad.renumberHelices(grid, helixPos, latticeType);
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

    private calculateScadnanoHelixPos(latticeType: ScadnanoGridType = 'square', wireframe = false): HelixPosMap {
        const { helixPos } = this.prepareScadnanoLayout(latticeType, false, wireframe);
        return this.cloneHelixPosMap(helixPos);
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

    window.scadnanoSelectHelixFromNucleotide = (nucleotideInput?: unknown) => {
        scadnanoManager.selectHelixFromNucleotide(nucleotideInput);
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
