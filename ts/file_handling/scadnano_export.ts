/// <reference path="../typescript_definitions/index.d.ts" />

// Scadnano export. All layout work belongs to toscad.layoutPipeline; this file
// only reads the dialog options, calls the pipeline, writes the .sc file, and
// hosts the grid pane so helix positions can be nudged before export.

type ScadnanoGridType = 'honeycomb' | 'square';
type ScadnanoRequestedGridType = ScadnanoGridType | 'automatic';
type HelixPosMap = Map<number, [number, number]>;

type ScadnanoLayout = ReturnType<typeof toscad.layoutPipeline>;

interface Window {
    view?: { longCalculation: (calc: () => void, message: string, callback?: () => void) => void };
    api?: { selectElements?: (elementsToSelect: any) => void };
    scadnano?: any;

    currentScadnanoHelixPos?: HelixPosMap;
    exportScadnanoFromGridView?: (helixPosInput?: unknown) => void;
    showScadnanoGridFromHelixPos?: (helixPosInput?: unknown, gridTypeInput?: unknown) => void;
    hideScadnanoGridPane?: () => void;
    scadnanoDialogExport?: () => void;
    toggleGridDropdown?: (checkboxElement: HTMLInputElement) => void;
    scadnanoSelectHelixFromNucleotide?: (nucleotideInput?: unknown, additive?: boolean) => void;
    scadnanoGetHelices?: () => Nucleotide[][] | null;
    scadnanoCombineSelectedHelices?: () => void;
    scadnanoSplitSelectedHelix?: () => void;
    scadnanoFocusOnHelixToggle?: (checkboxElement: HTMLInputElement) => void;
    scadnanoGridUndo?: () => void;
    scadnanoGridRedo?: () => void;
}

namespace scadnanoExport {
    const TOLERANCE = 3;
    const MAX_ITERATIONS = 10;

    let layout: ScadnanoLayout | null = null;
    let connections: Array<[number, number]> = [];
    let gridEditor: any = null;
    let gridEditorType: ScadnanoGridType | null = null;
    let suppressSelectionCallback = false;

    let lockedHelices = new Set<number>();
    const LOCKED_COLOR = 0x808080;
    const UNLOCKED_COLOR = 0xffd400;

    // ── pipeline + export ────────────────────────────────────────────────────

    function nucleotideMap(): Map<number, Nucleotide> {
        const out = new Map<number, Nucleotide>();
        elements.forEach((e, id) => { if (e instanceof Nucleotide) out.set(id, e); });
        return out;
    }

    // The one call sign. Everything else in this file consumes its output.
    function runPipeline(
        lattice: ScadnanoRequestedGridType,
        wireframe: boolean,
        pins: toscad.RelativePin[] = []
    ): ScadnanoLayout {
        const nucleotides = nucleotideMap();
        layout = toscad.layoutPipeline(nucleotides, {
            tolerance: TOLERANCE,
            lattice,
            maxIterations: MAX_ITERATIONS,
            wireframe,
            pins
        });

        const mapped = layout.helices.flat().length;
        if (mapped !== nucleotides.size) {
            notify(`Helix mapping error: ${mapped}/${nucleotides.size} nucleotides were mapped.`, 'alert', true);
        }

        connections = helixPairs(toscad.collectCrossovers(layout.grid).crossovers);
        return layout;
    }

    // Takes layoutPipeline output and writes the .sc file. helixPos overrides the
    // pipeline's positions (used after the grid pane has been edited).
    export function writeScadnanoFile(name: string, source: ScadnanoLayout, helixPos?: HelixPosMap): void {
        const { grid, helices, latticeType } = source;
        const sc = helixPos
            ? toscad.buildScadnano3(grid, helices, latticeType, helixPos)
            : toscad.buildScadnano3(grid, helices, latticeType);
        makeTextFile(name ? `${name}.sc` : 'output.sc', JSON.stringify(sc, null, 2));
    }

    // Undirected helix pairs for the grid view's connection lines.
    function helixPairs(
        crossovers: Map<number, Map<number, { sameWalk: number; diffWalk: number }>>
    ): Array<[number, number]> {
        const seen = new Set<string>();
        const pairs: Array<[number, number]> = [];
        crossovers.forEach((toMap, from) => {
            toMap.forEach((counts, to) => {
                if (from === to) return;
                if (Number(counts?.sameWalk ?? 0) + Number(counts?.diffWalk ?? 0) <= 0) return;
                const a = Math.min(from, to), b = Math.max(from, to);
                if (seen.has(`${a}:${b}`)) return;
                seen.add(`${a}:${b}`);
                pairs.push([a, b]);
            });
        });
        return pairs;
    }

    // ── dialog ───────────────────────────────────────────────────────────────

    function input(id: string): HTMLInputElement | null {
        return document.getElementById(id) as HTMLInputElement | null;
    }

    function requestedGridType(value?: string): ScadnanoRequestedGridType {
        return value === 'square' ? 'square' : value === 'honeycomb' ? 'honeycomb' : 'automatic';
    }

    function dialogOptions() {
        return {
            name: input('scadnanoFilename')?.value.trim() || 'output',
            lattice: requestedGridType(input('scadnanoGrid')?.value),
            includeHelixPos: Boolean(input('scadnanoIncludeHPos')?.checked),
            wireframe: Boolean(input('scadnanoWireframe')?.checked)
        };
    }

    function longCalculation(calc: () => void, done?: () => void): void {
        const runner = window.view?.longCalculation;
        if (typeof runner === 'function') {
            runner(calc, 'Preparing scadnano export, please be patient...', done);
            return;
        }
        calc();
        done?.();
    }

    function closeDialog(): void {
        try {
            (window as any)?.Metro?.dialog?.close('#scadnanoDialog');
        } catch (err) {
            console.warn('Failed to close scadnano dialog via Metro API:', err);
        }
        const dialogEl = document.getElementById('scadnanoDialog');
        if (dialogEl) {
            dialogEl.classList.remove('open');
            dialogEl.setAttribute('aria-hidden', 'true');
        }
    }

    // Ribbon dialog "Export". Without helix positions it writes the file directly;
    // with them it opens the grid pane on the pipeline's layout instead.
    export function handleDialogExport(): void {
        const { name, lattice, includeHelixPos, wireframe } = dialogOptions();
        closeDialog();

        requestAnimationFrame(() => requestAnimationFrame(() => {
            let result: ScadnanoLayout | null = null;
            longCalculation(
                () => {
                    try {
                        result = runPipeline(lattice, wireframe);
                        if (!includeHelixPos) writeScadnanoFile(name, result);
                    } catch (err) {
                        result = null;
                        notify(`Scadnano export failed: ${err}`, 'alert');
                    }
                },
                () => {
                    if (!includeHelixPos || !result) return;
                    window.currentScadnanoHelixPos = new Map(result.helixPos);
                    clearHistory(); // fresh pipeline run: start a new, empty history
                    showGrid(result.helixPos, result.latticeType);
                }
            );
        }));
    }

    // Grid pane "Export as scadnano": same file, with the edited positions.
    export function exportFromGridView(helixPosInput?: unknown): void {
        const helixPos = toHelixPosMap(helixPosInput ?? window.currentScadnanoHelixPos);
        if (!helixPos || helixPos.size === 0) {
            notify('No edited helix positions available to export.', 'warning');
            return;
        }

        const { name, lattice, wireframe } = dialogOptions();
        try {
            writeScadnanoFile(name, layout ?? runPipeline(lattice, wireframe), helixPos);
        } catch (err) {
            notify(`Scadnano export failed: ${err}`, 'alert');
        }
    }

    export function toggleGridDropdown(checkboxElement: HTMLInputElement): void {
        const dropdown = document.getElementById('scadnanoGrid') as HTMLSelectElement | null;
        if (dropdown) dropdown.disabled = !checkboxElement.checked;
    }

    // ── grid pane ────────────────────────────────────────────────────────────

    function toHelixPosMap(value: unknown): HelixPosMap | null {
        if (!value) return null;
        const entries: any[] = value instanceof Map ? [...value.entries()]
            : Array.isArray(value) ? value
            : typeof value === 'object' ? Object.entries(value) : [];

        const out: HelixPosMap = new Map();
        entries.forEach(entry => {
            const id = Number(entry?.[0]), pos = entry?.[1];
            if (!Array.isArray(pos) || pos.length < 2) return;
            const col = Number(pos[0]), row = Number(pos[1]);
            if (Number.isFinite(id) && Number.isFinite(col) && Number.isFinite(row)) out.set(id, [col, row]);
        });
        return out;
    }

    function canvas(): HTMLCanvasElement | null {
        return document.getElementById('scadnanoGridCanvas') as HTMLCanvasElement | null;
    }

    function resizeCanvas(): void {
        const pane = document.getElementById('scadnanoGridPane');
        const c = canvas();
        if (!pane || !c) return;
        c.width = pane.clientWidth;
        c.height = pane.clientHeight;
        if (typeof gridEditor?.resize === 'function') gridEditor.resize();
    }

    function publishHelixPos(): void {
        if (typeof gridEditor?.getNodes !== 'function') return;
        const out: HelixPosMap = new Map();
        (gridEditor.getNodes() as Array<{ id: number; col: number; row: number }>).forEach(n => {
            out.set(Number(n.id), [Number(n.col), Number(n.row)]);
        });
        window.currentScadnanoHelixPos = out;
        if (layout) layout.helixPos = out;
    }

    // Grid node ids are helix indices, so they index layout.helices directly.
    function selectHelicesInScene(ids: number[]): void {
        const helices = getHelices();
        const select = window.api?.selectElements;
        if (!helices || typeof select !== 'function') return;
        const nts: any[] = [];
        ids.forEach(id => { if (Array.isArray(helices[id])) nts.push(...helices[id]); });
        if (nts.length) select(nts);
    }

    function ensureEditor(gridType: ScadnanoGridType): any | null {
        if (gridEditor && gridEditorType === gridType) return gridEditor;
        if (gridEditor && typeof gridEditor.dispose === 'function') gridEditor.dispose();
        gridEditor = null;
        gridEditorType = null;

        const c = canvas();
        const Ctor = window.scadnano?.[gridType === 'square' ? 'SquareEditor' : 'HoneycombEditor'];
        if (!c || typeof Ctor !== 'function') return null;

        gridEditor = new Ctor(c);
        gridEditorType = gridType;
        gridEditor.onNodesChanged = () => publishHelixPos();
        // Genuine user drags only (suppressed for programmatic moves). Record one undo entry per
        // drag: snapshot the post-drag state, then rewind the moved node to its pre-drag cell so the
        // stored snapshot is the state as it was BEFORE the drag.
        gridEditor.onNodeMoved = (info: { id: number; from: [number, number]; to: [number, number] }) => {
            const before = snapshot();
            if (!before) return;
            const entry = before.helixPos.find(([id]) => id === Number(info.id));
            if (entry) entry[1] = [info.from[0], info.from[1]];
            pushHistory(before);
        };
        gridEditor.onNodeSelected = (node: any) => {
            if (suppressSelectionCallback) return;
            const id = Number(node?.id);
            if (Number.isFinite(id)) selectHelicesInScene([id]);
        };
        gridEditor.onSelectionChanged = (ids: number[]) => {
            if (!suppressSelectionCallback) selectHelicesInScene(ids);
        };
        return gridEditor;
    }

    export function showGrid(helixPosInput?: unknown, gridTypeInput?: unknown): void {
        if (!document.getElementById('scadnanoGridPane')) {
            notify('Scadnano grid pane is unavailable.', 'alert');
            return;
        }

        document.body.classList.add('scadnano-grid-open');
        resizeCanvas();

        const editor = ensureEditor(gridTypeInput === 'square' ? 'square' : 'honeycomb');
        if (!editor) {
            notify('Unable to open scadnano grid view.', 'alert');
            return;
        }

        const helixPos = toHelixPosMap(helixPosInput);
        if (!helixPos || helixPos.size === 0) {
            notify('No helix positions available for the grid view.', 'warning');
            return;
        }

        editor.loadFromHelixPos(helixPos);
        if (typeof editor.setConnections === 'function') editor.setConnections(connections);
        publishHelixPos();

        // loadFromHelixPos replaces every node, so lock colours and the editor's guard set
        // have to be reapplied. Deferred because loadFromHelixPos fits the view on a timeout.
        applyLockedColors();
    }

    export function hideGrid(): void {
        document.body.classList.remove('scadnano-grid-open');
        if (focusToggleEl()?.checked) focusOnHelixToggle(false); // restore view if closed while focused
        clearHistory(); // history lives only for the duration of the pane session
        lockedHelices.clear();
        if (gridEditor) gridEditor.lockedHelices = new Set();
    }

    // ── lock / unlock ────────────────────────────────────────────────────────

    // color the nodes based on their lock status.
    function applyLockedColors(): void {
        if (typeof gridEditor?.setNodeColor !== 'function') return;
        const nodes: Array<{ id: number }> = typeof gridEditor.getNodes === 'function' ? gridEditor.getNodes() : [];
        nodes.forEach(n => {
            const id = Number(n.id);
            gridEditor.setNodeColor(id, lockedHelices.has(id) ? LOCKED_COLOR : UNLOCKED_COLOR);
        });
        gridEditor.lockedHelices = new Set(lockedHelices);
    }

    // "Lock/Unlock Helices" toolbar button. Toggles every selected helix.
    export function lockSelectedHelices(): void {
        if (typeof gridEditor?.setNodeColor !== 'function') return;
        const ids: number[] = typeof gridEditor.getSelectedHelixIds === 'function'
            ? gridEditor.getSelectedHelixIds()
            : [];
        if (ids.length === 0) {
            notify('Select one or more helices in the grid first.', 'warning');
            return;
        }

        ids.forEach(id => {
            if (lockedHelices.has(id)) lockedHelices.delete(id);
            else lockedHelices.add(id);
        });
        applyLockedColors();
        notify(`${lockedHelices.size} helix/helices locked.`, 'success');
    }

    // ── history (undo / redo) ──────────────────────────────────────────────────

    // A snapshot is the whole grid-view state in cheap, cloneable form: helices as nucleotide-id
    // arrays (rebuilt via elements.get on restore), the per-nucleotide grid marks, helix positions,
    // the lock set, and the lattice. connections are derived from the grid, so they're recomputed on
    // restore rather than stored. One linear undo stack + redo stack, capped at HISTORY_CAP entries.
    interface GridSnapshot {
        helices: number[][];
        grid: Array<[number, { helixId: number; offset: number; direction: toscad.GridMark['direction'] }]>;
        helixPos: Array<[number, [number, number]]>;
        locked: number[];
        latticeType: ScadnanoGridType;
    }

    const HISTORY_CAP = 40;
    let undoStack: GridSnapshot[] = [];
    let redoStack: GridSnapshot[] = [];

    // Clone the current layout + lock state. Returns null when there's no layout to capture.
    function snapshot(): GridSnapshot | null {
        if (!layout) return null;
        return {
            helices: layout.helices.map(h => h.map(nt => nt.id)),
            grid: [...layout.grid.entries()].map(([id, m]) =>
                [id, { helixId: m.helixId, offset: m.offset, direction: m.direction }] as GridSnapshot['grid'][number]),
            helixPos: [...layout.helixPos.entries()].map(([id, p]) =>
                [id, [p[0], p[1]]] as GridSnapshot['helixPos'][number]),
            locked: [...lockedHelices],
            latticeType: layout.latticeType as ScadnanoGridType
        };
    }

    // Swap a snapshot back into the live layout and redraw. Nucleotide objects are looked up fresh
    // from the global elements map (they never leave it), so only their grouping is restored.
    function restore(snap: GridSnapshot): void {
        if (!layout) return;

        layout.helices = snap.helices.map(ids => {
            const arr: Nucleotide[] = [];
            ids.forEach(id => { const e = elements.get(id); if (e instanceof Nucleotide) arr.push(e); });
            return arr;
        });

        const grid: toscad.GridMap = new Map();
        snap.grid.forEach(([id, m]) => grid.set(id, { helixId: m.helixId, offset: m.offset, direction: m.direction }));
        layout.grid = grid;

        const helixPos: HelixPosMap = new Map();
        snap.helixPos.forEach(([id, p]) => helixPos.set(id, [p[0], p[1]]));
        layout.helixPos = helixPos;

        lockedHelices = new Set(snap.locked);
        (layout as any).latticeType = snap.latticeType;

        connections = helixPairs(toscad.collectCrossovers(layout.grid).crossovers);
        window.currentScadnanoHelixPos = new Map(layout.helixPos);
        showGrid(layout.helixPos, layout.latticeType);
    }

    function updateHistoryButtons(): void {
        const undoBtn = document.getElementById('scadnanoGridUndoBtn') as HTMLButtonElement | null;
        const redoBtn = document.getElementById('scadnanoGridRedoBtn') as HTMLButtonElement | null;
        if (undoBtn) undoBtn.disabled = undoStack.length === 0;
        if (redoBtn) redoBtn.disabled = redoStack.length === 0;
    }

    // Push a pre-operation snapshot. Any fresh operation invalidates the redo branch. Pass the state
    // captured BEFORE the mutation was applied.
    function pushHistory(snap: GridSnapshot | null): void {
        if (!snap) return;
        undoStack.push(snap);
        if (undoStack.length > HISTORY_CAP) undoStack.shift();
        redoStack = [];
        updateHistoryButtons();
    }

    function clearHistory(): void {
        undoStack = [];
        redoStack = [];
        updateHistoryButtons();
    }

    export function undo(): void {
        if (undoStack.length === 0) { notify('Nothing to undo.', 'warning'); return; }
        const current = snapshot();
        if (current) redoStack.push(current);
        restore(undoStack.pop()!);
        updateHistoryButtons();
    }

    export function redo(): void {
        if (redoStack.length === 0) { notify('Nothing to redo.', 'warning'); return; }
        const current = snapshot();
        if (current) {
            undoStack.push(current);
            if (undoStack.length > HISTORY_CAP) undoStack.shift();
        }
        restore(redoStack.pop()!);
        updateHistoryButtons();
    }

    // ── combine ────────────────────────────────────────────────────────────────

    // Remaps a helix-position map through a combineHelices result: merged-away helices drop out
    // (their nucleotides now live in the kept helix, which keeps its own cell) and every survivor
    // shifts down to its new index.
    function remapHelixPos(source: HelixPosMap, mergedIdx: number[], idRemap: Map<number, number>): HelixPosMap {
        const removed = new Set(mergedIdx);
        const out: HelixPosMap = new Map();
        source.forEach((pos, oldId) => {
            if (removed.has(oldId)) return;
            const next = idRemap.get(oldId);
            out.set(next !== undefined ? next : oldId, pos);
        });
        return out;
    }

    // "Combine Helices" toolbar button. Folds every selected helix into the lowest-numbered one via
    // helix.combineHelices, then rebuilds the grid view on the remapped layout. combineHelices takes
    // the whole index list and collapses it into a single helix in one pass, so there's nothing to
    // recurse over here — one call merges all of the selection.
    export function combineSelectedHelices(): void {
        const ids: number[] = typeof gridEditor?.getSelectedHelixIds === 'function'
            ? gridEditor.getSelectedHelixIds()
            : [];
        if (ids.length < 2) {
            notify('Select two or more helices in the grid to combine.', 'warning');
            return;
        }

        // Locked helices are pinned in place, so refuse to fold them into another helix.
        const lockedInSelection = ids.filter(id => lockedHelices.has(id));
        if (lockedInSelection.length > 0) {
            const s = lockedInSelection.length === 1;
            notify(`Cannot combine: helix ${lockedInSelection.join(', ')} ${s ? 'is' : 'are'} locked. Unlock ${s ? 'it' : 'them'} first.`, 'warning');
            return;
        }

        // Capture any pending drags first, so kept helices keep the positions the user is looking at.
        publishHelixPos();
        const source = window.currentScadnanoHelixPos ?? layout.helixPos;

        // Snapshot the pre-combine state for undo (pushed only once the merge actually happens).
        const before = snapshot();

        // combineHelices mutates layout.helices AND layout.grid in place, returning the id remap.
        const result = helix.combineHelices(layout.helices, ids, layout.grid, layout.latticeType);
        if (!result) {
            notify('Nothing to combine.', 'warning');
            return;
        }
        pushHistory(before);

        const { keptIdx, mergedIdx, idRemap } = result;

        // Positions: drop merged-away helices, shift survivors down through the remap.
        const newHelixPos = remapHelixPos(source, mergedIdx, idRemap);
        layout.helixPos = newHelixPos;
        window.currentScadnanoHelixPos = new Map(newHelixPos);

        // Connections come straight off the freshly-remapped grid.
        connections = helixPairs(toscad.collectCrossovers(layout.grid).crossovers);

        // Locks follow their helices through the same remap (survivor indices shifted too).
        const removed = new Set(mergedIdx);
        const newLocked = new Set<number>();
        lockedHelices.forEach(id => {
            if (removed.has(id)) { newLocked.add(keptIdx); return; }
            const next = idRemap.get(id);
            if (next !== undefined) newLocked.add(next);
        });
        lockedHelices = newLocked;

        // Rebuild the grid pane on the remapped layout (also reapplies connections + lock colours).
        showGrid(newHelixPos, layout.latticeType);
        notify(`Combined helix ${mergedIdx.join(', ')} into helix ${keptIdx}.`, 'success');
    }

    // ── split ────────────────────────────────────────────────────────────────

    // Which helix the currently selected nucleotides belong to, per the grid marks (the same
    // membership toscad.splitHelix trusts). Returns null when nothing usable is selected, or the
    // selection straddles more than one helix — split only makes sense within a single helix.
    function selectedHelixForSplit(): { helixId: number; nucleotides: Nucleotide[] } | null {
        const nucleotides: Nucleotide[] = [];
        const helixIds = new Set<number>();
        selectedBases.forEach(e => {
            if (!(e instanceof Nucleotide)) return;
            const mark = layout!.grid.get(e.id);
            if (!mark) return;
            nucleotides.push(e);
            helixIds.add(mark.helixId);
        });

        if (nucleotides.length === 0) return null;
        if (helixIds.size > 1) return null;
        return { helixId: [...helixIds][0], nucleotides };
    }

    // "Split from selected" toolbar button. Moves the selected nucleotides out of their helix into
    // a new one via toscad.splitHelix. splitHelix appends the new helix (newHelixId = helices.length)
    // and touches only the moved nts' grid marks, so no existing index shifts — the cascade is just
    // recompute-connections + redraw, with no remap of helixPos or locks.
    export function splitSelectedHelix(): void {
        if (!layout) {
            notify('Open the grid view first before splitting a helix.', 'warning');
            return;
        }

        const selection = selectedHelixForSplit();
        if (!selection) {
            notify('Select nucleotides within a single helix to split off.', 'warning');
            return;
        }

        // Snapshot the pre-split state for undo (pushed only once the split actually happens).
        const before = snapshot();

        // splitHelix mutates layout.helices, layout.grid and layout.helixPos in place, and does its
        // own guarding (nothing to keep, none belong to the helix, etc.), returning null on refusal.
        const result = toscad.splitHelix(
            layout.grid,
            layout.helixPos,
            selection.helixId,
            layout.helices,
            selection.nucleotides
        );
        if (!result) {
            notify(`Could not split helix ${selection.helixId} (see console for details).`, 'warning');
            return;
        }
        pushHistory(before);

        // New helix appended at the end, existing indices unchanged: locks stay valid as-is.
        connections = helixPairs(toscad.collectCrossovers(layout.grid).crossovers);
        window.currentScadnanoHelixPos = new Map(layout.helixPos);

        showGrid(layout.helixPos, layout.latticeType);
        const movedCount = layout.helices[result.newHelixId]?.length ?? 0;
        notify(`Split ${movedCount} nucleotide(s) off helix ${result.keptHelixId} into new helix ${result.newHelixId}.`, 'success');
    }

    // ── focus on helix ─────────────────────────────────────────────────────────

    const FOCUS_TOGGLE_ID = 'scadnanoGridSplitHelicesToggle';

    function focusToggleEl(): HTMLInputElement | null {
        return document.getElementById(FOCUS_TOGGLE_ID) as HTMLInputElement | null;
    }

    // Splitting only makes sense while focused on a helix, so the button tracks the toggle.
    function setSplitEnabled(on: boolean): void {
        const btn = document.getElementById('scadnanoGridSplitFromSelectedBtn') as HTMLButtonElement | null;
        if (btn) btn.disabled = !on;
    }

    // "Focus on helix" switch. On: invert the selection and reuse the existing "Hide Selection"
    // command (toggleVisArbitrary) to hide everything except the selected helix. Off: unhide all.
    export function focusOnHelixToggle(checked: boolean): void {
        if (checked) {
            if (selectedBases.size === 0) {
                notify('Select a helix first to focus on it.', 'warning');
                const chk = focusToggleEl(); if (chk) chk.checked = false;
                return;
            }
            invertSelection();      // select everything outside the helix
            toggleVisArbitrary();   // hide it (also clears the selection)
        } else {
            const chk = focusToggleEl(); if (chk) chk.checked = false;
            clearSelection();       // force the "unhide all" branch
            toggleVisArbitrary();   // restore the full view
        }
        setSplitEnabled(checked);
    }

    // ── recalculate ──────────────────────────────────────────────────────────

    // Locked helices become hard relative pins. Only differences are meaningful, so the
    // lowest locked id is the reference and every other locked helix contributes one pin
    // carrying its current offset from it. Fewer than two locks means nothing to constrain.
    function pinsFromLocked(): toscad.RelativePin[] {
        if (lockedHelices.size < 2) return [];
        const pos = window.currentScadnanoHelixPos;
        if (!pos) return [];

        const ids = [...lockedHelices].filter(id => pos.has(id)).sort((x, y) => x - y);
        if (ids.length < 2) return [];

        const [refId, ...rest] = ids;
        const ref = pos.get(refId)!;
        return rest.map(id => {
            const p = pos.get(id)!;
            return { a: refId, b: id, dCol: p[0] - ref[0], dRow: p[1] - ref[1] };
        });
    }

    // Absolute (col,row) cells for every locked helix. Pins (above) only carry
    // relative offsets; these hand kruskals the user's exact placement so the
    // locked component is anchored on its true cell/parity instead of drifting
    // to wherever the 2-coloring seed lands it.
    function anchorsFromLocked(): Map<number, [number, number]> {
        const anchors = new Map<number, [number, number]>();
        if (lockedHelices.size < 2) return anchors;
        const pos = window.currentScadnanoHelixPos;
        if (!pos) return anchors;
        for (const id of lockedHelices) {
            const p = pos.get(id);
            if (p) anchors.set(id, [p[0], p[1]]);
        }
        return anchors;
    }

    // "Recalculate Grid" toolbar button. Re-runs ONLY the positioning stage (angles → kruskals →
    // posCorr5) on the CURRENT layout.helices / layout.grid — never findHelices — so any manual
    // combines or splits the user has made are preserved. Locked helices feed in as pins to hold
    // their relative placement. Helix ids are untouched, so locks stay valid without re-finding.
    export function recalculateGrid(): void {
        if (!layout) {
            notify('Open the grid view first before recalculating.', 'warning');
            return;
        }

        publishHelixPos();
        const before = snapshot();
        const pins = pinsFromLocked();
        const anchors = anchorsFromLocked();
        const latticeType = layout.latticeType as ScadnanoGridType;
        let ok = false;

        longCalculation(
            () => {
                try {
                    const networkMap = toscad.getAngles(layout!.grid, layout!.helices, latticeType);
                    const votes = toscad.voteOrientations(networkMap, layout!.grid, latticeType);
                    const kr = toscad.kruskals(networkMap, layout!.grid, latticeType, votes, pins, anchors);
                    const pinnedIds = new Set<number>();
                    pins.forEach(p => { pinnedIds.add(p.a); pinnedIds.add(p.b); });
                    layout!.helixPos = toscad.posCorr5(kr, layout!.helices, pinnedIds);
                    connections = helixPairs(toscad.collectCrossovers(layout!.grid).crossovers);
                    ok = true;
                } catch (err) {
                    ok = false;
                    notify(`Recalculate grid failed: ${err}`, 'alert');
                }
            },
            () => {
                if (!ok) return;
                pushHistory(before);
                window.currentScadnanoHelixPos = new Map(layout!.helixPos);
                showGrid(layout!.helixPos, latticeType);
                applyLockedColors();
                notify(`Grid recalculated${pins.length ? ` with ${pins.length} pin(s)` : ''}.`, 'success');
            }
        );
    }

    // Used by "Helix" selection mode in base_selector.ts as well as the grid pane.
    // Falls back to plain helix detection when no layout has been computed yet.
    export function getHelices(): Nucleotide[][] | null {
        if (layout?.helices.length) return layout.helices;
        try {
            return (helix.findHelices(nucleotideMap(), TOLERANCE) as { helices: Nucleotide[][] }).helices ?? null;
        } catch (err) {
            notify(`Unable to map grid helix selection: ${err}`, 'warning');
            return null;
        }
    }

    export function selectHelixFromNucleotide(nucleotideInput?: unknown, additive = false): void {
        if (!document.body.classList.contains('scadnano-grid-open')) return;
        if (typeof gridEditor?.selectNodeById !== 'function') return;

        const nt = nucleotideInput instanceof Nucleotide ? nucleotideInput
            : typeof nucleotideInput === 'number' ? elements.get(nucleotideInput) : null;
        if (!(nt instanceof Nucleotide)) return;

        const helices = getHelices();
        if (!helices) return;
        const helixId = api.helix.findHelixID(nt.id, helices);
        if (helixId === null) return;

        suppressSelectionCallback = true;
        try {
            gridEditor.selectNodeById(helixId, additive);
        } finally {
            suppressSelectionCallback = false;
        }
    }

    // ── wiring ───────────────────────────────────────────────────────────────

    function initPaneControls(): void {
        // Toolbar tabs live in scadnano_gridview. Guard the namespace: this runs at
        // script-load time and a ReferenceError here would break the inline handlers.
        if (typeof scadnano !== 'undefined' && typeof (scadnano as any).initGridToolbar === 'function') {
            (scadnano as any).initGridToolbar();
        }

        document.getElementById('scdgridClose')?.addEventListener('click', hideGrid);
        document.getElementById('scadnanoGridExportBtn')?.addEventListener('click', () => {
            publishHelixPos();
            exportFromGridView(window.currentScadnanoHelixPos);
        });
        document.getElementById('scadnanoGridLockHelicesBtn')?.addEventListener('click', lockSelectedHelices);
        document.getElementById('scadnanoGridCombineBtn')?.addEventListener('click', combineSelectedHelices);
        document.getElementById('scadnanoGridSplitFromSelectedBtn')?.addEventListener('click', splitSelectedHelix);
        document.getElementById('scadnanoGridRecalcBtn')?.addEventListener('click', recalculateGrid);
        document.getElementById('scadnanoGridUndoBtn')?.addEventListener('click', undo);
        document.getElementById('scadnanoGridRedoBtn')?.addEventListener('click', redo);

        // While the grid pane is open, ctrl/cmd+z / ctrl/cmd+shift+z (and ctrl/cmd+y) drive the grid
        // history instead of the molecular editHistory. Captured on window so it runs before the 3D
        // canvas keydown handler, and stopped from propagating so undo doesn't fire in both places.
        window.addEventListener('keydown', (e: KeyboardEvent) => {
            if (!document.body.classList.contains('scadnano-grid-open')) return;
            if (!(e.ctrlKey || e.metaKey)) return;
            const key = e.key.toLowerCase();
            if (key !== 'z' && key !== 'y') return;
            e.preventDefault();
            e.stopImmediatePropagation();
            if (key === 'y' || e.shiftKey) redo();
            else undo();
        }, true);

        let resizing = false;
        document.getElementById('scadnanoGridResizeHandle')?.addEventListener('mousedown', (e: MouseEvent) => {
            resizing = true;
            document.body.classList.add('scadnano-grid-resizing');
            e.preventDefault();
        });
        window.addEventListener('mousemove', (e: MouseEvent) => {
            if (!resizing) return;
            const max = Math.max(240, Math.floor(window.innerWidth * 0.75));
            const width = Math.max(240, Math.min(max, Math.round(e.clientX)));
            document.documentElement.style.setProperty('--scadnano-pane-width', `${width}px`);
            resizeCanvas();
        });
        window.addEventListener('mouseup', () => {
            if (!resizing) return;
            resizing = false;
            document.body.classList.remove('scadnano-grid-resizing');
        });
        window.addEventListener('resize', resizeCanvas);
    }

    initPaneControls();

    window.scadnanoDialogExport = handleDialogExport;
    window.toggleGridDropdown = toggleGridDropdown;
    window.exportScadnanoFromGridView = exportFromGridView;
    window.showScadnanoGridFromHelixPos = showGrid;
    window.hideScadnanoGridPane = hideGrid;
    window.scadnanoGetHelices = getHelices;
    window.scadnanoCombineSelectedHelices = combineSelectedHelices;
    window.scadnanoSplitSelectedHelix = splitSelectedHelix;
    window.scadnanoGridUndo = undo;
    window.scadnanoGridRedo = redo;
    window.scadnanoFocusOnHelixToggle = (chk?: unknown) =>
        focusOnHelixToggle(Boolean((chk as HTMLInputElement)?.checked));
    window.scadnanoSelectHelixFromNucleotide = (nt?: unknown, additive?: boolean) =>
        selectHelixFromNucleotide(nt, additive === true);
}

// Named wrappers for the inline handlers in index.html.
function scadnanoDialogExport(): void {
    scadnanoExport.handleDialogExport();
}

function toggleGridDropdown(checkboxElement: HTMLInputElement): void {
    scadnanoExport.toggleGridDropdown(checkboxElement);
}

function scadnanoFocusOnHelixToggle(checkboxElement: HTMLInputElement): void {
    scadnanoExport.focusOnHelixToggle(Boolean(checkboxElement?.checked));
}
