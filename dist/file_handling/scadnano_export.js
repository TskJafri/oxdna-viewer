"use strict";
/// <reference path="../typescript_definitions/index.d.ts" />
var scadnanoExport;
(function (scadnanoExport) {
    const TOLERANCE = 3;
    let layout = null;
    let connections = [];
    let gridEditor = null;
    let gridEditorType = null;
    let suppressSelectionCallback = false;
    let lockedHelices = new Set();
    const LOCKED_COLOR = 0x808080;
    const UNLOCKED_COLOR = 0xffd400;
    // ── pipeline + export ────────────────────────────────────────────────────
    function nucleotideMap() {
        const out = new Map();
        elements.forEach((e, id) => { if (e instanceof Nucleotide)
            out.set(id, e); });
        return out;
    }
    // The one call sign. Everything else in this file consumes its output.
    function runPipeline(lattice, wireframe, pins = []) {
        const nucleotides = nucleotideMap();
        layout = toscad.layoutPipeline(nucleotides, {
            tolerance: TOLERANCE,
            lattice,
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
    function writeScadnanoFile(name, source, helixPos) {
        const { grid, helices, latticeType } = source;
        const sc = helixPos
            ? toscad.buildScadnano3(grid, helices, latticeType, helixPos)
            : toscad.buildScadnano3(grid, helices, latticeType);
        makeTextFile(name ? `${name}.sc` : 'output.sc', JSON.stringify(sc, null, 2));
    }
    scadnanoExport.writeScadnanoFile = writeScadnanoFile;
    // Undirected helix pairs for the grid view's connection lines.
    function helixPairs(crossovers) {
        const seen = new Set();
        const pairs = [];
        crossovers.forEach((toMap, from) => {
            toMap.forEach((counts, to) => {
                if (from === to)
                    return;
                if (Number(counts?.sameWalk ?? 0) + Number(counts?.diffWalk ?? 0) <= 0)
                    return;
                const a = Math.min(from, to), b = Math.max(from, to);
                if (seen.has(`${a}:${b}`))
                    return;
                seen.add(`${a}:${b}`);
                pairs.push([a, b]);
            });
        });
        return pairs;
    }
    // ── dialog ───────────────────────────────────────────────────────────────
    function input(id) {
        return document.getElementById(id);
    }
    function requestedGridType(value) {
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
    function longCalculation(calc, done) {
        const runner = window.view?.longCalculation;
        if (typeof runner === 'function') {
            runner(calc, 'Preparing scadnano export, please be patient...', done);
            return;
        }
        calc();
        done?.();
    }
    function closeDialog() {
        try {
            window?.Metro?.dialog?.close('#scadnanoDialog');
        }
        catch (err) {
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
    function handleDialogExport() {
        const { name, lattice, includeHelixPos, wireframe } = dialogOptions();
        closeDialog();
        requestAnimationFrame(() => requestAnimationFrame(() => {
            let result = null;
            longCalculation(() => {
                try {
                    result = runPipeline(lattice, wireframe);
                    if (!includeHelixPos)
                        writeScadnanoFile(name, result);
                }
                catch (err) {
                    result = null;
                    notify(`Scadnano export failed: ${err}`, 'alert');
                }
            }, () => {
                if (!includeHelixPos || !result)
                    return;
                window.currentScadnanoHelixPos = new Map(result.helixPos);
                clearHistory(); // fresh pipeline run: start a new, empty history
                showGrid(result.helixPos, result.latticeType);
            });
        }));
    }
    scadnanoExport.handleDialogExport = handleDialogExport;
    // Grid pane "Export as scadnano": same file, with the edited positions.
    function exportFromGridView(helixPosInput) {
        const helixPos = toHelixPosMap(helixPosInput ?? window.currentScadnanoHelixPos);
        if (!helixPos || helixPos.size === 0) {
            notify('No edited helix positions available to export.', 'warning');
            return;
        }
        const { name, lattice, wireframe } = dialogOptions();
        try {
            writeScadnanoFile(name, layout ?? runPipeline(lattice, wireframe), helixPos);
        }
        catch (err) {
            notify(`Scadnano export failed: ${err}`, 'alert');
        }
    }
    scadnanoExport.exportFromGridView = exportFromGridView;
    function toggleGridDropdown(checkboxElement) {
        const dropdown = document.getElementById('scadnanoGrid');
        if (dropdown)
            dropdown.disabled = !checkboxElement.checked;
    }
    scadnanoExport.toggleGridDropdown = toggleGridDropdown;
    // ── grid pane ────────────────────────────────────────────────────────────
    function toHelixPosMap(value) {
        if (!value)
            return null;
        const entries = value instanceof Map ? [...value.entries()]
            : Array.isArray(value) ? value
                : typeof value === 'object' ? Object.entries(value) : [];
        const out = new Map();
        entries.forEach(entry => {
            const id = Number(entry?.[0]), pos = entry?.[1];
            if (!Array.isArray(pos) || pos.length < 2)
                return;
            const col = Number(pos[0]), row = Number(pos[1]);
            if (Number.isFinite(id) && Number.isFinite(col) && Number.isFinite(row))
                out.set(id, [col, row]);
        });
        return out;
    }
    function canvas() {
        return document.getElementById('scadnanoGridCanvas');
    }
    function resizeCanvas() {
        const pane = document.getElementById('scadnanoGridPane');
        const c = canvas();
        if (!pane || !c)
            return;
        c.width = pane.clientWidth;
        c.height = pane.clientHeight;
        if (typeof gridEditor?.resize === 'function')
            gridEditor.resize();
    }
    function publishHelixPos() {
        if (typeof gridEditor?.getNodes !== 'function')
            return;
        const out = new Map();
        gridEditor.getNodes().forEach(n => {
            out.set(Number(n.id), [Number(n.col), Number(n.row)]);
        });
        window.currentScadnanoHelixPos = out;
        if (layout)
            layout.helixPos = out;
    }
    // Grid node ids are helix indices, so they index layout.helices directly.
    function selectHelicesInScene(ids) {
        const helices = getHelices();
        const select = window.api?.selectElements;
        if (!helices || typeof select !== 'function')
            return;
        const nts = [];
        ids.forEach(id => { if (Array.isArray(helices[id]))
            nts.push(...helices[id]); });
        if (nts.length)
            select(nts);
    }
    function ensureEditor(gridType) {
        if (gridEditor && gridEditorType === gridType)
            return gridEditor;
        if (gridEditor && typeof gridEditor.dispose === 'function')
            gridEditor.dispose();
        gridEditor = null;
        gridEditorType = null;
        const c = canvas();
        const Ctor = window.scadnano?.[gridType === 'square' ? 'SquareEditor' : 'HoneycombEditor'];
        if (!c || typeof Ctor !== 'function')
            return null;
        gridEditor = new Ctor(c);
        gridEditorType = gridType;
        gridEditor.onNodesChanged = () => publishHelixPos();
        // Genuine user drags only (suppressed for programmatic moves). Record one undo entry per
        // drag: snapshot the post-drag state, then rewind every moved node to its pre-drag cell so
        // the stored snapshot is the state as it was BEFORE the drag.
        gridEditor.onNodeMoved = (info) => {
            const before = snapshot();
            if (!before)
                return;
            (info.moves ?? [info]).forEach(move => {
                const entry = before.helixPos.find(([id]) => id === Number(move.id));
                if (entry)
                    entry[1] = [move.from[0], move.from[1]];
            });
            pushHistory(before);
        };
        gridEditor.onNodeSelected = (node) => {
            if (suppressSelectionCallback)
                return;
            const id = Number(node?.id);
            if (Number.isFinite(id))
                selectHelicesInScene([id]);
        };
        gridEditor.onSelectionChanged = (ids) => {
            if (!suppressSelectionCallback)
                selectHelicesInScene(ids);
        };
        return gridEditor;
    }
    function showGrid(helixPosInput, gridTypeInput) {
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
        if (typeof editor.setConnections === 'function')
            editor.setConnections(connections);
        publishHelixPos();
        // loadFromHelixPos replaces every node, so lock colours and the editor's guard set
        // have to be reapplied. Deferred because loadFromHelixPos fits the view on a timeout.
        applyLockedColors();
    }
    scadnanoExport.showGrid = showGrid;
    function hideGrid() {
        document.body.classList.remove('scadnano-grid-open');
        if (focusToggleEl()?.checked)
            focusOnHelixToggle(false); // restore view if closed while focused
        clearHistory(); // history lives only for the duration of the pane session
        lockedHelices.clear();
        if (gridEditor)
            gridEditor.lockedHelices = new Set();
    }
    scadnanoExport.hideGrid = hideGrid;
    // ── lock / unlock ────────────────────────────────────────────────────────
    // color the nodes based on their lock status.
    function applyLockedColors() {
        if (typeof gridEditor?.setNodeColor !== 'function')
            return;
        const nodes = typeof gridEditor.getNodes === 'function' ? gridEditor.getNodes() : [];
        nodes.forEach(n => {
            const id = Number(n.id);
            gridEditor.setNodeColor(id, lockedHelices.has(id) ? LOCKED_COLOR : UNLOCKED_COLOR);
        });
        gridEditor.lockedHelices = new Set(lockedHelices);
    }
    // "Lock/Unlock Helices" toolbar button. Toggles every selected helix.
    function lockSelectedHelices() {
        if (typeof gridEditor?.setNodeColor !== 'function')
            return;
        const ids = typeof gridEditor.getSelectedHelixIds === 'function'
            ? gridEditor.getSelectedHelixIds()
            : [];
        if (ids.length === 0) {
            notify('Select one or more helices in the grid first.', 'warning');
            return;
        }
        ids.forEach(id => {
            if (lockedHelices.has(id))
                lockedHelices.delete(id);
            else
                lockedHelices.add(id);
        });
        applyLockedColors();
        notify(`${lockedHelices.size} helix/helices locked.`, 'success');
    }
    scadnanoExport.lockSelectedHelices = lockSelectedHelices;
    const HISTORY_CAP = 40;
    let undoStack = [];
    let redoStack = [];
    // Clone the current layout + lock state. Returns null when there's no layout to capture.
    function snapshot() {
        if (!layout)
            return null;
        return {
            helices: layout.helices.map(h => h.map(nt => nt.id)),
            grid: [...layout.grid.entries()].map(([id, m]) => [id, { helixId: m.helixId, offset: m.offset, direction: m.direction }]),
            helixPos: [...layout.helixPos.entries()].map(([id, p]) => [id, [p[0], p[1]]]),
            locked: [...lockedHelices],
            latticeType: layout.latticeType
        };
    }
    // Swap a snapshot back into the live layout and redraw. Nucleotide objects are looked up fresh
    // from the global elements map (they never leave it), so only their grouping is restored.
    function restore(snap) {
        if (!layout)
            return;
        layout.helices = snap.helices.map(ids => {
            const arr = [];
            ids.forEach(id => { const e = elements.get(id); if (e instanceof Nucleotide)
                arr.push(e); });
            return arr;
        });
        const grid = new Map();
        snap.grid.forEach(([id, m]) => grid.set(id, { helixId: m.helixId, offset: m.offset, direction: m.direction }));
        layout.grid = grid;
        const helixPos = new Map();
        snap.helixPos.forEach(([id, p]) => helixPos.set(id, [p[0], p[1]]));
        layout.helixPos = helixPos;
        lockedHelices = new Set(snap.locked);
        layout.latticeType = snap.latticeType;
        connections = helixPairs(toscad.collectCrossovers(layout.grid).crossovers);
        window.currentScadnanoHelixPos = new Map(layout.helixPos);
        showGrid(layout.helixPos, layout.latticeType);
    }
    function updateHistoryButtons() {
        const undoBtn = document.getElementById('scadnanoGridUndoBtn');
        const redoBtn = document.getElementById('scadnanoGridRedoBtn');
        if (undoBtn)
            undoBtn.disabled = undoStack.length === 0;
        if (redoBtn)
            redoBtn.disabled = redoStack.length === 0;
    }
    // Push a pre-operation snapshot. Any fresh operation invalidates the redo branch. Pass the state
    // captured BEFORE the mutation was applied.
    function pushHistory(snap) {
        if (!snap)
            return;
        undoStack.push(snap);
        if (undoStack.length > HISTORY_CAP)
            undoStack.shift();
        redoStack = [];
        updateHistoryButtons();
    }
    function clearHistory() {
        undoStack = [];
        redoStack = [];
        updateHistoryButtons();
    }
    function undo() {
        if (undoStack.length === 0) {
            notify('Nothing to undo.', 'warning');
            return;
        }
        const current = snapshot();
        if (current)
            redoStack.push(current);
        restore(undoStack.pop());
        updateHistoryButtons();
    }
    scadnanoExport.undo = undo;
    function redo() {
        if (redoStack.length === 0) {
            notify('Nothing to redo.', 'warning');
            return;
        }
        const current = snapshot();
        if (current) {
            undoStack.push(current);
            if (undoStack.length > HISTORY_CAP)
                undoStack.shift();
        }
        restore(redoStack.pop());
        updateHistoryButtons();
    }
    scadnanoExport.redo = redo;
    // ── combine ────────────────────────────────────────────────────────────────
    // Remaps a helix-position map through a combineHelices result: merged-away helices drop out
    // (their nucleotides now live in the kept helix, which keeps its own cell) and every survivor
    // shifts down to its new index.
    function remapHelixPos(source, mergedIdx, idRemap) {
        const removed = new Set(mergedIdx);
        const out = new Map();
        source.forEach((pos, oldId) => {
            if (removed.has(oldId))
                return;
            const next = idRemap.get(oldId);
            out.set(next !== undefined ? next : oldId, pos);
        });
        return out;
    }
    // "Combine Helices" toolbar button. Folds every selected helix into the lowest-numbered one via
    // helix.combineHelices, then rebuilds the grid view on the remapped layout. combineHelices takes
    // the whole index list and collapses it into a single helix in one pass, so there's nothing to
    // recurse over here — one call merges all of the selection.
    function combineSelectedHelices() {
        const ids = typeof gridEditor?.getSelectedHelixIds === 'function'
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
        const newLocked = new Set();
        lockedHelices.forEach(id => {
            if (removed.has(id)) {
                newLocked.add(keptIdx);
                return;
            }
            const next = idRemap.get(id);
            if (next !== undefined)
                newLocked.add(next);
        });
        lockedHelices = newLocked;
        // Rebuild the grid pane on the remapped layout (also reapplies connections + lock colours).
        showGrid(newHelixPos, layout.latticeType);
        notify(`Combined helix ${mergedIdx.join(', ')} into helix ${keptIdx}.`, 'success');
    }
    scadnanoExport.combineSelectedHelices = combineSelectedHelices;
    // ── split ────────────────────────────────────────────────────────────────
    // Which helix the currently selected nucleotides belong to, per the grid marks (the same
    // membership toscad.splitHelix trusts). Returns null when nothing usable is selected, or the
    // selection straddles more than one helix — split only makes sense within a single helix.
    function selectedHelixForSplit() {
        const nucleotides = [];
        const helixIds = new Set();
        selectedBases.forEach(e => {
            if (!(e instanceof Nucleotide))
                return;
            const mark = layout.grid.get(e.id);
            if (!mark)
                return;
            nucleotides.push(e);
            helixIds.add(mark.helixId);
        });
        if (nucleotides.length === 0)
            return null;
        if (helixIds.size > 1)
            return null;
        return { helixId: [...helixIds][0], nucleotides };
    }
    // "Split from selected" toolbar button. Moves the selected nucleotides out of their helix into
    // a new one via toscad.splitHelix. splitHelix appends the new helix (newHelixId = helices.length)
    // and touches only the moved nts' grid marks, so no existing index shifts — the cascade is just
    // recompute-connections + redraw, with no remap of helixPos or locks.
    function splitSelectedHelix() {
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
        const result = toscad.splitHelix(layout.grid, layout.helixPos, selection.helixId, layout.helices, selection.nucleotides);
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
    scadnanoExport.splitSelectedHelix = splitSelectedHelix;
    // ── focus on helix ─────────────────────────────────────────────────────────
    const FOCUS_TOGGLE_ID = 'scadnanoGridSplitHelicesToggle';
    function focusToggleEl() {
        return document.getElementById(FOCUS_TOGGLE_ID);
    }
    // Splitting only makes sense while focused on a helix, so the button tracks the toggle.
    function setSplitEnabled(on) {
        const btn = document.getElementById('scadnanoGridSplitFromSelectedBtn');
        if (btn)
            btn.disabled = !on;
    }
    // "Focus on helix" switch. On: invert the selection and reuse the existing "Hide Selection"
    // command (toggleVisArbitrary) to hide everything except the selected helix. Off: unhide all.
    function focusOnHelixToggle(checked) {
        if (checked) {
            if (selectedBases.size === 0) {
                notify('Select a helix first to focus on it.', 'warning');
                const chk = focusToggleEl();
                if (chk)
                    chk.checked = false;
                return;
            }
            invertSelection(); // select everything outside the helix
            toggleVisArbitrary(); // hide it (also clears the selection)
        }
        else {
            const chk = focusToggleEl();
            if (chk)
                chk.checked = false;
            clearSelection(); // force the "unhide all" branch
            toggleVisArbitrary(); // restore the full view
        }
        setSplitEnabled(checked);
    }
    scadnanoExport.focusOnHelixToggle = focusOnHelixToggle;
    // ── recalculate ──────────────────────────────────────────────────────────
    // Locked helices become hard relative pins. Only differences are meaningful, so the
    // lowest locked id is the reference and every other locked helix contributes one pin
    // carrying its current offset from it. Fewer than two locks means nothing to constrain.
    function pinsFromLocked() {
        if (lockedHelices.size < 2)
            return [];
        const pos = window.currentScadnanoHelixPos;
        if (!pos)
            return [];
        const ids = [...lockedHelices].filter(id => pos.has(id)).sort((x, y) => x - y);
        if (ids.length < 2)
            return [];
        const [refId, ...rest] = ids;
        const ref = pos.get(refId);
        return rest.map(id => {
            const p = pos.get(id);
            return { a: refId, b: id, dCol: p[0] - ref[0], dRow: p[1] - ref[1] };
        });
    }
    // Absolute (col,row) cells for every locked helix. Pins (above) only carry
    // relative offsets; these hand kruskals the user's exact placement so the
    // locked component is anchored on its true cell/parity instead of drifting
    // to wherever the 2-coloring seed lands it.
    function anchorsFromLocked() {
        const anchors = new Map();
        if (lockedHelices.size < 2)
            return anchors;
        const pos = window.currentScadnanoHelixPos;
        if (!pos)
            return anchors;
        for (const id of lockedHelices) {
            const p = pos.get(id);
            if (p)
                anchors.set(id, [p[0], p[1]]);
        }
        return anchors;
    }
    // "Recalculate Grid" toolbar button. Re-runs ONLY the positioning stage (angles → kruskals →
    // posCorr5) on the CURRENT layout.helices / layout.grid — never findHelices — so any manual
    // combines or splits the user has made are preserved. Locked helices feed in as pins to hold
    // their relative placement. Helix ids are untouched, so locks stay valid without re-finding.
    function recalculateGrid() {
        if (!layout) {
            notify('Open the grid view first before recalculating.', 'warning');
            return;
        }
        publishHelixPos();
        const before = snapshot();
        const pins = pinsFromLocked();
        const anchors = anchorsFromLocked();
        const latticeType = layout.latticeType;
        let ok = false;
        longCalculation(() => {
            try {
                const networkMap = toscad.getAngles(layout.grid, layout.helices, latticeType);
                const votes = toscad.voteOrientations(networkMap, layout.grid, latticeType);
                const kr = toscad.kruskals(networkMap, layout.grid, latticeType, votes, pins, anchors);
                const pinnedIds = new Set();
                pins.forEach(p => { pinnedIds.add(p.a); pinnedIds.add(p.b); });
                layout.helixPos = toscad.posCorr5(kr, layout.helices, pinnedIds);
                connections = helixPairs(toscad.collectCrossovers(layout.grid).crossovers);
                ok = true;
            }
            catch (err) {
                ok = false;
                notify(`Recalculate grid failed: ${err}`, 'alert');
            }
        }, () => {
            if (!ok)
                return;
            pushHistory(before);
            window.currentScadnanoHelixPos = new Map(layout.helixPos);
            showGrid(layout.helixPos, latticeType);
            applyLockedColors();
            notify(`Grid recalculated${pins.length ? ` with ${pins.length} pin(s)` : ''}.`, 'success');
        });
    }
    scadnanoExport.recalculateGrid = recalculateGrid;
    // Used by "Helix" selection mode in base_selector.ts as well as the grid pane.
    // Falls back to plain helix detection when no layout has been computed yet.
    function getHelices() {
        if (layout?.helices.length)
            return layout.helices;
        try {
            return helix.findHelices(nucleotideMap(), TOLERANCE).helices ?? null;
        }
        catch (err) {
            notify(`Unable to map grid helix selection: ${err}`, 'warning');
            return null;
        }
    }
    scadnanoExport.getHelices = getHelices;
    function selectHelixFromNucleotide(nucleotideInput, additive = false) {
        if (!document.body.classList.contains('scadnano-grid-open'))
            return;
        if (typeof gridEditor?.selectNodeById !== 'function')
            return;
        const nt = nucleotideInput instanceof Nucleotide ? nucleotideInput
            : typeof nucleotideInput === 'number' ? elements.get(nucleotideInput) : null;
        if (!(nt instanceof Nucleotide))
            return;
        const helices = getHelices();
        if (!helices)
            return;
        const helixId = api.helix.findHelixID(nt.id, helices);
        if (helixId === null)
            return;
        suppressSelectionCallback = true;
        try {
            gridEditor.selectNodeById(helixId, additive);
        }
        finally {
            suppressSelectionCallback = false;
        }
    }
    scadnanoExport.selectHelixFromNucleotide = selectHelixFromNucleotide;
    // ── wiring ───────────────────────────────────────────────────────────────
    function initPaneControls() {
        // Toolbar tabs live in scadnano_gridview. Guard the namespace: this runs at
        // script-load time and a ReferenceError here would break the inline handlers.
        if (typeof scadnano !== 'undefined' && typeof scadnano.initGridToolbar === 'function') {
            scadnano.initGridToolbar();
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
        window.addEventListener('keydown', (e) => {
            if (!document.body.classList.contains('scadnano-grid-open'))
                return;
            if (!(e.ctrlKey || e.metaKey))
                return;
            const key = e.key.toLowerCase();
            if (key !== 'z' && key !== 'y')
                return;
            e.preventDefault();
            e.stopImmediatePropagation();
            if (key === 'y' || e.shiftKey)
                redo();
            else
                undo();
        }, true);
        let resizing = false;
        document.getElementById('scadnanoGridResizeHandle')?.addEventListener('mousedown', (e) => {
            resizing = true;
            document.body.classList.add('scadnano-grid-resizing');
            e.preventDefault();
        });
        window.addEventListener('mousemove', (e) => {
            if (!resizing)
                return;
            const max = Math.max(240, Math.floor(window.innerWidth * 0.75));
            const width = Math.max(240, Math.min(max, Math.round(e.clientX)));
            document.documentElement.style.setProperty('--scadnano-pane-width', `${width}px`);
            resizeCanvas();
        });
        window.addEventListener('mouseup', () => {
            if (!resizing)
                return;
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
    window.scadnanoFocusOnHelixToggle = (chk) => focusOnHelixToggle(Boolean(chk?.checked));
    window.scadnanoSelectHelixFromNucleotide = (nt, additive) => selectHelixFromNucleotide(nt, additive === true);
})(scadnanoExport || (scadnanoExport = {}));
// Named wrappers for the inline handlers in index.html.
function scadnanoDialogExport() {
    scadnanoExport.handleDialogExport();
}
function toggleGridDropdown(checkboxElement) {
    scadnanoExport.toggleGridDropdown(checkboxElement);
}
function scadnanoFocusOnHelixToggle(checkboxElement) {
    scadnanoExport.focusOnHelixToggle(Boolean(checkboxElement?.checked));
}
