"use strict";
/// <reference path="../typescript_definitions/index.d.ts" />
class ScadnanoExportManager {
    currentScadnanoHelices = null;
    currentScadnanoConnections = [];
    currentScadnanoLayout = null;
    scadnanoGridEditor = null;
    scadnanoGridEditorType = null;
    suppressNodeSelectedCallback = false;
    // Edit-history journal. Lives only while the grid pane is open; cleared on export / close /
    // reload. Holds combine and move operations in the readable JSON form discussed.
    history = this.createEmptyHistory();
    createEmptyHistory() {
        return { entries: {}, order: [], cursor: 0, nextMove: 0, nextCombine: 0 };
    }
    clearHistory() {
        this.history = this.createEmptyHistory();
        this.refreshHistoryButtons();
    }
    // Drop any redo-tail entries before a new operation extends the journal.
    truncateRedoTail() {
        if (this.history.cursor >= this.history.order.length)
            return;
        const tail = this.history.order.splice(this.history.cursor);
        tail.forEach(key => { delete this.history.entries[key]; });
    }
    pushHistoryEntry(entry) {
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
    canUndo() { return this.history.cursor > 0; }
    canRedo() { return this.history.cursor < this.history.order.length; }
    refreshHistoryButtons() {
        const undoBtn = document.getElementById('scadnanoGridUndoBtn');
        const redoBtn = document.getElementById('scadnanoGridRedoBtn');
        if (undoBtn)
            undoBtn.disabled = !this.canUndo();
        if (redoBtn)
            redoBtn.disabled = !this.canRedo();
    }
    constructor() {
        this.initScadnanoGridPaneControls();
    }
    handleDialogExport() {
        const options = this.readDialogOptions();
        if (!options)
            return;
        this.closeScadnanoDialog();
        requestAnimationFrame(() => {
            requestAnimationFrame(() => {
                this.runDialogExport(options);
            });
        });
    }
    exportFromGridView(helixPosInput) {
        const { name, gridType, wireframe } = this.readCurrentExportTarget();
        const map = this.normalizeHelixPosMap(helixPosInput ?? window.currentScadnanoHelixPos);
        if (!map || map.size === 0) {
            notify('No edited helix positions available to export.', 'warning');
            return;
        }
        try {
            this.exportToScadnano(name, gridType, map, wireframe);
        }
        catch (err) {
            notify(`Scadnano export failed: ${err}`, 'alert');
        }
    }
    // Triggered by the "Combine" button in the grid view.
    // Merges the helices currently selected in the grid editor into the lowest-numbered one,
    combineSelectedHelicesFromGridView() {
        const editor = this.scadnanoGridEditor;
        if (!editor) {
            notify('Open the scadnano grid view before combining helices.', 'warning');
            return;
        }
        const ids = typeof editor.getSelectedHelixIds === 'function'
            ? editor.getSelectedHelixIds()
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
                        const inverse = new Map();
                        preCombineShifts.forEach((delta, hid) => inverse.set(hid, -delta));
                        toscad.applyCombineShifts(grid, inverse);
                    }
                    notify(`Cannot combine: helices ${h1} and ${h2} are not disjoint (overlapping nucleotides on the grid).`, 'alert');
                    return;
                }
            }
        }
        const willRemoveOld = sortedIds.slice(1);
        const editorNodesBefore = typeof editor.getNodes === 'function'
            ? editor.getNodes()
            : [];
        const cellByOldId = new Map();
        editorNodesBefore.forEach(n => cellByOldId.set(Number(n.id), [Number(n.col), Number(n.row)]));
        const removedSnapshots = willRemoveOld.map(oldIdx => {
            const cell = cellByOldId.get(oldIdx) ?? [0, 0];
            const slot = helices[oldIdx] || [];
            return {
                oldIdx,
                col: cell[0],
                row: cell[1],
                ntIds: slot.map(nt => nt.id)
            };
        });
        const connectionsBefore = this.currentScadnanoConnections.map(([a, b]) => [a, b]);
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
        const entry = {
            op: 'combine',
            kept: keptIdx,
            indices: sortedIds,
            idRemap: Array.from(idRemap.entries()),
            removed: removedSnapshots,
            connections: {
                before: connectionsBefore,
                after: this.currentScadnanoConnections.map(([a, b]) => [a, b])
            },
            preCombineShifts: Array.from(preCombineShifts.entries())
        };
        this.pushHistoryEntry(entry);
        if (preCombineShifts.size > 0) {
            const shiftSummary = Array.from(preCombineShifts.entries())
                .map(([hid, d]) => `${hid}:${d > 0 ? '+' : ''}${d}`)
                .join(', ');
            notify(`Combined helix ${mergedIdx.join(', ')} into helix ${keptIdx}. Shifted overlapping helices to make them disjoint: ${shiftSummary}.`);
        }
        else {
            notify(`Combined helix ${mergedIdx.join(', ')} into helix ${keptIdx}.`);
        }
    }
    // ── Undo / Redo ─────────────────────────────────────────────────────────
    undoFromGridView() {
        if (!this.canUndo())
            return;
        const key = this.history.order[this.history.cursor - 1];
        const entry = this.history.entries[key];
        if (!entry)
            return;
        if (entry.op === 'move')
            this.undoMove(entry);
        else if (entry.op === 'combine')
            this.undoCombine(entry);
        this.history.cursor -= 1;
        this.refreshHistoryButtons();
    }
    redoFromGridView() {
        if (!this.canRedo())
            return;
        const key = this.history.order[this.history.cursor];
        const entry = this.history.entries[key];
        if (!entry)
            return;
        if (entry.op === 'move')
            this.redoMove(entry);
        else if (entry.op === 'combine')
            this.redoCombine(entry);
        this.history.cursor += 1;
        this.refreshHistoryButtons();
    }
    getHistorySnapshot() {
        // Returned by reference for inspection / debugging. Don't mutate from outside.
        return this.history;
    }
    redoMove(entry) {
        const editor = this.scadnanoGridEditor;
        if (!editor || typeof editor.moveNodeById !== 'function')
            return;
        editor.moveNodeById(entry.id, entry.to[0], entry.to[1]);
        this.publishCurrentHelixPosFromEditor();
        this.syncLayoutHelixPosFromEditor();
    }
    undoMove(entry) {
        const editor = this.scadnanoGridEditor;
        if (!editor || typeof editor.moveNodeById !== 'function')
            return;
        editor.moveNodeById(entry.id, entry.from[0], entry.from[1]);
        this.publishCurrentHelixPosFromEditor();
        this.syncLayoutHelixPosFromEditor();
    }
    // Replay a combine forward by re-running combineHelices on the saved indices, then handing
    // off to the same cascade helper the live combine path uses. The only difference is that
    // redo applies the cached post-merge connection list directly (passed via override), instead
    // of re-deriving it through the remap.
    redoCombine(entry) {
        const editor = this.scadnanoGridEditor;
        if (!editor)
            return;
        const helices = this.ensureScadnanoHelicesCache();
        if (!helices)
            return;
        if (!this.currentScadnanoLayout)
            return;
        // Reapply the pre-combine offset shifts BEFORE combining, mirroring the live path.
        // Keyed by pre-merge helix id, which is what the grid marks still carry at this point.
        if (Array.isArray(entry.preCombineShifts) && entry.preCombineShifts.length > 0) {
            const shifts = new Map(entry.preCombineShifts);
            toscad.applyCombineShifts(this.currentScadnanoLayout.grid, shifts);
        }
        const result = helix.combineHelices(helices, entry.indices, this.currentScadnanoLayout.grid);
        if (!result)
            return;
        this.applyCombineCascade(result, entry.connections.after);
    }
    // Reverse a combine: re-insert the removed helix slots, pull their nucleotides back out of the
    // kept helix, then renumber everything in the GridMap and editor through the inverse remap.
    undoCombine(entry) {
        const editor = this.scadnanoGridEditor;
        if (!editor)
            return;
        const helices = this.ensureScadnanoHelicesCache();
        if (!helices)
            return;
        if (!this.currentScadnanoLayout)
            return;
        // splitHelices mutates `helices` AND `currentScadnanoLayout.grid` in place, mirroring
        // what combineHelices did on the forward path. It hands back the inverse remap the
        // cascade needs to renumber editor nodes that are still in post-merge id space.
        const result = helix.splitHelices(helices, this.currentScadnanoLayout.grid, entry);
        if (!result)
            return;
        // AFTER splitting, grid marks are back in their pre-combine helix ids, but any offset
        // shifts applied by the shift-solver are still present. Reverse them with negated deltas.
        if (Array.isArray(entry.preCombineShifts) && entry.preCombineShifts.length > 0) {
            const inverse = new Map();
            entry.preCombineShifts.forEach(([hid, delta]) => inverse.set(hid, -delta));
            toscad.applyCombineShifts(this.currentScadnanoLayout.grid, inverse);
        }
        this.applyCombineCascadeInverse(entry, result.inverseRemap);
    }
    remapCachedGridIds(remap) {
        if (!this.currentScadnanoLayout)
            return;
        this.currentScadnanoLayout.grid.forEach(mark => {
            const newId = remap(mark.helixId);
            if (newId === null)
                return;
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
    applyCombineCascade(result, connectionsOverride) {
        const editor = this.scadnanoGridEditor;
        if (!editor)
            return;
        const { keptIdx, mergedIdx, idRemap } = result;
        const removed = new Set(mergedIdx);
        const remapId = (oldId) => {
            if (removed.has(oldId))
                return keptIdx;
            const next = idRemap.get(oldId);
            return next === undefined ? null : next;
        };
        // Renumber editor nodes through the remap, drop merged-away duplicates. Survivors keep
        // their col/row — the kept helix stays in its original cell.
        if (typeof editor.getNodes === 'function') {
            const nodes = editor.getNodes();
            const remapped = [];
            const seenNew = new Set();
            nodes.forEach(node => {
                const oldId = Number(node.id);
                const newId = remapId(oldId);
                if (newId === null)
                    return;
                if (seenNew.has(newId))
                    return; // skip duplicates (e.g. merged-away nodes mapping to keptIdx)
                seenNew.add(newId);
                remapped.push({ ...node, id: newId, label: String(newId) });
            });
            if (typeof editor.setNodes === 'function')
                editor.setNodes(remapped);
        }
        // Connections: redo applies the cached post-merge list verbatim; live combines remap the
        // current list (drop self-loops + dedup).
        if (connectionsOverride) {
            this.currentScadnanoConnections = connectionsOverride.map(([a, b]) => [a, b]);
        }
        else {
            const remappedConnKeys = new Set();
            const updatedConnections = [];
            this.currentScadnanoConnections.forEach(([from, to]) => {
                const a = remapId(from);
                const b = remapId(to);
                if (a === null || b === null)
                    return;
                if (a === b)
                    return;
                const key = a < b ? `${a}:${b}` : `${b}:${a}`;
                if (remappedConnKeys.has(key))
                    return;
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
        if (typeof editor.clearSelection === 'function')
            editor.clearSelection();
        // Refresh the published helix-pos map from the editor (it just lost a few nodes).
        this.publishCurrentHelixPosFromEditor();
    }
    // Inverse of applyCombineCascade. Run after splitHelices has already restored helices + grid
    // back to their pre-merge state. Renumbers existing editor nodes through the inverse remap,
    // re-adds the removed nodes at their saved cells, and restores the pre-merge connection list
    // from the journal entry.
    applyCombineCascadeInverse(entry, inverseRemap) {
        const editor = this.scadnanoGridEditor;
        if (!editor)
            return;
        // Renumber existing editor nodes (still in post-merge id space) back to old ids, then
        // re-insert the removed nodes at their captured cells.
        const currentNodes = typeof editor.getNodes === 'function'
            ? editor.getNodes()
            : [];
        const restoredNodes = [];
        currentNodes.forEach(node => {
            const oldId = inverseRemap(Number(node.id));
            restoredNodes.push({ ...node, id: oldId, label: String(oldId) });
        });
        entry.removed.forEach(slot => {
            restoredNodes.push({ id: slot.oldIdx, col: slot.col, row: slot.row, label: String(slot.oldIdx) });
        });
        if (typeof editor.setNodes === 'function')
            editor.setNodes(restoredNodes);
        // Restore the pre-merge connection list verbatim.
        this.currentScadnanoConnections = entry.connections.before.map(([a, b]) => [a, b]);
        if (typeof editor.setConnections === 'function') {
            editor.setConnections(this.currentScadnanoConnections);
        }
        this.syncLayoutHelixPosFromEditor();
        this.publishCurrentHelixPosFromEditor();
        if (typeof editor.clearSelection === 'function')
            editor.clearSelection();
    }
    syncLayoutHelixPosFromEditor() {
        if (!this.currentScadnanoLayout)
            return;
        const editor = this.scadnanoGridEditor;
        if (!editor || typeof editor.getNodes !== 'function')
            return;
        const refreshed = new Map();
        editor.getNodes().forEach(node => {
            refreshed.set(Number(node.id), [Number(node.col), Number(node.row)]);
        });
        this.currentScadnanoLayout.helixPos = refreshed;
    }
    showGridFromHelixPos(helixPosInput, gridTypeInput) {
        const pane = this.getScadnanoGridPane();
        if (!pane) {
            notify('Scadnano grid pane is unavailable.', 'alert');
            return;
        }
        const gridType = gridTypeInput === 'square' ? 'square' : 'honeycomb';
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
    hideScadnanoGridPane() {
        document.body.classList.remove('scadnano-grid-open');
        this.clearHistory();
    }
    toggleGridDropdown(checkboxElement) {
        const gridDropdown = document.getElementById('scadnanoGrid');
        if (gridDropdown) {
            gridDropdown.disabled = !checkboxElement.checked;
        }
    }
    getHelices() {
        return this.ensureScadnanoHelicesCache();
    }
    selectHelixFromNucleotide(nucleotideInput) {
        if (!document.body.classList.contains('scadnano-grid-open'))
            return;
        if (!this.scadnanoGridEditor || typeof this.scadnanoGridEditor.selectNodeById !== 'function')
            return;
        const nucleotide = nucleotideInput instanceof Nucleotide
            ? nucleotideInput
            : (typeof nucleotideInput === 'number' ? elements.get(nucleotideInput) : null);
        if (!nucleotide || !(nucleotide instanceof Nucleotide))
            return;
        const helices = this.ensureScadnanoHelicesCache();
        if (!helices)
            return;
        const helixId = toscad.findHelixID(nucleotide.id, helices);
        if (helixId === null)
            return;
        this.suppressNodeSelectedCallback = true;
        try {
            this.scadnanoGridEditor.selectNodeById(helixId);
        }
        finally {
            this.suppressNodeSelectedCallback = false;
        }
    }
    runDialogExport(options) {
        if (!options.includeHelixPos) {
            this.runScadnanoLongCalculation(() => {
                try {
                    this.exportToScadnano(options.name, options.gridType, undefined, options.wireframe);
                }
                catch (err) {
                    notify(`Scadnano export failed: ${err}`, 'alert');
                }
            });
            return;
        }
        let helixPos = null;
        let resolvedGridType = 'honeycomb';
        let failed = false;
        this.runScadnanoLongCalculation(() => {
            try {
                const result = this.calculateScadnanoHelixPos(options.gridType, options.wireframe);
                helixPos = result.helixPos;
                resolvedGridType = result.latticeType;
                if (helixPos) {
                    window.currentScadnanoHelixPos = this.cloneHelixPosMap(helixPos);
                }
            }
            catch (err) {
                failed = true;
                notify(`Scadnano export failed: ${err}`, 'alert');
            }
        }, () => {
            if (failed || !helixPos)
                return;
            this.showGridFromHelixPos(helixPos, resolvedGridType);
        });
    }
    readDialogOptions() {
        const nameInput = document.getElementById('scadnanoFilename');
        const helixPosCheckbox = document.getElementById('scadnanoIncludeHPos');
        const scadnanoGrid = document.getElementById('scadnanoGrid');
        const wireframeCheckbox = document.getElementById('scadnanoWireframe');
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
    readCurrentExportTarget() {
        const nameInput = document.getElementById('scadnanoFilename');
        const scadnanoGrid = document.getElementById('scadnanoGrid');
        const wireframeCheckbox = document.getElementById('scadnanoWireframe');
        return {
            name: nameInput?.value.trim() || 'output',
            gridType: this.normalizeRequestedGridType(scadnanoGrid?.value),
            wireframe: Boolean(wireframeCheckbox?.checked),
        };
    }
    // Dialog/dropdown value — accepts 'automatic' for downstream resolution.
    normalizeRequestedGridType(value) {
        if (value === 'square')
            return 'square';
        if (value === 'honeycomb')
            return 'honeycomb';
        return 'automatic';
    }
    runScadnanoLongCalculation(calc, callback) {
        const longCalculation = window.view?.longCalculation;
        if (typeof longCalculation === 'function') {
            longCalculation(calc, 'Preparing scadnano export, please be patient...', callback);
            return;
        }
        calc();
        if (callback)
            callback();
    }
    closeScadnanoDialog() {
        let closedByMetro = false;
        const metroDialog = window?.Metro?.dialog;
        if (metroDialog && typeof metroDialog.close === 'function') {
            try {
                metroDialog.close('#scadnanoDialog');
                closedByMetro = true;
            }
            catch (err) {
                console.warn('Failed to close scadnano dialog via Metro API:', err);
            }
        }
        if (!closedByMetro) {
            const closeBtn = document.querySelector('#scadnanoDialog .js-dialog-close');
            if (closeBtn)
                closeBtn.click();
        }
        const dialogEl = document.getElementById('scadnanoDialog');
        if (dialogEl) {
            dialogEl.classList.remove('open');
            dialogEl.setAttribute('aria-hidden', 'true');
        }
    }
    exportToScadnano(name, gridType, helixPos, wireframe = false) {
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
    getCurrentNucleotideCount() {
        let count = 0;
        elements.forEach((element) => {
            if (element instanceof Nucleotide)
                count += 1;
        });
        return count;
    }
    cloneHelixPosMap(input) {
        const out = new Map();
        input.forEach((value, key) => {
            out.set(Number(key), [Number(value[0]), Number(value[1])]);
        });
        return out;
    }
    calculateScadnanoHelices() {
        const nucleotideElements = new Map();
        elements.forEach((element, id) => {
            if (element instanceof Nucleotide) {
                nucleotideElements.set(id, element);
            }
        });
        const result = helix.findHelices(nucleotideElements, 3);
        const helices = result?.helices ?? [];
        this.notifyHelixCoverageMismatch(helices, nucleotideElements);
        return helices;
    }
    prepareScadnanoLayout(requestedLatticeType, forceRecompute = false, wireframe = false) {
        const nucleotideCount = this.getCurrentNucleotideCount();
        // Cache hit:
        //   - For concrete kinds: cached latticeType must match.
        //   - For 'automatic': any cached latticeType is acceptable (it was
        //     either detected the same way last time or explicitly chosen).
        if (!forceRecompute &&
            this.currentScadnanoLayout &&
            this.currentScadnanoLayout.nucleotideCount === nucleotideCount &&
            this.currentScadnanoLayout.wireframe === wireframe &&
            (requestedLatticeType === 'automatic' ||
                this.currentScadnanoLayout.latticeType === requestedLatticeType)) {
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
        let latticeType;
        if (requestedLatticeType === 'automatic') {
            latticeType = toscad.detectLatticeKind(grid, binderHelices);
            notify(`Auto-detected lattice: ${latticeType}`, 'success');
        }
        else {
            latticeType = requestedLatticeType;
        }
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
    calculateScadnanoHelixPos(latticeType = 'automatic', wireframe = false) {
        const layout = this.prepareScadnanoLayout(latticeType, false, wireframe);
        return {
            helixPos: this.cloneHelixPosMap(layout.helixPos),
            latticeType: layout.latticeType
        };
    }
    notifyHelixCoverageMismatch(helices, inputMap) {
        const helixCount = helices.flat().length;
        const totalCount = inputMap.size;
        if (helixCount === totalCount)
            return;
        const missingCount = totalCount - helixCount;
        notify(`Helix mapping error: ${helixCount}/${totalCount} nucleotides were mapped. Missing ${missingCount} nucleotides.`, 'alert', true);
    }
    buildScadnanoConnections(crossovers) {
        const uniquePairs = new Set();
        const pairs = [];
        for (const [fromHelix, toMap] of crossovers.entries()) {
            for (const [toHelix, counts] of toMap.entries()) {
                const totalConnections = Number(counts?.sameWalk ?? 0) + Number(counts?.diffWalk ?? 0);
                if (totalConnections <= 0)
                    continue;
                const a = Math.min(fromHelix, toHelix);
                const b = Math.max(fromHelix, toHelix);
                if (a === b)
                    continue;
                const key = `${a}:${b}`;
                if (uniquePairs.has(key))
                    continue;
                uniquePairs.add(key);
                pairs.push([a, b]);
            }
        }
        return pairs;
    }
    normalizeHelixPosMap(input) {
        if (!input)
            return null;
        if (input instanceof Map) {
            const out = new Map();
            input.forEach((value, key) => {
                if (!Array.isArray(value) || value.length < 2)
                    return;
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
            const out = new Map();
            input.forEach((entry) => {
                if (!Array.isArray(entry) || entry.length < 2)
                    return;
                const helixId = Number(entry[0]);
                const value = entry[1];
                if (!Array.isArray(value) || value.length < 2)
                    return;
                const col = Number(value[0]);
                const row = Number(value[1]);
                if (Number.isFinite(helixId) && Number.isFinite(col) && Number.isFinite(row)) {
                    out.set(helixId, [col, row]);
                }
            });
            return out;
        }
        if (typeof input === 'object') {
            const out = new Map();
            Object.keys(input).forEach((k) => {
                const value = input[k];
                if (!Array.isArray(value) || value.length < 2)
                    return;
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
    getScadnanoGridPane() {
        return document.getElementById('scadnanoGridPane');
    }
    getScadnanoGridCanvas() {
        return document.getElementById('scadnanoGridCanvas');
    }
    setScadnanoPaneWidth(widthPx) {
        const minW = 240;
        const maxW = Math.max(minW, Math.floor(window.innerWidth * 0.75));
        const clamped = Math.max(minW, Math.min(maxW, Math.round(widthPx)));
        document.documentElement.style.setProperty('--scadnano-pane-width', `${clamped}px`);
    }
    resizeScadnanoGridCanvas() {
        const pane = this.getScadnanoGridPane();
        const canvas = this.getScadnanoGridCanvas();
        if (!pane || !canvas)
            return;
        canvas.width = pane.clientWidth;
        canvas.height = pane.clientHeight;
    }
    mapFromEditorNodes(editor) {
        const out = new Map();
        const nodes = typeof editor.getNodes === 'function' ? editor.getNodes() : [];
        nodes.forEach((node) => {
            out.set(Number(node.id), [Number(node.col), Number(node.row)]);
        });
        return out;
    }
    publishCurrentHelixPosFromEditor() {
        if (!this.scadnanoGridEditor)
            return;
        window.currentScadnanoHelixPos = this.mapFromEditorNodes(this.scadnanoGridEditor);
    }
    ensureScadnanoHelicesCache() {
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
        }
        catch (err) {
            notify(`Unable to map grid helix selection: ${err}`, 'warning');
            return null;
        }
    }
    selectHelixFromGridNode(helixId) {
        const helices = this.ensureScadnanoHelicesCache();
        if (!helices)
            return;
        const helix = helices[helixId];
        if (!Array.isArray(helix) || helix.length === 0)
            return;
        const selectElements = window.api?.selectElements;
        if (typeof selectElements !== 'function')
            return;
        selectElements(helix);
    }
    // Called when the grid editor's multi-select set changes (cmd/ctrl+click).
    // Selects all nucleotides belonging to every currently-selected helix in the 3D scene.
    selectHelicesFromMultiSelect(ids) {
        const helices = this.ensureScadnanoHelicesCache();
        if (!helices)
            return;
        const selectElements = window.api?.selectElements;
        if (typeof selectElements !== 'function')
            return;
        const allNucleos = [];
        ids.forEach(id => {
            const helix = helices[id];
            if (Array.isArray(helix))
                allNucleos.push(...helix);
        });
        // Pass the full combined set in one call so the scene clears its old
        // selection and replaces it with exactly the highlighted helices.
        selectElements(allNucleos);
    }
    ensureScadnanoGridEditor(gridType) {
        if (this.scadnanoGridEditor && this.scadnanoGridEditorType === gridType)
            return this.scadnanoGridEditor;
        if (this.scadnanoGridEditor && this.scadnanoGridEditorType !== gridType) {
            if (typeof this.scadnanoGridEditor.dispose === 'function') {
                this.scadnanoGridEditor.dispose();
            }
            this.scadnanoGridEditor = null;
            this.scadnanoGridEditorType = null;
        }
        const canvas = this.getScadnanoGridCanvas();
        if (!canvas)
            return null;
        const scadnanoNs = window.scadnano;
        if (!scadnanoNs)
            return null;
        const editorCtorName = gridType === 'square' ? 'SquareEditor' : 'HoneycombEditor';
        const EditorCtor = scadnanoNs[editorCtorName];
        if (typeof EditorCtor !== 'function')
            return null;
        this.scadnanoGridEditor = new EditorCtor(canvas);
        this.scadnanoGridEditorType = gridType;
        this.scadnanoGridEditor.onNodesChanged = () => {
            this.publishCurrentHelixPosFromEditor();
        };
        this.scadnanoGridEditor.onNodeSelected = (node) => {
            if (this.suppressNodeSelectedCallback)
                return;
            const helixId = Number(node?.id);
            if (!Number.isFinite(helixId))
                return;
            this.selectHelixFromGridNode(helixId);
        };
        // Multi-select (cmd/ctrl+click): sync all highlighted grid nodes to the 3D scene.
        this.scadnanoGridEditor.onSelectionChanged = (ids) => {
            if (this.suppressNodeSelectedCallback)
                return;
            this.selectHelicesFromMultiSelect(ids);
        };
        // Genuine user-initiated drags push a move entry onto the history journal.
        // Programmatic moves (undo/redo) suppress this callback inside the editor.
        this.scadnanoGridEditor.onNodeMoved = (info) => {
            const id = Number(info?.id);
            if (!Number.isFinite(id))
                return;
            const fromCol = Number(info.from?.[0]);
            const fromRow = Number(info.from?.[1]);
            const toCol = Number(info.to?.[0]);
            const toRow = Number(info.to?.[1]);
            if (!Number.isFinite(fromCol) || !Number.isFinite(fromRow))
                return;
            if (!Number.isFinite(toCol) || !Number.isFinite(toRow))
                return;
            if (fromCol === toCol && fromRow === toRow)
                return;
            this.pushHistoryEntry({
                op: 'move',
                id,
                from: [fromCol, fromRow],
                to: [toCol, toRow]
            });
            this.syncLayoutHelixPosFromEditor();
        };
        return this.scadnanoGridEditor;
    }
    initScadnanoGridPaneControls() {
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
            gridCanvas.addEventListener('keydown', (e) => {
                const cmd = e.ctrlKey || e.metaKey;
                if (!cmd)
                    return;
                if (e.key.toLowerCase() === 'z') {
                    e.preventDefault();
                    if (e.shiftKey)
                        this.redoFromGridView();
                    else
                        this.undoFromGridView();
                }
                else if (e.key.toLowerCase() === 'y') {
                    e.preventDefault();
                    this.redoFromGridView();
                }
            });
        }
        this.refreshHistoryButtons();
        const resizeHandle = document.getElementById('scadnanoGridResizeHandle');
        let resizing = false;
        if (resizeHandle) {
            resizeHandle.addEventListener('mousedown', (e) => {
                resizing = true;
                document.body.classList.add('scadnano-grid-resizing');
                e.preventDefault();
            });
        }
        window.addEventListener('mousemove', (e) => {
            if (!resizing)
                return;
            this.setScadnanoPaneWidth(e.clientX);
            this.resizeScadnanoGridCanvas();
            if (this.scadnanoGridEditor && typeof this.scadnanoGridEditor.resize === 'function') {
                this.scadnanoGridEditor.resize();
            }
        });
        window.addEventListener('mouseup', () => {
            if (!resizing)
                return;
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
function registerScadnanoWindowApi() {
    window.exportScadnanoFromGridView = (helixPosInput) => {
        scadnanoManager.exportFromGridView(helixPosInput);
    };
    window.showScadnanoGridFromHelixPos = (helixPosInput, gridTypeInput) => {
        scadnanoManager.showGridFromHelixPos(helixPosInput, gridTypeInput);
    };
    window.hideScadnanoGridPane = () => {
        scadnanoManager.hideScadnanoGridPane();
    };
    window.scadnanoDialogExport = () => {
        scadnanoManager.handleDialogExport();
    };
    window.toggleGridDropdown = (checkboxElement) => {
        scadnanoManager.toggleGridDropdown(checkboxElement);
    };
    window.scadnanoSelectHelixFromNucleotide = (nucleotideInput) => {
        scadnanoManager.selectHelixFromNucleotide(nucleotideInput);
    };
    window.scadnanoGetHelices = () => scadnanoManager.getHelices();
    window.scadnanoGridUndo = () => scadnanoManager.undoFromGridView();
    window.scadnanoGridRedo = () => scadnanoManager.redoFromGridView();
    window.scadnanoGridGetHistory = () => scadnanoManager.getHistorySnapshot();
}
registerScadnanoWindowApi();
// Keep these named wrappers for inline HTML handlers and backwards compatibility.
function scadnanoDialogExport() {
    scadnanoManager.handleDialogExport();
}
function toggleGridDropdown(checkboxElement) {
    scadnanoManager.toggleGridDropdown(checkboxElement);
}
