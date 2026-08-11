"use strict";
/// <reference path="../typescript_definitions/index.d.ts" />
class ScadnanoExportManager {
    currentScadnanoHelices = null;
    currentScadnanoConnections = [];
    currentScadnanoLayout = null;
    // Original helices from calculateScadnanoHelices(), before any merges.
    // Used by recalculateGridFromScratch to ensure the pipeline starts from
    // the same initial state as the first export, making it truly idempotent.
    originalScadnanoHelices = null;
    // Cached findHelices() output. Populated lazily by calculateScadnanoHelices and
    // reused by the export pipeline so the partials + usedSides ledger computed by
    // generateHelix is available to applyAxisOverlapMerge without re-walking the strands.
    currentScadnanoPartials = null;
    currentScadnanoUsedSides = null;
    scadnanoGridEditor = null;
    scadnanoGridEditorType = null;
    suppressNodeSelectedCallback = false;
    // Edit-history journal. Lives only while the grid pane is open; cleared on export / close /
    // reload. Holds combine and move operations in the readable JSON form discussed.
    history = this.createEmptyHistory();
    createEmptyHistory() {
        return { entries: {}, order: [], cursor: 0, nextMove: 0, nextCombine: 0 };
    }
    // Helix ids that have been locked by the user. Persists for the lifetime of the grid pane
    // session; cleared when the pane is closed or the grid is reloaded.
    lockedHelices = new Set();
    // "Split mode" state. When splitHelicesFromGridView successfully hides everything outside a
    // single helix, the target helix id is captured here so the follow-up "Split from selected"
    // button knows which helix to break apart. null while not in split mode.
    splitTargetHelixId = null;
    // Toggle state for the "Focus on helix" button. Set when focus mode is on
    // (hides everything outside the selected nucleotides), null when off. Lets
    // the button behave as a toggle: first click focuses, second click restores
    // the full view.
    //
    // Multi-helix support: when the user selects nucleotides spanning N helices
    // and toggles focus on, this set holds all N ids so the toggle-off path
    // knows what was hidden. splitTargetHelixId, however, is only populated
    // when N === 1 — that's the only case "Split from selected" can act on.
    focusedHelixIds = null;
    static LOCKED_COLOR = 0x808080;
    static UNLOCKED_COLOR = 0xffd400;
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
    // Triggered by the "Recalculate Grid" button. Takes the existing
    // helices[][] array from currentScadnanoLayout, rebuilds the grid from
    // scratch via toscad.setGrid, and then reruns the full layout pipeline
    // (directionAlign2 → alignGridPrim → getAngles → filter binders →
    // anglecomb → anglecorr → runAxisOverlapMerge → reintroduce binders →
    // getAngles → anglecorr → calculateGlobalPositions → renumberHelicesGNN →
    // applyHelixRenumber) to a fingerprint fixed point via convergeLayout.
    // Rebuilding the grid ensures the layout reflects the helix membership
    // after any combines / splits the user has done; otherwise the old grid
    // marks (with their previous helix ids and offsets) would leak into the
    // new pipeline. Locked helix positions are tracked by nucleotide
    // membership so they survive renumbering.
    recalculateGridFromScratch() {
        if (!this.currentScadnanoLayout) {
            notify('Open the grid view first before recalculating.', 'warning');
            return;
        }
        // Save one nucleotide id per locked helix as an anchor. Nucleotide ids are stable
        // across any renumbering, so we can find the new helix index after the pipeline runs.
        const lockedNtAnchors = [];
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
        this.runScadnanoLongCalculation(() => {
            try {
                // Take the existing helices[][] array as the source of truth.
                // Shallow-clone helices so anglecomb's in-place splices don't
                // mutate currentScadnanoLayout.helices, keeping recalculate idempotent.
                const helices = this.currentScadnanoLayout.helices
                    .map(slot => slot.slice());
                // Pass the stored grid so convergeLayout's first iteration can
                // preserve grid marks from it. This makes recalculateGridFromScratch
                // behave the same as the loop's subsequent iterations.
                const storedGrid = this.currentScadnanoLayout.grid;
                const dummyBinders = [];
                const { helices: finalHelices, grid: finalGrid, helixPos, latticeType } = this.convergeLayout(helices, storedGrid, dummyBinders, gridType, wireframe);
                const { crossovers } = toscad.collectCrossovers(finalGrid);
                this.currentScadnanoConnections = this.buildScadnanoConnections(crossovers);
                this.currentScadnanoHelices = finalHelices;
                this.currentScadnanoLayout = {
                    latticeType,
                    nucleotideCount: this.currentScadnanoLayout.nucleotideCount,
                    helices: finalHelices,
                    grid: finalGrid,
                    helixPos,
                    wireframe
                };
            }
            catch (err) {
                failed = true;
                notify(`Recalculate grid failed: ${err}`, 'alert');
            }
        }, () => {
            if (failed || !this.currentScadnanoLayout)
                return;
            const newLayout = this.currentScadnanoLayout;
            // Remap locked helix ids through the new numbering using nucleotide anchors.
            const newLocked = new Set();
            lockedNtAnchors.forEach(ntId => {
                const newHelixId = toscad.findHelixID(ntId, newLayout.helices);
                if (newHelixId !== null)
                    newLocked.add(newHelixId);
            });
            const helixPos = this.cloneHelixPosMap(newLayout.helixPos);
            window.currentScadnanoHelixPos = helixPos;
            // Reload the editor. showGridFromHelixPos clears history and lockedHelices,
            // so we restore the remapped set immediately after.
            this.showGridFromHelixPos(helixPos, newLayout.latticeType);
            this.lockedHelices = newLocked;
            this.applyLockedColors();
            notify('Grid recalculated.', 'success');
        });
    }
    // Triggered by the "Lock/Unlock Helices" button in the grid view.
    // Toggles the locked state of every currently-selected helix:
    //   - Unlocked → grey fill, added to lockedHelices.
    //   - Already locked → blue (unlocked) fill restored, removed from lockedHelices.
    lockSelectedHelices() {
        const editor = this.scadnanoGridEditor;
        if (!editor || typeof editor.setNodeColor !== 'function')
            return;
        const ids = typeof editor.getSelectedHelixIds === 'function'
            ? editor.getSelectedHelixIds()
            : [];
        if (ids.length === 0)
            return;
        ids.forEach(id => {
            if (this.lockedHelices.has(id)) {
                this.lockedHelices.delete(id);
                editor.setNodeColor(id, ScadnanoExportManager.UNLOCKED_COLOR);
            }
            else {
                this.lockedHelices.add(id);
                editor.setNodeColor(id, ScadnanoExportManager.LOCKED_COLOR);
            }
        });
        // Keep the editor's own set in sync so its drag/move guards work.
        editor.lockedHelices = new Set(this.lockedHelices);
    }
    // Triggered by the "Split helices" button in the grid view.
    // Hides everything that is NOT currently selected in the 3D scene, while preserving the
    // visibility state of nucleotides that were already hidden (so previously-hidden bases
    // don't get flipped back on). The user's original selection is restored afterwards.
    //
    // Toggles the "Focus on helix" mode based on a Metro 4 switch checkbox's
    // checked state. Wired to onchange="scadnanoFocusOnHelixToggle(this)" in
    // index.html, exactly like the Arrows/Box/Fog/3' markers switches in the
    // View toolbar. When `checked === true`, focuses on the selected nucleotides
    // (hiding everything else — works for 1 OR multiple helices). When
    // `checked === false`, restores the full view. If focus fails (e.g.
    // nothing selected), flips the switch back off so the UI doesn't lie.
    focusOnHelixToggleFromGridView(checked) {
        const chk = document.getElementById('scadnanoGridSplitHelicesToggle');
        if (!checked) {
            // Unfocus path: restore everything, exit split mode, leave the
            // switch in its off state.
            if (this.focusedHelixIds !== null) {
                this.clearFocusMode();
                this.splitTargetHelixId = null;
                this.setSplitFromSelectedEnabled(false);
                notify('Restored full view.', 'success');
                render();
            }
            return;
        }
        // Focus path: bail if nothing's selected. Revert the switch so the UI
        // doesn't show "focus mode on" when nothing happened.
        if (selectedBases.size === 0) {
            notify('Select a helix (or any nucleotides) before focusing.', 'warning');
            if (chk)
                chk.checked = false;
            return;
        }
        // Collect every helix id represented in the current selection. This
        // only works when the grid pane is open and its layout has been
        // prepared — otherwise we still hide, but can't track which helices
        // are focused.
        const helixIds = new Set();
        if (this.scadnanoGridEditor && this.currentScadnanoLayout) {
            const grid = this.currentScadnanoLayout.grid;
            selectedBases.forEach(e => {
                const mark = grid.get(e.id);
                if (mark)
                    helixIds.add(mark.helixId);
            });
        }
        // Invert the selection so `selectedBases` now holds every element outside
        // the focused set — those are the ones we want to hide. This works
        // whether the user picked 1 helix or N: any non-selected nucleotide
        // is hidden, regardless of which helix it belongs to.
        invertSelection();
        const affectedSystems = new Set();
        selectedBases.forEach(e => {
            // Only hide currently-visible elements. Toggling an already-hidden nucleotide
            // would make it visible again, which we explicitly want to avoid.
            if (e.getInstanceParameter3('visibility').x !== 0) {
                e.toggleVisibility();
            }
            affectedSystems.add(e.getSystem());
        });
        affectedSystems.forEach(sys => sys.callUpdates(['instanceVisibility']));
        if (tmpSystems.length > 0) {
            tmpSystems.forEach(sys => sys.callUpdates(['instanceVisibility']));
        }
        // Restore the user's original selection.
        invertSelection();
        // "Split from selected" only works on a single helix — gate the button
        // and the splitTargetHelixId pointer on `helixIds.size === 1`. With 0
        // (e.g. selection had no grid marks) or 2+ helices, the button stays
        // disabled and the user gets a warning when they try to click it.
        const singleHelix = helixIds.size === 1;
        this.splitTargetHelixId = singleHelix ? helixIds.values().next().value : null;
        this.setSplitFromSelectedEnabled(singleHelix);
        this.focusedHelixIds = helixIds.size > 0 ? helixIds : null;
        if (this.focusedHelixIds !== null) {
            if (singleHelix) {
                notify(`Focus mode: helix ${this.splitTargetHelixId}. Toggle off "Focus on helix" to restore the view.`, 'success');
            }
            else {
                notify(`Focus mode: ${this.focusedHelixIds.size} helices visible. Toggle off "Focus on helix" to restore the view.`, 'success');
            }
        }
        render();
    }
    // Enable/disable the "Split from selected" button by id. Safe to call when the DOM isn't
    // ready — the lookup just no-ops.
    setSplitFromSelectedEnabled(enabled) {
        const btn = document.getElementById('scadnanoGridSplitFromSelectedBtn');
        if (btn)
            btn.disabled = !enabled;
    }
    // Clear focus-mode state: restore visibility, reset the focused-helix set,
    // and flip the Metro switch back to off. Used by any code path that exits
    // the focus session (unfocus, close pane, after a split, etc.).
    clearFocusMode() {
        if (this.focusedHelixIds !== null) {
            this.restoreAllVisibility();
            this.focusedHelixIds = null;
        }
        const chk = document.getElementById('scadnanoGridSplitHelicesToggle');
        if (chk)
            chk.checked = false;
    }
    // Make every nucleotide across all systems visible. Used after a split completes so the user
    // sees the full structure again. Iterates like toggleVisArbitrary: only toggle bases whose
    // visibility.x is currently 0.
    restoreAllVisibility() {
        const affected = new Set();
        systems.forEach(sys => {
            sys.strands.forEach(strand => {
                strand.forEach((mono) => {
                    if (mono.getInstanceParameter3('visibility').x === 0) {
                        mono.toggleVisibility();
                        affected.add(sys);
                    }
                });
            });
        });
        affected.forEach(sys => sys.callUpdates(['instanceVisibility']));
        if (tmpSystems.length > 0) {
            tmpSystems.forEach(sys => sys.callUpdates(['instanceVisibility']));
        }
    }
    // Triggered by the "Split from selected" button. Splits `splitTargetHelixId` into two helices
    // along the boundary defined by the currently-selected nucleotides: the selected bases move
    // into a brand-new helix, the rest of the original helix stays put. Cascades the split
    // through the grid editor, helixPos, connection list, and locked-helix colors, then restores
    // visibility so the user sees the whole scene again.
    splitFromSelectedGridView() {
        if (this.splitTargetHelixId === null) {
            notify('Click "Split helices" on a selected helix first.', 'warning');
            return;
        }
        // Defensive guard: the button should be disabled when more than one
        // helix is focused, but if the user manages to invoke this with a
        // multi-helix focus set, refuse and warn rather than silently picking
        // a wrong target helix.
        if (this.focusedHelixIds !== null && this.focusedHelixIds.size !== 1) {
            notify(`"Split from selected" requires exactly one focused helix — you have ${this.focusedHelixIds.size}. Focus on a single helix to enable splitting.`, 'warning');
            return;
        }
        const editor = this.scadnanoGridEditor;
        if (!editor || !this.currentScadnanoLayout) {
            notify('Open the scadnano grid view before splitting.', 'warning');
            return;
        }
        if (selectedBases.size === 0) {
            notify('Select the nucleotides to peel off first.', 'warning');
            return;
        }
        if (this.lockedHelices.has(this.splitTargetHelixId)) {
            notify(`Cannot split: helix ${this.splitTargetHelixId} is locked. Unlock it first.`, 'warning');
            return;
        }
        const helices = this.ensureScadnanoHelicesCache();
        if (!helices) {
            notify('Helix data is not yet available.', 'alert');
            return;
        }
        const targetHelixId = this.splitTargetHelixId;
        const grid = this.currentScadnanoLayout.grid;
        const helixPos = this.currentScadnanoLayout.helixPos;
        // Only keep selected bases that (a) are Nucleotides and (b) currently belong to the
        // target helix. Bases hidden or from other helices are silently ignored.
        const nucleotidesToMove = [];
        selectedBases.forEach(e => {
            if (!(e instanceof Nucleotide))
                return;
            const mark = grid.get(e.id);
            if (!mark || mark.helixId !== targetHelixId)
                return;
            nucleotidesToMove.push(e);
        });
        if (nucleotidesToMove.length === 0) {
            notify(`None of the selected nucleotides belong to helix ${targetHelixId}.`, 'warning');
            return;
        }
        const result = toscad.splitHelix(grid, helixPos, targetHelixId, helices, nucleotidesToMove);
        if (!result) {
            notify(`Split failed for helix ${targetHelixId} (check console).`, 'alert');
            return;
        }
        const { keptHelixId, newHelixId, helixApos, helixBpos } = result;
        // Cascade the split into the grid editor. Existing node ids/positions are unchanged;
        // just add a new node for the freshly-created helix at helixBpos.
        if (typeof editor.addNode === 'function') {
            editor.addNode({ id: newHelixId, col: helixBpos[0], row: helixBpos[1], label: String(newHelixId) });
        }
        // The two halves used to share physical crossovers, so add a visual connection between
        // them. Dedupes against the existing list.
        const connKey = keptHelixId < newHelixId
            ? `${keptHelixId}:${newHelixId}`
            : `${newHelixId}:${keptHelixId}`;
        const existing = new Set(this.currentScadnanoConnections.map(([a, b]) => a < b ? `${a}:${b}` : `${b}:${a}`));
        if (!existing.has(connKey)) {
            this.currentScadnanoConnections.push([keptHelixId, newHelixId]);
        }
        if (typeof editor.setConnections === 'function') {
            editor.setConnections(this.currentScadnanoConnections);
        }
        // Refresh downstream caches from the editor (adds newHelixId to helixPos).
        this.syncLayoutHelixPosFromEditor();
        this.publishCurrentHelixPosFromEditor();
        // Keep lock colors correct — the new node inherits the default (unlocked) color, but
        // applyLockedColors is idempotent so this is a cheap way to guarantee consistency.
        this.applyLockedColors();
        if (typeof editor.clearSelection === 'function')
            editor.clearSelection();
        // Bring every hidden nucleotide back into view so the user sees the whole structure.
        this.restoreAllVisibility();
        // Exit split mode.
        this.splitTargetHelixId = null;
        this.setSplitFromSelectedEnabled(false);
        notify(`Split helix ${keptHelixId}: moved ${nucleotidesToMove.length} nt into new helix ${newHelixId} ` +
            `at [${helixBpos[0]},${helixBpos[1]}] (kept helix ${keptHelixId} at [${helixApos[0]},${helixApos[1]}]).`, 'success');
        render();
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
        this.currentScadnanoLayout.grid.forEach((mark) => {
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
        // Remap lockedHelices through the same id transformation so locks follow their nodes.
        this.remapLockedHelicesForward(remapId);
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
        // Restore the locked set back to pre-merge ids and reapply colors.
        this.remapLockedHelicesInverse(inverseRemap, entry);
    }
    // Remap lockedHelices forward after a combine: each old id is mapped through remapId to its
    // new post-merge id. If a merged-away helix was locked, its lock carries over to the kept
    // helix (keptIdx). Reapplies grey fill to all still-locked nodes.
    remapLockedHelicesForward(remapId) {
        if (this.lockedHelices.size === 0)
            return;
        const next = new Set();
        this.lockedHelices.forEach(id => {
            const newId = remapId(id);
            if (newId !== null)
                next.add(newId);
        });
        this.lockedHelices = next;
        this.applyLockedColors();
    }
    // Remap lockedHelices back to pre-merge ids after an undo-combine. inverseRemap maps the
    // current post-merge id → the original pre-merge id. Reapplies grey fill to locked nodes
    // and restores the default blue on the re-split nodes that are no longer locked.
    remapLockedHelicesInverse(inverseRemap, entry) {
        if (this.lockedHelices.size === 0 && entry.removed.length === 0)
            return;
        const next = new Set();
        this.lockedHelices.forEach(id => {
            next.add(inverseRemap(id));
        });
        this.lockedHelices = next;
        this.applyLockedColors();
    }
    // Reapply the locked (grey) color to every node currently in lockedHelices, and restore
    // the default blue to every other node. Called after any operation that rebuilds the
    // editor node set (combine, undo, redo).
    applyLockedColors() {
        const editor = this.scadnanoGridEditor;
        if (!editor || typeof editor.setNodeColor !== 'function')
            return;
        const nodes = typeof editor.getNodes === 'function'
            ? editor.getNodes()
            : [];
        nodes.forEach(node => {
            const id = Number(node.id);
            editor.setNodeColor(id, this.lockedHelices.has(id)
                ? ScadnanoExportManager.LOCKED_COLOR
                : ScadnanoExportManager.UNLOCKED_COLOR);
        });
        // Keep the editor's drag/move guard in sync.
        editor.lockedHelices = new Set(this.lockedHelices);
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
        // Fresh grid view = fresh history and no locks.
        this.lockedHelices.clear();
        if (this.scadnanoGridEditor)
            this.scadnanoGridEditor.lockedHelices = new Set();
        this.clearHistory();
        this.splitTargetHelixId = null;
        this.setSplitFromSelectedEnabled(false);
    }
    hideScadnanoGridPane() {
        document.body.classList.remove('scadnano-grid-open');
        this.lockedHelices.clear();
        if (this.scadnanoGridEditor)
            this.scadnanoGridEditor.lockedHelices = new Set();
        this.clearHistory();
        this.splitTargetHelixId = null;
        this.setSplitFromSelectedEnabled(false);
        this.clearFocusMode();
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
    // Exports the helices[][] as a JSON file containing just the nucleotide
    // ids grouped by helix (e.g. [[1, 2, 3], [4, 5], ...]). Intentionally does
    // NOT include any grid positions, lattice info, or the grid map — only
    // the helix grouping. Reuses the existing cache when possible; otherwise
    // runs the cheap helix-detection path only (no full layout pipeline,
    // since grid positions are deliberately excluded from the output).
    //
    // Selection behaviour: if the grid view is open AND the user has at least
    // one helix (or any nucleotide within a helix) selected, only the
    // selected helices are exported — preserving the helices[][] shape. With
    // no selection, every helix is exported.
    exportHelicesAsJson(name) {
        let helices = this.ensureScadnanoHelicesCache();
        if (!helices || helices.length === 0) {
            try {
                helices = this.calculateScadnanoHelices();
                this.currentScadnanoHelices = helices;
            }
            catch (err) {
                notify(`Helices export failed: ${err}`, 'alert');
                return;
            }
        }
        if (!helices || helices.length === 0) {
            notify('No helices available to export.', 'warning');
            return;
        }
        // If a selection exists, restrict the export to the selected helices.
        // The grid editor's getSelectedHelixIds() returns the renumbered helix
        // indices of every selected node (each grid node represents one helix,
        // and its `id` is the renumbered helix index — see
        // scadnano_gridview.loadFromHelixPos, which sets node.id = helixId
        // from the helixPos map). Those indices match the helices[][] array
        // indices produced by applyHelixRenumber, so we filter by INDEX here,
        // not by nucleotide id.
        let helicesToExport = helices;
        const editor = this.scadnanoGridEditor;
        if (editor && typeof editor.getSelectedHelixIds === 'function') {
            const selectedIds = editor.getSelectedHelixIds();
            if (Array.isArray(selectedIds) && selectedIds.length > 0) {
                const selectedSet = new Set(selectedIds);
                helicesToExport = helices.filter((_, idx) => selectedSet.has(idx));
                if (helicesToExport.length === 0) {
                    notify('Selected helices are no longer present in the layout.', 'warning');
                    return;
                }
            }
        }
        const idsOnly = helicesToExport.map(helix => helix.map(n => n.id));
        const fileName = name && name.trim() ? `${name.trim()}.json` : 'helices.json';
        makeTextFile(fileName, JSON.stringify(idsOnly, null, 2));
    }
    selectHelixFromNucleotide(nucleotideInput, additive = false) {
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
            this.scadnanoGridEditor.selectNodeById(helixId, additive);
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
        this.currentScadnanoPartials = result?.partials ?? null;
        this.currentScadnanoUsedSides = result?.usedSides ?? null;
        this.notifyHelixCoverageMismatch(helices, nucleotideElements);
        // Store a deep copy as the original (pre-merge) helices for recalculateGridFromScratch.
        this.originalScadnanoHelices = helices.map(h => h.slice());
        return helices;
    }
    // Deterministic fingerprint of every mutable field the pipeline touches on
    // the grid, plus the current helix count. `convergeLayout` compares
    // fingerprints across iterations to decide when the pipeline has reached a
    // fixed point — this is more reliable than checking return counts from
    // individual stages because it catches ALL mutations (directionAlign2 flips,
    // alignGridPrim offset shifts, anglecomb/anglecorr/axisOverlap changes)
    // in a single signal.
    gridFingerprint(grid, helicesLength) {
        const ntIds = Array.from(grid.keys()).sort((a, b) => a - b);
        const parts = [`h=${helicesLength}`];
        for (const ntId of ntIds) {
            const m = grid.get(ntId);
            parts.push(`${ntId}:${m.helixId}:${m.offset}:${m.direction === 'forward' ? 'f' : 'b'}`);
        }
        return parts.join('|');
    }
    convergeLayout(helicesIn, grid, binderHelices, requestedLatticeType, wireframe) {
        let helices = helicesIn;
        let latticeType;
        let networkMap;
        let helixPos;
        // ── Wireframe: EXACTLY one pass, no loop ─────────────────────────────
        if (wireframe) {
            toscad.directionAlign2(grid);
            toscad.alignGridPrim(grid, binderHelices);
            if (requestedLatticeType === 'automatic') {
                latticeType = toscad.detectLatticeKind(grid, binderHelices);
                notify(`Auto-detected lattice: ${latticeType}`, 'success');
            }
            else {
                latticeType = requestedLatticeType;
            }
            networkMap = toscad.getAngles(grid, helices, latticeType);
            helixPos = toscad.calculateGlobalPositions(networkMap, undefined, undefined, latticeType);
            const renumber = toscad.renumberHelicesGNN(grid, helixPos, latticeType, binderHelices);
            const renumbered = toscad.applyHelixRenumber(helices, grid, helixPos, renumber.remap);
            helices = renumbered.helices;
            helixPos = renumbered.helixPos;
            console.log(`[scadnano] convergeLayout (wireframe) — single pass, no iteration`);
            return { helices, grid, helixPos, latticeType, networkMap };
        }
        // ── Non-wireframe: fingerprint fixed-point loop ─────────────────────
        // Each iteration rebuilds the grid from scratch via setGrid on the
        // previous iteration's helices, then runs the full pipeline. This
        // ensures that the grid is always consistent with the helices array,
        // even after merges/renumbering. The lattice type is detected once
        // on the first iteration and frozen for subsequent iterations.
        // 
        // After the first iteration, we pass the previous grid to setGrid as
        // `preserveGrid`. This carries forward grid marks for nucleotides in
        // helices that haven't been merged, while only reassigning positions
        // for nucleotides in newly-merged or altered helices.
        let latticeTypeSet = null;
        networkMap = new Map();
        helixPos = new Map();
        let prevFp = '';
        // Initialize prevGrid from the input grid parameter. This allows
        // callers (like recalculateGridFromScratch) to pass a stored grid
        // so the first iteration can preserve marks from it.
        let prevGrid = grid;
        // Initialize prevHelicesSnapshot from the input helices. This allows
        // the first iteration to detect which helix slots are "stable" vs
        // merged/altered compared to the previous state.
        let prevHelicesSnapshot = helices.map(h => h.slice());
        // Merge-provenance products from the previous iteration. Each entry is
        // one merged helix, expressed as a list of origin-groups (nt-id arrays,
        // one per pre-merge helix). Rebuilt at the end of every iteration from
        // that iteration's explicit merge log (anglecomb2 + axis-overlap) and
        // mapped onto the current helix slots at the start of the next.
        let prevMergeProducts = [];
        // Iteration cap for the fixed-point loop.
        const MAX_ITER = 7;
        for (let iter = 1; iter <= MAX_ITER; iter++) {
            // Merge events recorded during this pass, folded into provenance
            // products at the end of the pass.
            const iterMerges = [];
            // Rebuild grid from scratch on the current helices array.
            // After the first iteration, preserve grid marks from the previous
            // iteration for nucleotides in helices that haven't been merged.
            // We detect merges by comparing nucleotide IDs: if a nucleotide
            // was in helix[i] before and is still in helix[i] now, preserve it.
            // If nucleotides moved (e.g. helix was merged into another), don't
            // preserve — let setGrid reassign their grid positions.
            let preservedNtIds;
            if (prevGrid && prevHelicesSnapshot) {
                preservedNtIds = new Set();
                const prevNtsBySlot = prevHelicesSnapshot.map(h => new Set(h.map(nt => nt.id)));
                for (let i = 0; i < helices.length; i++) {
                    const currNts = helices[i];
                    if (!currNts)
                        continue;
                    // A helix is "stable" if its nucleotide set exactly matches
                    // a helix from the previous iteration. This means no merge
                    // touched it.
                    for (let j = 0; j < prevNtsBySlot.length; j++) {
                        const prevSlot = prevNtsBySlot[j];
                        if (currNts.length !== prevSlot.size)
                            continue;
                        let match = true;
                        for (const nt of currNts) {
                            if (!prevSlot.has(nt.id)) {
                                match = false;
                                break;
                            }
                        }
                        if (match) {
                            // This helix slot is identical to a previous slot.
                            // Preserve all its nucleotides' grid marks.
                            for (const nt of currNts)
                                preservedNtIds.add(nt.id);
                            break;
                        }
                    }
                }
            }
            // Map the previous iteration's merge products onto the current
            // helix slots. Nucleotide ids are stable across renumbering, so a
            // product is located via any of its member nucleotides (no splits
            // exist in the pipeline, so all members land in one slot).
            let mergedGroups;
            if (prevGrid && prevMergeProducts.length > 0) {
                const productByNt = new Map();
                for (const groups of prevMergeProducts) {
                    for (const group of groups) {
                        for (const ntId of group)
                            productByNt.set(ntId, groups);
                    }
                }
                mergedGroups = new Map();
                for (let i = 0; i < helices.length; i++) {
                    const slot = helices[i];
                    if (!slot || slot.length === 0)
                        continue;
                    const groups = productByNt.get(slot[0].id);
                    if (groups)
                        mergedGroups.set(i, groups);
                }
                if (mergedGroups.size === 0)
                    mergedGroups = undefined;
            }
            const { grid: freshGrid, binderHelices: freshBinders } = toscad.setGrid(helices, prevGrid, preservedNtIds, mergedGroups);
            grid = freshGrid;
            binderHelices = freshBinders ?? [];
            // setGrid only places nucleotides; a merged helix's origin-groups
            // land at their collision-free pre-merge offsets, unaligned. Align
            // them to the lattice here, BEFORE directionAlign2 — that pass
            // flips helices and rewrites offsets, which would change what the
            // shift observations see. No-op when nothing was merged.
            toscad.alignMergedGroups(grid, mergedGroups, binderHelices);
            toscad.directionAlign2(grid);
            toscad.alignGridPrim(grid, binderHelices);
            // Resolve 'automatic' once, after the first alignment (detection
            // reads helixId / offset / direction off the aligned grid). Freeze
            // the choice for the rest of the loop.
            if (latticeTypeSet === null) {
                if (requestedLatticeType === 'automatic') {
                    latticeTypeSet = toscad.detectLatticeKind(grid, binderHelices);
                    notify(`Auto-detected lattice: ${latticeTypeSet}`, 'success');
                }
                else {
                    latticeTypeSet = requestedLatticeType;
                }
            }
            // ── Phase 1: main pipeline WITHOUT binder helices ────────────────
            // Binder helices can interfere with angle resolution for the main
            // structure. We track binder helices by their nucleotide IDs (which
            // are immutable across merges/renumbering) so that the filter
            // stays correct even after anglecomb2's mergeHelixInto remaps
            // helixIds and splices the helices array.
            networkMap = toscad.getAngles(grid, helices, latticeTypeSet);
            // Snapshot binder nucleotide IDs ONCE. These survive any merge/
            // renumber because nucleotide ids are assigned at creation and
            // never change — only helixId indices shift.
            const binderNtIds = new Set();
            for (const bHid of binderHelices) {
                const slot = helices[bHid];
                if (Array.isArray(slot)) {
                    for (const nt of slot)
                        binderNtIds.add(nt.id);
                }
            }
            // Re-derive the set of binder helixIds from the grid, using
            // the immutable nucleotide IDs as the anchor. This stays
            // correct after any helixId remapping.
            const deriveBinderHelixSet = () => {
                const s = new Set();
                for (const [ntId, m] of grid.entries()) {
                    if (binderNtIds.has(ntId))
                        s.add(m.helixId);
                }
                return s;
            };
            const filterBinderEntries = (nm) => {
                const bSet = deriveBinderHelixSet();
                const out = new Map();
                for (const [hid, neighbors] of nm.entries()) {
                    if (bSet.has(hid))
                        continue;
                    const filteredNeighbors = new Map();
                    for (const [nid, angle] of neighbors.entries()) {
                        if (bSet.has(nid))
                            continue;
                        filteredNeighbors.set(nid, angle);
                    }
                    out.set(hid, filteredNeighbors);
                }
                return out;
            };
            const filteredNetworkMap = filterBinderEntries(networkMap);
            // ── anglecomb3 replaces anglecomb2 ─────────────────────────────
            // New rule: merge iff hashAxisOverlap picked the pair AND their
            // crossover-count cosine similarity > 0.5. Precompute the hash
            // pair list here (same inputs runAxisOverlapMerge would use) and
            // hand it to anglecomb3, which delegates the actual merges to
            // applyAxisOverlapMerge (identical mechanics to runAxisOverlapMerge).
            //
            // Note: this REPLACES anglecomb2's collision-bucket signal. Pairs
            // that anglecomb2 would have merged via same-cell colocation but
            // that hashAxisOverlap didn't nominate will no longer merge here.
            let hashMergePairs = [];
            {
                const partials = this.currentScadnanoPartials;
                const usedSides = this.currentScadnanoUsedSides;
                if (Array.isArray(partials) && usedSides && partials.length > 0) {
                    const partialEnds = helix.mapPartialEnds(partials);
                    if (partialEnds.size > 0) {
                        const partialAxes = helix.partialAxesTowardFreeSide(partials, partialEnds, usedSides);
                        if (partialAxes.size > 0) {
                            hashMergePairs = helix.hashAxisOverlap(partials, partialEnds, usedSides, partialAxes);
                        }
                    }
                }
            }
            const combResult = toscad.anglecomb3(grid, helices, latticeTypeSet, this.currentScadnanoPartials ?? [], hashMergePairs, filteredNetworkMap);
            for (const mp of combResult.mergedPairs) {
                iterMerges.push({ keepNtIds: mp.keepNtIds, mergedNtIds: mp.mergedNtIds });
            }
            // anglecomb3 internally re-derives networkMap via getAngles after
            // each merge, which would include the (now-remapped) binder
            // helixIds. Re-filter to keep binders out for the subsequent
            // anglecorr2 call.
            networkMap = filterBinderEntries(combResult.networkMap);
            const corrResult = toscad.anglecorr2(grid, helices, latticeTypeSet, networkMap);
            networkMap = corrResult.networkMap;
            // ── Phase 2: reintroduce binder helices ──────────────────────────
            // Grid and helices array have stayed self-consistent through the
            // merges above (binder helixIds shifted along with the rest, and
            // our nucleotide-ID tracking knows which helices they are now).
            // Run getAngles on the full grid to bring binders back into the
            // network map, then anglecorr2 to resolve their angles.
            networkMap = toscad.getAngles(grid, helices, latticeTypeSet);
            const corrResult2 = toscad.anglecorr2(grid, helices, latticeTypeSet, networkMap);
            networkMap = corrResult2.networkMap;
            helixPos = toscad.calculateGlobalPositions(networkMap, undefined, undefined, latticeTypeSet);
            // Renumber inside the loop. Each pass' renumber rewrites helix ids
            // on the grid; the next pass' directionAlign2 / alignGridPrim then
            // anchor on the newly-designated helix 0. Once the numbering
            // matches the GNN canonical order and mutations have settled, the
            // remap is identity and the fingerprint stops changing.
            //
            // Pass the CURRENT binder helix set (derived from immutable
            // nucleotide ids, so it's post-merge/post-remap correct) so the
            // renumber's anchor selection never assigns helix 0 to a binder.
            // Helix 0 seeds calculateGlobalPositions and directionAlign2 on
            // the next iteration; a binder anchor there misroots the whole
            // lattice.
            const currentBinderHelices = Array.from(deriveBinderHelixSet());
            const renumber = toscad.renumberHelicesGNN(grid, helixPos, latticeTypeSet, currentBinderHelices);
            const renumbered = toscad.applyHelixRenumber(helices, grid, helixPos, renumber.remap);
            helices = renumbered.helices;
            helixPos = renumbered.helixPos;
            const fp = this.gridFingerprint(grid, helices.length);
            if (fp === prevFp) {
                console.log(`[scadnano] convergeLayout converged in ${iter} pass${iter === 1 ? '' : 'es'}`);
                return { helices, grid, helixPos, latticeType: latticeTypeSet, networkMap };
            }
            prevFp = fp;
            // Save the current grid for preservation in the next iteration.
            // Helices that survived this iteration's merges keep their grid marks;
            // only nucleotides in newly-merged/altered helices get reassigned.
            prevGrid = grid;
            // Snapshot helices so the next iteration can detect which slots
            // are stable (identical nt sets) vs merged/altered.
            prevHelicesSnapshot = helices.map(h => h.slice());
            // Fold this iteration's merge events into provenance products so
            // the next iteration's setGrid can preserve each origin-group's
            // internal marks and re-align the groups, instead of re-walking
            // the merged helix as one disconnected component.
            prevMergeProducts = this.buildMergeProducts(iterMerges);
        }
        console.warn(`[scadnano] convergeLayout hit iteration cap ${MAX_ITER}; using last state.`);
        return { helices, grid, helixPos, latticeType: latticeTypeSet ?? 'honeycomb', networkMap };
    }
    // Fold an iteration's ordered merge events into "products": for every helix
    // that received merges, the list of origin-groups (nt-id arrays) — one
    // group per pre-merge helix. Grouping is one level deep: a helix that was
    // itself merged in a previous iteration contributes all of its nucleotides
    // as a single group.
    buildMergeProducts(mergeLog) {
        const products = [];
        const productIndexByNt = new Map();
        const retarget = (fromIdx, toIdx) => {
            if (fromIdx === toIdx)
                return;
            for (const group of products[fromIdx]) {
                for (const ntId of group)
                    productIndexByNt.set(ntId, toIdx);
            }
            products[fromIdx] = [];
        };
        for (const { keepNtIds, mergedNtIds } of mergeLog) {
            if (!keepNtIds.length || !mergedNtIds.length)
                continue;
            const keepIdx = productIndexByNt.get(keepNtIds[0]);
            const mergedIdx = productIndexByNt.get(mergedNtIds[0]);
            if (keepIdx !== undefined && keepIdx === mergedIdx)
                continue;
            const keepGroups = keepIdx !== undefined ? products[keepIdx] : [keepNtIds.slice()];
            const mergedGroups = mergedIdx !== undefined ? products[mergedIdx] : [mergedNtIds.slice()];
            let target;
            if (keepIdx !== undefined) {
                target = keepIdx;
                products[target] = [...keepGroups, ...mergedGroups];
            }
            else {
                target = products.length;
                products.push([...keepGroups, ...mergedGroups]);
                for (const ntId of keepNtIds)
                    productIndexByNt.set(ntId, target);
            }
            if (mergedIdx !== undefined) {
                retarget(mergedIdx, target);
            }
            else {
                for (const ntId of mergedNtIds)
                    productIndexByNt.set(ntId, target);
            }
        }
        return products.filter(p => p.length > 1);
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
        // Iterate the whole pipeline — directionAlign2 → alignGridPrim →
        // getAngles → (filter out binders) → anglecomb → anglecorr →
        // runAxisOverlapMerge → (reintroduce binders) → getAngles → anglecorr →
        // calculateGlobalPositions → renumberHelicesGNN → applyHelixRenumber
        // — to a fingerprint fixed point. Binder helices are excluded during
        // the main angle-resolution pipeline so they don't interfere, then
        // reintroduced for a final anglecorr pass before global positioning.
        // Renumber is INSIDE the loop because it rewrites helix ids on the grid,
        // which changes the anchor for the next iteration's directionAlign2 /
        // alignGridPrim; see convergeLayout.
        const { helices: finalHelices, grid: finalGrid, helixPos, latticeType } = this.convergeLayout(helices, grid, binderHelices, requestedLatticeType, wireframe);
        this.currentScadnanoHelices = finalHelices;
        // Build connections AFTER convergence so they reference the final IDs.
        const { crossovers } = toscad.collectCrossovers(finalGrid);
        this.currentScadnanoConnections = this.buildScadnanoConnections(crossovers);
        this.currentScadnanoLayout = {
            latticeType,
            nucleotideCount,
            helices: finalHelices,
            grid: finalGrid,
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
        // Wire the top toolbar (Edit/Export tabs → content sections). The wiring lives in
        // scadnano_gridview so it can be exported/imported cleanly. Safe to call multiple
        // times — listeners are deduped inside.
        // Guard against `scadnano` being undeclared: this constructor runs at script load
        // time, before scadnano_gridview.js has declared its namespace. Without the
        // `typeof scadnano !== 'undefined'` check, dereferencing `scadnano.initGridToolbar`
        // throws a ReferenceError, which aborts the top-level `const scadnanoManager = ...`
        // initialization and leaves it in permanent TDZ. Every subsequent inline onclick
        // (e.g. scadnanoDialogExport) then fails with "Cannot access 'scadnanoManager'
        // before initialization".
        if (typeof scadnano !== 'undefined' && typeof scadnano.initGridToolbar === 'function') {
            scadnano.initGridToolbar();
        }
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
        const exportHelicesBtn = document.getElementById('scadnanoGridExportHelicesBtn');
        if (exportHelicesBtn) {
            exportHelicesBtn.addEventListener('click', () => {
                this.exportHelicesAsJson();
            });
        }
        const combineBtn = document.getElementById('scadnanoGridCombineBtn');
        if (combineBtn) {
            combineBtn.addEventListener('click', () => {
                this.combineSelectedHelicesFromGridView();
            });
        }
        // The "Focus on helix" toggle is a Metro 4 data-role="switch" checkbox
        // (see index.html). It binds directly via inline onchange in the HTML,
        // so no addEventListener wiring is needed here.
        const splitFromSelectedBtn = document.getElementById('scadnanoGridSplitFromSelectedBtn');
        if (splitFromSelectedBtn) {
            splitFromSelectedBtn.addEventListener('click', () => {
                this.splitFromSelectedGridView();
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
    window.scadnanoSelectHelixFromNucleotide = (nucleotideInput, additive) => {
        scadnanoManager.selectHelixFromNucleotide(nucleotideInput, additive === true);
    };
    window.scadnanoGetHelices = () => scadnanoManager.getHelices();
    window.scadnanoExportHelices = (name) => {
        scadnanoManager.exportHelicesAsJson(name);
    };
    window.scadnanoFocusOnHelixToggle = (chkBox) => {
        scadnanoManager.focusOnHelixToggleFromGridView(Boolean(chkBox?.checked));
    };
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
