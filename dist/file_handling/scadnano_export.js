/// <reference path="../typescript_definitions/index.d.ts" />
async function makeScadnanoJsonFile(name, gridType = 'square') {
    try {
        const nucleotideElements = new Map();
        elements.forEach((element, id) => {
            if (element instanceof Nucleotide) {
                nucleotideElements.set(id, element);
            }
        });
        const helices = await honda.findHelices(nucleotideElements, 3);
        const { grid, binderHelices } = toscad.setGrid(helices);
        toscad.directionAlign2(grid);
        toscad.alignGridPrim(grid, binderHelices);
        toscad.combinedHelices(15, grid, helices, binderHelices);
        const scadnano = toscad.buildScadnano2(grid, helices, gridType);
        const fileName = name ? `${name}.sc` : "output.sc";
        makeTextFile(fileName, JSON.stringify(scadnano, null, 2));
    }
    catch (err) {
        notify(`Scadnano export failed: ${err}`, "alert");
    }
}
async function makeScadnanowHexPos(name, gridType = 'honeycomb') {
    try {
        const nucleotideElements = new Map();
        elements.forEach((element, id) => {
            if (element instanceof Nucleotide) {
                nucleotideElements.set(id, element);
            }
        });
        const helices = await honda.findHelices(nucleotideElements, 3);
        const { grid, binderHelices } = toscad.setGrid(helices);
        toscad.directionAlign2(grid);
        toscad.alignGridPrim(grid, binderHelices);
        toscad.combinedHelices(15, grid, helices, binderHelices);
        const helixpos = toscad.HelixPos(grid, helices);
        const scadnano = toscad.buildScadnano2(grid, helices, gridType, helixpos);
        const fileName = name ? `${name}.sc` : "output.sc";
        makeTextFile(fileName, JSON.stringify(scadnano, null, 2));
    }
    catch (err) {
        notify(`Scadnano export failed: ${err}`, "alert");
    }
}
async function calculateScadnanoHelixPos() {
    const nucleotideElements = new Map();
    elements.forEach((element, id) => {
        if (element instanceof Nucleotide) {
            nucleotideElements.set(id, element);
        }
    });
    const helices = await honda.findHelices(nucleotideElements, 3);
    const { grid, binderHelices } = toscad.setGrid(helices);
    toscad.directionAlign2(grid);
    toscad.alignGridPrim(grid, binderHelices);
    toscad.combinedHelices(15, grid, helices, binderHelices);
    return toscad.HelixPos(grid, helices);
}
function normalizeHelixPosMap(input) {
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
async function exportScadnanoWithHelixPos(name, gridType, helixPos) {
    const nucleotideElements = new Map();
    elements.forEach((element, id) => {
        if (element instanceof Nucleotide) {
            nucleotideElements.set(id, element);
        }
    });
    const helices = await honda.findHelices(nucleotideElements, 3);
    const { grid, binderHelices } = toscad.setGrid(helices);
    toscad.directionAlign2(grid);
    toscad.alignGridPrim(grid, binderHelices);
    toscad.combinedHelices(15, grid, helices, binderHelices);
    const scadnano = toscad.buildScadnano2(grid, helices, gridType, helixPos);
    const fileName = name ? `${name}.sc` : 'output.sc';
    makeTextFile(fileName, JSON.stringify(scadnano, null, 2));
}
async function scadnanoDialogExport() {
    const nameInput = document.getElementById('scadnanoFilename');
    const helixPosCheckbox = document.getElementById('scadnanoIncludeHPos');
    const scadnanoGrid = document.getElementById('scadnanoGrid');
    if (!nameInput || !helixPosCheckbox || !scadnanoGrid) {
        console.warn('scadnano export dialog missing inputs');
        return;
    }
    try {
        const helixPos = await calculateScadnanoHelixPos();
        window.currentScadnanoHelixPos = helixPos;
        const showGrid = window.showScadnanoGridFromHelixPos;
        if (typeof showGrid !== 'function') {
            notify('Scadnano grid view is unavailable in this page.', 'alert');
            return;
        }
        showGrid(helixPos);
    }
    catch (err) {
        notify(`Scadnano export failed: ${err}`, 'alert');
    }
}
window.exportScadnanoFromGridView = async function (helixPosInput) {
    const nameInput = document.getElementById('scadnanoFilename');
    const scadnanoGrid = document.getElementById('scadnanoGrid');
    const name = nameInput?.value.trim() || 'output';
    const gridType = scadnanoGrid?.value === 'honeycomb' ? 'honeycomb' : 'square';
    const map = normalizeHelixPosMap(helixPosInput ?? window.currentScadnanoHelixPos);
    if (!map || map.size === 0) {
        notify('No edited helix positions available to export.', 'warning');
        return;
    }
    try {
        await exportScadnanoWithHelixPos(name, gridType, map);
    }
    catch (err) {
        notify(`Scadnano export failed: ${err}`, 'alert');
    }
};
let scadnanoGridEditor = null;
function getScadnanoGridPane() {
    return document.getElementById('scadnanoGridPane');
}
function getScadnanoGridCanvas() {
    return document.getElementById('scadnanoGridCanvas');
}
function resizeScadnanoGridCanvas() {
    const pane = getScadnanoGridPane();
    const canvas = getScadnanoGridCanvas();
    if (!pane || !canvas)
        return;
    canvas.width = pane.clientWidth;
    canvas.height = pane.clientHeight;
}
function mapFromEditorNodes(editor) {
    const out = new Map();
    const nodes = typeof editor.getNodes === 'function' ? editor.getNodes() : [];
    nodes.forEach((node) => {
        out.set(Number(node.id), [Number(node.col), Number(node.row)]);
    });
    return out;
}
function publishCurrentHelixPosFromEditor() {
    if (!scadnanoGridEditor)
        return;
    window.currentScadnanoHelixPos = mapFromEditorNodes(scadnanoGridEditor);
}
function ensureScadnanoGridEditor() {
    if (scadnanoGridEditor)
        return scadnanoGridEditor;
    const canvas = getScadnanoGridCanvas();
    if (!canvas)
        return null;
    const scadnanoNs = window.scadnano;
    if (!scadnanoNs || typeof scadnanoNs.HoneycombEditor !== 'function')
        return null;
    scadnanoGridEditor = new scadnanoNs.HoneycombEditor(canvas);
    scadnanoGridEditor.onNodesChanged = publishCurrentHelixPosFromEditor;
    return scadnanoGridEditor;
}
window.showScadnanoGridFromHelixPos = function (helixPosInput) {
    const pane = getScadnanoGridPane();
    if (!pane) {
        notify('Scadnano grid pane is unavailable.', 'alert');
        return;
    }
    document.body.classList.add('scadnano-grid-open');
    resizeScadnanoGridCanvas();
    const editor = ensureScadnanoGridEditor();
    if (!editor) {
        notify('Unable to open scadnano grid view.', 'alert');
        return;
    }
    if (typeof editor.resize === 'function') {
        editor.resize();
    }
    const map = normalizeHelixPosMap(helixPosInput);
    if (!map || map.size === 0) {
        notify('No helix positions available for the grid view.', 'warning');
        return;
    }
    editor.loadFromHelixPos(map);
    publishCurrentHelixPosFromEditor();
};
window.hideScadnanoGridPane = function () {
    document.body.classList.remove('scadnano-grid-open');
};
function initScadnanoGridPaneControls() {
    const exportBtn = document.getElementById('scadnanoGridExportBtn');
    if (exportBtn) {
        exportBtn.addEventListener('click', () => {
            publishCurrentHelixPosFromEditor();
            const doExport = window.exportScadnanoFromGridView;
            if (typeof doExport !== 'function') {
                notify('Scadnano export is unavailable.', 'alert');
                return;
            }
            doExport(window.currentScadnanoHelixPos);
        });
    }
    window.addEventListener('resize', () => {
        resizeScadnanoGridCanvas();
        if (scadnanoGridEditor && typeof scadnanoGridEditor.resize === 'function') {
            scadnanoGridEditor.resize();
        }
    });
}
initScadnanoGridPaneControls();
/**
 * Toggles the disabled state of the grid dropdown based on the checkbox.
 * @param checkboxElement - The HTMLInputElement (checkbox) that was clicked.
 */
function toggleGridDropdown(checkboxElement) {
    // We cast to HTMLSelectElement so TypeScript knows it has a 'disabled' property
    const gridDropdown = document.getElementById('scadnanoGrid');
    if (gridDropdown) {
        gridDropdown.disabled = !checkboxElement.checked;
    }
}
