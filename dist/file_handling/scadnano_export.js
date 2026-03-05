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
function scadnanoDialogExport() {
    const nameInput = document.getElementById('scadnanoFilename');
    const helixPosCheckbox = document.getElementById('scadnanoIncludeHPos');
    const scadnanoGrid = document.getElementById('scadnanoGrid');
    if (!nameInput || !helixPosCheckbox || !scadnanoGrid) {
        console.warn('scadnano export dialog missing inputs');
        return;
    }
    const fileName = nameInput.value.trim() || 'output';
    if (helixPosCheckbox.checked) {
        const grid = scadnanoGrid.value === 'honeycomb' ? 'honeycomb' : 'square';
        makeScadnanowHexPos(fileName, grid);
    }
    else {
        makeScadnanoJsonFile(fileName, 'square');
    }
}
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
