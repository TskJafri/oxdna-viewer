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
}

namespace scadnanoExport {
    const TOLERANCE = 3;
    const MAX_ITERATIONS = 10;

    let layout: ScadnanoLayout | null = null;
    let connections: Array<[number, number]> = [];
    let gridEditor: any = null;
    let gridEditorType: ScadnanoGridType | null = null;
    let suppressSelectionCallback = false;

    // ── pipeline + export ────────────────────────────────────────────────────

    function nucleotideMap(): Map<number, Nucleotide> {
        const out = new Map<number, Nucleotide>();
        elements.forEach((e, id) => { if (e instanceof Nucleotide) out.set(id, e); });
        return out;
    }

    // The one call sign. Everything else in this file consumes its output.
    function runPipeline(lattice: ScadnanoRequestedGridType, wireframe: boolean): ScadnanoLayout {
        const nucleotides = nucleotideMap();
        layout = toscad.layoutPipeline(nucleotides, {
            tolerance: TOLERANCE,
            lattice,
            maxIterations: MAX_ITERATIONS,
            wireframe
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
        resizeCanvas();

        const helixPos = toHelixPosMap(helixPosInput);
        if (!helixPos || helixPos.size === 0) {
            notify('No helix positions available for the grid view.', 'warning');
            return;
        }

        editor.loadFromHelixPos(helixPos);
        if (typeof editor.setConnections === 'function') editor.setConnections(connections);
        publishHelixPos();
    }

    export function hideGrid(): void {
        document.body.classList.remove('scadnano-grid-open');
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
