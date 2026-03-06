/**
 * scadnano.ts  –  Standalone 2D Helix Position Editor
 *
 * Renders editable helix lattices using Three.js.
 *
 * Honeycomb layout is intentionally NOT standard odd-q offset. It is driven by parity:
 *   if (col + row) is even => shift UP by D/4
 *   if (col + row) is odd  => shift DOWN by D/4
 *
 * Screen-space formula (y grows downward):
 *   Px = col * (D * 0.866025)
 *   Py = row * (1.5 * D) + Offset
 *
 * Three.js world uses y-up, so we store:
 *   world_y = -Py
 *
 * Usage:
 *   const editor = new scadnano.HoneycombEditor(canvas);
 *   editor.addNode({ col: 0, row: 0, id: 0 });
 *   editor.loadFromHelixPos(toscad.HelixPos(grid, helices));
 */

namespace scadnano {

    export type GridLayout = 'honeycomb' | 'square';

    // ── Grid geometry constants ────────────────────────────────────────────────
    //
    // D is the center-to-center distance between touching circles.
    // Horizontal spacing = D * 0.866025
    // Row baseline spacing = 1.5 * D
    // Parity offset = ±D/4 based on (col + row) parity
    export const COL_SPACING: number = 2.5;
    export const ROW_SPACING: number = 1.5 * COL_SPACING;
    const HONEYCOMB_X_SCALE = 0.866025 * COL_SPACING;
    const HONEYCOMB_PARITY_OFFSET = COL_SPACING / 4;
    const SQUARE_X_SCALE = COL_SPACING;
    const SQUARE_ROW_SPACING = COL_SPACING;

    // Node circle radius in Three.js world units
    const NODE_RADIUS  = 0.55;
    // Ghost dot radius (background grid marker)
    const GHOST_RADIUS = 0.14;
    // All dots (grid ghosts + nodes) use one shared yellow.
    const DOT_COLOR = 0xffd400;
    const RING_DEFAULT_COLOR = 0x000000;
    const RING_SELECTED_COLOR = 0xff4da6;

    // ── Types ─────────────────────────────────────────────────────────────────

    export interface HelixNode {
        col:    number;
        row:    number;
        id:     number;
        color?: number;
        label?: string;
    }

    export interface HelixConnection {
        from: number;
        to: number;
    }

    interface NodeRecord {
        node: HelixNode;
        mesh: THREE.Mesh;
        labelSprite: THREE.Sprite;
    }

    interface EditorOptions {
        gridCols?: [number, number];
        gridRows?: [number, number];
        initialNodes?: HelixNode[];
        layout?: GridLayout;
    }

    // ── Coordinate helpers ────────────────────────────────────────────────────

    /** Convert a grid coordinate (col, row) to Three.js world coords. */
    export function oddQToWorld(col: number, row: number): THREE.Vector2 {
        return honeycombToWorld(col, row);
    }

    export function honeycombToWorld(col: number, row: number): THREE.Vector2 {
        const isEvenParity = (((col + row) & 1) === 0);
        const offset = isEvenParity ? -HONEYCOMB_PARITY_OFFSET : HONEYCOMB_PARITY_OFFSET;
        const x = col * HONEYCOMB_X_SCALE;
        const screenY = row * ROW_SPACING + offset;
        const y = -screenY;
        return new THREE.Vector2(x, y);
    }

    export function squareToWorld(col: number, row: number): THREE.Vector2 {
        const x = col * SQUARE_X_SCALE;
        const screenY = row * SQUARE_ROW_SPACING;
        return new THREE.Vector2(x, -screenY);
    }

    function gridToWorld(col: number, row: number, layout: GridLayout): THREE.Vector2 {
        return layout === 'square' ? squareToWorld(col, row) : honeycombToWorld(col, row);
    }

    function honeycombWorldToNearestCell(wx: number, wy: number): { col: number; row: number } | null {
        const colEst = Math.round(wx / HONEYCOMB_X_SCALE);
        const rowEst = Math.round((-wy) / ROW_SPACING);

        let best: { col: number; row: number } | null = null;
        let bestDist = Infinity;
        for (let dc = -2; dc <= 2; dc++) {
            for (let dr = -2; dr <= 2; dr++) {
                const c = colEst + dc;
                const r = rowEst + dr;
                const p = honeycombToWorld(c, r);
                const d = Math.hypot(p.x - wx, p.y - wy);
                if (d < bestDist) { bestDist = d; best = { col: c, row: r }; }
            }
        }

        if (!best) return null;
        return bestDist <= COL_SPACING * 0.8 ? best : null;
    }

    function squareWorldToNearestCell(wx: number, wy: number): { col: number; row: number } {
        return {
            col: Math.round(wx / SQUARE_X_SCALE),
            row: Math.round((-wy) / SQUARE_ROW_SPACING),
        };
    }

    /** Find the nearest grid cell to a world position. */
    export function worldToNearestCell(wx: number, wy: number, layout: GridLayout = 'honeycomb'): { col: number; row: number } | null {
        if (layout === 'square') return squareWorldToNearestCell(wx, wy);
        return honeycombWorldToNearestCell(wx, wy);
    }

    function estimateVisibleGridBounds(camera: THREE.OrthographicCamera, layout: GridLayout) {
        const worldLeft = camera.left + camera.position.x;
        const worldRight = camera.right + camera.position.x;
        const worldTop = camera.top + camera.position.y;
        const worldBottom = camera.bottom + camera.position.y;

        if (layout === 'square') {
            const colMin = Math.floor(worldLeft / SQUARE_X_SCALE) - 3;
            const colMax = Math.ceil(worldRight / SQUARE_X_SCALE) + 3;
            const screenYMin = -worldTop;
            const screenYMax = -worldBottom;
            const rowMin = Math.floor(screenYMin / SQUARE_ROW_SPACING) - 3;
            const rowMax = Math.ceil(screenYMax / SQUARE_ROW_SPACING) + 3;
            return { colMin, colMax, rowMin, rowMax };
        }

        const colMin = Math.floor(worldLeft / HONEYCOMB_X_SCALE) - 3;
        const colMax = Math.ceil(worldRight / HONEYCOMB_X_SCALE) + 3;
        const screenYMin = -worldTop;
        const screenYMax = -worldBottom;
        const rowMin = Math.floor(screenYMin / ROW_SPACING) - 3;
        const rowMax = Math.ceil(screenYMax / ROW_SPACING) + 3;
        return { colMin, colMax, rowMin, rowMax };
    }

    // ── Main editor class ──────────────────────────────────────────────────────

    export class HoneycombEditor {

        // Three.js core
        private scene:    THREE.Scene;
        private camera:   THREE.OrthographicCamera;
        private renderer: THREE.WebGLRenderer;

        // Node state
        private records: Map<string, NodeRecord> = new Map();
        protected readonly layout: GridLayout;

        // Ghost (background grid) meshes – reused geometry
        private ghostGeo: THREE.CircleGeometry;
        private ghostMat: THREE.MeshBasicMaterial;
        private ghostMeshes: Map<string, THREE.Mesh> = new Map();

        // Shared geometry/material for user-placed nodes
        private nodeGeo: THREE.CircleGeometry;

        // Shared ring geometry drawn around each node
        private ringGeo: THREE.RingGeometry;
        private ringMat: THREE.MeshBasicMaterial;

        // Interaction state
        private isPanning  = false;
        private hasDragged = false;
        private panLast    = new THREE.Vector2();
        private selectedKey: string | null = null;
        private draggingNodeKey: string | null = null;

        // Grid extent (inclusive)
        private minCol: number;
        private maxCol: number;
        private minRow: number;
        private maxRow: number;

        // Frustum size (half-height in world units) – modified by zoom
        private frustumHalfH: number = 14;

        // Callback fired when the node set changes
        public onNodesChanged: (() => void) | null = null;
        // Callback fired when the selected node changes
        public onNodeSelected: ((node: HelixNode | null) => void) | null = null;

        // Crossover connection visualization state
        private connections: HelixConnection[] = [];
        private connectionLines: THREE.LineSegments | null = null;
        private connectionLineMaterial: THREE.LineBasicMaterial;

        constructor(private canvas: HTMLCanvasElement, options: EditorOptions = {}) {
            const { gridCols = [-2, 8], gridRows = [-2, 10], initialNodes = [] } = options;
            [this.minCol, this.maxCol] = gridCols;
            [this.minRow, this.maxRow] = gridRows;
            this.layout = options.layout === 'square' ? 'square' : 'honeycomb';

            // ── Scene ──────────────────────────────────────────────────────────
            this.scene = new THREE.Scene();
            this.scene.background = new THREE.Color(0xfafafa);   // light mode

            // ── Camera ─────────────────────────────────────────────────────────
            const aspect = canvas.clientWidth / canvas.clientHeight || 1;
            this.frustumHalfH = 14;
            this.camera = new THREE.OrthographicCamera(
                -this.frustumHalfH * aspect,  this.frustumHalfH * aspect,
                 this.frustumHalfH,           -this.frustumHalfH,
                -500, 500
            );
            this.camera.position.z = 10;

            // ── Renderer ───────────────────────────────────────────────────────
            this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
            this.renderer.setPixelRatio(window.devicePixelRatio || 1);
            this.renderer.setSize(canvas.clientWidth, canvas.clientHeight, false);

            // ── Shared geometries ──────────────────────────────────────────────
            this.ghostGeo = new THREE.CircleGeometry(GHOST_RADIUS, 12);
            this.ghostMat = new THREE.MeshBasicMaterial({ color: DOT_COLOR });
            this.nodeGeo  = new THREE.CircleGeometry(NODE_RADIUS, 36);
            this.ringGeo  = new THREE.RingGeometry(NODE_RADIUS + 0.04, NODE_RADIUS + 0.20, 36);
            this.ringMat  = new THREE.MeshBasicMaterial({ color: RING_DEFAULT_COLOR, opacity: 0.18, transparent: true, side: THREE.DoubleSide });
            this.connectionLineMaterial = new THREE.LineBasicMaterial({ color: 0x8c8c8c });

            // ── Build background grid ──────────────────────────────────────────
            this._buildGrid(this.minCol, this.maxCol, this.minRow, this.maxRow);
            this._ensureGridCoverage();

            // ── Grid axis labels (HTML overlay handled in CSS/HTML) ────────────
            this._buildAxisLabels();

            // ── Initial nodes ──────────────────────────────────────────────────
            initialNodes.forEach(n => this.addNode(n));

            // ── Events ─────────────────────────────────────────────────────────
            this._bindEvents();

            // ── Render loop ────────────────────────────────────────────────────
            this._animate();
        }

        // ── Public API ─────────────────────────────────────────────────────────

        /** Add a helix node to the grid. No-op if that cell is already occupied. */
        addNode(node: HelixNode): void {
            const key = this._key(node.col, node.row);
            if (this.records.has(key)) return;

            const pos   = gridToWorld(node.col, node.row, this.layout);
            const color = DOT_COLOR;

            const mat  = new THREE.MeshBasicMaterial({ color });
            const mesh = new THREE.Mesh(this.nodeGeo, mat);
            mesh.position.set(pos.x, pos.y, 0);
            mesh.userData = { col: node.col, row: node.row, key };

            // Ring highlight
            const ringMaterial = this.ringMat.clone();
            const ring = new THREE.Mesh(this.ringGeo, ringMaterial);
            ring.position.set(pos.x, pos.y, 0.5);
            this.scene.add(ring);
            mesh.userData.ring = ring;

            // Always show helix number for each node.
            const labelSprite = this._createNodeLabelSprite(String(node.id));
            labelSprite.position.set(pos.x, pos.y, 1.2);
            this.scene.add(labelSprite);
            mesh.userData.labelSprite = labelSprite;

            this.scene.add(mesh);
            this.records.set(key, { node: { ...node, color }, mesh, labelSprite });

            this._rebuildConnectionLines();

            this.onNodesChanged?.();
        }

        /** Remove a node by grid coordinate. No-op if cell is empty. */
        removeNode(col: number, row: number): void {
            const key = this._key(col, row);
            const rec = this.records.get(key);
            if (!rec) return;

            this.scene.remove(rec.mesh);
            const ring = rec.mesh.userData.ring as THREE.Mesh | undefined;
            if (ring) this.scene.remove(ring);
            this.scene.remove(rec.labelSprite);
            const labelMaterial = rec.labelSprite.material as THREE.SpriteMaterial;
            const labelTexture = labelMaterial.map;
            if (labelTexture) labelTexture.dispose();
            labelMaterial.dispose();

            this.records.delete(key);
            if (this.selectedKey === key) this.selectedKey = null;
            if (this.draggingNodeKey === key) this.draggingNodeKey = null;
            this._rebuildConnectionLines();
            this.onNodesChanged?.();
        }

        /** Replace the entire node set. */
        setNodes(nodes: HelixNode[]): void {
            this._setSelectedKey(null);
            this.draggingNodeKey = null;
            const keys = Array.from(this.records.keys());
            keys.forEach(k => {
                const rec = this.records.get(k)!;
                const ring = rec.mesh.userData.ring as THREE.Mesh | undefined;
                if (ring) this.scene.remove(ring);
                this.scene.remove(rec.labelSprite);
                const labelMaterial = rec.labelSprite.material as THREE.SpriteMaterial;
                const labelTexture = labelMaterial.map;
                if (labelTexture) labelTexture.dispose();
                labelMaterial.dispose();
                this.scene.remove(rec.mesh);
            });
            this.records.clear();
            nodes.forEach(n => this.addNode(n));
            this._rebuildConnectionLines();
        }

        /** Return a snapshot of all current nodes. */
        getNodes(): HelixNode[] {
            return Array.from(this.records.values()).map(r => ({ ...r.node }));
        }

        /** Set crossover connections to render as thin, constant gray lines. */
        setConnections(connections: Array<[number, number]> | HelixConnection[]): void {
            const normalized: HelixConnection[] = [];
            connections.forEach((conn: any) => {
                const from = Number(Array.isArray(conn) ? conn[0] : conn?.from);
                const to = Number(Array.isArray(conn) ? conn[1] : conn?.to);
                if (!Number.isFinite(from) || !Number.isFinite(to) || from === to) return;
                normalized.push({ from, to });
            });
            this.connections = normalized;
            this._rebuildConnectionLines();
        }

        /** Expand the pre-rendered ghost extent. */
        expandGrid(minCol: number, maxCol: number, minRow: number, maxRow: number): void {
            this.minCol = minCol; this.maxCol = maxCol;
            this.minRow = minRow; this.maxRow = maxRow;
            this._buildGrid(minCol, maxCol, minRow, maxRow);
        }

        /**
         * Load helix positions produced by toscad.HelixPos().
         *
         *   const posMap = toscad.HelixPos(grid, helices);
         *   honeyEditor.loadFromHelixPos(posMap);
         *
         * @param posMap   Map<helixIndex, [col, row]> as returned by HelixPos.
         * @param colors   Optional Map<helixIndex, cssHexNumber> override.
         */
        loadFromHelixPos(
            posMap: Map<number, [number, number]>,
            colors?: Map<number, number>
        ): void {
            const nodes: HelixNode[] = [];
            posMap.forEach(([col, row], helixId) => {
                nodes.push({
                    id:    helixId,
                    col,
                    row,
                    label: String(helixId),
                    color: colors?.get(helixId),
                });
            });
            this.setNodes(nodes);
            // Fit the view after a short delay so the renderer has time to lay out.
            setTimeout(() => this.fitView(), 60);
        }

        /** Programmatically zoom/pan to fit all current nodes. */
        fitView(): void {
            if (!this.records.size) return;
            let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
            this.records.forEach(({ mesh }) => {
                minX = Math.min(minX, mesh.position.x);
                maxX = Math.max(maxX, mesh.position.x);
                minY = Math.min(minY, mesh.position.y);
                maxY = Math.max(maxY, mesh.position.y);
            });
            const cx = (minX + maxX) / 2;
            const cy = (minY + maxY) / 2;
            this.camera.position.x = cx;
            this.camera.position.y = cy;

            const span = Math.max(maxX - minX, maxY - minY) / 2 + COL_SPACING * 2;
            const aspect = this.canvas.clientWidth / this.canvas.clientHeight || 1;
            this.frustumHalfH = span;
            this.camera.left   = -span * aspect;
            this.camera.right  =  span * aspect;
            this.camera.top    =  span;
            this.camera.bottom = -span;
            this.camera.updateProjectionMatrix();
        }

        /** Recompute renderer/camera sizing from current canvas dimensions. */
        resize(): void {
            this._onResize();
        }

        dispose(): void {
            this.renderer.dispose();
        }

        // ── Private helpers ────────────────────────────────────────────────────

        private _key(col: number, row: number): string { return `${col},${row}`; }

        private _setRecordSelected(rec: NodeRecord, selected: boolean): void {
            const ring = rec.mesh.userData.ring as THREE.Mesh | undefined;
            if (!ring) return;
            const material = ring.material as THREE.MeshBasicMaterial;
            material.color.setHex(selected ? RING_SELECTED_COLOR : RING_DEFAULT_COLOR);
            material.opacity = selected ? 0.9 : 0.18;
            material.needsUpdate = true;
        }

        private _setSelectedKey(key: string | null): void {
            const nextKey = key && this.records.has(key) ? key : null;
            if (this.selectedKey === nextKey) return;

            if (this.selectedKey) {
                const prev = this.records.get(this.selectedKey);
                if (prev) this._setRecordSelected(prev, false);
            }

            this.selectedKey = nextKey;

            if (this.selectedKey) {
                const current = this.records.get(this.selectedKey);
                if (current) this._setRecordSelected(current, true);
            }

            const selectedNode = this.selectedKey ? this.records.get(this.selectedKey)?.node ?? null : null;
            this.onNodeSelected?.(selectedNode ? { ...selectedNode } : null);
        }

        private _cellFromMouseEvent(e: MouseEvent): { col: number; row: number } | null {
            const ndc = this._screenToNDC(e.clientX, e.clientY);
            const world = this._ndcToWorld(ndc);
            return worldToNearestCell(world.x, world.y, this.layout);
        }

        private _moveNode(fromKey: string, toCol: number, toRow: number): string {
            const rec = this.records.get(fromKey);
            if (!rec) return fromKey;

            const toKey = this._key(toCol, toRow);
            if (toKey === fromKey) return fromKey;
            if (this.records.has(toKey)) return fromKey;

            const pos = gridToWorld(toCol, toRow, this.layout);
            rec.node.col = toCol;
            rec.node.row = toRow;
            rec.mesh.position.set(pos.x, pos.y, 0);
            rec.mesh.userData.col = toCol;
            rec.mesh.userData.row = toRow;
            rec.mesh.userData.key = toKey;

            const ring = rec.mesh.userData.ring as THREE.Mesh | undefined;
            if (ring) ring.position.set(pos.x, pos.y, 0.5);
            rec.labelSprite.position.set(pos.x, pos.y, 1.2);

            this.records.delete(fromKey);
            this.records.set(toKey, rec);

            if (this.selectedKey === fromKey) this.selectedKey = toKey;
            this._rebuildConnectionLines();
            this.onNodesChanged?.();
            return toKey;
        }

        private _rebuildConnectionLines(): void {
            if (this.connectionLines) {
                this.scene.remove(this.connectionLines);
                this.connectionLines.geometry.dispose();
                this.connectionLines = null;
            }

            if (!this.connections.length) return;

            const idToPos = new Map<number, THREE.Vector3>();
            this.records.forEach((rec) => {
                idToPos.set(Number(rec.node.id), rec.mesh.position);
            });

            const positions: number[] = [];
            this.connections.forEach((conn) => {
                const fromPos = idToPos.get(conn.from);
                const toPos = idToPos.get(conn.to);
                if (!fromPos || !toPos) return;

                positions.push(fromPos.x, fromPos.y, 0.2);
                positions.push(toPos.x, toPos.y, 0.2);
            });

            if (!positions.length) return;

            const geometry = new THREE.BufferGeometry();
            const positionAttr = new THREE.Float32BufferAttribute(positions, 3);
            const geometryAny = geometry as any;
            // some older version issue, thus the fallback check. As of 2026-06, setAttribute does not exist, but addAttribute does. In the future, we can remove the addAttribute fallback.
            if (typeof geometryAny.setAttribute === 'function') {
                geometryAny.setAttribute('position', positionAttr);
            } else if (typeof geometryAny.addAttribute === 'function') {
                geometryAny.addAttribute('position', positionAttr);
            }
            this.connectionLines = new THREE.LineSegments(geometry, this.connectionLineMaterial);
            this.scene.add(this.connectionLines);
        }

        private _buildGrid(minCol: number, maxCol: number, minRow: number, maxRow: number): void {
            for (let col = minCol; col <= maxCol; col++) {
                for (let row = minRow; row <= maxRow; row++) {
                    const key = this._key(col, row);
                    if (this.ghostMeshes.has(key)) continue;
                    const p    = gridToWorld(col, row, this.layout);
                    const mesh = new THREE.Mesh(this.ghostGeo, this.ghostMat);
                    mesh.position.set(p.x, p.y, -1);
                    mesh.userData.isGhost = true;
                    mesh.userData.key = key;
                    this.scene.add(mesh);
                    this.ghostMeshes.set(key, mesh);
                }
            }
        }

        private _ensureGridCoverage(): void {
            const bounds = estimateVisibleGridBounds(this.camera, this.layout);
            this._buildGrid(bounds.colMin, bounds.colMax, bounds.rowMin, bounds.rowMax);
        }

        private _buildAxisLabels(): void {
            // We rely on HTML overlay labels in editor-dev.html; nothing to do here.
            // Stub kept for future canvas-based label rendering.
        }

        private _createNodeLabelSprite(text: string): THREE.Sprite {
            const canvas = document.createElement('canvas');
            canvas.width = 128;
            canvas.height = 64;
            const ctx = canvas.getContext('2d');
            if (!ctx) {
                const fallbackMat = new THREE.SpriteMaterial({ color: 0x111111 });
                return new THREE.Sprite(fallbackMat);
            }

            ctx.clearRect(0, 0, canvas.width, canvas.height);
            ctx.font = 'bold 30px Arial';
            ctx.fillStyle = '#111111';
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            ctx.fillText(text, canvas.width / 2, canvas.height / 2);

            const texture = new THREE.CanvasTexture(canvas);
            texture.needsUpdate = true;

            const material = new THREE.SpriteMaterial({
                map: texture,
                transparent: true,
                depthTest: false,
            });
            const sprite = new THREE.Sprite(material);
            sprite.scale.set(1.5, 0.75, 1);
            return sprite;
        }

        private _screenToNDC(clientX: number, clientY: number): THREE.Vector2 {
            const rect = this.canvas.getBoundingClientRect();
            return new THREE.Vector2(
                 ((clientX - rect.left) / rect.width)  * 2 - 1,
                -((clientY - rect.top)  / rect.height) * 2 + 1
            );
        }

        private _ndcToWorld(ndc: THREE.Vector2): THREE.Vector2 {
            const cam = this.camera;
            const x   = ndc.x * (cam.right - cam.left) / 2 + (cam.right + cam.left) / 2 + cam.position.x;
            const y   = ndc.y * (cam.top - cam.bottom) / 2 + (cam.top + cam.bottom) / 2 + cam.position.y;
            return new THREE.Vector2(x, y);
        }

        private _bindEvents(): void {
            this.canvas.addEventListener('click',       this._onClick.bind(this));
            this.canvas.addEventListener('mousedown',   this._onMouseDown.bind(this));
            this.canvas.addEventListener('mousemove',   this._onMouseMove.bind(this));
            this.canvas.addEventListener('mouseup',     this._onMouseUp.bind(this));
            this.canvas.addEventListener('wheel',       this._onWheel.bind(this), { passive: false });
            this.canvas.addEventListener('contextmenu', e => e.preventDefault());
            window.addEventListener('resize', this._onResize.bind(this));
        }

        private _onClick(e: MouseEvent): void {
            if (this.hasDragged) return;         // suppress toggle after a pan
            if (e.button !== 0) return;

            const cell  = this._cellFromMouseEvent(e);
            if (!cell) return;

            const key = this._key(cell.col, cell.row);
            if (this.records.has(key)) {
                this._setSelectedKey(key);
            }
        }

        private _onMouseDown(e: MouseEvent): void {
            // Middle (button 1) or right-click (button 2) → pan
            if (e.button === 1 || e.button === 2) {
                this.isPanning  = true;
                this.hasDragged = false;
                this.panLast.set(e.clientX, e.clientY);
                e.preventDefault();
                return;
            }

            if (e.button === 0) {
                const cell = this._cellFromMouseEvent(e);
                if (!cell) return;
                const key = this._key(cell.col, cell.row);
                if (key === this.selectedKey && this.records.has(key)) {
                    this.draggingNodeKey = key;
                    this.hasDragged = false;
                    e.preventDefault();
                }
            }
        }

        private _onMouseMove(e: MouseEvent): void {
            if (this.draggingNodeKey) {
                const cell = this._cellFromMouseEvent(e);
                if (!cell) return;
                const nextKey = this._moveNode(this.draggingNodeKey, cell.col, cell.row);
                if (nextKey !== this.draggingNodeKey) {
                    this.draggingNodeKey = nextKey;
                    this.hasDragged = true;
                }
                return;
            }

            if (!this.isPanning) return;
            const dx = e.clientX - this.panLast.x;
            const dy = e.clientY - this.panLast.y;
            this.panLast.set(e.clientX, e.clientY);

            const cam = this.camera;
            const scaleX = (cam.right - cam.left) / this.canvas.clientWidth;
            const scaleY = (cam.top - cam.bottom)  / this.canvas.clientHeight;
            cam.position.x -= dx * scaleX;
            cam.position.y += dy * scaleY;

            if (Math.abs(dx) > 2 || Math.abs(dy) > 2) this.hasDragged = true;
        }

        private _onMouseUp(_e: MouseEvent): void {
            const wasDraggingNode = this.draggingNodeKey !== null;
            if (this.draggingNodeKey) this.draggingNodeKey = null;

            if (this.isPanning) {
                this.isPanning = false;
                // Reset hasDragged on the next tick so click handler can check it first
                setTimeout(() => { this.hasDragged = false; }, 0);
                return;
            }

            if (wasDraggingNode) {
                // Reset hasDragged on the next tick so click handler can check it first
                setTimeout(() => { this.hasDragged = false; }, 0);
            }
        }

        private _onWheel(e: WheelEvent): void {
            e.preventDefault();
            const factor = e.deltaY > 0 ? 1.12 : 1 / 1.12;
            const cam    = this.camera;

            // Zoom around the mouse cursor
            const ndc   = this._screenToNDC(e.clientX, e.clientY);
            const wx    = this._ndcToWorld(ndc).x;
            const wy    = this._ndcToWorld(ndc).y;

            cam.left   = (cam.left   - wx) * factor + wx;
            cam.right  = (cam.right  - wx) * factor + wx;
            cam.top    = (cam.top    - wy) * factor + wy;
            cam.bottom = (cam.bottom - wy) * factor + wy;
            cam.updateProjectionMatrix();
        }

        private _onResize(): void {
            const w = this.canvas.clientWidth;
            const h = this.canvas.clientHeight;
            if (w === 0 || h === 0) return;

            const aspect  = w / h;
            const halfH   = (this.camera.top - this.camera.bottom) / 2;

            // Keep the vertical extent, adjust horizontal to match new aspect ratio
            const cx = (this.camera.left + this.camera.right) / 2;
            this.camera.left   = cx - halfH * aspect;
            this.camera.right  = cx + halfH * aspect;
            this.camera.updateProjectionMatrix();

            this.renderer.setSize(w, h, false);
            this._ensureGridCoverage();
        }

        private _animate(): void {
            requestAnimationFrame(this._animate.bind(this));
            this._ensureGridCoverage();
            this.renderer.render(this.scene, this.camera);
        }
    }

    export class SquareEditor extends HoneycombEditor {
        constructor(canvas: HTMLCanvasElement, options: Omit<EditorOptions, 'layout'> = {}) {
            super(canvas, { ...options, layout: 'square' });
        }
    }
}
