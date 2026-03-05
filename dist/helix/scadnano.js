/**
 * scadnano.ts  –  Standalone 2D Helix Position Editor
 *
 * Renders a parity-staggered honeycomb-like lattice using Three.js.
 *
 * This is intentionally NOT standard odd-q offset. Layout is driven by parity:
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
var scadnano;
(function (scadnano) {
    // ── Grid geometry constants ────────────────────────────────────────────────
    //
    // D is the center-to-center distance between touching circles.
    // Horizontal spacing = D * 0.866025
    // Row baseline spacing = 1.5 * D
    // Parity offset = ±D/4 based on (col + row) parity
    scadnano.COL_SPACING = 2.5;
    scadnano.ROW_SPACING = 1.5 * scadnano.COL_SPACING;
    const X_SCALE = 0.866025 * scadnano.COL_SPACING;
    const PARITY_OFFSET = scadnano.COL_SPACING / 4;
    // Node circle radius in Three.js world units
    const NODE_RADIUS = 0.55;
    // Ghost dot radius (background grid marker)
    const GHOST_RADIUS = 0.14;
    const RING_DEFAULT_COLOR = 0x000000;
    const RING_SELECTED_COLOR = 0xff4da6;
    // ── Coordinate helpers ────────────────────────────────────────────────────
    /** Convert a grid coordinate (col, row) to Three.js world coords. */
    function oddQToWorld(col, row) {
        const isEvenParity = (((col + row) & 1) === 0);
        const offset = isEvenParity ? -PARITY_OFFSET : PARITY_OFFSET;
        const x = col * X_SCALE;
        const screenY = row * scadnano.ROW_SPACING + offset;
        const y = -screenY;
        return new THREE.Vector2(x, y);
    }
    scadnano.oddQToWorld = oddQToWorld;
    /** Find the nearest grid cell to a world position. */
    function worldToNearestCell(wx, wy) {
        const colEst = Math.round(wx / X_SCALE);
        const rowEst = Math.round((-wy) / scadnano.ROW_SPACING);
        let best = null;
        let bestDist = Infinity;
        for (let dc = -2; dc <= 2; dc++) {
            for (let dr = -2; dr <= 2; dr++) {
                const c = colEst + dc;
                const r = rowEst + dr;
                const p = oddQToWorld(c, r);
                const d = Math.hypot(p.x - wx, p.y - wy);
                if (d < bestDist) {
                    bestDist = d;
                    best = { col: c, row: r };
                }
            }
        }
        if (!best)
            return null;
        return bestDist <= scadnano.COL_SPACING * 0.8 ? best : null;
    }
    scadnano.worldToNearestCell = worldToNearestCell;
    // ── Colour palette for auto-assigning colours to new nodes ────────────────
    // Saturated, clearly distinguishable colours that all look good on white.
    const PALETTE = [
        0x1565c0, 0x2e7d32, 0xc62828, 0x6a1b9a,
        0xe65100, 0x00695c, 0xad1457, 0x4527a0,
        0x0277bd, 0x558b2f, 0x4e342e, 0x00838f,
    ];
    let _paletteIdx = 0;
    function nextColor() { return PALETTE[_paletteIdx++ % PALETTE.length]; }
    // ── Main editor class ──────────────────────────────────────────────────────
    class HoneycombEditor {
        canvas;
        // Three.js core
        scene;
        camera;
        renderer;
        // Node state
        records = new Map();
        // Ghost (background grid) meshes – reused geometry
        ghostGeo;
        ghostMat;
        ghostMeshes = new Map();
        // Shared geometry/material for user-placed nodes
        nodeGeo;
        // Shared ring geometry drawn around each node
        ringGeo;
        ringMat;
        // Interaction state
        isPanning = false;
        hasDragged = false;
        panLast = new THREE.Vector2();
        selectedKey = null;
        draggingNodeKey = null;
        // Grid extent (inclusive)
        minCol;
        maxCol;
        minRow;
        maxRow;
        // Frustum size (half-height in world units) – modified by zoom
        frustumHalfH = 14;
        // Callback fired when the node set changes
        onNodesChanged = null;
        constructor(canvas, options = {}) {
            this.canvas = canvas;
            const { gridCols = [-2, 8], gridRows = [-2, 10], initialNodes = [] } = options;
            [this.minCol, this.maxCol] = gridCols;
            [this.minRow, this.maxRow] = gridRows;
            // ── Scene ──────────────────────────────────────────────────────────
            this.scene = new THREE.Scene();
            this.scene.background = new THREE.Color(0xfafafa); // light mode
            // ── Camera ─────────────────────────────────────────────────────────
            const aspect = canvas.clientWidth / canvas.clientHeight || 1;
            this.frustumHalfH = 14;
            this.camera = new THREE.OrthographicCamera(-this.frustumHalfH * aspect, this.frustumHalfH * aspect, this.frustumHalfH, -this.frustumHalfH, -500, 500);
            this.camera.position.z = 10;
            // ── Renderer ───────────────────────────────────────────────────────
            this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
            this.renderer.setPixelRatio(window.devicePixelRatio || 1);
            this.renderer.setSize(canvas.clientWidth, canvas.clientHeight, false);
            // ── Shared geometries ──────────────────────────────────────────────
            this.ghostGeo = new THREE.CircleGeometry(GHOST_RADIUS, 12);
            this.ghostMat = new THREE.MeshBasicMaterial({ color: 0xc0c8d8 }); // light-mode dot
            this.nodeGeo = new THREE.CircleGeometry(NODE_RADIUS, 36);
            this.ringGeo = new THREE.RingGeometry(NODE_RADIUS + 0.04, NODE_RADIUS + 0.20, 36);
            this.ringMat = new THREE.MeshBasicMaterial({ color: RING_DEFAULT_COLOR, opacity: 0.18, transparent: true, side: THREE.DoubleSide });
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
        addNode(node) {
            const key = this._key(node.col, node.row);
            if (this.records.has(key))
                return;
            const pos = oddQToWorld(node.col, node.row);
            const color = node.color ?? nextColor();
            const mat = new THREE.MeshBasicMaterial({ color });
            const mesh = new THREE.Mesh(this.nodeGeo, mat);
            mesh.position.set(pos.x, pos.y, 0);
            mesh.userData = { col: node.col, row: node.row, key };
            // Ring highlight
            const ringMaterial = this.ringMat.clone();
            const ring = new THREE.Mesh(this.ringGeo, ringMaterial);
            ring.position.set(pos.x, pos.y, 0.5);
            this.scene.add(ring);
            mesh.userData.ring = ring;
            this.scene.add(mesh);
            this.records.set(key, { node: { ...node, color }, mesh });
            this.onNodesChanged?.();
        }
        /** Remove a node by grid coordinate. No-op if cell is empty. */
        removeNode(col, row) {
            const key = this._key(col, row);
            const rec = this.records.get(key);
            if (!rec)
                return;
            this.scene.remove(rec.mesh);
            const ring = rec.mesh.userData.ring;
            if (ring)
                this.scene.remove(ring);
            this.records.delete(key);
            if (this.selectedKey === key)
                this.selectedKey = null;
            if (this.draggingNodeKey === key)
                this.draggingNodeKey = null;
            this.onNodesChanged?.();
        }
        /** Replace the entire node set. */
        setNodes(nodes) {
            this._setSelectedKey(null);
            this.draggingNodeKey = null;
            const keys = Array.from(this.records.keys());
            keys.forEach(k => {
                const rec = this.records.get(k);
                const ring = rec.mesh.userData.ring;
                if (ring)
                    this.scene.remove(ring);
                this.scene.remove(rec.mesh);
            });
            this.records.clear();
            nodes.forEach(n => this.addNode(n));
        }
        /** Return a snapshot of all current nodes. */
        getNodes() {
            return Array.from(this.records.values()).map(r => ({ ...r.node }));
        }
        /** Expand the pre-rendered ghost extent. */
        expandGrid(minCol, maxCol, minRow, maxRow) {
            this.minCol = minCol;
            this.maxCol = maxCol;
            this.minRow = minRow;
            this.maxRow = maxRow;
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
        loadFromHelixPos(posMap, colors) {
            const nodes = [];
            posMap.forEach(([col, row], helixId) => {
                nodes.push({
                    id: helixId,
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
        fitView() {
            if (!this.records.size)
                return;
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
            const span = Math.max(maxX - minX, maxY - minY) / 2 + scadnano.COL_SPACING * 2;
            const aspect = this.canvas.clientWidth / this.canvas.clientHeight || 1;
            this.frustumHalfH = span;
            this.camera.left = -span * aspect;
            this.camera.right = span * aspect;
            this.camera.top = span;
            this.camera.bottom = -span;
            this.camera.updateProjectionMatrix();
        }
        dispose() {
            this.renderer.dispose();
        }
        // ── Private helpers ────────────────────────────────────────────────────
        _key(col, row) { return `${col},${row}`; }
        _setRecordSelected(rec, selected) {
            const ring = rec.mesh.userData.ring;
            if (!ring)
                return;
            const material = ring.material;
            material.color.setHex(selected ? RING_SELECTED_COLOR : RING_DEFAULT_COLOR);
            material.opacity = selected ? 0.9 : 0.18;
            material.needsUpdate = true;
        }
        _setSelectedKey(key) {
            const nextKey = key && this.records.has(key) ? key : null;
            if (this.selectedKey === nextKey)
                return;
            if (this.selectedKey) {
                const prev = this.records.get(this.selectedKey);
                if (prev)
                    this._setRecordSelected(prev, false);
            }
            this.selectedKey = nextKey;
            if (this.selectedKey) {
                const current = this.records.get(this.selectedKey);
                if (current)
                    this._setRecordSelected(current, true);
            }
        }
        _cellFromMouseEvent(e) {
            const ndc = this._screenToNDC(e.clientX, e.clientY);
            const world = this._ndcToWorld(ndc);
            return worldToNearestCell(world.x, world.y);
        }
        _moveNode(fromKey, toCol, toRow) {
            const rec = this.records.get(fromKey);
            if (!rec)
                return fromKey;
            const toKey = this._key(toCol, toRow);
            if (toKey === fromKey)
                return fromKey;
            if (this.records.has(toKey))
                return fromKey;
            const pos = oddQToWorld(toCol, toRow);
            rec.node.col = toCol;
            rec.node.row = toRow;
            rec.mesh.position.set(pos.x, pos.y, 0);
            rec.mesh.userData.col = toCol;
            rec.mesh.userData.row = toRow;
            rec.mesh.userData.key = toKey;
            const ring = rec.mesh.userData.ring;
            if (ring)
                ring.position.set(pos.x, pos.y, 0.5);
            this.records.delete(fromKey);
            this.records.set(toKey, rec);
            if (this.selectedKey === fromKey)
                this.selectedKey = toKey;
            this.onNodesChanged?.();
            return toKey;
        }
        _buildGrid(minCol, maxCol, minRow, maxRow) {
            for (let col = minCol; col <= maxCol; col++) {
                for (let row = minRow; row <= maxRow; row++) {
                    const key = this._key(col, row);
                    if (this.ghostMeshes.has(key))
                        continue;
                    const p = oddQToWorld(col, row);
                    const mesh = new THREE.Mesh(this.ghostGeo, this.ghostMat);
                    mesh.position.set(p.x, p.y, -1);
                    mesh.userData.isGhost = true;
                    mesh.userData.key = key;
                    this.scene.add(mesh);
                    this.ghostMeshes.set(key, mesh);
                }
            }
        }
        _ensureGridCoverage() {
            const cam = this.camera;
            const worldLeft = cam.left + cam.position.x;
            const worldRight = cam.right + cam.position.x;
            const worldTop = cam.top + cam.position.y;
            const worldBottom = cam.bottom + cam.position.y;
            const colMin = Math.floor(worldLeft / X_SCALE) - 3;
            const colMax = Math.ceil(worldRight / X_SCALE) + 3;
            const screenYMin = -worldTop;
            const screenYMax = -worldBottom;
            const rowMin = Math.floor(screenYMin / scadnano.ROW_SPACING) - 3;
            const rowMax = Math.ceil(screenYMax / scadnano.ROW_SPACING) + 3;
            this._buildGrid(colMin, colMax, rowMin, rowMax);
        }
        _buildAxisLabels() {
            // We rely on HTML overlay labels in editor-dev.html; nothing to do here.
            // Stub kept for future canvas-based label rendering.
        }
        _screenToNDC(clientX, clientY) {
            const rect = this.canvas.getBoundingClientRect();
            return new THREE.Vector2(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1);
        }
        _ndcToWorld(ndc) {
            const cam = this.camera;
            const x = ndc.x * (cam.right - cam.left) / 2 + (cam.right + cam.left) / 2 + cam.position.x;
            const y = ndc.y * (cam.top - cam.bottom) / 2 + (cam.top + cam.bottom) / 2 + cam.position.y;
            return new THREE.Vector2(x, y);
        }
        _bindEvents() {
            this.canvas.addEventListener('click', this._onClick.bind(this));
            this.canvas.addEventListener('mousedown', this._onMouseDown.bind(this));
            this.canvas.addEventListener('mousemove', this._onMouseMove.bind(this));
            this.canvas.addEventListener('mouseup', this._onMouseUp.bind(this));
            this.canvas.addEventListener('wheel', this._onWheel.bind(this), { passive: false });
            this.canvas.addEventListener('contextmenu', e => e.preventDefault());
            window.addEventListener('resize', this._onResize.bind(this));
        }
        _onClick(e) {
            if (this.hasDragged)
                return; // suppress toggle after a pan
            if (e.button !== 0)
                return;
            const cell = this._cellFromMouseEvent(e);
            if (!cell)
                return;
            const key = this._key(cell.col, cell.row);
            if (this.records.has(key)) {
                this._setSelectedKey(key);
            }
            else {
                // Debug mode behavior: allow adding nodes on empty cells.
                this.addNode({ col: cell.col, row: cell.row, id: Date.now() });
                this._setSelectedKey(this._key(cell.col, cell.row));
            }
        }
        _onMouseDown(e) {
            // Middle (button 1) or right-click (button 2) → pan
            if (e.button === 1 || e.button === 2) {
                this.isPanning = true;
                this.hasDragged = false;
                this.panLast.set(e.clientX, e.clientY);
                e.preventDefault();
                return;
            }
            if (e.button === 0) {
                const cell = this._cellFromMouseEvent(e);
                if (!cell)
                    return;
                const key = this._key(cell.col, cell.row);
                if (key === this.selectedKey && this.records.has(key)) {
                    this.draggingNodeKey = key;
                    this.hasDragged = false;
                    e.preventDefault();
                }
            }
        }
        _onMouseMove(e) {
            if (this.draggingNodeKey) {
                const cell = this._cellFromMouseEvent(e);
                if (!cell)
                    return;
                const nextKey = this._moveNode(this.draggingNodeKey, cell.col, cell.row);
                if (nextKey !== this.draggingNodeKey) {
                    this.draggingNodeKey = nextKey;
                    this.hasDragged = true;
                }
                return;
            }
            if (!this.isPanning)
                return;
            const dx = e.clientX - this.panLast.x;
            const dy = e.clientY - this.panLast.y;
            this.panLast.set(e.clientX, e.clientY);
            const cam = this.camera;
            const scaleX = (cam.right - cam.left) / this.canvas.clientWidth;
            const scaleY = (cam.top - cam.bottom) / this.canvas.clientHeight;
            cam.position.x -= dx * scaleX;
            cam.position.y += dy * scaleY;
            if (Math.abs(dx) > 2 || Math.abs(dy) > 2)
                this.hasDragged = true;
        }
        _onMouseUp(_e) {
            const wasDraggingNode = this.draggingNodeKey !== null;
            if (this.draggingNodeKey)
                this.draggingNodeKey = null;
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
        _onWheel(e) {
            e.preventDefault();
            const factor = e.deltaY > 0 ? 1.12 : 1 / 1.12;
            const cam = this.camera;
            // Zoom around the mouse cursor
            const ndc = this._screenToNDC(e.clientX, e.clientY);
            const wx = this._ndcToWorld(ndc).x;
            const wy = this._ndcToWorld(ndc).y;
            cam.left = (cam.left - wx) * factor + wx;
            cam.right = (cam.right - wx) * factor + wx;
            cam.top = (cam.top - wy) * factor + wy;
            cam.bottom = (cam.bottom - wy) * factor + wy;
            cam.updateProjectionMatrix();
        }
        _onResize() {
            const w = this.canvas.clientWidth;
            const h = this.canvas.clientHeight;
            if (w === 0 || h === 0)
                return;
            const aspect = w / h;
            const halfH = (this.camera.top - this.camera.bottom) / 2;
            // Keep the vertical extent, adjust horizontal to match new aspect ratio
            const cx = (this.camera.left + this.camera.right) / 2;
            this.camera.left = cx - halfH * aspect;
            this.camera.right = cx + halfH * aspect;
            this.camera.updateProjectionMatrix();
            this.renderer.setSize(w, h, false);
            this._ensureGridCoverage();
        }
        _animate() {
            requestAnimationFrame(this._animate.bind(this));
            this._ensureGridCoverage();
            this.renderer.render(this.scene, this.camera);
        }
    }
    scadnano.HoneycombEditor = HoneycombEditor;
})(scadnano || (scadnano = {}));
