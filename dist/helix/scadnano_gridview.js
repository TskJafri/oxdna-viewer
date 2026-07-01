"use strict";
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
    const HONEYCOMB_X_SCALE = 0.866025 * scadnano.COL_SPACING;
    const HONEYCOMB_PARITY_OFFSET = scadnano.COL_SPACING / 4;
    const SQUARE_X_SCALE = scadnano.COL_SPACING;
    const SQUARE_ROW_SPACING = scadnano.COL_SPACING;
    const EDGE_PAN_MARGIN_PX = 48;
    const EDGE_PAN_MAX_SPEED_PX = 12;
    const GRID_EXPANSION_PAD = 4;
    // Node circle radius in Three.js world units
    const NODE_RADIUS = 0.55;
    // Ghost dot radius (background grid marker)
    const GHOST_RADIUS = 0.14;
    // All dots (grid ghosts + nodes) use one shared yellow.
    const DOT_COLOR = 0xffd400;
    const RING_DEFAULT_COLOR = 0x000000;
    const RING_SELECTED_COLOR = 0xff4da6;
    // ── Coordinate helpers ────────────────────────────────────────────────────
    /** Convert a grid coordinate (col, row) to Three.js world coords. */
    function oddQToWorld(col, row) {
        return honeycombToWorld(col, row);
    }
    scadnano.oddQToWorld = oddQToWorld;
    function honeycombToWorld(col, row) {
        const isEvenParity = (((col + row) & 1) === 0);
        const offset = isEvenParity ? -HONEYCOMB_PARITY_OFFSET : HONEYCOMB_PARITY_OFFSET;
        const x = col * HONEYCOMB_X_SCALE;
        const screenY = row * scadnano.ROW_SPACING + offset;
        const y = -screenY;
        return new THREE.Vector2(x, y);
    }
    scadnano.honeycombToWorld = honeycombToWorld;
    function squareToWorld(col, row) {
        const x = col * SQUARE_X_SCALE;
        const screenY = row * SQUARE_ROW_SPACING;
        return new THREE.Vector2(x, -screenY);
    }
    scadnano.squareToWorld = squareToWorld;
    function gridToWorld(col, row, layout) {
        return layout === 'square' ? squareToWorld(col, row) : honeycombToWorld(col, row);
    }
    function honeycombWorldToNearestCell(wx, wy) {
        const colEst = Math.round(wx / HONEYCOMB_X_SCALE);
        const rowEst = Math.round((-wy) / scadnano.ROW_SPACING);
        let best = null;
        let bestDist = Infinity;
        for (let dc = -2; dc <= 2; dc++) {
            for (let dr = -2; dr <= 2; dr++) {
                const c = colEst + dc;
                const r = rowEst + dr;
                const p = honeycombToWorld(c, r);
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
    function squareWorldToNearestCell(wx, wy) {
        return {
            col: Math.round(wx / SQUARE_X_SCALE),
            row: Math.round((-wy) / SQUARE_ROW_SPACING),
        };
    }
    /** Find the nearest grid cell to a world position. */
    function worldToNearestCell(wx, wy, layout = 'honeycomb') {
        if (layout === 'square')
            return squareWorldToNearestCell(wx, wy);
        return honeycombWorldToNearestCell(wx, wy);
    }
    scadnano.worldToNearestCell = worldToNearestCell;
    function estimateVisibleGridBounds(camera, layout) {
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
        const rowMin = Math.floor(screenYMin / scadnano.ROW_SPACING) - 3;
        const rowMax = Math.ceil(screenYMax / scadnano.ROW_SPACING) + 3;
        return { colMin, colMax, rowMin, rowMax };
    }
    // ── Main editor class ──────────────────────────────────────────────────────
    class HoneycombEditor {
        canvas;
        // Three.js core
        scene;
        camera;
        renderer;
        // Node state
        records = new Map();
        layout;
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
        selectedKeys = new Set();
        draggingNodeKey = null;
        draggingPointer = null;
        // Grid extent (inclusive)
        minCol;
        maxCol;
        minRow;
        maxRow;
        // Frustum size (half-height in world units) – modified by zoom
        frustumHalfH = 14;
        // Callback fired when the node set changes
        onNodesChanged = null;
        // Callback fired when the selected node changes
        onNodeSelected = null;
        // Callback fired when a user drag finishes and the node landed on a new cell.
        // Suppressed for programmatic moves so undo/redo don't pollute the history.
        onNodeMoved = null;
        // Records the cell a drag started in so we can report the full transition on drag-end.
        dragStartCell = null;
        // When true, _moveNode skips firing onNodeMoved (used by undo/redo paths).
        suppressMoveCallback = false;
        // Crossover connection visualization state
        connections = [];
        connectionLines = null;
        connectionLineMaterial;
        constructor(canvas, options = {}) {
            this.canvas = canvas;
            const { gridCols = [-2, 8], gridRows = [-2, 10], initialNodes = [] } = options;
            [this.minCol, this.maxCol] = gridCols;
            [this.minRow, this.maxRow] = gridRows;
            this.layout = options.layout === 'square' ? 'square' : 'honeycomb';
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
            this.ghostMat = new THREE.MeshBasicMaterial({ color: DOT_COLOR });
            this.nodeGeo = new THREE.CircleGeometry(NODE_RADIUS, 36);
            this.ringGeo = new THREE.RingGeometry(NODE_RADIUS + 0.04, NODE_RADIUS + 0.20, 36);
            this.ringMat = new THREE.MeshBasicMaterial({ color: RING_DEFAULT_COLOR, opacity: 0.18, transparent: true, side: THREE.DoubleSide });
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
        addNode(node) {
            const key = this._key(node.col, node.row);
            if (this.records.has(key))
                return;
            const pos = gridToWorld(node.col, node.row, this.layout);
            const color = DOT_COLOR;
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
        removeNode(col, row) {
            const key = this._key(col, row);
            const rec = this.records.get(key);
            if (!rec)
                return;
            this.scene.remove(rec.mesh);
            const ring = rec.mesh.userData.ring;
            if (ring)
                this.scene.remove(ring);
            this.scene.remove(rec.labelSprite);
            const labelMaterial = rec.labelSprite.material;
            const labelTexture = labelMaterial.map;
            if (labelTexture)
                labelTexture.dispose();
            labelMaterial.dispose();
            this.records.delete(key);
            if (this.selectedKey === key)
                this.selectedKey = null;
            if (this.selectedKeys.has(key))
                this.selectedKeys.delete(key);
            if (this.draggingNodeKey === key)
                this.draggingNodeKey = null;
            this._rebuildConnectionLines();
            this.onNodesChanged?.();
        }
        /** Replace the entire node set. */
        setNodes(nodes) {
            this._setSelectedKey(null);
            this.selectedKeys.clear();
            this.draggingNodeKey = null;
            const keys = Array.from(this.records.keys());
            keys.forEach(k => {
                const rec = this.records.get(k);
                const ring = rec.mesh.userData.ring;
                if (ring)
                    this.scene.remove(ring);
                this.scene.remove(rec.labelSprite);
                const labelMaterial = rec.labelSprite.material;
                const labelTexture = labelMaterial.map;
                if (labelTexture)
                    labelTexture.dispose();
                labelMaterial.dispose();
                this.scene.remove(rec.mesh);
            });
            this.records.clear();
            nodes.forEach(n => this.addNode(n));
            this._rebuildConnectionLines();
        }
        /** Return a snapshot of all current nodes. */
        getNodes() {
            return Array.from(this.records.values()).map(r => ({ ...r.node }));
        }
        /** Set crossover connections to render as thin, constant gray lines. */
        setConnections(connections) {
            const normalized = [];
            connections.forEach((conn) => {
                const from = Number(Array.isArray(conn) ? conn[0] : conn?.from);
                const to = Number(Array.isArray(conn) ? conn[1] : conn?.to);
                if (!Number.isFinite(from) || !Number.isFinite(to) || from === to)
                    return;
                normalized.push({ from, to });
            });
            this.connections = normalized;
            this._rebuildConnectionLines();
        }
        /** Select a node by helix id (no-op if not present). */
        selectNodeById(helixId) {
            for (const [key, record] of this.records.entries()) {
                if (record.node.id === helixId) {
                    this._setSelectedKey(key);
                    return;
                }
            }
        }
        /**
         * Move a node identified by helix id to a new (col, row). Used by undo/redo.
         * Does NOT fire onNodeMoved (that callback is reserved for genuine user drags).
         */
        moveNodeById(helixId, toCol, toRow) {
            for (const [key, record] of this.records.entries()) {
                if (record.node.id !== helixId)
                    continue;
                this.suppressMoveCallback = true;
                try {
                    this._moveNode(key, toCol, toRow);
                }
                finally {
                    this.suppressMoveCallback = false;
                }
                return true;
            }
            return false;
        }
        /** Get the helix ids of all currently selected nodes (primary + multi-selected). */
        getSelectedHelixIds() {
            const ids = [];
            const seen = new Set();
            const collect = (key) => {
                if (!key)
                    return;
                const rec = this.records.get(key);
                if (!rec)
                    return;
                const id = Number(rec.node.id);
                if (!Number.isFinite(id) || seen.has(id))
                    return;
                seen.add(id);
                ids.push(id);
            };
            collect(this.selectedKey);
            this.selectedKeys.forEach(k => collect(k));
            return ids;
        }
        /** Clear all selection state (primary + multi). */
        clearSelection() {
            // Walk multi-selection and turn off ring highlight on each.
            this.selectedKeys.forEach(k => {
                const rec = this.records.get(k);
                if (rec)
                    this._setRecordSelected(rec, false);
            });
            this.selectedKeys.clear();
            this._setSelectedKey(null);
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
        /** Recompute renderer/camera sizing from current canvas dimensions. */
        resize() {
            this._onResize();
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
                // Only clear the highlight if the previous primary isn't held by multi-select.
                if (prev && !this.selectedKeys.has(this.selectedKey)) {
                    this._setRecordSelected(prev, false);
                }
            }
            this.selectedKey = nextKey;
            if (this.selectedKey) {
                const current = this.records.get(this.selectedKey);
                if (current)
                    this._setRecordSelected(current, true);
            }
            const selectedNode = this.selectedKey ? this.records.get(this.selectedKey)?.node ?? null : null;
            this.onNodeSelected?.(selectedNode ? { ...selectedNode } : null);
        }
        // Toggle a key's membership in the multi-select set. Used by ctrl/cmd+click.
        _toggleMultiSelectKey(key) {
            const rec = this.records.get(key);
            if (!rec)
                return;
            if (this.selectedKeys.has(key)) {
                this.selectedKeys.delete(key);
                // Don't dim the ring if this key is also the primary selection.
                if (this.selectedKey !== key) {
                    this._setRecordSelected(rec, false);
                }
                return;
            }
            this.selectedKeys.add(key);
            this._setRecordSelected(rec, true);
        }
        _cellFromMouseEvent(e) {
            return this._cellFromClientPoint(e.clientX, e.clientY);
        }
        _cellFromClientPoint(clientX, clientY) {
            const ndc = this._screenToNDC(clientX, clientY);
            const world = this._ndcToWorld(ndc);
            return worldToNearestCell(world.x, world.y, this.layout);
        }
        _expandGridToIncludeCell(col, row) {
            let minCol = this.minCol;
            let maxCol = this.maxCol;
            let minRow = this.minRow;
            let maxRow = this.maxRow;
            if (col <= this.minCol + 1)
                minCol = Math.min(minCol, col - GRID_EXPANSION_PAD);
            if (col >= this.maxCol - 1)
                maxCol = Math.max(maxCol, col + GRID_EXPANSION_PAD);
            if (row <= this.minRow + 1)
                minRow = Math.min(minRow, row - GRID_EXPANSION_PAD);
            if (row >= this.maxRow - 1)
                maxRow = Math.max(maxRow, row + GRID_EXPANSION_PAD);
            if (minCol !== this.minCol || maxCol !== this.maxCol || minRow !== this.minRow || maxRow !== this.maxRow) {
                this.expandGrid(minCol, maxCol, minRow, maxRow);
            }
        }
        _edgePanFromPointer(clientX, clientY) {
            const rect = this.canvas.getBoundingClientRect();
            if (rect.width === 0 || rect.height === 0)
                return false;
            const leftDist = clientX - rect.left;
            const rightDist = rect.right - clientX;
            const topDist = clientY - rect.top;
            const bottomDist = rect.bottom - clientY;
            const panX = leftDist < EDGE_PAN_MARGIN_PX ? -(EDGE_PAN_MARGIN_PX - Math.max(0, leftDist)) / EDGE_PAN_MARGIN_PX
                : rightDist < EDGE_PAN_MARGIN_PX ? (EDGE_PAN_MARGIN_PX - Math.max(0, rightDist)) / EDGE_PAN_MARGIN_PX
                    : 0;
            const panY = topDist < EDGE_PAN_MARGIN_PX ? (EDGE_PAN_MARGIN_PX - Math.max(0, topDist)) / EDGE_PAN_MARGIN_PX
                : bottomDist < EDGE_PAN_MARGIN_PX ? -(EDGE_PAN_MARGIN_PX - Math.max(0, bottomDist)) / EDGE_PAN_MARGIN_PX
                    : 0;
            if (!panX && !panY)
                return false;
            const cam = this.camera;
            const scaleX = (cam.right - cam.left) / rect.width;
            const scaleY = (cam.top - cam.bottom) / rect.height;
            cam.position.x += panX * EDGE_PAN_MAX_SPEED_PX * scaleX;
            cam.position.y += panY * EDGE_PAN_MAX_SPEED_PX * scaleY;
            return true;
        }
        _syncDraggingNodeToPointer(clientX, clientY, allowEdgePan) {
            if (!this.draggingNodeKey)
                return;
            this.draggingPointer = { x: clientX, y: clientY };
            if (allowEdgePan) {
                this._edgePanFromPointer(clientX, clientY);
            }
            const cell = this._cellFromClientPoint(clientX, clientY);
            if (!cell)
                return;
            const nextKey = this._moveNode(this.draggingNodeKey, cell.col, cell.row);
            if (nextKey !== this.draggingNodeKey) {
                this.draggingNodeKey = nextKey;
            }
            this.hasDragged = true;
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
            const pos = gridToWorld(toCol, toRow, this.layout);
            rec.node.col = toCol;
            rec.node.row = toRow;
            rec.mesh.position.set(pos.x, pos.y, 0);
            rec.mesh.userData.col = toCol;
            rec.mesh.userData.row = toRow;
            rec.mesh.userData.key = toKey;
            const ring = rec.mesh.userData.ring;
            if (ring)
                ring.position.set(pos.x, pos.y, 0.5);
            rec.labelSprite.position.set(pos.x, pos.y, 1.2);
            this.records.delete(fromKey);
            this.records.set(toKey, rec);
            this._expandGridToIncludeCell(toCol, toRow);
            if (this.selectedKey === fromKey)
                this.selectedKey = toKey;
            if (this.selectedKeys.has(fromKey)) {
                this.selectedKeys.delete(fromKey);
                this.selectedKeys.add(toKey);
            }
            this._rebuildConnectionLines();
            this.onNodesChanged?.();
            return toKey;
        }
        _rebuildConnectionLines() {
            if (this.connectionLines) {
                this.scene.remove(this.connectionLines);
                this.connectionLines.geometry.dispose();
                this.connectionLines = null;
            }
            if (!this.connections.length)
                return;
            const idToPos = new Map();
            this.records.forEach((rec) => {
                idToPos.set(Number(rec.node.id), rec.mesh.position);
            });
            const positions = [];
            this.connections.forEach((conn) => {
                const fromPos = idToPos.get(conn.from);
                const toPos = idToPos.get(conn.to);
                if (!fromPos || !toPos)
                    return;
                positions.push(fromPos.x, fromPos.y, 0.2);
                positions.push(toPos.x, toPos.y, 0.2);
            });
            if (!positions.length)
                return;
            const geometry = new THREE.BufferGeometry();
            const positionAttr = new THREE.Float32BufferAttribute(positions, 3);
            const geometryAny = geometry;
            // some older version issue, thus the fallback check. As of 2026-06, setAttribute does not exist, but addAttribute does. In the future, we can remove the addAttribute fallback.
            if (typeof geometryAny.setAttribute === 'function') {
                geometryAny.setAttribute('position', positionAttr);
            }
            else if (typeof geometryAny.addAttribute === 'function') {
                geometryAny.addAttribute('position', positionAttr);
            }
            this.connectionLines = new THREE.LineSegments(geometry, this.connectionLineMaterial);
            this.scene.add(this.connectionLines);
        }
        _buildGrid(minCol, maxCol, minRow, maxRow) {
            for (let col = minCol; col <= maxCol; col++) {
                for (let row = minRow; row <= maxRow; row++) {
                    const key = this._key(col, row);
                    if (this.ghostMeshes.has(key))
                        continue;
                    const p = gridToWorld(col, row, this.layout);
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
            const bounds = estimateVisibleGridBounds(this.camera, this.layout);
            const minCol = Math.min(this.minCol, bounds.colMin);
            const maxCol = Math.max(this.maxCol, bounds.colMax);
            const minRow = Math.min(this.minRow, bounds.rowMin);
            const maxRow = Math.max(this.maxRow, bounds.rowMax);
            if (minCol !== this.minCol || maxCol !== this.maxCol || minRow !== this.minRow || maxRow !== this.maxRow) {
                this.expandGrid(minCol, maxCol, minRow, maxRow);
                return;
            }
            this._buildGrid(bounds.colMin, bounds.colMax, bounds.rowMin, bounds.rowMax);
        }
        _buildAxisLabels() {
            // We rely on HTML overlay labels in editor-dev.html; nothing to do here.
            // Stub kept for future canvas-based label rendering.
        }
        _createNodeLabelSprite(text) {
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
            window.addEventListener('mousemove', this._onMouseMove.bind(this));
            window.addEventListener('mouseup', this._onMouseUp.bind(this));
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
            if (!this.records.has(key))
                return;
            // Ctrl/Cmd+click: toggle multi-selection without disturbing the primary selection.
            if (e.ctrlKey || e.metaKey) {
                this._toggleMultiSelectKey(key);
                return;
            }
            // Plain click: clear any multi-selection extras and select this node.
            if (this.selectedKeys.size) {
                this.selectedKeys.forEach(k => {
                    if (k === key)
                        return;
                    const rec = this.records.get(k);
                    if (rec)
                        this._setRecordSelected(rec, false);
                });
                this.selectedKeys.clear();
            }
            this._setSelectedKey(key);
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
                // Skip drag init when ctrl/cmd is held — that gesture is for multi-select toggle.
                if (e.ctrlKey || e.metaKey)
                    return;
                if (key === this.selectedKey && this.records.has(key)) {
                    this.draggingNodeKey = key;
                    this.draggingPointer = { x: e.clientX, y: e.clientY };
                    this.hasDragged = false;
                    const startRec = this.records.get(key);
                    this.dragStartCell = startRec
                        ? { col: cell.col, row: cell.row, id: Number(startRec.node.id) }
                        : null;
                    e.preventDefault();
                }
            }
        }
        _onMouseMove(e) {
            if (this.draggingNodeKey) {
                this._syncDraggingNodeToPointer(e.clientX, e.clientY, false);
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
            const draggedKey = this.draggingNodeKey;
            const startCell = this.dragStartCell;
            if (this.draggingNodeKey)
                this.draggingNodeKey = null;
            this.draggingPointer = null;
            this.dragStartCell = null;
            if (this.isPanning) {
                this.isPanning = false;
                // Reset hasDragged on the next tick so click handler can check it first
                setTimeout(() => { this.hasDragged = false; }, 0);
                return;
            }
            if (wasDraggingNode) {
                // Emit a move-completed event only if the cell actually changed (skip no-op clicks).
                if (draggedKey && startCell && !this.suppressMoveCallback) {
                    const rec = this.records.get(draggedKey);
                    if (rec) {
                        const toCol = rec.node.col;
                        const toRow = rec.node.row;
                        if (toCol !== startCell.col || toRow !== startCell.row) {
                            this.onNodeMoved?.({
                                id: startCell.id,
                                from: [startCell.col, startCell.row],
                                to: [toCol, toRow]
                            });
                        }
                    }
                }
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
            if (this.draggingNodeKey && this.draggingPointer) {
                this._syncDraggingNodeToPointer(this.draggingPointer.x, this.draggingPointer.y, true);
            }
            this._ensureGridCoverage();
            this.renderer.render(this.scene, this.camera);
        }
    }
    scadnano.HoneycombEditor = HoneycombEditor;
    class SquareEditor extends HoneycombEditor {
        constructor(canvas, options = {}) {
            super(canvas, { ...options, layout: 'square' });
        }
    }
    scadnano.SquareEditor = SquareEditor;
})(scadnano || (scadnano = {}));
