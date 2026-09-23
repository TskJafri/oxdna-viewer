"use strict";
/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
var helix;
(function (helix_1) {
    const MAX_ANGULAR_SEPARATION = 40 * Math.PI / 180;
    const NOISE = -1;
    function hasInternalBasePair(helix) {
        const ids = new Set(helix.map(nt => nt.id));
        return helix.some(nt => nt.pair instanceof Nucleotide && ids.has(nt.pair.id));
    }
    function axisPoint(nt, helixIds) {
        const point = nt.getInstanceParameter3('bbOffsets').clone();
        if (nt.pair instanceof Nucleotide && helixIds.has(nt.pair.id)) {
            point.add(nt.pair.getInstanceParameter3('bbOffsets')).multiplyScalar(0.5);
        }
        return point;
    }
    function getHelixAxis(hx) {
        if (!Array.isArray(hx) || hx.length < 2)
            return null;
        const endpoints = toscad.helixEndpoints(hx);
        if (!endpoints || endpoints.end1.id === endpoints.end2.id)
            return null;
        const ids = new Set(hx.map(nt => nt.id));
        const start = axisPoint(endpoints.end1, ids);
        const end = axisPoint(endpoints.end2, ids);
        const axis = end.sub(start);
        const length = axis.length();
        if (!isFinite(length) || length === 0)
            return null;
        return axis.divideScalar(length);
    }
    function angularSeparation(a, b) {
        // A helix axis has no preferred endpoint order, so parallel and antiparallel
        // vectors describe the same orientation.
        const dot = Math.max(-1, Math.min(1, Math.abs(a.dot(b))));
        return Math.acos(dot);
    }
    function appendUnique(target, source) {
        const seen = new Set(target.map(nt => nt.id));
        for (const nt of source) {
            if (seen.has(nt.id))
                continue;
            seen.add(nt.id);
            target.push(nt);
        }
    }
    /**
     * Cluster helices by axis orientation using DBSCAN.
     *
     * Neighbours must be separated by at most 20 degrees. Binder helices are
     * excluded from the angular search and appended to their first backbone-
     * connected duplex parent. Noise and helices without a usable axis are
     * retained as singleton output helices.
     */
    function dbscan(helices, minPts = 1) {
        if (!Array.isArray(helices) || helices.length === 0)
            return [];
        minPts = Math.max(1, Math.floor(minPts));
        const binderIds = new Set();
        const axes = new Map();
        const nucleotideOwner = new Map();
        helices.forEach((hx, helixId) => {
            for (const nt of hx) {
                if (!nucleotideOwner.has(nt.id))
                    nucleotideOwner.set(nt.id, helixId);
            }
            if (!hasInternalBasePair(hx)) {
                binderIds.add(helixId);
                return;
            }
            const axis = getHelixAxis(hx);
            if (axis)
                axes.set(helixId, axis);
        });
        const candidates = Array.from(axes.keys());
        const labels = new Map();
        const visited = new Set();
        let clusterId = 0;
        const regionQuery = (helixId) => {
            const axis = axes.get(helixId);
            if (!axis)
                return [];
            return candidates.filter(otherId => {
                if (otherId === helixId)
                    return false;
                const otherAxis = axes.get(otherId);
                return !!otherAxis && angularSeparation(axis, otherAxis) <= MAX_ANGULAR_SEPARATION;
            });
        };
        for (const helixId of candidates) {
            if (visited.has(helixId))
                continue;
            visited.add(helixId);
            let neighbours = regionQuery(helixId);
            if (neighbours.length < minPts) {
                labels.set(helixId, NOISE);
                continue;
            }
            clusterId++;
            labels.set(helixId, clusterId);
            const queued = new Set(neighbours);
            for (let i = 0; i < neighbours.length; i++) {
                const neighbourId = neighbours[i];
                if (!visited.has(neighbourId)) {
                    visited.add(neighbourId);
                    const expanded = regionQuery(neighbourId);
                    if (expanded.length >= minPts) {
                        for (const expandedId of expanded) {
                            if (queued.has(expandedId))
                                continue;
                            queued.add(expandedId);
                            neighbours.push(expandedId);
                        }
                    }
                }
                if (!labels.has(neighbourId) || labels.get(neighbourId) === NOISE) {
                    labels.set(neighbourId, clusterId);
                }
            }
        }
        // Build merged helices only after DBSCAN finishes. This keeps every region
        // query based on the original endpoint-derived axes.
        const clustered = [];
        const clusterOutput = new Map();
        const helixOutput = new Map();
        helices.forEach((hx, helixId) => {
            if (binderIds.has(helixId))
                return;
            const label = labels.get(helixId);
            let outputId;
            if (label !== undefined && label !== NOISE && clusterOutput.has(label)) {
                outputId = clusterOutput.get(label);
                appendUnique(clustered[outputId], hx);
            }
            else {
                outputId = clustered.length;
                clustered.push([]);
                appendUnique(clustered[outputId], hx);
                if (label !== undefined && label !== NOISE)
                    clusterOutput.set(label, outputId);
            }
            helixOutput.set(helixId, outputId);
        });
        for (const binderId of binderIds) {
            let parentOutput;
            for (const nt of helices[binderId]) {
                for (const neighbour of [nt.n5, nt.n3]) {
                    if (!(neighbour instanceof Nucleotide))
                        continue;
                    const parentId = nucleotideOwner.get(neighbour.id);
                    if (parentId === undefined || parentId === binderId || binderIds.has(parentId))
                        continue;
                    parentOutput = helixOutput.get(parentId);
                    break;
                }
                if (parentOutput !== undefined)
                    break;
            }
            if (parentOutput === undefined) {
                parentOutput = clustered.length;
                clustered.push([]);
            }
            appendUnique(clustered[parentOutput], helices[binderId]);
        }
        return clustered;
    }
    helix_1.dbscan = dbscan;
})(helix || (helix = {}));
