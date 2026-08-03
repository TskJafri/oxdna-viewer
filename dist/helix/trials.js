"use strict";
/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
/*
Developer-console utilities for ad-hoc experiments / measurements on a loaded structure.

Usage (paste into the browser dev console):

    trials.crossovers();             // log + download crossovers.csv
    trials.crossovers({ download: false });   // just return the CSV string

The crossover counts are directional in the sense that we walk every strand
5'->3' and credit the helix we leave (parent) -> the helix we enter (target).
A symmetric crossover therefore contributes to BOTH (A,B) and (B,A).
*/
var trials;
(function (trials) {
    // ── Internal helpers ───────────────────────────────────────────────────
    function nucleotideOnly() {
        const out = new Map();
        elements.forEach((element, id) => {
            if (element instanceof Nucleotide)
                out.set(id, element);
        });
        return out;
    }
    function binCrossovers(grid, helixCount) {
        const bins = new Map();
        const helixIds = new Set();
        for (let i = 0; i < helixCount; i++)
            helixIds.add(i);
        const bump = (from, to) => {
            helixIds.add(from);
            helixIds.add(to);
            let inner = bins.get(from);
            if (!inner) {
                inner = new Map();
                bins.set(from, inner);
            }
            inner.set(to, (inner.get(to) ?? 0) + 1);
        };
        for (const c of toscad.crossoverNts(grid)) {
            if (c.fromHelix === c.toHelix)
                continue;
            bump(c.fromHelix, c.toHelix);
        }
        return { bins, helixIds };
    }
    function emitCsv(bins, helixIds, opts) {
        const sortedIds = Array.from(helixIds).sort((a, b) => a - b);
        const rows = ['parent_helix_id,target_helix_id,crossovers'];
        for (const parent of sortedIds) {
            const inner = bins.get(parent);
            if (!inner)
                continue;
            const targets = Array.from(inner.keys()).sort((a, b) => a - b);
            for (const target of targets) {
                const count = inner.get(target) ?? 0;
                if (!count)
                    continue;
                rows.push(`${parent},${target},${count}`);
            }
        }
        const csv = rows.join('\n');
        if (opts.log) {
            console.log(`[${opts.tag}] ${opts.helixCount} helices, ${rows.length - 1} non-zero bins`);
            console.log(csv);
        }
        if (opts.download) {
            try {
                makeTextFile(opts.filename, csv);
            }
            catch (e) {
                console.warn(`[${opts.tag}] Could not trigger download:`, e);
            }
        }
        return csv;
    }
    // ── Public API ─────────────────────────────────────────────────────────
    function crossovers(options = {}) {
        const { tolerance = 3, download = true, filename = 'crossovers.csv', log = true, } = options;
        const result = helix.findHelices(nucleotideOnly(), tolerance);
        const helices = result?.helices ?? [];
        if (!helices.length) {
            console.warn('[trials.crossovers] No helices found.');
            return '';
        }
        const { grid } = toscad.setGrid(helices);
        const { bins, helixIds } = binCrossovers(grid, helices.length);
        return emitCsv(bins, helixIds, {
            filename, download, log,
            tag: 'trials.crossovers',
            helixCount: helices.length,
        });
    }
    trials.crossovers = crossovers;
})(trials || (trials = {}));
