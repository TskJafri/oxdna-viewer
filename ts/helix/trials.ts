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

namespace trials {

    export interface CrossoverOptions {
        tolerance?: number;   // passed straight to helix.findHelices()
        download?: boolean;   // trigger a CSV download
        filename?: string;    // CSV filename when download is true
        log?: boolean;        // log the CSV to the console
    }

    // ── Internal helpers ───────────────────────────────────────────────────

    function nucleotideOnly(): Map<number, Nucleotide> {
        const out = new Map<number, Nucleotide>();
        elements.forEach((element, id) => {
            if (element instanceof Nucleotide) out.set(id, element);
        });
        return out;
    }

    function binCrossovers(grid: toscad.GridMap, helixCount: number): {
        bins: Map<number, Map<number, number>>;
        helixIds: Set<number>;
    } {
        const bins = new Map<number, Map<number, number>>();
        const helixIds = new Set<number>();
        for (let i = 0; i < helixCount; i++) helixIds.add(i);

        const bump = (from: number, to: number) => {
            helixIds.add(from);
            helixIds.add(to);
            let inner = bins.get(from);
            if (!inner) {
                inner = new Map<number, number>();
                bins.set(from, inner);
            }
            inner.set(to, (inner.get(to) ?? 0) + 1);
        };

        for (const c of toscad.crossoverNts(grid)) {
            if (c.fromHelix === c.toHelix) continue;
            bump(c.fromHelix, c.toHelix);
        }
        return { bins, helixIds };
    }

    function emitCsv(
        bins: Map<number, Map<number, number>>,
        helixIds: Set<number>,
        opts: { filename: string; download: boolean; log: boolean; tag: string; helixCount: number }
    ): string {
        const sortedIds = Array.from(helixIds).sort((a, b) => a - b);
        const rows: string[] = ['parent_helix_id,target_helix_id,crossovers'];
        for (const parent of sortedIds) {
            const inner = bins.get(parent);
            if (!inner) continue;
            const targets = Array.from(inner.keys()).sort((a, b) => a - b);
            for (const target of targets) {
                const count = inner.get(target) ?? 0;
                if (!count) continue;
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
            } catch (e) {
                console.warn(`[${opts.tag}] Could not trigger download:`, e);
            }
        }
        return csv;
    }

    // ── Public API ─────────────────────────────────────────────────────────

    export function crossovers(options: CrossoverOptions = {}): string {
        const {
            tolerance = 3,
            download = true,
            filename = 'crossovers.csv',
            log = true,
        } = options;

        const result = helix.findHelices(nucleotideOnly(), tolerance) as { helices: Nucleotide[][] };
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
}
