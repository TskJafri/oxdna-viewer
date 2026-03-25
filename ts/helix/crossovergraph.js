async function crossoverToCSV() {
    const helices = await honda.findHelices(elements, 2);

    const { grid, binderHelices } = toscad.setGrid(helices);
    toscad.directionAlign2(grid);
    toscad.alignGridPrim(grid, binderHelices);
    toscad.combinedHelices(15, grid, helices, binderHelices);

    let scaffoldId = honda.getScaffoldStrand();

    const cross = toscad.crossoverNts(grid);
    // Make a set because we don't want to deal with duplicates (it does not matter if you use arrays).
    const cross_set = new Set(cross
        .flatMap(item => [item.fromNt, item.toNt])
        .filter(nt => nt.strand !== scaffoldId));
    const cross_set_noScaff = new Set(cross
        .flatMap(item => [item.fromNt, item.toNt])
        .filter(nt => nt.strand !== scaffoldId));
    const cross_set_onlyScaff = new Set(cross
        .flatMap(item => [item.fromNt, item.toNt])
        .filter(nt => nt.strand === scaffoldId));
    // const cross_ids = new Set(cross.flatMap(item => [item.fromNt.id, item.toNt.id]));

    const ntToNeighbor = new Map();
    cross.forEach(c => {
        // For the nucleotide exiting a helix, the neighbor is the destination
        ntToNeighbor.set(c.fromNt.id, c.toHelix);
        // For the nucleotide entering a helix, the neighbor is the source
        ntToNeighbor.set(c.toNt.id, c.fromHelix);
    });

    // const visited = new Set();
    const ntToHelix = new Map();
    helices.forEach((helix, hId) => {
        helix.forEach(nt => ntToHelix.set(nt.id, hId));
    });

    const getCrossoverNtOnHelix = (nt, hid) => {
        if (nt && cross_set.has(nt) && ntToHelix.get(nt.id) === hid) {
            return nt;
        }
        if (nt?.pair && cross_set.has(nt.pair) && ntToHelix.get(nt.pair.id) === hid) {
            return nt.pair;
        }
        return null;
    };

    // constants for building the csv later
    const helix_id = [];
    const connected_nt_1 = [];
    const connected_nt_2 = [];
    const distance_nt = [];
    const direction = [];

    // fill the columns first.
    cross_set.forEach(startNt => {
        let cur = startNt;
        let hid = ntToHelix.get(cur.id);
        let requiredNeighbor = ntToNeighbor.get(startNt.id);
        let dist = 1;
        let dir = 'n3' || 'n5';
        let foundmatch = false;
        let matchNt = null;

        while (cur && foundmatch === false) {
            if (!cur.n3) {
                let jump = cur.pair?.n5?.pair;
                if (jump && ntToHelix.get(jump.id) === hid) {
                    cur = jump;
                    dist++;
                }
                else break;
            }
            else {
                if (ntToHelix.get(cur?.n3?.id) !== hid) break;
            }
            cur = cur.n3;
            dist++;
            dir = 'n3';

            const crossoverNt = getCrossoverNtOnHelix(cur, hid);
            if (crossoverNt) {
                if (ntToNeighbor.get(crossoverNt.id) === requiredNeighbor) {
                    foundmatch = true;
                    matchNt = crossoverNt;
                }
                else {
                    const skip = crossoverNt.pair?.n5?.n5?.pair;
                    if (!skip || ntToHelix.get(skip.id) !== hid) {
                        break;
                    }
                    cur = skip;
                    dist++;
                }
            }
        }

        if (foundmatch === false) {
            cur = startNt;
            dist = 1;
        }

        while (cur && foundmatch === false) {
            if (!cur.n5) {
                let jump = cur.pair?.n3?.pair;
                if (jump && ntToHelix.get(jump.id) === hid) {
                    cur = jump;
                    dist++;
                }
                else break;
            }
            else {
                if (ntToHelix.get(cur?.n5?.id) !== hid) break;
            };
            cur = cur.n5;
            dist++;
            dir = 'n5';

            const crossoverNt = getCrossoverNtOnHelix(cur, hid);
            if (crossoverNt) {
                if (ntToNeighbor.get(crossoverNt.id) === requiredNeighbor) {
                    foundmatch = true;
                    matchNt = crossoverNt;
                }
                else {
                    const skip = crossoverNt.pair?.n3?.n3?.pair;
                    if (!skip || ntToHelix.get(skip.id) !== hid) {
                        break;
                    }
                    cur = skip;
                    dist++;
                }
            }
        }
        if (foundmatch) {
            // if (ntToNeighbor.get(startNt.id) === ntToNeighbor.get(matchNt.id)) {
                    helix_id.push(hid);
                    connected_nt_1.push(startNt.id);
                    connected_nt_2.push(matchNt.id);
                    distance_nt.push(dist);
                    direction.push(dir);
                // }
            }
    });

    const rows = [];
    let i = 0;

    while (i < helix_id.length) {
        const rowString = `${helix_id[i]},${connected_nt_1[i]},${connected_nt_2[i]},${distance_nt[i]}`;
        rows.push(rowString);
        i++;
    }

    rows.sort((rowA, rowB) => {
        let colA = rowA.split(',');
        let val1a = colA[0];
        let val2a = colA[1];
        let val3a = colA[2];
        let val4a = colA[3];
        let colB = rowB.split(',');
        let val1b = colB[0];
        let val2b = colB[1];
        let val3b = colB[2];
        let val4b = colB[3];
        return Number(val1a) - Number(val1b) || Number(val4a) - Number(val4b) || Number(val2a) - Number(val2b) || Number(val3a) - Number(val3b);
    });
    rows.unshift("helix_id,connected_nt_1,connected_nt_2,distance_nt");

    const finalCsvString = rows.join('\n');
    console.log(finalCsvString);
}