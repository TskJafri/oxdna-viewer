/// <reference path="./typescript_definitions/index.d.ts" />
/// <reference path="./typescript_definitions/oxView.d.ts" />
/*
Hello my dear snooper, welcome to the source code for oxView!
I assume if you're reading this file, you're probably either a developer or looking for something in the code.
The main file here isn't super helpful at explaining what's going on, so I will try to give you a bit of a roadmap here.
main contains a few definitions of functions and data structures that the rest of the code uses.  Also a few functions that don't have a home yet.
The canvas, camera, renderer, periodic boundary condition handler and shader code can be found in the scene folder.
The file reading code, including the event listeners that handle drag/drop events are in file_handling
There you can also find the output options, including oxDNA files and videos.
Controls contains slightly modified stock Three.js control schemes.  These handle moving the camera and dragging objects.
Api and editing have most of the functions that let you control how things look and edit the actual structure.
Everything in the api can be called through the browser console by typing <apiName>.function(arguments) if you want to script some edits or visuals.
UI has all the functions relating to things that happen when you press buttons or hit the keyboard.
lib contains three.js and associated files.
typescript_definitions contains references between files to keep typescript editors happy.

If you add new files to your own copy of the viewer, you need to add it to tsconfig.json so the compiler knows to compile it.
The .js file will then appear in dist and you must add it to the script list at the bottom of index.html before it will take effect.

If you have any questions, feel free to open an issue on the GitHub page.
*/
// The ElementMap provies a mapping between particle ID in the simulation and JS objects here
class ElementMap extends Map {
    idCounter;
    constructor() {
        super();
        this.idCounter = 0;
    }
    // Avoid using this unless you really need to set
    // a specific id.
    set(id, element) {
        if (this.idCounter < id + 1) {
            this.idCounter = id + 1;
        }
        // Reading oxDNA files we set elements as undefined for
        // concurrency issues
        if (element) {
            element.id = id;
        }
        return super.set(id, element);
    }
    /**
     * Add an element, keeping track of
     * global id
     * @param element
     * @returns id
     */
    push(e) {
        e.id = this.idCounter++;
        super.set(e.id, e);
        return e.id;
    }
    /**
     * Remove element
     * @param id
     */
    delete(id) {
        // If we delete the last added, we can decrease the id counter.
        if (this.idCounter == id + 1) {
            this.idCounter = id;
        }
        return super.delete(id);
    }
    getNextId() {
        return this.idCounter;
    }
    reset() {
        this.clear();
        this.idCounter = 0;
    }
}
class smartSet extends Set {
    last;
    constructor() {
        super();
        this.last = undefined;
    }
    /**
     * Add an element to the set and remember the last added element
     * @param value
     * @returns void
     */
    add(value) {
        this.last = value;
        return super.add(value);
    }
    /**
     * Delete an element from the set and forget it if it was the last element added.
     * @param value
     * @returns
     */
    delete(value) {
        this.last = value;
        return super.delete(value);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////                  oxView's global variables                 ////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// Particle indexing stuff
const elements = new ElementMap(); //contains references to all BasicElements
const systems = []; // contains references to all systems
const selectedBases = new smartSet(); // contains the set of currently selected BasicElements
var clusterCounter = 0; //idk about this one...
// File reading stuff
var pdbtemp = []; // stores output from worker, so worker can terminate
const pdbFileInfo = []; //Stores all PDB Info (Necessary for future Protein Models)
const unfFileInfo = []; // Stores UNF file info (Necessary for writing out UNF files)
var confNum = 0; // Current configuration number in a trajectory
var box = new THREE.Vector3(); // Box size of the current scene
// BaseSelector stuff
var selectionMode = 'Monomer';
// ANM stuff
const networks = []; // Only used for networks, replaced anms
var selectednetwork = 0; // Only used for networks
const graphDatasets = []; // Only used for fluctuation graph
// Forces stuff
var forceHandler = new ForceHandler();
// color overlay stuff
var defaultColormap = "cooltowarm";
var lut, devs;
// Editing stuff
const editHistory = new EditHistory(); // Track do/undo
var tmpSystems = []; // Track memory for newly created systems
var topologyEdited = false; // to keep track of if the topology was edited at any point.
function resetScene(resetCamera = true) {
    elements.reset();
    while (systems.length > 0) {
        systems[systems.length - 1].backbone.parent.remove(systems[systems.length - 1].backbone);
        systems[systems.length - 1].nucleoside.parent.remove(systems[systems.length - 1].nucleoside);
        systems[systems.length - 1].connector.parent.remove(systems[systems.length - 1].connector);
        systems[systems.length - 1].bbconnector.parent.remove(systems[systems.length - 1].bbconnector);
        systems[systems.length - 1].dummyBackbone.parent.remove(systems[systems.length - 1].dummyBackbone);
        systems.pop();
    }
    while (tmpSystems.length > 0) {
        tmpSystems[tmpSystems.length - 1].backbone.parent.remove(tmpSystems[tmpSystems.length - 1].backbone);
        tmpSystems[tmpSystems.length - 1].nucleoside.parent.remove(tmpSystems[tmpSystems.length - 1].nucleoside);
        tmpSystems[tmpSystems.length - 1].connector.parent.remove(tmpSystems[tmpSystems.length - 1].connector);
        tmpSystems[tmpSystems.length - 1].bbconnector.parent.remove(tmpSystems[tmpSystems.length - 1].bbconnector);
        tmpSystems[tmpSystems.length - 1].dummyBackbone.parent.remove(tmpSystems[tmpSystems.length - 1].dummyBackbone);
        tmpSystems.pop();
    }
    selectedBases.clear();
    clusterCounter = 0; //idk about this one...
    // File reading stuff
    pdbtemp = []; // stores output from worker, so worker can terminate
    while (pdbFileInfo.length > 0) {
        pdbFileInfo.pop();
    }
    while (unfFileInfo.length > 0) {
        unfFileInfo.pop();
    }
    confNum = 0; // Current configuration number in a trajectory
    box = new THREE.Vector3(); // Box size of the current scene
    // BaseSelector stuff
    selectionMode = 'Monomer';
    // ANM stuff
    while (networks.length > 0) {
        networks.pop();
    }
    while (graphDatasets.length > 0) {
        graphDatasets.pop();
    }
    selectednetwork = 0; // Only used for networks
    // Forces stuff
    forceHandler.clearForcesFromScene();
    forceHandler = new ForceHandler();
    if (document.getElementById("forces")) {
        listForces();
    }
    // color overlay stuff
    defaultColormap = "cooltowarm";
    lut = [];
    devs = [];
    // Editing stuff
    editHistory.clear();
    tmpSystems = [];
    topologyEdited = false;
    if (resetCamera) {
        controls.reset();
    }
    render();
}
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////                       File input                           ////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//Check if there are files provided in the url (and load them if that is the case)
readFilesFromURLParams();
// define the drag and drop behavior of the scene
const target = renderer.domElement;
target.addEventListener("dragover", function (event) {
    event.preventDefault();
    target.classList.add('dragging');
}, false);
target.addEventListener("dragenter", function (event) {
    event.preventDefault();
    target.classList.add('dragging');
}, false);
target.addEventListener("dragexit", function (event) {
    event.preventDefault();
    target.classList.remove('dragging');
}, false);
// What to do if a file is dropped
target.addEventListener("drop", function (event) { event.preventDefault(); });
target.addEventListener("drop", handleDrop, false);
// Define message passing behavior
window.addEventListener("message", (event) => {
    if (event.data.message) { // do we have a message ?
        handleMessage(event.data);
    }
}, false);
render();
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////                      Random functions                      ////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
// These should probably be moved somewhere else...
function findBasepairs(min_length = 0) {
    systems.forEach(system => {
        if (!system.checkedForBasepairs) {
            system.strands.forEach(strand => {
                if (strand.getLength() >= min_length)
                    strand.forEach(e => {
                        if (e instanceof Nucleotide) {
                            if (!e.pair) {
                                e.pair = e.findPair();
                                if (e.pair) {
                                    e.pair.pair = e;
                                }
                            }
                        }
                    });
            });
        }
        system.checkedForBasepairs = true;
    });
}
;
// This one includes the A2 angle check. More complicated, but if there's cross-pairing or multi-pairing, then the code also checks A2 vectors, whose dot must be >0.85.
function findBasepairsOptimized(min_length = 0) {
    const CUTOFF_DIST = 0.65; // From your code
    const CELL_SIZE = 0.7; // Slightly larger than cutoff
    const CUTOFF_A1 = -0.85;
    const CUTOFF_A2 = 0.85;
    systems.forEach(system => {
        if (system.checkedForBasepairs)
            return;
        // 1. Flatten all valid nucleotides into a single list
        // This avoids nested looping through strands later
        let allNucs = [];
        system.strands.forEach(strand => {
            if (strand.getLength() >= min_length && strand.isNucleicAcid()) {
                strand.forEach(e => {
                    if (e instanceof Nucleotide)
                        allNucs.push(e);
                });
            }
        });
        // Recompute pairing from scratch for this system to avoid stale/asymmetric links.
        allNucs.forEach(n => {
            n.pair = null;
        });
        // 2. Build the Spatial Grid
        // Map Key: "x,y,z" coordinate of the cell
        // Map Value: Array of nucleotides in that cell
        const grid = new Map();
        // Helper to get grid key from position
        const getGridKey = (pos) => {
            const x = Math.floor(pos.x / CELL_SIZE);
            const y = Math.floor(pos.y / CELL_SIZE);
            const z = Math.floor(pos.z / CELL_SIZE);
            return `${x},${y},${z}`;
        };
        allNucs.forEach(n => {
            const pos = n.getInstanceParameter3("nsOffsets");
            n._cachedPos = pos; // Cache position to avoid re-fetching
            n._cachedKey = getGridKey(pos);
            if (!grid.has(n._cachedKey)) {
                grid.set(n._cachedKey, []);
            }
            grid.get(n._cachedKey).push(n);
        });
        // 3. Find Pairs using the Grid
        allNucs.forEach(curr => {
            if (curr.pair)
                return; // Already paired? Skip.
            let bestCandidate = null;
            let bestDist = CUTOFF_DIST;
            let bestOrient = 1;
            const EPS = 1e-9;
            const currPos = curr._cachedPos;
            // Calculate current cell coordinates
            const cx = Math.floor(currPos.x / CELL_SIZE);
            const cy = Math.floor(currPos.y / CELL_SIZE);
            const cz = Math.floor(currPos.z / CELL_SIZE);
            // 4. Iterate ONLY through neighbor cells (3x3x3 area)
            // This reduces checks from ~10,000 to ~20-50 per nucleotide
            for (let x = -1; x <= 1; x++) {
                for (let y = -1; y <= 1; y++) {
                    for (let z = -1; z <= 1; z++) {
                        const neighborKey = `${cx + x},${cy + y},${cz + z}`;
                        const cellNucs = grid.get(neighborKey);
                        if (!cellNucs)
                            continue;
                        // Check candidates in this cell
                        for (let other of cellNucs) {
                            if (curr === other)
                                continue; // Don't pair with self
                            // --- Original Logic from findPair() starts here ---
                            // 1. Topology Check (No neighbors)
                            if (curr.n3 === other || curr.n5 === other)
                                continue;
                            // 2. Complementary Rule Check
                            // (Combined your boolean logic for readability)
                            const typeSum = curr.getTypeNumber() + other.getTypeNumber();
                            const isWatsonCrick = (typeSum % 3 == 0) && (curr.getTypeNumber() !== other.getTypeNumber());
                            let isWobble = false;
                            if (curr.isRNA || other.isRNA) { // Assuming isRNA is on the nuc
                                const t1 = curr.type;
                                const t2 = other.type;
                                isWobble = (t1 == 'G' && t2 == 'U') || (t1 == 'U' && t2 == 'G');
                            }
                            if (isWatsonCrick || isWobble) {
                                // 3. Distance Check
                                const dist = other._cachedPos.distanceTo(currPos);
                                if (!(dist < CUTOFF_DIST))
                                    continue;
                                // 4. A1 Orientation gate
                                const orient = other.getA1().dot(curr.getA1());
                                if (orient < CUTOFF_A1) {
                                    // Base selection (no A2 competition unless this target is contested)
                                    const isBetterDist = dist < bestDist - EPS;
                                    const isSameDist = Math.abs(dist - bestDist) <= EPS;
                                    const isBetterOrientTieBreak = isSameDist && orient < bestOrient;
                                    if (isBetterDist || isBetterOrientTieBreak) {
                                        bestCandidate = other;
                                        bestDist = dist;
                                        bestOrient = orient;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // Apply the pair if found
            if (bestCandidate) {
                const incumbent = bestCandidate.pair;
                // No competition on target: pair directly.
                if (!incumbent || incumbent === curr) {
                    curr.pair = bestCandidate;
                    bestCandidate.pair = curr;
                }
                else {
                    // Competition case only:
                    // if multiple nucleotides target the same partner, require A2 cutoff
                    // and select the nucleotide with higher A2 dot product.
                    const currA2 = bestCandidate.getA2().dot(curr.getA2());
                    const incumbentA2 = bestCandidate.getA2().dot(incumbent.getA2());
                    const currPassA2 = currA2 > CUTOFF_A2;
                    const incumbentPassA2 = incumbentA2 > CUTOFF_A2;
                    let currWins = false;
                    if (currPassA2 && !incumbentPassA2) {
                        currWins = true;
                    }
                    else if (currPassA2 && incumbentPassA2) {
                        currWins = currA2 > incumbentA2 + EPS;
                    }
                    if (currWins) {
                        incumbent.pair = null;
                        curr.pair = bestCandidate;
                        bestCandidate.pair = curr;
                    }
                }
            }
        });
        // Cleanup temporary cache props if you want strictly clean objects
        allNucs.forEach(n => {
            delete n._cachedPos;
            delete n._cachedKey;
        });
        system.checkedForBasepairs = true;
    });
}
// Variant of findBasepairsOptimized:
// in conflict cases (multiple nucleotides competing for the same target),
// choose the challenger with the LOWEST RMSD score.
function findBasepairsOptim2(min_length = 0) {
    const CUTOFF_DIST = 0.65;
    const CELL_SIZE = 0.7;
    const CUTOFF_A1 = -0.85;
    // Set to false (or comment out) to exclude A2 term from RMSD scoring.
    // const USE_A2_IN_RMSD = false;
    systems.forEach(system => {
        if (system.checkedForBasepairs)
            return;
        let allNucs = [];
        system.strands.forEach(strand => {
            if (strand.getLength() >= min_length && strand.isNucleicAcid()) {
                strand.forEach(e => {
                    if (e instanceof Nucleotide)
                        allNucs.push(e);
                });
            }
        });
        allNucs.forEach(n => {
            n.pair = null;
        });
        const grid = new Map();
        const getGridKey = (pos) => {
            const x = Math.floor(pos.x / CELL_SIZE);
            const y = Math.floor(pos.y / CELL_SIZE);
            const z = Math.floor(pos.z / CELL_SIZE);
            return `${x},${y},${z}`;
        };
        allNucs.forEach(n => {
            const pos = n.getInstanceParameter3("nsOffsets");
            n._cachedPos = pos;
            n._cachedKey = getGridKey(pos);
            if (!grid.has(n._cachedKey)) {
                grid.set(n._cachedKey, []);
            }
            grid.get(n._cachedKey).push(n);
        });
        const computeCompetitionRmsd = (target, candidate) => {
            const dist = target.getInstanceParameter3("nsOffsets").distanceTo(candidate.getInstanceParameter3("nsOffsets"));
            const a1Dot = target.getA1().dot(candidate.getA1());
            // const a2Dot = target.getA2().dot(candidate.getA2());
            // Convert each metric into an error-like term where lower is better.
            const distTerm = dist / CUTOFF_DIST;
            const a1Term = 1 + a1Dot; // ideal antiparallel a1Dot ~ -1 => 0
            const terms = [distTerm, a1Term];
            // if (USE_A2_IN_RMSD) {
            //     const a2Term = 1 - a2Dot; // ideal parallel a2Dot ~ 1 => 0
            //     terms.push(a2Term);
            // }
            const sumSq = terms.reduce((acc, t) => acc + t * t, 0);
            return Math.sqrt(sumSq / terms.length);
        };
        allNucs.forEach(curr => {
            if (curr.pair)
                return;
            let bestCandidate = null;
            let bestDist = CUTOFF_DIST;
            let bestOrient = 1;
            const EPS = 1e-9;
            const currPos = curr._cachedPos;
            const cx = Math.floor(currPos.x / CELL_SIZE);
            const cy = Math.floor(currPos.y / CELL_SIZE);
            const cz = Math.floor(currPos.z / CELL_SIZE);
            for (let x = -1; x <= 1; x++) {
                for (let y = -1; y <= 1; y++) {
                    for (let z = -1; z <= 1; z++) {
                        const neighborKey = `${cx + x},${cy + y},${cz + z}`;
                        const cellNucs = grid.get(neighborKey);
                        if (!cellNucs)
                            continue;
                        for (let other of cellNucs) {
                            if (curr === other)
                                continue;
                            if (curr.n3 === other || curr.n5 === other)
                                continue;
                            const typeSum = curr.getTypeNumber() + other.getTypeNumber();
                            const isWatsonCrick = (typeSum % 3 == 0) && (curr.getTypeNumber() !== other.getTypeNumber());
                            let isWobble = false;
                            if (curr.isRNA || other.isRNA) {
                                const t1 = curr.type;
                                const t2 = other.type;
                                isWobble = (t1 == 'G' && t2 == 'U') || (t1 == 'U' && t2 == 'G');
                            }
                            if (isWatsonCrick || isWobble) {
                                const dist = other._cachedPos.distanceTo(currPos);
                                if (!(dist < CUTOFF_DIST))
                                    continue;
                                const orient = other.getA1().dot(curr.getA1());
                                if (orient < CUTOFF_A1) {
                                    const isBetterDist = dist < bestDist - EPS;
                                    const isSameDist = Math.abs(dist - bestDist) <= EPS;
                                    const isBetterOrientTieBreak = isSameDist && orient < bestOrient;
                                    if (isBetterDist || isBetterOrientTieBreak) {
                                        bestCandidate = other;
                                        bestDist = dist;
                                        bestOrient = orient;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (bestCandidate) {
                const incumbent = bestCandidate.pair;
                if (!incumbent || incumbent === curr) {
                    curr.pair = bestCandidate;
                    bestCandidate.pair = curr;
                }
                else {
                    // Competition case: lower RMSD wins.
                    const currRmsd = computeCompetitionRmsd(bestCandidate, curr);
                    const incumbentRmsd = computeCompetitionRmsd(bestCandidate, incumbent);
                    if (currRmsd + EPS < incumbentRmsd) {
                        incumbent.pair = null;
                        curr.pair = bestCandidate;
                        bestCandidate.pair = curr;
                    }
                }
            }
        });
        allNucs.forEach(n => {
            delete n._cachedPos;
            delete n._cachedKey;
        });
        system.checkedForBasepairs = true;
    });
}
function findBasepairsOptimwA2(min_length = 0) {
    const CUTOFF_DIST = 0.65;
    const CELL_SIZE = 0.7;
    const CUTOFF_A1 = -0.85;
    // Set to false (or comment out) to exclude A2 term from RMSD scoring.
    const USE_A2_IN_RMSD = false;
    systems.forEach(system => {
        if (system.checkedForBasepairs)
            return;
        let allNucs = [];
        system.strands.forEach(strand => {
            if (strand.getLength() >= min_length && strand.isNucleicAcid()) {
                strand.forEach(e => {
                    if (e instanceof Nucleotide)
                        allNucs.push(e);
                });
            }
        });
        allNucs.forEach(n => {
            n.pair = null;
        });
        const grid = new Map();
        const getGridKey = (pos) => {
            const x = Math.floor(pos.x / CELL_SIZE);
            const y = Math.floor(pos.y / CELL_SIZE);
            const z = Math.floor(pos.z / CELL_SIZE);
            return `${x},${y},${z}`;
        };
        allNucs.forEach(n => {
            const pos = n.getInstanceParameter3("nsOffsets");
            n._cachedPos = pos;
            n._cachedKey = getGridKey(pos);
            if (!grid.has(n._cachedKey)) {
                grid.set(n._cachedKey, []);
            }
            grid.get(n._cachedKey).push(n);
        });
        const computeCompetitionRmsd = (target, candidate) => {
            const dist = target.getInstanceParameter3("nsOffsets").distanceTo(candidate.getInstanceParameter3("nsOffsets"));
            const a1Dot = target.getA1().dot(candidate.getA1());
            const a2Dot = target.getA2().dot(candidate.getA2());
            // Convert each metric into an error-like term where lower is better.
            // const distTerm = dist / CUTOFF_DIST;
            const distTerm = 1;
            // const a1Term = 1 + a1Dot; // ideal antiparallel a1Dot ~ -1 => 0
            const a1Term = 1;
            const terms = [distTerm, a1Term];
            if (USE_A2_IN_RMSD) {
                const a2Term = 1 - a2Dot; // ideal parallel a2Dot ~ 1 => 0
                terms.push(a2Term);
            }
            const sumSq = terms.reduce((acc, t) => acc + t * t, 0);
            return Math.sqrt(sumSq / terms.length);
        };
        allNucs.forEach(curr => {
            if (curr.pair)
                return;
            let bestCandidate = null;
            let bestDist = CUTOFF_DIST;
            let bestOrient = 1;
            const EPS = 1e-9;
            const currPos = curr._cachedPos;
            const cx = Math.floor(currPos.x / CELL_SIZE);
            const cy = Math.floor(currPos.y / CELL_SIZE);
            const cz = Math.floor(currPos.z / CELL_SIZE);
            for (let x = -1; x <= 1; x++) {
                for (let y = -1; y <= 1; y++) {
                    for (let z = -1; z <= 1; z++) {
                        const neighborKey = `${cx + x},${cy + y},${cz + z}`;
                        const cellNucs = grid.get(neighborKey);
                        if (!cellNucs)
                            continue;
                        for (let other of cellNucs) {
                            if (curr === other)
                                continue;
                            if (curr.n3 === other || curr.n5 === other)
                                continue;
                            const typeSum = curr.getTypeNumber() + other.getTypeNumber();
                            const isWatsonCrick = (typeSum % 3 == 0) && (curr.getTypeNumber() !== other.getTypeNumber());
                            let isWobble = false;
                            if (curr.isRNA || other.isRNA) {
                                const t1 = curr.type;
                                const t2 = other.type;
                                isWobble = (t1 == 'G' && t2 == 'U') || (t1 == 'U' && t2 == 'G');
                            }
                            if (isWatsonCrick || isWobble) {
                                const dist = other._cachedPos.distanceTo(currPos);
                                if (!(dist < CUTOFF_DIST))
                                    continue;
                                const orient = other.getA1().dot(curr.getA1());
                                if (orient < CUTOFF_A1) {
                                    const isBetterDist = dist < bestDist - EPS;
                                    const isSameDist = Math.abs(dist - bestDist) <= EPS;
                                    const isBetterOrientTieBreak = isSameDist && orient < bestOrient;
                                    if (isBetterDist || isBetterOrientTieBreak) {
                                        bestCandidate = other;
                                        bestDist = dist;
                                        bestOrient = orient;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (bestCandidate) {
                const incumbent = bestCandidate.pair;
                if (!incumbent || incumbent === curr) {
                    curr.pair = bestCandidate;
                    bestCandidate.pair = curr;
                }
                else {
                    // Competition case: lower RMSD wins.
                    const currRmsd = computeCompetitionRmsd(bestCandidate, curr);
                    const incumbentRmsd = computeCompetitionRmsd(bestCandidate, incumbent);
                    if (currRmsd + EPS < incumbentRmsd) {
                        incumbent.pair = null;
                        curr.pair = bestCandidate;
                        bestCandidate.pair = curr;
                    }
                }
            }
        });
        allNucs.forEach(n => {
            delete n._cachedPos;
            delete n._cachedKey;
        });
        system.checkedForBasepairs = true;
    });
}
// Utility function to pick a random element from list
function randomChoice(l) {
    return l[Math.floor(Math.random() * l.length)];
}
function findBasepairsOrigami(min_length = 1000) {
    findBasepairs(min_length);
}
// Ugly hacks for testing
function getElements() {
    return elements;
}
function getSystems() {
    return systems;
}
//Temporary solution to adding configuration storage
//This section sets interface values from the storage 
if (window.sessionStorage.centerOption) {
    view.centeringMode.set(window.sessionStorage.centerOption);
}
if (window.sessionStorage.inboxingOption) {
    view.inboxingMode.set(window.sessionStorage.inboxingOption);
}
//https://stackoverflow.com/questions/326069/how-to-identify-if-a-webpage-is-being-loaded-inside-an-iframe-or-directly-into-t
function inIframe() {
    try {
        return window.self !== window.top;
    }
    catch (e) {
        return true;
    }
}
