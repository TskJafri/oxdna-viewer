// If you see this, ignore it.
// Its for fun only.
// Look in funny.html for the graph

let helices = await honda.findHelices(elements, 2);
let { grid, binderHelices } = toscad.setGrid(helices);
let { crossovers, helixIds } = toscad.collectCrossovers(grid);

// 1. Generate Nodes from helixIds
const nodes = [];
helixIds.forEach(helixId => {
    nodes.push({
        data: { id: helixId.toString(), label: `Helix ${helixId}` }
    });
});

// 2. Generate Edges from crossovers
const edges = [];
// Assuming crossovers is a Map<number, Map<number, { sameWalk: number; diffWalk: number }>> 
// sourceHelix -> destinationHelix -> connection stats
crossovers.forEach((destinations, sourceId) => {
    destinations.forEach((stats, destId) => {
        // Use the total connections as the weight of the edge
        const totalConnections = stats.sameWalk + stats.diffWalk;
        if (totalConnections > 0) {
            edges.push({
                data: {
                    id: `e_${sourceId}_${destId}`,
                    source: sourceId.toString(),
                    target: destId.toString(),
                    weight: totalConnections
                }
            });
        }
    });
});

// 3. Assemble and log the format
const graphData = {
    elements: {
        nodes: nodes,
        edges: edges
    }
};

console.log("Graph Data JSON: ");
console.log(JSON.stringify(graphData, null, 2));


/*
// Simply because I dont want to type this every time

let helix1_center1 = elements.get(12730).getPos(); let helix1_center2 = elements.get(1136).getPos();
let h1_center = new THREE.Vector3();
h1_center.addVectors(helix1_center1, helix1_center2).divideScalar(2);
console.log(h1_center);

let helix2_center1 = elements.get(7152).getPos(); let helix2_center2 = elements.get(11530).getPos();
let h2_center = new THREE.Vector3();
h2_center.addVectors(helix2_center1, helix2_center2).divideScalar(2);
console.log(h2_center);

let geometry = new THREE.SphereGeometry(0.5, 5,5);
let material = new THREE.MeshBasicMaterial( { color: 0xffff00 } );
let sphere = new THREE.Mesh( geometry, material );
sphere.position.copy(h3_center)
scene.add(sphere)

let a1 = new THREE.Vector3().subVectors(h2_center, h1_center);
let length = a1.length();
a1.normalize();
let arrow = new THREE.ArrowHelper(a1, h1_center, length);
scene.add(arrow);
*/

/*
let helices = await honda.findHelices(elements, 2);
let {grid, binderHelices} = toscad.setGrid(helices);
console.log(grid);
console.log(toscad.validateGrid(grid));
toscad.directionAlign2(grid);
// toscad.gridFlip(grid,33);
// toscad.gridFlip(grid,35);
// toscad.gridFlip(grid,36);
// toscad.gridFlip(grid,37);
// toscad.gridFlip(grid,38);
// toscad.gridFlip(grid,39);
// toscad.gridFlip(grid,40);
// toscad.gridFlip(grid,41);
// toscad.gridFlip(grid,42);
// toscad.gridFlip(grid,43);
toscad.alignGridPrim(grid, binderHelices);
let scadnano = toscad.buildScadnano2(grid,helices);
// copy(scadnano);

honda.findHelices(elements, 2).then(helices => {
    let positionsMap = toscad.getRelativePositions(helices);
    let h0_nodes = helices[0];
    let anchorZ = 0;
    h0_nodes.forEach(n => anchorZ += n.getPos().z);
    anchorZ /= h0_nodes.length;
    console.log("--- 2D Relative Positions ---");
    for (const [id, pos] of positionsMap.entries()) {
        console.log(`Helix ${id}: { x: ${pos.x.toFixed(2)}, y: ${pos.y.toFixed(2)} } nm`);

        let geometry = new THREE.SphereGeometry(0.5, 32, 16);
        let material = new THREE.MeshBasicMaterial({ color: id === 0 ? 0x00ff00 : 0xffff00 });
        let sphere = new THREE.Mesh(geometry, material);
        sphere.position.set(pos.x + 20, pos.y, anchorZ);
        scene.add(sphere);
    }

    render();
    console.log(`Drew ${positionsMap.size} flat 2D helix positions into the scene!`);
});
*/