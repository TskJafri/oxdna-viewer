// First analysis: minimum distances between the centers of every helices.
const helixCOM = (helix) => {
    if (!helix || helix.length === 0) return null;
    const sum = new THREE.Vector3();
    let count = 0;
    for (const nt of helix) {
        if (!(nt instanceof Nucleotide)) continue;
        sum.add(nt.getPos());
        count++;
    }
    if (count === 0) return null;
    return sum.divideScalar(count);
};

// Creates a coordinate system based on 
// input: helices[][]
const estimateHexBasis3D = (allHelices, grid) => {
    const centers = [];
    const axisSamples = [];

    // "PCA" is basically using endpoints 1 and 2 as a rough axis direction, and then averaging the centroids to find the origin.
    for (const helix of allHelices) {
        if (!helix || helix.length === 0) continue;
        const center = helixCOM(helix);
        if (center) centers.push(center);

        const ep = toscad.helixEndpoints(helix);
        if (!ep) continue;
        const axis = ep.end2.getPos().clone().sub(ep.end1.getPos());
        if (axis.lengthSq() > 1e-8) axisSamples.push(axis.normalize());
    }

    if (centers.length === 0) return null;

    // find the average center to use as the origin
    const origin = new THREE.Vector3();
    for (const c of centers) origin.add(c);
    origin.divideScalar(centers.length);

    // Average helix direction to get a rough axis. Not really for any rigorous reason, it works fine for this and thus, we will use it.
    let axisVec = new THREE.Vector3(0, 0, 1);
    if (axisSamples.length > 0) {
        axisVec.set(0, 0, 0);
        for (const s of axisSamples) axisVec.add(s);
        if (axisVec.lengthSq() < 1e-8) axisVec.set(0, 0, 1);
        else axisVec.normalize();
    }

    // smart and concise way of finding the perpendicular direction of the spread...
    // removes the entire axisVec vector from the axes, so only the perpendicular spread remains.
    const projectedCenters = centers.map((c) => {
        const rel = c.clone().sub(origin);
        // normalized already, so no need to divide by axisVec.lengthSq()
        return rel.sub(axisVec.clone().multiplyScalar(rel.dot(axisVec)));
    });

    const {crossovers, helixIds} = toscad.collectCrossovers(grid); 

    // Here's what i need now that the coordinate grid has technically been established:
    // For every helices' connections, find the distance between those centers, and then plot the output as a continuous graph. See what exists.
};



//     let e1 = new THREE.Vector3(1, 0, 0);
//     let bestD2 = Infinity;
//     for (let i = 0; i < projectedCenters.length; i++) {
//         for (let j = i + 1; j < projectedCenters.length; j++) {
//             const d = projectedCenters[j].clone().sub(projectedCenters[i]);
//             const d2 = d.lengthSq();
//             if (d2 > 1e-8 && d2 < bestD2) {
//                 bestD2 = d2;
//                 e1 = d.normalize();
//             }
//         }
//     }

//     let e2 = axisVec.clone().cross(e1);
//     if (e2.lengthSq() < 1e-8) {
//         const fallback = Math.abs(axisVec.x) < 0.9
//             ? new THREE.Vector3(1, 0, 0)
//             : new THREE.Vector3(0, 1, 0);
//         e2 = axisVec.clone().cross(fallback);
//     }
//     e2.normalize();
//     e1 = e2.clone().cross(axisVec).normalize();

//     const nearestDistances = [];
//     for (let i = 0; i < projectedCenters.length; i++) {
//         let nearest = Infinity;
//         for (let j = 0; j < projectedCenters.length; j++) {
//             if (i === j) continue;
//             const d = projectedCenters[j].clone().sub(projectedCenters[i]).length();
//             if (d > 1e-6 && d < nearest) nearest = d;
//         }
//         if (nearest < Infinity) nearestDistances.push(nearest);
//     }

//     let spacing = 1;
//     if (nearestDistances.length > 0) {
//         nearestDistances.sort((a, b) => a - b);
//         spacing = nearestDistances[Math.floor(nearestDistances.length / 2)] || 1;
//         if (spacing <= 1e-6) spacing = 1;
//     }

//     const qVec = e1.clone().multiplyScalar(spacing);
//     const rVec = e1.clone().multiplyScalar(0.5 * spacing)
//         .add(e2.clone().multiplyScalar((Math.sqrt(3) / 2) * spacing));

//     return { origin, qVec, rVec, axisVec };
// };
