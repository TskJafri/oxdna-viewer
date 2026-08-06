/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
import THREE from "../typescript_definitions";
// Helper function to show the endpoints of all helices in the console. Useful for debugging.
export function showHelixEndpoints(helices) {
    // const helices = await helix.findHelices(elements, 2);
    const endpoints = helices.map((helix, i) => {
        const res = toscad.helixEndpoints(helix);
        return {
            helixIndex: i,
            endpointA: res?.end1?.id,
            endpointB: res?.end2?.id,
            diameter: res?.diameter
        };
    });
    console.log(endpoints);
    return endpoints;
}
;
// Just for a visualization and good only for debugging...
export function addPartialAxisToScene(d) {
    const { planeVector } = helix.getPartialAxis(d);
    const origin = d.start1.getInstanceParameter3('bbOffsets')
        .add(d.end2.getInstanceParameter3('bbOffsets'))
        .multiplyScalar(0.5);
    if (typeof THREE !== 'undefined' && typeof scene !== 'undefined' && scene?.add) {
        const arrow = new THREE.ArrowHelper(planeVector.clone().normalize(), origin, 10);
        scene.add(arrow);
    }
}
// helper function for adding the average a3 vector in canvas. Really should not be here.
export function averageA3a(list) {
    if (!list.length)
        return new THREE.Vector3(0, 0, 0);
    // Align all A3 vectors so they point in a consistent direction before averaging.
    const ref = list[0].getA3().clone().normalize();
    const acc = ref.clone();
    for (let i = 1; i < list.length; i++) {
        const v = list[i].getA3().clone().normalize();
        acc.add(v.dot(ref) < 0 ? v.multiplyScalar(-1) : v);
    }
    acc.divideScalar(list.length);
    const avg = acc.normalize();
    // Visualize the averaged orientation from the first nucleotide origin when possible.
    const origin = list[0]?.getPos();
    if (origin && typeof THREE !== 'undefined' && typeof scene !== 'undefined' && scene?.add) {
        const helper = new THREE.ArrowHelper(avg.clone(), origin, 5);
        scene.add(helper);
    }
    return avg;
}
;
