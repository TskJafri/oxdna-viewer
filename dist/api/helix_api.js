"use strict";
/// <reference path="../typescript_definitions/index.d.ts" />
/// <reference path="../typescript_definitions/oxView.d.ts" />
/// <reference path="../main.ts" />
// Alias the global `helix` namespace before entering `api.helix`, where the
// bare name would shadow it.
const helixNs = helix;
var api;
(function (api) {
    var helix;
    (function (helix_1) {
        // Helper function to show the endpoints of all helices in the console. Useful for debugging.
        function showHelixEndpoints(helices) {
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
        helix_1.showHelixEndpoints = showHelixEndpoints;
        ;
        // Just for a visualization and good only for debugging...
        function addPartialAxisToScene(d) {
            const { planeVector } = helixNs.getPartialAxis(d);
            const origin = d.start1.getInstanceParameter3('bbOffsets')
                .add(d.end2.getInstanceParameter3('bbOffsets'))
                .multiplyScalar(0.5);
            if (typeof THREE !== 'undefined' && typeof scene !== 'undefined' && scene?.add) {
                const arrow = new THREE.ArrowHelper(planeVector.clone().normalize(), origin, 10);
                scene.add(arrow);
            }
        }
        helix_1.addPartialAxisToScene = addPartialAxisToScene;
        // helper function for adding the average a3 vector in canvas. Really should not be here.
        function averageA3a(list) {
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
        helix_1.averageA3a = averageA3a;
        ;
        function findHelixID(targetId, helices) {
            for (let i = 0; i < helices.length; i++) {
                const helix = helices[i];
                for (const nt of helix) {
                    if (nt.id === targetId)
                        return i;
                }
            }
            return null;
        }
        helix_1.findHelixID = findHelixID;
        function positionsToJSON(positions) {
            return JSON.stringify([...positions.entries()].sort((a, b) => a[0] - b[0]));
        }
        helix_1.positionsToJSON = positionsToJSON;
    })(helix = api.helix || (api.helix = {}));
})(api || (api = {}));
