/// <reference path="../typescript_definitions/index.d.ts" />

function _handleScadnanoDrop(data) {
    if (data.action === 'initialize') {
        camera.up.multiplyScalar(-1);
        scadnanoOrigin = data.scadnano_origin;
        console.log('Scadnano origin: ' + scadnanoOrigin);
        oxviewOrigin = window.location.origin;
        console.log('OxView origin: ' + oxviewOrigin);
    }
}