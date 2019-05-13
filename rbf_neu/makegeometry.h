#include "../Library/tri2d.h"   // 

pclustergeometry newClgeom_mesh2d(uint numRef) {
	
	/* BUILD MESH */
	ptri2d   *meshHierarchy = (ptri2d *) allocmem((size_t) sizeof(ptri2d) * (numRef + 1));
	meshHierarchy[0] = new_unitsquare_tri2d();	/* Set domain */
	//meshHierarchy[0] = new_unitcircle_tri2d();
	//meshHierarchy[0] = new_lshape_tri2d();
	for (uint idx = 0; idx < numRef; idx++) {	/* Mesh refinements */
		meshHierarchy[idx + 1] = refine_tri2d(meshHierarchy[idx], NULL);
	}
	uint numPts = meshHierarchy[numRef]->vertices;
	
	uint dim = 2;
	pclustergeometry clGeom = new_clustergeometry(dim, numPts);
	
	/* BUILD CLUSTERGEOMETRY */
	for(uint i = 0; i < numPts; i++){
			(clGeom->x)[i] = meshHierarchy[numRef]->x[i];
			(clGeom->smax)[i] = meshHierarchy[numRef]->x[i];
			(clGeom->smin)[i] = meshHierarchy[numRef]->x[i];
	}

	/* DELETE MEMORY */
	uint i, j;
    for (i = 0; i <= numRef; i++) {
      j = numRef - i;
      del_tri2d(meshHierarchy[j]);
    }
    freemem(meshHierarchy);
    printf("fine");
    return clGeom;
}



pclustergeometry newClgeom_mesh1d(uint numPts) {
	uint dim = 1;
	pclustergeometry clGeom = new_clustergeometry(dim, numPts);

	/* BUILD CLUSTERGEOMETRY */
	for(uint i = 0; i < numPts; i++){
			(clGeom->x)[i][0] = (float)(i + 0.5)/numPts;
			(clGeom->smax)[i][0] = (float)(i + 0.5)/numPts;
			(clGeom->smin)[i][0] = (float)(i + 0.5)/numPts;
	}
	return clGeom;
}



pclustergeometry newClgeom_rand2d(int numPts) {
	uint dim = 2;
	pclustergeometry clGeom = new_clustergeometry(dim, numPts);

	/* BUILD CLUSTERGEOMETRY */
	for(uint i = 0; i < numPts; i++){
		for(uint d = 0; d < dim; d++){
			(clGeom->x)[i][d] = (float)rand()/(float)(RAND_MAX);
			(clGeom->smax)[i][d] = (clGeom->x)[i][d];
			(clGeom->smin)[i][d] = (clGeom->x)[i][d];
		}
	}
	return clGeom;
}



pclustergeometry newClgeom_rand1d(int numPts) {
	uint dim = 1;
	pclustergeometry clGeom = new_clustergeometry(dim, numPts);

	/* BUILD CLUSTERGEOMETRY */
	for(uint i = 0; i < numPts; i++){
		for(uint d = 0; d < dim; d++){
			(clGeom->x)[i][d] = (float)rand()/(float)(RAND_MAX);
			(clGeom->smax)[i][d] = (clGeom->x)[i][d];
			(clGeom->smin)[i][d] = (clGeom->x)[i][d];
		}
	}
	return clGeom;
}

