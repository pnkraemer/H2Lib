/*
NAME: 'selection_hyperparam.c'

PURPOSE: Simulations for Section 7.5.3

AUTHOR: kraemer(at)ins.uni-bonn.de
*/


#include <stdio.h>
#include <stdlib.h>   	// für rand
#include <math.h>   	// für sqrt und so (kernfunktionen)

#include "../Library/parameters.h"   // --> für askforint
#include "../Library/cluster.h"   // 
#include "../Library/block.h"   //
#include "../Library/h2matrix.h"   // 
#include "../Library/h2arith.h"   // 
#include "../Library/harith.h"   // 
#include "../Library/hmatrix.h"   // 
#include "../Library/tri2d.h"   // 
#include "../Library/krylovsolvers.h"   // 


#include "riley.h"
#include "kernelmat.h"
#include "makegeometry.h"
#include "helpfcts.h"



real interpolant(real x){

    return exp(x) * sin(x);

}


int main(int argc, char *argv[]){


	/* DETERMINE DIMENSIONALITIES */
	uint numPts = atoi(argv[1]);
	HPar *h2MatPar = malloc(sizeof(HPar));
	h2MatPar->admPar = 1.0;
	h2MatPar->leafSize = 12;
	h2MatPar->truncAcc = atof(argv[2]);
    real gmresAcc = atof(argv[3]);

	/* BUILD CLUSTERGEOMETRY */
    pclustergeometry clGeom = newClgeom_mesh1d(numPts);
	numPts = clGeom->nidx;

	/* BUILD CLUSTERTREE(S) */
	uint* idxSet = allocuint(numPts);// = allocuint(N);
	for(uint i = 0; i < numPts; i++){
		idxSet[i] = i;
	}
	pcluster clust = build_cluster(clGeom, numPts, idxSet, h2MatPar->leafSize, H2_REGULAR);

	/* BUILD BLOCK TREE */
	pblock blockTree = build_strict_block(clust, clust, &(h2MatPar->admPar), admissible_2_cluster);


	/* MAKE RHS */
	pavector rhsVec = new_zero_avector(numPts);
	
    for(uint i = 0; i < numPts; i++){
    	rhsVec->v[i] = interpolant(clGeom->x[i][0]);
    }


    /* Compute solution */
    pavector startVec = clone_avector(rhsVec);

    pamatrix kern2 = new_kernelamatrix(clGeom, clGeom, 0.0);
	ph2matrix h2kernelMtrx = makeH2Matfrob(kern2, blockTree, h2MatPar);
	uint iter = solve_gmres_h2matrix_avector(h2kernelMtrx, rhsVec, startVec, gmresAcc, numPts, numPts);
    del_amatrix(kern2);


    real h_nu = 1.0/((real)numPts);
    h_nu = pow(h_nu, 2.0);
    printf("\nExpected convergence rate: %.1e\n\n", h_nu);


    /* COMPUTE RMSE OF INTERPOLATION */
    uint numEvalPts = 100000;
	pclustergeometry clGeomEval = newClgeom_rand1d(numEvalPts);

	pamatrix kernelMtrxEval = new_kernelamatrix(clGeomEval, clGeom, 0.0);
	pavector matVecProd2 = new_zero_avector(numEvalPts);
	addeval_amatrix_avector(1.0, kernelMtrxEval, startVec, matVecProd2);
    

	pavector trueEval = new_zero_avector(numEvalPts);
	for(uint i = 0; i < numEvalPts; i++){
		trueEval->v[i] = interpolant((clGeomEval->x[i][0]));
	}

	add_avector(-1.0, trueEval, matVecProd2);
	real currentRelError = norm2_avector(matVecProd2) / sqrt((real)numEvalPts);
	printf("\nNormalised interpolation error: %.1e\n\n", currentRelError);

    del_clustergeometry(clGeomEval);
    del_amatrix(kernelMtrxEval);
    del_avector(matVecProd2);
    del_avector(trueEval);
	del_clustergeometry(clGeom);
	del_cluster(clust);
	del_block(blockTree);
	del_h2matrix(h2kernelMtrx);
	del_avector(rhsVec);
	del_avector(startVec);
	free(h2MatPar);
	freemem(idxSet);

	return 0;

}






