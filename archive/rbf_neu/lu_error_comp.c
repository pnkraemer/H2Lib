/*
NAME: 'lu_error_comp.c'

PURPOSE: Simulations for Example 7.2

AUTHOR: kraemer@ins.uni-bonn.de
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

#include "kernelmat.h"
#include "makegeometry.h"
#include "helpfcts.h"


int main(int argc, char *argv[]){



	/* DETERMINE DIMENSIONALITIES */
	uint numPts, dim, leafSize;
	real truncAcc, admPar, shiftPar;
	numPts = atoi(argv[1]);
	leafSize = 12;
	truncAcc = 1e-30;
	shiftPar = atof(argv[2]);
	dim = 1;
	admPar = 1.0;
	real luAcc = 1e-16;



	/* BUILD CLUSTERGEOMETRY */
	pclustergeometry clGeom = newClgeom_mesh1d(numPts);

	/* BUILD CLUSTERTREE(S) */
	uint* idxSet = allocuint(numPts);// = allocuint(N);
	for(uint i = 0; i < numPts; i++){
		idxSet[i] = i;
	}

	pcluster clust = build_cluster(clGeom, numPts, idxSet, leafSize, H2_REGULAR);

	/* BUILD BLOCK TREE */
	pblock blockTree = build_strict_block(clust, clust, &admPar, admissible_2_cluster);

	/* BUILD FULLMATRIX (KERNELMATRIX) */
	pamatrix kernelMtrx = new_kernelamatrix_maternsquare(clGeom, clGeom, shiftPar);

	/* COMPRESS FULL INTO H2- AND H-MATRIX */
	ptruncmode truncMode = new_truncmode();
	ph2matrix h2KernelMtrx = compress_amatrix_h2matrix(kernelMtrx, blockTree, truncMode, truncAcc);
	phmatrix hKernelMtrx = convert_h2matrix_hmatrix(h2KernelMtrx);



	/* COMPUTE LU-DECOMP. OF FULL MATRIX */
	pavector rhsVecFull = new_zero_avector(numPts);
	rhsVecFull->v[0] = 1.0;
	lrdecomp_amatrix(kernelMtrx);
	lrsolve_amatrix_avector(0, kernelMtrx, rhsVecFull);


	/* COMPUTE LU-DECOMP. OF H-MATRIX */
	pavector rhsVecH = new_zero_avector(numPts);
	ptruncmode truncModeLU = new_truncmode();
	rhsVecH->v[0] = 1.0;
	lrdecomp_hmatrix(hKernelMtrx, 0, luAcc);
	lrsolve_hmatrix_avector(0, hKernelMtrx, rhsVecH);

	/* CHECK ERROR OF SOLUTIONS */
	pavector diffOfSols = new_zero_avector(numPts);
	copy_avector(rhsVecH, diffOfSols);

	add_avector(-1.0, rhsVecFull, diffOfSols);
	real comprError = norm2_avector(diffOfSols)/norm2_avector(rhsVecFull);
	printf("\nRelative Diskrepanz in der l2-norm:\
\n\t||VOLL - H||_2/||VOLL||_2 \n\t= %e\n", comprError);



    /* FREE POINTERS */
	del_avector(diffOfSols);
	del_avector(rhsVecH);
	del_avector(rhsVecFull);
	del_hmatrix(hKernelMtrx);
	del_h2matrix(h2KernelMtrx);
	del_truncmode(truncModeLU);
	del_truncmode(truncMode);
	del_block(blockTree);
	del_amatrix(kernelMtrx);
	del_cluster(clust);
	del_clustergeometry(clGeom);
	freemem(idxSet);

	return 0;

}
