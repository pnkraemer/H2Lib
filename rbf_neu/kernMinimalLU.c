/*
NAME: 'kernMinimalLU.c'

PURPOSE: Provide a minimal example for LU decomposition 
of kernel matrices using H2 compression

DESCRIPTION: We construct a full kernel matrix, 
compress it into an H-matrix, compute the LU 
decompositions of the H and the full matrix 
and check the discrepancy of the solutions

AUTHOR: NK
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


real gaussKernel(real coord){
	return exp(-0.5 * coord);
}

real maternKernel(real coord){
	return (1 + sqrt(3)*coord/2) * exp(-0.5 *sqrt(3)* coord);
}

void 
lookatcluster(pccluster t)
{
	printf("\n\nSize of the cluster: %d indices", t->size);
	printf("\nIndex set of cluster: ");
	for (uint i = 0; i < t->size; i++){
		printf("%d ", t->idx[i]);
	}
	printf("\nNumber of Sons: %d", t->sons);
	printf("\nDimensionality of bdg box: %d", t->dim);
	printf("\nMinimal coordinate(s) of bdg box:");
	for (uint i = 0; i < t->dim; i++){
		printf("%f ", t->bmin[i]);
	}
	printf("\nMaximal coordinate(s) of bdg box:");
	for (uint i = 0; i < t->dim; i++){
		printf("%f ", t->bmax[i]);
	}
	printf("\nNumber of descendants %d", t->desc);
	printf("\nType of the cluster: %d", t->type);
	printf("\nHow deep is the cluster?: %d level\n\n", getmindepth_cluster(t));
}





int main(){
	printf("\n\n\n\n\n\n\n\n\n-----------------------\n");

	printf("\nNAME:\n\
	'kernMinimalLU.c'\n");

	printf("\nPURPOSE:\n\
	Provide a minimal example for LU decomposition\n\
	of kernel matrices using H2 compression\n");

	printf("\nDESCRIPTION:\n\
	We construct a full kernel matrix,\n\
	compress it into an H-matrix, compute the LU\n\
	decompositions of the H and the full matrix\n\
	and check the discrepancy of the solutions\n");

	printf("\nAUTHOR:\n\
	NK\n");
	printf("\n-----------------------\n");


	/* DETERMINE DIMENSIONALITIES */
	uint numPts, dim, leafSize;
	real truncAcc, admPar, shiftPar;
	numPts = askforint("\nMit wie vielen Punkten sollen wir arbeiten?", "zahl", 150);
	//dim = askforint("In wie vielen Dimensionen?", "zahl", 1);
	//leafSize = askforint("Wie groß sollen die klein(st)en Blätter sein?", "zahl", 8);
	leafSize = 12;
	truncAcc = 1e-14;
	shiftPar = 1e-15;
	dim = 1;
	//truncAcc = askforreal("Welche (relative) Komprimierungsgenauigkeit?", "zahl", 1e-16);
	//admPar  = askforreal("Welchen Zulässigkeitsparameter?", "zahl", 1.0);
	admPar = 1.0;
	//shiftPar  = askforreal("Welchen Shiftparameter sigma für (K+sigma I)?", "zahl", 0.1);
	real luAcc = 1e-16;

	printf("-----------------------\n");


	/* BUILD CLUSTERGEOMETRY */
	pclustergeometry clGeom = new_clustergeometry(dim, numPts);
	for(uint i = 0; i < numPts; i++){
		for(uint j = 0; j < dim; j++){
			const real randomPt = (double)(i)/numPts;
			(clGeom->x)[i][j] = randomPt;
			(clGeom->smax)[i][j] = randomPt;
			(clGeom->smin)[i][j] = randomPt;
		}
	}	

	/* BUILD CLUSTERTREE(S) */
	uint* idxSet = allocuint(numPts);// = allocuint(N);
	for(uint i = 0; i < numPts; i++){
		idxSet[i] = i;
	}

	pcluster clust = build_cluster(clGeom, numPts, idxSet, leafSize, H2_REGULAR);

	/* BUILD BLOCK TREE */
	pblock blockTree = build_strict_block(clust, clust, &admPar, admissible_2_cluster);

	/* BUILD FULLMATRIX (KERNELMATRIX) */
	pamatrix kernelMtrx = new_amatrix(numPts,numPts);
	real norm;
	for(uint i = 0; i < numPts; i++){
		real* X = allocreal(dim);
		for(uint d = 0; d < dim; d++){
			X[d] = (clGeom->x)[i][d];
		}
		for(uint j = 0; j < numPts; j++){
			real* Y = allocreal(dim);
			for(uint d = 0; d < dim; d++){
				Y[d] = (clGeom->x)[j][d];
			}
			norm = 0;
			for(uint d = 0; d < dim; d++){
				norm += (X[d] - Y[d])*(X[d] - Y[d]);
			}
			kernelMtrx->a[i * numPts + j] = maternKernel(sqrt(norm));
			freemem(Y);
			if(i==j){
				kernelMtrx->a[i * numPts + j] += shiftPar * 1.0;
			}
		}
		freemem(X);
	}

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
	printf("(LU-Zerlegung mit Genauigkeit %e)\n\n", luAcc);



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
