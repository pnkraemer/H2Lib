/*
NAME: 'sparseGrid.c'

PURPOSE: Find the bottleneck in the computations


1d

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


#include "riley.h"
#include "kernelmat.h"
#include "makegeometry.h"
#include "helpfcts.h"


/*
typedef struct {
	uint leafSize;
	real admPar;
	real truncAcc;
	real luAcc;
}HPar;


typedef struct  {
	uint maxIt;
	real acc;
	real shift;
}RileyPar;


typedef struct {
	uint numit;
	real lastRelError;
	pavector relError;
	pavector solution;
}RileyOutput;
*/
// typedef struct _RileyOutput RileyOutput;
// typedef RileyOutput* pRileyOutput;

// void del_RileyOutput(pRileyOutput RilOut) {
// 	del_avector(RilOut->relError);
// 	del_avector(RilOut->solution);
// 	free(RilOut);
// }

// pRileyOutPut new_RileyOutput() {
// 	pRileyOutput rilout;
// 	rilout = (pRileyOutput) allocmen(sizeof(pRileyOutput));
// 	return rilout;
// }





int main(int argc, char *argv[]){


	/* DETERMINE DIMENSIONALITIES */
	uint numPts = atoi(argv[1]);
	HPar *h2MatPar = malloc(sizeof(HPar));
	h2MatPar->admPar = 1.0;
	h2MatPar->leafSize = 12;
	h2MatPar->truncAcc = 1e-14;
	h2MatPar->luAcc = 1e-14;

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


	/* BUILD FULLMATRIX (KERNELMATRIX) */
	real shift = atof(argv[2]);
	pamatrix kernelMtrx = new_kernelamatrix_maternsquare(clGeom, clGeom, shift);
	ph2matrix h2kernelMtrx = makeH2Mat(kernelMtrx, blockTree, h2MatPar);
/*	phmatrix hkernelMtrx = makeHMat(kernelMtrx, blockTree, h2MatPar);
*/

 //   print_amatrix(kernelMtrx);
   // print_amatrix(kernelMtrxShift);



	//print_amatrix(kernelMtrxShift);

	/* MAKE RHS */
	pavector rhsVec = new_zero_avector(numPts);
	rhsVec->v[0] = 1.0;
	//pavector rhsVec2 = clone_avector(rhsVec);
	printf("\n--------\nN = %i, shift = %.1e\n--------\n", numPts, shift);



	/* LR DECOMPOSITION OF FULL MATRIX */
	pstopwatch luFull = new_stopwatch();
	start_stopwatch(luFull);
	lrdecomp_amatrix(kernelMtrx);	

	pavector trueSol = clone_avector(rhsVec);
	lrsolve_amatrix_avector(0, kernelMtrx, trueSol);
	real secLuFull = stop_stopwatch(luFull);
	printf("\nTime taken for LU of full matrix:\n\t%.1f sec\n", secLuFull);
	del_stopwatch(luFull);



	pclusterbasis rblow = build_from_cluster_clusterbasis(clust);
  	pclusterbasis cblow = build_from_cluster_clusterbasis(clust);
  	ph2matrix L = build_from_block_lower_h2matrix(blockTree, rblow, cblow);

  	pclusterbasis rbup = build_from_cluster_clusterbasis(clust);
  	pclusterbasis cbup = build_from_cluster_clusterbasis(clust);
  	ph2matrix R = build_from_block_upper_h2matrix(blockTree, rbup, cbup);

  	ptruncmode tm = new_releucl_truncmode();
  	pclusteroperator rwf = prepare_row_clusteroperator(h2kernelMtrx->rb, h2kernelMtrx->cb, tm);
  	pclusteroperator cwf = prepare_col_clusteroperator(h2kernelMtrx->rb, h2kernelMtrx->cb, tm);
  	pclusteroperator rwflow = prepare_row_clusteroperator(L->rb, L->cb, tm);
  	pclusteroperator cwflow = prepare_col_clusteroperator(L->rb, L->cb, tm);
  	pclusteroperator rwfup = prepare_row_clusteroperator(R->rb, R->cb, tm);
  	pclusteroperator cwfup = prepare_col_clusteroperator(R->rb, R->cb, tm);


	/* COMPUTE LU-DECOMP. OF H-MATRIX */
	pstopwatch luH = new_stopwatch();
	start_stopwatch(luH);



  	lrdecomp_h2matrix(h2kernelMtrx, rwf, cwf, L, rwflow, cwflow, R, rwfup, cwfup, tm,
		    h2MatPar->luAcc);
	lrsolve_h2matrix_avector(L, R, rhsVec);


	/* STARTING VECTOR */
	
//	startVec->v[0] = 1.0;
	
	//RileyOutput *rilOut = rileyAlgo2(rileyPar, h2kernelMtrx, L, R,  rhsVec, rhsVec);



	real secLuH = stop_stopwatch(luH);
	del_stopwatch(luH);
	printf("\nTime taken for LU of H2 matrix:\n\t%.1f sec\n", secLuH);
//	printf("\nNumber of iterations:\n\tnumIt = %i\n", rilOut->numit);



	/* CHECK ERROR OF SOLUTIONS */
	real normm = norm2_avector(trueSol);
	add_avector(-1.0, rhsVec, trueSol);
	real currentRelError = norm2_avector(trueSol) / normm;
	printf("\nRelative error of LU decomposition with H2:\n\t r = %.2e\n\n", currentRelError);

//	print_avector(rilOut->solution);
//	for(uint i = 0; i < rilOut->numit; i++){
//		printf("( %i, %.2e )\n", i, rilOut->relError->v[i]);
//	}

	/* FREE POINTERS */
	del_clusteroperator(rwf);
	del_clusteroperator(cwf);
	del_clusteroperator(rwflow);
	del_clusteroperator(cwflow);
	del_clusteroperator(rwfup);
	del_clusteroperator(cwfup);

/*	del_clusterbasis(rbup);
	del_clusterbasis(rblow);
	del_clusterbasis(cbup);
	del_clusterbasis(cblow);
*/
	del_clustergeometry(clGeom);
	del_cluster(clust);
	del_block(blockTree);
	del_amatrix(kernelMtrx);
	del_h2matrix(h2kernelMtrx);
	del_h2matrix(L);
	del_h2matrix(R);
	del_avector(rhsVec);
	del_avector(trueSol);
//	del_avector(startVec);
	del_truncmode(tm);
//	del_avector(rilOut->solution);
//	del_avector(rilOut->relError);
//	free(rilOut);

	free(h2MatPar);

	freemem(idxSet);

	return 0;

}









// void 
// lookatcluster(pccluster t)
// {
// 	printf("\n\nSize of the cluster: %d indices", t->size);
// 	printf("\nIndex set of cluster: ");
// 	for (uint i = 0; i < t->size; i++){
// 		printf("%d ", t->idx[i]);
// 	}
// 	printf("\nNumber of Sons: %d", t->sons);
// 	printf("\nDimensionality of bdg box: %d", t->dim);
// 	printf("\nMinimal coordinate(s) of bdg box:");
// 	for (uint i = 0; i < t->dim; i++){
// 		printf("%f ", t->bmin[i]);
// 	}
// 	printf("\nMaximal coordinate(s) of bdg box:");
// 	for (uint i = 0; i < t->dim; i++){
// 		printf("%f ", t->bmax[i]);
// 	}
// 	printf("\nNumber of descendants %d", t->desc);
// 	printf("\nType of the cluster: %d", t->type);
// 	printf("\nHow deep is the cluster?: %d level\n\n", getmindepth_cluster(t));
// }

