
#include <stdio.h>
#include <stdlib.h>   	// für rand
#include <math.h>   	// für sqrt und so (kernfunktionen)
#include "parameters.h"   // --> für askforint
#include "clustergeometry.h"   // selbsterklärend
#include "cluster.h"   // --> selbsterklärend
#include "block.h"   // --> selbsterklärend
#include "hmatrix.h"   // --> selbsterklärend
#include "h2matrix.h"   // --> selbsterklärend
#include "h2arith.h"   // --> selbsterklärend
#include "h2compression.h"   // --> selbsterklärend
#include "amatrix.h"   // --> selbsterklärend
#include "matrixnorms.h"   // --> selbsterklärend

#include "riley.h"
#include "kernelmat.h"
#include "makegeometry.h"
#include "helpfcts.h"


void prepare_lrdecomp_h2matrix(ph2matrix h2kernelMtrx, ph2matrix L, ph2matrix R, HPar *h2MatPar){

	pblock blockk = build_from_h2matrix_block(h2kernelMtrx);
 	L = build_from_block_lower_h2matrix(blockk, h2kernelMtrx->rb, h2kernelMtrx->cb);
	R = build_from_block_upper_h2matrix(blockk, h2kernelMtrx->rb, h2kernelMtrx->cb);
	
	
  	ptruncmode tm = new_releucl_truncmode();	
	pclusteroperator rwf = prepare_row_clusteroperator(h2kernelMtrx->rb, h2kernelMtrx->cb, tm);
  	pclusteroperator cwf = prepare_col_clusteroperator(h2kernelMtrx->rb, h2kernelMtrx->cb, tm);
  	pclusteroperator rwflow = prepare_row_clusteroperator(L->rb, L->cb, tm);
  	pclusteroperator cwflow = prepare_col_clusteroperator(L->rb, L->cb, tm);
  	pclusteroperator rwfup = prepare_row_clusteroperator(R->rb, R->cb, tm);
  	pclusteroperator cwfup = prepare_col_clusteroperator(R->rb, R->cb, tm);

  	lrdecomp_h2matrix(h2kernelMtrx, rwf, cwf, L, rwflow, cwflow, R, rwfup, cwfup, tm,
		    h2MatPar->luAcc);


	del_clusteroperator(rwf);
	del_clusteroperator(cwf);
	del_clusteroperator(rwflow);
	del_clusteroperator(cwflow);
	del_clusteroperator(rwfup);
	del_clusteroperator(cwfup);
	del_block(blockk);

}


/*
Run as:
./rbf/storagecomp 100
*/
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

	/* BUILD CLUSTERTREE(S) */
	uint* idxSet = allocuint(numPts);// = allocuint(N);
	for(uint i = 0; i < numPts; i++){
		idxSet[i] = i;
	}
	pcluster clust = build_cluster(clGeom, numPts, idxSet, h2MatPar->leafSize, H2_REGULAR);


	/* BUILD BLOCK TREE */
	pblock blockTree = build_strict_block(clust, clust, &(h2MatPar->admPar), admissible_2_cluster);


	/* BUILD FULLMATRIX (KERNELMATRIX) */
	real shift = 0.1;
	pamatrix kernelMtrx = new_kernelamatrix_maternsquare(clGeom, clGeom, shift);
	ph2matrix h2kernelMtrx = makeH2Mat(kernelMtrx, blockTree, h2MatPar);
/*	phmatrix hkernelMtrx = makeHMat(kernelMtrx, blockTree, h2MatPar);
*/

//	print_amatrix(kernelMtrx);

	/* MAKE RHS */
	pavector rhsVec = new_zero_avector(numPts);
	rhsVec->v[0] = 1.0;
	pavector rhsVec2 = clone_avector(rhsVec);
	printf("\n--------\nN = %i\n--------\n", numPts);


	lrdecomp_amatrix(kernelMtrx);	

	lrsolve_amatrix_avector(0, kernelMtrx, rhsVec2);
    real norm123 = norm2_avector(rhsVec2);
	/* HOW GOOD WAS THE COMPRESSION? */
/*	real relDiffCompression = norm2diff_amatrix_h2matrix(h2kernelMtrx, kernelMtrx)/norm2_amatrix(kernelMtrx);
	printf("\nKompressionsgenauigkeit in der Spektralnorm:\n\t||VOLL - H2||_2/||VOLL||_2 \n\t= %e\n", relDiffCompression);
	size_t sizefull = getsize_amatrix(kernelMtrx);
	printf("\nGröße der VOLLEN MATRIX:\n\t%f MB", sizefull/1024.0/1024.0);
	size_t sizeh = getsize_hmatrix(hkernelMtrx);
	printf("\nGröße der HMATRIX:\n\t%f MB", sizeh/1024.0/1024.0);
	size_t sizeh2 = getsize_h2matrix(h2kernelMtrx);
	printf("\nGröße der H2MATRIX:\n\t%f MB\n\n", sizeh2/1024.0/1024.0);
*/


/*	print_avector(rhsVec);
*/



    /* DRAW H2MATRIX */
/*	cairo_surface_t *surface_h2mat = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 1536, 1024);
    cairo_t *cr_h2mat = cairo_create (surface_h2mat);
    draw_cairo_h2matrix(cr_h2mat, h2kernelMtrx, true,0);
    cairo_surface_write_to_png (surface_h2mat, "h2mat.png");
    cairo_surface_destroy (surface_h2mat);
*/





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

  	lrdecomp_h2matrix(h2kernelMtrx, rwf, cwf, L, rwflow, cwflow, R, rwfup, cwfup, tm,
		    h2MatPar->luAcc);

	lrsolve_h2matrix_avector(L, R, rhsVec);

 //   	print_avector(rhsVec);

   	add_avector(-1.0, rhsVec, rhsVec2);
	real currentRelError = norm2_avector(rhsVec2) / norm123;
	printf("\nRelative error of LU decomposition with H2:\n\t r = %.2e\n\n", currentRelError);
    	
/* 	del_clusterbasis(rblow);
  	del_clusterbasis(cblow);
  	del_clusterbasis(rbup);
  	del_clusterbasis(cbup);
*/ 	del_h2matrix(L);
  	del_h2matrix(R);
  	del_truncmode(tm);
  	del_clusteroperator(rwf);
  	del_clusteroperator(cwf);
  	del_clusteroperator(rwflow);
  	del_clusteroperator(cwflow);
  	del_clusteroperator(rwfup);
  	del_clusteroperator(cwfup);
	/* LR DECOMP */
/*	pclusterbasis clBasis = build_from_cluster_clusterbasis(clust);
	ph2matrix h2Lower = build_from_block_lower_h2matrix(blockTree, clBasis, clBasis);
	ph2matrix h2Upper = build_from_block_upper_h2matrix(blockTree, clBasis, clBasis);

	ptruncmode tm = new_truncmode();
	pclusteroperator clOp = prepare_col_clusteroperator(clBasis, clBasis, tm);

//	lrdecomp_h2matrix(h2kernelMtrx, clOp, clOp, h2Lower, clOp, clOp, h2Upper, clOp, clOp, tm, h2MatPar->luAcc);



	lrsolve_h2matrix_avector(h2Lower, h2Upper, rhsVec);
	print_avector(rhsVec);
*/






    /* FREE POINTERS */
	del_h2matrix(h2kernelMtrx);
	del_avector(rhsVec);
	del_avector(rhsVec2);
/*	del_hmatrix(hkernelMtrx);
*/	del_block(blockTree);
	del_amatrix(kernelMtrx);
	del_cluster(clust);
	del_clustergeometry(clGeom);
	freemem(idxSet);
	freemem(h2MatPar);

/*	del_truncmode(tm);
	del_clusteroperator(clOp);
//	del_clusteroperator(rowOp);
	del_h2matrix(h2Lower);
	del_h2matrix(h2Upper);
//	unref_clusterbasis(clBasis);
//	del_clusterbasis(clBasis);
*/

	return 0;

}
