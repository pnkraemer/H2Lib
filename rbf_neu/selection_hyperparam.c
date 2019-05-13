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
#include "../Library/krylovsolvers.h"   // 


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






static void
precondition(void *pdata, pavector r){
    
    pamatrix A = (pamatrix) pdata;
    pavector rcopy = clone_avector(r);
    clear_avector(r);
    addeval_amatrix_avector(1.0, A, rcopy, r);
    del_avector(rcopy);

}



static void
loadfromtxt_matrix(pamatrix A){

    uint numPts = A->rows;
    char string[50];
    sprintf(string, "./rbf/matrix_files/%d_kernelMtrx.txt", numPts);

    FILE *myFile;
    myFile = fopen(string, "r");
    if(myFile == NULL){
        printf("Error reading file\n");
        exit(0);
    }
    
    for(uint i = 0; i < numPts; i++){
        for(uint j = 0; j < numPts; j++){
            fscanf(myFile, "%lf ", &(A->a[i + numPts * j]));      
        }
    }
    fclose(myFile);

}




static void
loadfromtxt_precon(pamatrix A){

    uint numPts = A->rows;
    char string[50];
    sprintf(string, "./rbf/matrix_files/%d_preconditioner.txt", numPts);

    FILE *myFile;
    myFile = fopen(string, "r");
    if(myFile == NULL){
        printf("Error reading file\n");
        exit(0);
    }
    
    for(uint i = 0; i < numPts; i++){
        for(uint j = 0; j < numPts; j++){
            fscanf(myFile, "%lf ", &(A->a[j + numPts * i]));      
        }
    }
    fclose(myFile);

}



static void
loadfromtxt_pts(pclustergeometry A){

    uint numPts = A->nidx;
    char string[50];
    sprintf(string, "./rbf/matrix_files/%d_ptSet.txt", numPts);

    FILE *myFile;
    myFile = fopen(string, "r");

    if(myFile == NULL){
        printf("Error reading file\n");
        exit(0);
    }
    
    for(uint i = 0; i < numPts; i++){
        for(uint j = 0; j < 1; j++){
            fscanf(myFile, "%lf ", &((A->x)[i][j]));
            (A->x)[i][j] = (A->x)[i][j];    
            (A->smax)[i][j] = (A->x)[i][j];      
            (A->smin)[i][j] = (A->x)[i][j];
        }
    }
 

 /*    for(uint i = 0; i < numPts; i++){
        for(uint j = 0; j < 3; j++){
            printf("%.2lf ",  (A->x)[i][j]);
        }
        printf("\n");
    }
*/    fclose(myFile);
}


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
//	h2MatPar->luAcc = 1e-12;
//	RileyPar *rileyPar = malloc(sizeof(RileyPar));
//	rileyPar->shift = atof(argv[2]);
//	rileyPar->acc = 1e-7;
//	rileyPar->maxIt = 1e4;

	/* BUILD CLUSTERGEOMETRY */
//	pclustergeometry clGeom = new_clustergeometry(1, numPts);
    pclustergeometry clGeom = newClgeom_mesh1d(numPts);
  //  	loadfromtxt_pts(clGeom);
	numPts = clGeom->nidx;

	/* BUILD CLUSTERTREE(S) */
	uint* idxSet = allocuint(numPts);// = allocuint(N);
	for(uint i = 0; i < numPts; i++){
		idxSet[i] = i;
	}

	pcluster clust = build_cluster(clGeom, numPts, idxSet, h2MatPar->leafSize, H2_REGULAR);
 //   printf("\n\nalrighty\n\n");
   // lookatcluster(clust);

	/* BUILD BLOCK TREE */
	pblock blockTree = build_strict_block(clust, clust, &(h2MatPar->admPar), admissible_2_cluster);


	/* BUILD FULLMATRIX (KERNELMATRIX) */

	//pamatrix kernelMtrx = new_kernelamatrix(clGeom, clGeom, shift);

//	real shift2 = shift + rileyPar->shift;
//	pamatrix kernelMtrxShift = new_kernelamatrix_exp(clGeom, clGeom, shift2);
//	ph2matrix h2kernelMtrxShift = makeH2Mat(kernelMtrxShift, blockTree, h2MatPar);
/*	phmatrix hkernelMtrx = makeHMat(kernelMtrx, blockTree, h2MatPar);
*/


	/* MAKE RHS */
	pavector rhsVec = new_zero_avector(numPts);
	
    for(uint i = 0; i < numPts; i++){
    	rhsVec->v[i] = interpolant(clGeom->x[i][0]);
    }
    //print_avector(rhsVec);
    //pamatrix kern2 = new_kernelamatrix(clGeom, clGeom, shift);
 //   pamatrix kernelMtrx = new_amatrix(numPts, numPts);

 // loadfromtxt_matrix(kernelMtrx);
//    print_amatrix(kernelMtrx);


    

/*
	cairo_surface_t *surface_h2mat = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 1536, 1024);
    cairo_t *cr_h2mat = cairo_create (surface_h2mat);
    draw_cairo_h2matrix(cr_h2mat, h2kernelMtrx, true,0);
    cairo_surface_write_to_png (surface_h2mat, "./bilder/h2mat.png");
    cairo_surface_destroy (surface_h2mat);
*/
    /* Compute SVD */
//    pamatrix U = new_amatrix(numPts, numPts);
  //  pamatrix Vt = new_amatrix(numPts, numPts);
//    prealavector sigma = new_realavector(numPts);

 /*   uint checkInt = svd_amatrix(kern2, sigma, U, Vt);
    if(checkInt == 0){
        printf("\n\nsvd was fine\n\n");
    }else{
        printf("\n\nSVD WAS NOT NOT NOT NOT FINE\n\n");
    }
//    print_amatrix(U);
//    print_amatrix(Vt);
//    print_realavector(sigma);
*/
    /* Compute Preconditioner */
//   print_realavector(sigma);
  //  pamatrix precon = new_amatrix(numPts, numPts);
    //loadfromtxt_precon(precon);







  //  loadfromtxt_matrix(kernelMtrx);




    /* Compute CG */
    pavector startVec = clone_avector(rhsVec);
    //random_avector(startVec);

    pamatrix kern2 = new_kernelamatrix(clGeom, clGeom, 0.0);
	ph2matrix h2kernelMtrx = makeH2Matfrob(kern2, blockTree, h2MatPar);
//    uint iter = solve_pgmres_amatrix_avector(kern2, precondition, precon, rhsVec, startVec, 1e-15, numPts, numPts)   ;
uint iter = solve_gmres_h2matrix_avector(h2kernelMtrx, rhsVec, startVec, gmresAcc, numPts, numPts);
//    print_avector(startVec);
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
//	random_real_avector(trueEval);
	for(uint i = 0; i < numEvalPts; i++){
		trueEval->v[i] = interpolant((clGeomEval->x[i][0]));
//		printf("%.2f ", trueEval->v[i]);
	}
//	printf("\n\n");
    //print_avector(trueEval);
//	print_avector(matVecProd2);
//	print_avector(trueEval);

	add_avector(-1.0, trueEval, matVecProd2);
	real currentRelError = norm2_avector(matVecProd2) / sqrt((real)numEvalPts);
	printf("\nNormalised interpolation error: %.1e\n\n", currentRelError);

    del_clustergeometry(clGeomEval);
    del_amatrix(kernelMtrxEval);
    del_avector(matVecProd2);
    del_avector(trueEval);















 //   del_amatrix(precon);
	del_clustergeometry(clGeom);
	del_cluster(clust);
	del_block(blockTree);
//	del_amatrix(kernelMtrx);
//	del_amatrix(kernelMtrxShift);
	del_h2matrix(h2kernelMtrx);
//	del_h2matrix(h2kernelMtrxShift);
//	del_h2matrix(L);
//	del_h2matrix(R);
	del_avector(rhsVec);
	del_avector(startVec);
//	del_avector(trueSol);
//	del_avector(startVec);
//	del_truncmode(tm);
//	del_avector(rilOut->solution);
//	del_avector(rilOut->relError);
//	free(rilOut);

//	free(rileyPar);
	free(h2MatPar);

	freemem(idxSet);

	return 0;

}






