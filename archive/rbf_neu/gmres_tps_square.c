/*
NAME: 'gmres_tps_square.c'

PURPOSE: Simulations for Example 7.4

AUTHOR: kraemer(at)ins.uni-bonn.de
*/

#include <stdio.h>
#include <stdlib.h>   	// für rand
#include <math.h>   	// für sqrt und so (kernfunktionen)
#include <time.h>

#include "../Library/parameters.h"   // 
#include "../Library/cluster.h"   // 
#include "../Library/block.h"   //
#include "../Library/h2matrix.h"   // 
#include "../Library/h2arith.h"   // 
#include "../Library/harith.h"   // 
#include "../Library/hmatrix.h"   // 
#include "../Library/tri2d.h"   // 
#include "../Library/krylovsolvers.h"   // 

#include "kernelmat.h"
#include "makegeometry.h"
#include "helpfcts.h"






real interpolant(real x, real y){
    return 2.0;
}



void make_latticepoints(pclustergeometry A, bool printt){

    uint numPts = A->nidx - 3;
    real genvec0 = 1.0;
    real genvec1 = 433461.0;
    
    for(uint i = 0; i < numPts; i++){
        A->x[i][0] = fmod(genvec0 * i / (numPts), 1.0);
        A->smin[i][0] = A->x[i][0];
        A->smin[i][0] = A->x[i][0];
        A->x[i][1] = fmod(genvec1 * i / (numPts), 1.0);
        A->smin[i][1] = A->x[i][1];
        A->smin[i][1] = A->x[i][1];
        if(printt == 1){
        printf("(%.3f, %.3f), ", A->x[i][0], A->x[i][1] );
        }
        
    }
    pavector randvec = new_avector(2);
    for(uint i = numPts; i < numPts + 3; i++){
        random_avector(randvec);
        A->x[i][0] = randvec->v[0] ;
        A->smax[i][0] = randvec->v[0] ;
        A->smin[i][0] = randvec->v[0] ;
        A->x[i][1] = randvec->v[1] ;
        A->smax[i][1] = randvec->v[1] ;
        A->smin[i][1] = randvec->v[1] ;
    }
    del_avector(randvec);

}





static void
loadfromtxt_precon(pavector preconVals, pavector preconRowIdx, pavector preconColIdx, uint N, uint n){

    uint numPts = preconVals->dim;
    char strVals[150], strRowIdx[150], strColIdx[150];
    sprintf(strVals, "/home/kraemer/Documents/rbf-tools-master/largeSphere/precon_txt_tps_square/preconVals_N%d_n%d.txt", N, n);
    sprintf(strRowIdx, "/home/kraemer/Documents/rbf-tools-master/largeSphere/precon_txt_tps_square/preconRowIdx_N%d_n%d.txt", N, n);
    sprintf(strColIdx, "/home/kraemer/Documents/rbf-tools-master/largeSphere/precon_txt_tps_square/preconColIdx_N%d_n%d.txt", N, n);
    FILE *myFileVals, *myFileColIdx, *myFileRowIdx;
    myFileVals = fopen(strVals, "r");
    myFileColIdx = fopen(strColIdx, "r");
    myFileRowIdx = fopen(strRowIdx, "r");
    if(myFileVals == NULL){
        printf("Error reading file\n");
        exit(0);
    }
    for(uint i = 0; i < numPts; i++){
            fscanf(myFileVals, "%lf ", &(preconVals->v[i]));      
            fscanf(myFileColIdx, "%lf ", &(preconColIdx->v[i]));      
            fscanf(myFileRowIdx, "%lf ", &(preconRowIdx->v[i]));      
    }
    fclose(myFileVals);
    fclose(myFileColIdx);
    fclose(myFileRowIdx);
}

psparsematrix make_preconditioner_sparse(pavector preconVals, pavector preconRowIdx, pavector preconColIdx, uint N, uint n){

    /* SPARSITY PATTERN */
    psparsepattern sp = new_sparsepattern(N + 3, N + 3);
    for(uint i = 0; i < n*N; i++){
        addnz_sparsepattern(sp, preconRowIdx->v[i], preconColIdx->v[i]);
    }
    for(uint i = 0; i < 3; i++){
        addnz_sparsepattern(sp, N + i, N + i);
    }

    /* SPARSE MATRIX */
    psparsematrix spm = new_zero_sparsematrix(sp);
    for(uint i = 0; i < n*N; i++){
        setentry_sparsematrix(spm, preconRowIdx->v[i], preconColIdx->v[i], preconVals->v[i]);
    }
        for(uint i = 0; i < 3; i++){
        setentry_sparsematrix(spm, N+i, N+i, 1.0);
    }
    del_sparsepattern(sp);
    return spm;
}



pamatrix make_preconditioner(pavector preconVals, pavector preconRowIdx, pavector preconColIdx, uint N, uint n){

    pamatrix spm = new_zero_amatrix(N + 3, N+3);
    for(uint i = 0; i < n*N; i++){
        setentry_amatrix(spm, preconRowIdx->v[i], preconColIdx->v[i], preconVals->v[i]);
    }
    return spm;
}




typedef struct{

    ph2matrix h;
	psparsematrix p;
}h2matrixandamatrixandprecon;



void
addeval_h2matrix_avector_withprecon(field alpha, h2matrixandamatrixandprecon *H2P, pcavector src, pavector trg)
{
    pavector rdash = new_zero_avector(src->dim);
    addeval_sparsematrix_avector(1.0, H2P->p, src, rdash);
    addeval_h2matrix_avector(alpha, H2P->h, rdash, trg);
    del_avector(rdash);
}



typedef struct{

    pamatrix am;
	psparsematrix p;
}amatrixandamatrixandprecon;



void
addeval_amatrix_avector_withprecon(field alpha, amatrixandamatrixandprecon *H2P, pcavector src, pavector trg)
{
    pavector rdash = new_zero_avector(src->dim);
    addeval_sparsematrix_avector(1.0, H2P->p, src, rdash);
    addeval_amatrix_avector(alpha, H2P->am, rdash, trg);
    del_avector(rdash);
}







int main(int argc, char *argv[]){
    
    
    /* METHOD PARAMETERS*/
    uint N = atoi(argv[1]);
    uint n = atoi(argv[2]);
    real gmresAcc = 1e-8;
    HPar *h2MatPar = malloc(sizeof(HPar));
    h2MatPar->admPar = 1.0;
    h2MatPar->leafSize = 12;
    h2MatPar->truncAcc = 1e-12;
    longindex sz_h, sz_full, sz_sparse;
    h2matrixandamatrixandprecon *H2P = malloc(sizeof(h2matrixandamatrixandprecon));
    amatrixandamatrixandprecon *AmP = malloc(sizeof(amatrixandamatrixandprecon));
    time_t start_kernel, end_kernel, start_sparse, end_sparse, start_h2, end_h2, start_gmres_a, end_gmres_a,start_gmres_a_noprecon,  end_gmres_a_noprecon, start_gmres_h2, end_gmres_h2;
    time_t start_all, end_all;
    start_all = time(NULL);

    /* LOAD FIBONACCIPOINTS*/
    pclustergeometry clGeom = new_clustergeometry(2, N + 3);
    make_latticepoints(clGeom, 0);

    /* MAKE RHS */
    pavector rhsVec = new_zero_avector(N + 3);
    for(uint i = 0; i < N; i++){
        rhsVec->v[i] = interpolant(clGeom->x[i][0], clGeom->x[i][1]);
    }

    /* LOAD H2-MATRIX*/
    start_kernel = time(NULL);
    AmP->am = new_kernelamatrix_tpssquare(clGeom, clGeom, 0.0);
    end_kernel = time(NULL) - start_kernel;
    
    start_h2 = time(NULL);
    uint* idxSet = allocuint(N + 3);// = allocuint(N);
    for(uint i = 0; i < N + 3; i++){
        idxSet[i] = i;
    }
    pcluster clust = build_cluster(clGeom, N + 3, idxSet, h2MatPar->leafSize, H2_REGULAR);
    pblock blockTree = build_strict_block(clust, clust, &(h2MatPar->admPar), admissible_2_cluster);
    H2P->h = makeH2Mat(AmP->am, blockTree, h2MatPar);
    end_h2 = time(NULL) - start_h2;

    sz_h = getsize_h2matrix(H2P->h);
    sz_full = getsize_amatrix(AmP->am);

    /* CHECK APPROXIMATION ERROR OF H2MATRIX*/
    pavector testvec = new_avector(N+3);
    random_avector(testvec);
    pavector residual_a = new_zero_avector(N+3);
    pavector residual_h2 = new_zero_avector(N+3);
    addeval_amatrix_avector(1.0, AmP->am, testvec, residual_a);
    addeval_h2matrix_avector(1.0, H2P->h, testvec, residual_h2);
    add_avector(-1.0, residual_a, residual_h2);
    real diff = norm2_avector(residual_h2) / norm2_avector(testvec);
    printf("\nH2-MVM discrepancy:\n\tdiscr = %.1e (truncacc: %.1e)\n", diff, h2MatPar->truncAcc);
    del_avector(testvec);
    del_avector(residual_a);
    del_avector(residual_h2);
    
    
    /* LOAD PRECONDITIONER*/
    start_sparse = time(NULL);
    pavector preconVals = new_avector(N * n);
    pavector preconRowIdx = new_avector(N * n);
    pavector preconColIdx = new_avector(N * n);
    loadfromtxt_precon(preconVals, preconRowIdx, preconColIdx, N, n);
    H2P->p = make_preconditioner_sparse(preconVals, preconRowIdx, preconColIdx, N, n);
    sz_sparse = getsize_sparsematrix(H2P->p);
    end_sparse = time(NULL) - start_sparse;
    AmP->p = make_preconditioner_sparse(preconVals, preconRowIdx, preconColIdx, N, n);

    /* SOLVE LINEAR SYSTEM */
    pavector startVec = new_zero_avector(rhsVec->dim);
    uint iter_lang_noprecon = 0;
    start_gmres_a_noprecon = time(NULL);
    if(N + 3 <= 9999){
        copy_avector(rhsVec, startVec);
        iter_lang_noprecon = solve_gmres_amatrix_avector(AmP->am, rhsVec, startVec, gmresAcc, 10000, 20);
    }
    end_gmres_a_noprecon = time(NULL) - start_gmres_a_noprecon;

    uint iter_lang = 0;
    start_gmres_a = time(NULL);
    if(N + 3 <= 99999){
        copy_avector(rhsVec, startVec);
        iter_lang = solve_gmres_avector(AmP, addeval_amatrix_avector_withprecon, rhsVec, startVec, gmresAcc, 10000, 20);
    }
    end_gmres_a = time(NULL) - start_gmres_a;


    start_gmres_h2 = time(NULL);
    copy_avector(rhsVec, startVec);
    uint iter = solve_gmres_avector(H2P, addeval_h2matrix_avector_withprecon, rhsVec, startVec, gmresAcc, 10000, 20);
    end_gmres_h2 = time(NULL) - start_gmres_h2;

    if(N + 3 <= 9999){
        printf("\nNumber of GMRES iterations with H2 matrix and sparse preconditioner:\n\t#It = %i (long: %i, maxit 10000)(full matrix and precon: %i)\n\n", iter, iter_lang_noprecon, iter_lang);
    }else{
        printf("\nNumber of GMRES iterations with H2 matrix and sparse preconditioner:\n\t#It = %i (full matrix and precon: %i)(discard the timing for full gmres below!)\n\n", iter, iter_lang);
    }
    pavector solVec = new_zero_avector(N+3);
    addeval_sparsematrix_avector(1.0, H2P->p, startVec, solVec);

    /* COMPUTE RMSE OF INTERPOLATION */
    uint numEvalPts = 9999;
    pclustergeometry clGeomEval = new_clustergeometry(2, numEvalPts);
    make_latticepoints(clGeomEval, 0);
    pamatrix kernelMtrxEval = new_kernelamatrix_tpssquare(clGeomEval, clGeom, 0.0);
    pavector matVecProd2 = new_zero_avector(numEvalPts);
    addeval_amatrix_avector(1.0, kernelMtrxEval, solVec, matVecProd2);
    pavector trueApprox = new_zero_avector(numEvalPts - 3);
    for(uint i = 0; i < numEvalPts - 3; i++){
        trueApprox->v[i] = matVecProd2->v[i];
    }
    pavector trueEval = new_zero_avector(numEvalPts - 3);
    for(uint i = 0; i < numEvalPts - 3; i++){
        trueEval->v[i] = interpolant(clGeomEval->x[i][0], clGeomEval->x[i][1]);
    }
    add_avector(-1.0, trueEval, trueApprox);
    real currentRelError = norm2_avector(trueApprox) / sqrt((real)(trueApprox->dim));
    printf("\nRMSE of interpolation with H2 matrix and sparse preconditioner:\n\tRMSE =: %.1e\n\n", currentRelError);

    printf("\nMemory footprints:\n\tfull \t= %.2f MB\n\tH2 \t= %.2f MB\n\tsparse \t= %.2f MB\n", sz_full/1024.0/1024.0, sz_h/1024.0/1024.0,  sz_sparse/1024.0/1024.0);

    
    printf("\nStopwatch allocation:\n\tfull \t= %.2f s\n\tH2 \t= %.2f s\n\tsparse \t= %.2f s\n", (double)end_kernel, (double)end_h2, (double)end_sparse);
    printf("\nStopwatch gmres:\n\tfull \t= %.2f s\n\tH2 \t= %.2f s\n\tnoprec\t= %.2f s (full)\n", (double)end_gmres_a, (double)end_gmres_h2, (double)end_gmres_a_noprecon);

    /* CLEAN UP */
    del_clustergeometry(clGeomEval);
    del_amatrix(kernelMtrxEval);
    del_avector(matVecProd2);
    del_avector(trueEval);
    del_avector(trueApprox);
    del_amatrix(AmP->am);
    del_clustergeometry(clGeom);
    del_avector(rhsVec);
    del_avector(startVec);
    del_avector(solVec);
    free(h2MatPar);
    del_avector(preconVals);
    del_avector(preconRowIdx);
    del_avector(preconColIdx);
    del_sparsematrix(H2P->p);
    del_sparsematrix(AmP->p);
    del_cluster(clust);
    del_block(blockTree);
    del_h2matrix(H2P->h);
    freemem(idxSet);
    free(H2P);
    free(AmP);
    end_all = time(NULL) - start_all;
    printf("\nStopwatch all the rest:\n\tall \t= %.2f s\n\tergo: \t~ %.2f s for all the other stuff\n\n", (double)end_all, (double)(end_all - end_kernel - end_h2 - end_sparse - end_gmres_a - end_gmres_h2));

    return 0;

}






