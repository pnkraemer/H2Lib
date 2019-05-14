/*
NAME: 'gmres_tps_sphere.c'

PURPOSE: Simulations for Example 7.3

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






real interpolant(real x, real y, real z){
    return x + y + z;
}


void make_fibonaccipoints(pclustergeometry A){
    uint numPts = A->nidx - 4;
    real offset = 2.0/(1.0*numPts);
    real increment = M_PI * (3.0 - sqrt(5.0));
    for(uint i = 0; i < numPts; i++){
        const real y = (((i*1.0) * offset)- 1 + offset * 0.5);
        const real r = sqrt(1 - y*y);
        const real phi = ((i + 1) % numPts) * increment;
        const real x = cos(phi) * r;
        const real z = sin(phi) * r;
        A->x[i][0] = x;
        A->smax[i][0] = x;
        A->smin[i][0] = x;
        A->x[i][1] = y;
        A->smax[i][1] = y;
        A->smin[i][1] = y;
        A->x[i][2] = z;
        A->smax[i][2] = z;
        A->smin[i][2] = z;
}
    pavector randvec = new_avector(3);
    for(uint i = numPts; i < numPts + 4; i++){
        random_avector(randvec);
        const real normofrandvec = norm2_avector(randvec);
        A->x[i][0] = randvec->v[0] / normofrandvec;
        A->smax[i][0] = randvec->v[0] / normofrandvec;
        A->smin[i][0] = randvec->v[0] / normofrandvec;
        A->x[i][1] = randvec->v[1] / normofrandvec;
        A->smax[i][1] = randvec->v[1] / normofrandvec;
        A->smin[i][1] = randvec->v[1] / normofrandvec;
        A->x[i][2] = randvec->v[2] / normofrandvec;
        A->smax[i][2] = randvec->v[2] / normofrandvec;
        A->smin[i][2] = randvec->v[2] / normofrandvec;
    }
    del_avector(randvec);
}








static void
loadfromtxt_precon(pavector preconVals, pavector preconRowIdx, pavector preconColIdx, uint N, uint n){

    uint numPts = preconVals->dim;
    char strVals[150], strRowIdx[150], strColIdx[150];
    sprintf(strVals, "/home/kraemer/Documents/rbf-tools-master/largeSphere/precon_txt_tps_sphere/preconVals_N%d_n%d.txt", N, n);
    sprintf(strRowIdx, "/home/kraemer/Documents/rbf-tools-master/largeSphere/precon_txt_tps_sphere/preconRowIdx_N%d_n%d.txt", N, n);
    sprintf(strColIdx, "/home/kraemer/Documents/rbf-tools-master/largeSphere/precon_txt_tps_sphere/preconColIdx_N%d_n%d.txt", N, n);
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
    psparsepattern sp = new_sparsepattern(N + 4, N + 4);
    for(uint i = 0; i < n*N; i++){
        addnz_sparsepattern(sp, preconRowIdx->v[i], preconColIdx->v[i]);
    }
    for(uint i = 0; i < 4; i++){
        addnz_sparsepattern(sp, N + i, N + i);
    }

    /* SPARSE MATRIX */
    psparsematrix spm = new_zero_sparsematrix(sp);
    for(uint i = 0; i < n*N; i++){
        setentry_sparsematrix(spm, preconRowIdx->v[i], preconColIdx->v[i], preconVals->v[i]);
    }
        for(uint i = 0; i < 4; i++){
        setentry_sparsematrix(spm, N+i, N+i, 1.0);
    }
    del_sparsepattern(sp);
    return spm;
}



pamatrix make_preconditioner(pavector preconVals, pavector preconRowIdx, pavector preconColIdx, uint N, uint n){

    pamatrix spm = new_zero_amatrix(N + 4, N+4);
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
    pclustergeometry clGeom = new_clustergeometry(3, N + 4);
    make_fibonaccipoints(clGeom);

    /* MAKE RHS */
    pavector rhsVec = new_zero_avector(N + 4);
    for(uint i = 0; i < N; i++){
        rhsVec->v[i] = interpolant(clGeom->x[i][0], clGeom->x[i][1], clGeom->x[i][2]);
    }

    /* LOAD H2-MATRIX*/
    start_kernel = time(NULL);
    AmP->am = new_kernelamatrix_tpssphere(clGeom, clGeom, 0.0);
    end_kernel = time(NULL) - start_kernel;
    
    start_h2 = time(NULL);
    uint* idxSet = allocuint(N + 4);// = allocuint(N);
    for(uint i = 0; i < N + 4; i++){
        idxSet[i] = i;
    }
    pcluster clust = build_cluster(clGeom, N + 4, idxSet, h2MatPar->leafSize, H2_REGULAR);
    pblock blockTree = build_strict_block(clust, clust, &(h2MatPar->admPar), admissible_2_cluster);
    H2P->h = makeH2Mat(AmP->am, blockTree, h2MatPar);
    end_h2 = time(NULL) - start_h2;

    sz_h = getsize_h2matrix(H2P->h);
    sz_full = getsize_amatrix(AmP->am);

    /* CHECK APPROXIMATION ERROR OF H2MATRIX*/
    pavector testvec = new_avector(N+4);
    random_avector(testvec);
    pavector residual_a = new_zero_avector(N+4);
    pavector residual_h2 = new_zero_avector(N+4);
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
    if(N + 4 <= 9999){
        copy_avector(rhsVec, startVec);
        iter_lang_noprecon = solve_gmres_amatrix_avector(AmP->am, rhsVec, startVec, gmresAcc, 1000, 20);
    }
    end_gmres_a_noprecon = time(NULL) - start_gmres_a_noprecon;

    uint iter_lang = 0;
    start_gmres_a = time(NULL);
    if(N + 4 <= 100001){
        copy_avector(rhsVec, startVec);
        iter_lang = solve_gmres_avector(AmP, addeval_amatrix_avector_withprecon, rhsVec, startVec, gmresAcc, 1000, 20);
    }
    end_gmres_a = time(NULL) - start_gmres_a;


    start_gmres_h2 = time(NULL);
    copy_avector(rhsVec, startVec);
    uint iter = solve_gmres_avector(H2P, addeval_h2matrix_avector_withprecon, rhsVec, startVec, gmresAcc, 1000, 20);
    end_gmres_h2 = time(NULL) - start_gmres_h2;

    if(N + 4 <= 9999){
        printf("\nNumber of GMRES iterations with H2 matrix and sparse preconditioner:\n\t#It = %i (long: %i, maxit 1000)(full matrix and precon: %i)\n\n", iter, iter_lang_noprecon, iter_lang);
    }else{
        printf("\nNumber of GMRES iterations with H2 matrix and sparse preconditioner:\n\t#It = %i (full matrix and precon: %i)(discard the timing for nonprecon gmres below!)\n\n", iter, iter_lang);
    }
    pavector solVec = new_zero_avector(N+4);
    addeval_sparsematrix_avector(1.0, H2P->p, startVec, solVec);

    /* COMPUTE RMSE OF INTERPOLATION */
    uint numEvalPts = 9999;
    pclustergeometry clGeomEval = new_clustergeometry(3, numEvalPts);
    make_fibonaccipoints(clGeomEval);
    pamatrix kernelMtrxEval = new_kernelamatrix_tpssphere(clGeomEval, clGeom, 0.0);
    pavector matVecProd2 = new_zero_avector(numEvalPts);
    addeval_amatrix_avector(1.0, kernelMtrxEval, solVec, matVecProd2);
    pavector trueApprox = new_zero_avector(numEvalPts - 4);
    for(uint i = 0; i < numEvalPts - 4; i++){
        trueApprox->v[i] = matVecProd2->v[i];
    }
    pavector trueEval = new_zero_avector(numEvalPts - 4);
    for(uint i = 0; i < numEvalPts - 4; i++){
        trueEval->v[i] = interpolant(clGeomEval->x[i][0], clGeomEval->x[i][1], clGeomEval->x[i][2]);
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






