/*
* A few add-ons for the krylovsolvers.h module,
* mostly regarding conditionally positive definite kernels
*
* Author: Nicholas Krämer
*/

#ifndef ADDON_KRYLOVSOLVERS
#define ADDON_KRYLOVSOLVERS

#include "settings.h"
#include "krylovsolvers.h"
#include "krylov.h"
#include "kernelmatrix.h"
#include "addon_kernelmatrix.h"
#include "avector.h"


HEADER_PREFIX void 
addeval_cond_kernelh2matrix(field alpha, ph2matrix h2km, pamatrix pb, pavector src, pavector trg);
HEADER_PREFIX void 
addeval_cond_kernelh2matrix_precon(field alpha, ph2matrix h2km, pamatrix pb, psparsematrix spm, pavector src, pavector trg);

HEADER_PREFIX uint 
solve_gmres_h2precond_avector(ph2matrix h2m, pamatrix polblock, psparsematrix spm, pcavector b, pavector x, real eps, uint maxiter, uint kmax);

HEADER_PREFIX void 
loadfromtxt_precon(pavector preconVals, pavector preconRowIdx, pavector preconColIdx, uint N, uint n);

HEADER_PREFIX psparsematrix 
make_precon_sparse(pavector preconVals, pavector preconRowIdx, pavector preconColIdx, uint N, uint n);



HEADER_PREFIX pamatrix 
make_precon_full(pavector preconVals, pavector preconRowIdx, pavector preconColIdx, uint N, uint n);

HEADER_PREFIX void
init_gmres_special(ph2matrix h2m, pamatrix polblock, psparsematrix spm, pcavector b, pavector x, pavector rhat, pavector q, uint *kk, pamatrix qr, pavector tau);






HEADER_PREFIX void
step_gmres_special(ph2matrix h2m, pamatrix polblock, psparsematrix spm, pcavector b, pavector x, pavector rhat, pavector q, uint *kk, pamatrix qr, pavector tau);



HEADER_PREFIX void
finish_gmres_special(ph2matrix h2m, pamatrix polblock, psparsematrix spm, pcavector b, pavector x, pavector rhat, pavector q, uint *kk, pamatrix qr, pavector tau);



#endif