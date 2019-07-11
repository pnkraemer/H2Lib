#include "addon_krylovsolvers.h"
#include "settings.h"
#include "krylovsolvers.h"
#include "krylov.h"
#include "kernelmatrix.h"
#include "addon_kernelmatrix.h"
#include "avector.h"


static field
findapply_givens(pfield a, pfield b)
{
  field     aa = *a;
  field     bb = *b;
  field     s, c, t;
  field     rho;

  if (ABSSQR(bb) < ABSSQR(aa)) {
    t = bb / aa;
    c = REAL_RSQRT(ABSSQR(t) + 1.0);
    s = t * c;

    rho = s;
  }
  else if (ABSSQR(aa) < ABSSQR(bb)) {
    t = aa / bb;
    s = REAL_RSQRT(ABSSQR(t) + 1.0);
    c = t * s;

    rho = 1.0 / c;
  }
  else {
    c = 1.0;
    s = 0.0;

    rho = 1.0;
  }

  *a = c * aa + s * bb;
  *b = 0.0;

  return rho;
}

static void
apply_givens(field rho, pfield a, pfield b)
{
  field     aa = *a;
  field     bb = *b;
  field     s, c;

  if (ABSSQR(rho) < 1.0) {
    s = rho;
    c = REAL_SQRT(1.0 - ABSSQR(s));
  }
  else if (ABSSQR(rho) > 1.0) {
    c = 1.0 / rho;
    s = REAL_SQRT(1.0 - ABSSQR(c));
  }
  else {
    c = 1.0;
    s = 0.0;
  }

  *a = c * aa + s * bb;
  *b = c * bb - s * aa;
}


HEADER_PREFIX void
init_gmres_special(ph2matrix h2m, pamatrix polblock, psparsematrix spm, pcavector b, pavector x, pavector rhat, pavector q, uint *kk, pamatrix qr, pavector tau)
{

  avector   tmp1;
  amatrix   tmp2;
  pavector  r;
  pamatrix  qr_k;
  uint      kmax = qr->cols;

  assert(b->dim == x->dim);
  assert(b->dim == q->dim);
  assert(b->dim == qr->rows);
  assert(kmax <= tau->dim);

  if (kmax < 1)
    return;

  /* Residual r in the first column of qr */
  r = init_column_avector(&tmp1, qr, 0);
  copy_avector(b, r);
//  addeval(-1.0, matrix, x, r);
  addeval_cond_kernelh2matrix_precon(-1.0, h2m, polblock, spm, x, r);
  uninit_avector(r);

  /* Compute factorization */
  qr_k = init_sub_amatrix(&tmp2, qr, qr->rows, 0, 1, 0);
  qrdecomp_amatrix(qr_k, tau);

  /* Construct first orthogonal direction */
  clear_avector(q);
  q->v[0] = 1.0;
  qreval_amatrix_avector(false, qr_k, tau, q);
  uninit_amatrix(qr_k);

  /* Set up transformed residual */
  clear_avector(rhat);
  rhat->v[0] = qr->a[0];

  /* Set dimension */
  *kk = 0;
}







HEADER_PREFIX void
step_gmres_special(ph2matrix h2m, pamatrix polblock, psparsematrix spm, pcavector b, pavector x, pavector rhat, pavector q, uint *kk, pamatrix qr, pavector tau)
{
  avector   tmp1, tmp2;
  amatrix   tmp3;
  pavector  a, tau_k;
  pamatrix  qr_k;
  field     rho;
  uint      k = *kk;
  uint      kmax = qr->cols;
  uint      i;

  (void) b;
  (void) x;

  if (k + 1 >= kmax)
    return;

  /* (k+1)-th Krylov vector A q in the (k+1)-th column of qr */
  a = init_column_avector(&tmp1, qr, k + 1);
  clear_avector(a);
//    addeval(1.0, matrix, q, a);

  addeval_cond_kernelh2matrix_precon(1.0, h2m, polblock, spm, q, a);
  uninit_avector(a);

  /* Apply previous reflections */
  qr_k = init_sub_amatrix(&tmp3, qr, qr->rows, 0, k + 1, 0);
  qreval_amatrix_avector(true, qr_k, tau, a);
  uninit_amatrix(qr_k);

  /* Compute next reflection */
  qr_k = init_sub_amatrix(&tmp3, qr, qr->rows - (k + 1), k + 1, 1, k + 1);
  tau_k = init_sub_avector(&tmp2, tau, 1, k + 1);
  qrdecomp_amatrix(qr_k, tau_k);
  uninit_avector(tau_k);
  uninit_amatrix(qr_k);

  /* Construct next orthogonal direction */
  qr_k = init_sub_amatrix(&tmp3, qr, qr->rows, 0, k + 2, 0);
  clear_avector(q);
  q->v[k + 1] = 1.0;
  qreval_amatrix_avector(false, qr_k, tau, q);
  uninit_amatrix(qr_k);

  /* Apply preceding Givens rotations */
  for (i = 0; i < k; i++) {
    rho = qr->a[(i + 1) + (i + 1) * qr->ld];
    apply_givens(rho, qr->a + i + (k + 1) * qr->ld,
         qr->a + (i + 1) + (k + 1) * qr->ld);
  }

  /* Eliminate subdiagonal */
  rho = findapply_givens(qr->a + k + (k + 1) * qr->ld,
             qr->a + (k + 1) + (k + 1) * qr->ld);
  qr->a[(k + 1) + (k + 1) * qr->ld] = rho;
  apply_givens(rho, rhat->v + k, rhat->v + (k + 1));

  /* Increase dimension */
  *kk = k + 1;
}





HEADER_PREFIX void
finish_gmres_special(ph2matrix h2m, pamatrix polblock, psparsematrix spm, pcavector b, pavector x, pavector rhat, pavector q, uint *kk, pamatrix qr, pavector tau)
{

  avector   tmp1;
  amatrix   tmp2;
  pamatrix  qr_k;
  pavector  rhat_k;
  uint      k = *kk;

  (void) b;
  (void) q;

  rhat_k = init_sub_avector(&tmp1, rhat, k, 0);
  qr_k = init_sub_amatrix(&tmp2, qr, k, 0, k, 1);

  //On entry to DTRSM bla fehler ist hier in triangularsolve
  // und kommt, wenn bei solve_gmres das kmax zu klein ist
  triangularsolve_amatrix_avector(false, false, false, qr_k, rhat_k);

  uninit_amatrix(qr_k);
  uninit_avector(rhat_k);

  rhat_k = init_sub_avector(&tmp1, rhat, rhat->dim - k, k);
  clear_avector(rhat_k);
  uninit_avector(rhat_k);

  qr_k = init_sub_amatrix(&tmp2, qr, qr->rows, 0, k, 0);
  qreval_amatrix_avector(false, qr_k, tau, rhat);
  uninit_amatrix(qr_k);

  add_avector(1.0, rhat, x);

  init_gmres_special(h2m, polblock, spm, b, x, rhat, q, kk, qr, tau);
//  init_gmres(addeval, matrix, b, x, rhat, q, kk, qr, tau);

}





HEADER_PREFIX uint 
solve_gmres_h2precond_avector(ph2matrix h2m, pamatrix polblock, psparsematrix spm, pcavector b, pavector x, real eps, uint maxiter, uint kmax)
{
  pavector  rhat, r, q, tau;
  pamatrix  qr;
  real      norm, error;
  uint      n, iter, k;

  n = x->dim;

  assert(b->dim == n);

  rhat = new_avector(n);
  r = new_avector(n);
  q = new_avector(n);
  qr = new_amatrix(n, kmax);
  tau = new_avector(kmax);

  norm = norm2_avector(b);
  init_gmres_special(h2m, polblock, spm, b, x, rhat, q, &k, qr, tau);
  error = residualnorm_gmres(rhat, k);
  iter = 0;

  while (error > eps * norm && iter + 1 != maxiter) {
    if (k + 1 >= kmax) {
      finish_gmres_special(h2m, polblock, spm, b, x, rhat, q, &k, qr, tau);
    }

//    step_gmres(addeval_A, A, b, x, rhat, q, &k, qr, tau);
    step_gmres_special(h2m, polblock, spm, b, x, rhat, q, &k, qr, tau);
    error = residualnorm_gmres(rhat, k);

    iter++;
  }
//  finish_gmres(addeval_A, A, b, x, rhat, q, &k, qr, tau);
  finish_gmres_special(h2m, polblock, spm, b, x, rhat, q, &k, qr, tau);

  del_avector(tau);
  del_amatrix(qr);
  del_avector(q);
  del_avector(r);
  del_avector(rhat);

  return iter;
}

HEADER_PREFIX void 
addeval_cond_kernelh2matrix(field alpha, ph2matrix h2km, pamatrix pb, pavector src, pavector trg)
{

    pavector    top_subvec_src, bottom_subvec_src, top_subvec_trg, bottom_subvec_trg;
    uint        pts, pbsize;

    pts = pb->rows;
    pbsize = pb->cols;
    top_subvec_src = new_sub_avector(src, pts, 0);
    bottom_subvec_src = new_sub_avector(src, pbsize, pts);
    top_subvec_trg = new_sub_avector(trg, pts, 0);
    bottom_subvec_trg = new_sub_avector(trg, pbsize, pts);

    addeval_h2matrix_avector(alpha, h2km, top_subvec_src, top_subvec_trg);
    addeval_amatrix_avector(alpha, pb, bottom_subvec_src, top_subvec_trg);
    addevaltrans_amatrix_avector(alpha, pb, top_subvec_src, bottom_subvec_trg);

    del_avector(top_subvec_src);
    del_avector(bottom_subvec_src);
    del_avector(top_subvec_trg);
    del_avector(bottom_subvec_trg);
}


HEADER_PREFIX void 
addeval_cond_kernelh2matrix_precon(field alpha, ph2matrix h2km, pamatrix pb, psparsematrix spm, pavector src, pavector trg)
{

    pavector    copy, top_subvec_src, bottom_subvec_src, top_subvec_trg, bottom_subvec_trg;
    uint        pts, pbsize;

    pts = pb->rows;
    pbsize = pb->cols;

    assert(src->dim == trg->dim);
    copy = new_zero_avector(src->dim);
    addeval_sparsematrix_avector(1.0, spm, src, copy);


    top_subvec_src = new_sub_avector(copy, pts, 0);
    bottom_subvec_src = new_sub_avector(copy, pbsize, pts);
    top_subvec_trg = new_sub_avector(trg, pts, 0);
    bottom_subvec_trg = new_sub_avector(trg, pbsize, pts);

    addeval_h2matrix_avector(alpha, h2km, top_subvec_src, top_subvec_trg);
    addeval_amatrix_avector(alpha, pb, bottom_subvec_src, top_subvec_trg);
    addevaltrans_amatrix_avector(alpha, pb, top_subvec_src, bottom_subvec_trg);

    del_avector(copy);
    del_avector(top_subvec_src);
    del_avector(bottom_subvec_src);
    del_avector(top_subvec_trg);
    del_avector(bottom_subvec_trg);

}

HEADER_PREFIX void 
loadfromtxt_precon(pavector preconVals, pavector preconRowIdx, pavector preconColIdx, uint N, uint n)
{

    uint numPts = preconVals->dim;
    char strVals[150], strRowIdx[150], strColIdx[150];
    sprintf(strVals, "/home/kraemer/Programmieren/gp-emulators/lshape_simulation/precon_txt_tps_square/preconVals_N%d_n%d.txt", N, n);
    sprintf(strRowIdx, "/home/kraemer/Programmieren/gp-emulators/lshape_simulation/precon_txt_tps_square/preconRowIdx_N%d_n%d.txt", N, n);
    sprintf(strColIdx, "/home/kraemer/Programmieren/gp-emulators/lshape_simulation/precon_txt_tps_square/preconColIdx_N%d_n%d.txt", N, n);
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

HEADER_PREFIX psparsematrix 
make_precon_sparse(pavector preconVals, pavector preconRowIdx, pavector preconColIdx, uint N, uint n)
{

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



HEADER_PREFIX pamatrix 
make_precon_full(pavector preconVals, pavector preconRowIdx, pavector preconColIdx, uint N, uint n)
{

    pamatrix spm = new_zero_amatrix(N + 3, N+3);
    for(uint i = 0; i < n*N; i++){
        setentry_amatrix(spm, preconRowIdx->v[i], preconColIdx->v[i], preconVals->v[i]);
    }

    for(uint i = 0; i < 3; i++){
        setentry_amatrix(spm, N + i, N + i, 1.0);
    }
    return spm;
}




