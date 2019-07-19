/*
 * TODO: blockkernelmatrix support
*/
#include <stdio.h>
#include "amatrix.h"
#include "krylovsolvers.h"

static void
jacobi(void *pdata, pavector r)
{
  pamatrix  A = (pamatrix) pdata;

  diagsolve_amatrix_avector(false, A, r);
}

static void
gauss_seidel(void *pdata, pavector r)
{
  pamatrix  A = (pamatrix) pdata;

  triangularsolve_amatrix_avector(true, false, false, A, r);
}

static void
multi(void *pdata, pavector r)
{
    pamatrix  A = (pamatrix) pdata;

    pavector rcopy = new_zero_avector(r->dim);
    assert(A->cols == rcopy->dim);
    addeval_amatrix_avector(1.0, A, r, rcopy);
    copy_avector(rcopy, r);
    del_avector(rcopy);

}



int
main()
{
  pamatrix  A;
  pavector  b, x;
  pavector  r;
  pamatrix  Id;

  real      eps, norm, error;
  uint      n, kmax;
  uint      iter;
  uint      problems;

  problems = 0;

  n = 47;
  kmax = 5;
  eps = 1e-6;

  A = new_amatrix(n, n);
  b = new_avector(n);
  x = new_avector(n);

  r = new_avector(n);

  (void) printf("Testing conjugate gradient method\n");
  random_spd_amatrix(A, 1.0);
  random_avector(b);
  norm = norm2_avector(b);

  clear_avector(x);
  iter = solve_cg_amatrix_avector(A, b, x, eps, 0);
  copy_avector(b, r);
  addeval_amatrix_avector(-1.0, A, x, r);
  error = norm2_avector(r);
  (void) printf("  %u steps\n"
		"  Residual %.2e (%.2e)", iter, error, error / norm);

  if (iter <= n && error <= eps * norm)
    printf("    Okay\n");
  else {
    printf("    NOT Okay\n");
    problems++;
  }

  (void) printf("Testing preconditioned conjugate gradient method\n");
  random_spd_amatrix(A, 1.0);
  random_avector(b);
  norm = norm2_avector(b);

  clear_avector(x);
  iter = solve_pcg_amatrix_avector(A, jacobi, A, b, x, eps, 0);
  copy_avector(b, r);
  addeval_amatrix_avector(-1.0, A, x, r);
  error = norm2_avector(r);
  (void) printf("  %u steps\n"
		"  Residual %.2e (%.2e)", iter, error, error / norm);

  if (iter <= n && error <= eps * norm)
    printf("    Okay\n");
  else {
    printf("    NOT Okay\n");
    problems++;
  }

  (void) printf("Testing GMRES method\n");
  random_invertible_amatrix(A, 1.0);
  random_avector(b);
  norm = norm2_avector(b);

  clear_avector(x);
  iter = solve_gmres_amatrix_avector(A, b, x, eps, 0, kmax);
  copy_avector(b, r);
  addeval_amatrix_avector(-1.0, A, x, r);
  error = norm2_avector(r);
  (void) printf("  %u steps\n"
		"  Residual %.2e (%.2e)", iter, error, error / norm);

  if (iter <= n && error <= eps * norm)
    printf("    Okay\n");
  else {
    printf("    NOT Okay\n");
    problems++;
  }

  (void) printf("Testing (left) preconditioned GMRES method\n");
  random_invertible_amatrix(A, 1.0);
  random_avector(b);
  copy_avector(b, x);
  gauss_seidel(A, x);
  norm = norm2_avector(x);

  clear_avector(x);
  iter = solve_lpgmres_amatrix_avector(A, gauss_seidel, A, b, x, eps, 0, kmax);
  copy_avector(b, r);
  addeval_amatrix_avector(-1.0, A, x, r);
  gauss_seidel(A, r);
  error = norm2_avector(r);
  (void) printf("  %u steps\n"
		"  Preconditioned residual %.2e (%.2e)",
		iter, error, error / norm);

  if (iter <= n && error <= eps * norm)
    printf("    Okay\n");
  else {
    printf("    NOT Okay\n");
    problems++;
  }

  (void) printf("Testing (right) preconditioned GMRES method\n");
  random_invertible_amatrix(A, 1.0);
  Id = new_identity_amatrix(A->rows, A->cols);
//  random_invertible_amatrix(Id, 1.0);
  random_avector(b);
  norm = norm2_avector(b);

  clear_avector(x);
  iter = solve_rpgmres_amatrix_avector(A, multi, Id, b, x, eps, 0, kmax);
  print_avector(x);
//  gauss_seidel(A, x);
  copy_avector(b, r);
  addeval_amatrix_avector(-1.0, A, x, r);
  error = norm2_avector(r);
  (void) printf("  %u steps\n"
		"  Preconditioned residual %.2e (%.2e)",
		iter, error, error / norm);

  if (iter <= n && error <= eps * norm)
    printf("    Okay (but not the greatest of tests...)\n");
  else {
    printf("    NOT Okay (but not the greatest of tests...)\n");
    problems++;
  }

  del_avector(r);
  del_avector(x);
  del_avector(b);
  del_amatrix(A);
  del_amatrix(Id);

  printf("----------------------------------------\n"
	 "  %u matrices and\n"
	 "  %u vectors still active\n"
	 "  %u errors found\n",
	 getactives_amatrix(), getactives_avector(), problems);

  return (int) problems;
}
