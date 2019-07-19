

#include "kernelmatrix.h"
#include "kernelfcts.h"
#include "blockkernelmatrix.h"
#include "krylovsolvers.h"


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
main(int argc, char **argv)
{
    init_h2lib(&argc, &argv);

    pkernelmatrix km;
    pblockkernelmatrix bkm_a, bkm_h2;
    pclustergeometry cg;
    pcluster root;
    pblock broot;
    pclusterbasis cb;
    uint *idx;
    pavector tvec, res_a, res_h2;
    pamatrix bmat;
    pamatrix precon;

    bool use_h2;
    uint dim, points, m, cpos;
    uint i, j;
    uint leafsize;
    real eta;
    uint gmres_maxit;
    uint gmres_k;
    real gmres_tol;
    uint iter;
    real length;
    real discrep;

    dim = 3;
    points = 100;
    m = 10;
    cpos = 0;
    leafsize = 2 * m * m;
    eta = 2.0;

    /* Define kernelmatrix object */
    km = new_kernelmatrix(dim, points, m, cpos);
    for(i=0; i<points; i++) {
    length = 0.0;
        for(j=0; j<dim; j++){
            km->x[i][j] = FIELD_RAND();	/* Random points in [-1,1]^2 */
            length += km->x[i][j]*km->x[i][j];
        }
        length = REAL_SQRT(length);
        for(j=0; j<dim; j++){
            km->x[i][j] = km->x[i][j] / length; /* Get to sphere */
        }
    }
    km->kernel = tps_kernel_s2;


    /* Construct block and clusterbasis */
    cg = creategeometry_kernelmatrix(km);
    idx = (uint *) allocmem(sizeof(uint) * points);
    for(i=0; i<points; i++){
        idx[i] = i;
    }
    root = build_adaptive_cluster(cg, points, idx, leafsize);
    broot = build_strict_block(root, root, &eta, admissible_2_cluster);
    cb = build_from_cluster_clusterbasis(root);
    fill_clusterbasis_kernelmatrix(km, cb);

    /* Define blockkernelmatrix object */
    use_h2 = true;
    bkm_h2 = build_from_kernelmatrix_blockkernelmatrix(km, use_h2, broot, cb);
    use_h2 = false;
    bkm_a = build_from_kernelmatrix_blockkernelmatrix(km, use_h2, 0, 0);
    bmat = new_zero_amatrix(bkm_a->dof, bkm_a->dof);
    convert_blockkernelmatrix_amatrix(bkm_a, bmat);

    /* Check MVM */
    tvec = new_zero_avector(bkm_a->dof);
    res_a = new_zero_avector(bkm_a->dof);
    res_h2 = new_zero_avector(bkm_a->dof);
    random_avector(tvec);
    addeval_blockkernelmatrix_avector(1.0, bkm_a, tvec, res_a);
    addeval_blockkernelmatrix_avector(1.0, bkm_h2, tvec, res_h2);
    add_avector(-1.0, res_a, res_h2);
    discrep = norm2_avector(res_h2) / norm2_avector(res_a);
    printf("\nDiscrepancy between bkms: %.1f\n", discrep);

    clear_avector(res_a);
    clear_avector(res_h2);
    addeval_amatrix_avector(1.0, bmat, tvec, res_a);
    addeval_blockkernelmatrix_avector(1.0, bkm_h2, tvec, res_h2);
    add_avector(-1.0, res_a, res_h2);
    discrep = norm2_avector(res_h2) / norm2_avector(res_a);
    printf("\nDiscrepancy bkm & amatrix: %.1f\n", discrep);

    /* Check GMRES */
    gmres_maxit = 1000;
    gmres_k = 25;
    gmres_tol = 1e-5;
    precon = new_identity_amatrix(bkm_a->dof, bkm_a->dof);
    clear_avector(res_a);
    iter = solve_gmres_amatrix_avector(bmat, tvec, res_a, gmres_tol, gmres_maxit, gmres_k);
    printf("Number of iterations (amatrix, no precon):\n\titer = %u\n", iter);
    clear_avector(res_h2);
    iter = solve_gmres_blockkernelmatrix_avector(bkm_a, tvec, res_h2, gmres_tol, gmres_maxit, gmres_k);
    printf("Number of iterations (blockkernelmatrix, no precon):\n\titer = %u\n", iter);
    clear_avector(res_a);
    iter = solve_lpgmres_blockkernelmatrix_avector(bkm_a, multi, precon, tvec, res_a, gmres_tol, gmres_maxit, gmres_k);
    printf("Number of iterations (blockkernelmatrix, left identity precon):\n\titer = %u\n", iter);
    clear_avector(res_h2);
    iter = solve_rpgmres_blockkernelmatrix_avector(bkm_a, multi, precon, tvec, res_h2, gmres_tol, gmres_maxit, gmres_k);
    printf("Number of iterations (blockkernelmatrix, right identity precon):\n\titer = %u\n", iter);
    add_avector(-1.0, res_a, res_h2);
    discrep = norm2_avector(res_h2) / norm2_avector(res_a);
    printf("\nDiscrepancy between last two:\n\tdisc = %.1f\n", discrep);




    del_amatrix(precon);
    del_amatrix(bmat);
    del_avector(res_a);
    del_avector(res_h2);
    del_avector(tvec);
    del_kernelmatrix(km);
    del_clustergeometry(cg);
    del_cluster(root);
    del_block(broot);
//    if(bkm_h2->h2kmat && bkm_a->h2kmat){
//        del_clusterbasis(cb);
//    }
    del_blockkernelmatrix(bkm_h2);
    del_blockkernelmatrix(bkm_a);
    freemem(idx);



    uninit_h2lib();

    return 0;
}
