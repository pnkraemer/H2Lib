

#include "kernelmatrix.h"
#include "kernelfcts.h"
#include "blockkernelmatrix.h"

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

    bool use_h2;
    uint dim, points, m, cpos;
    uint i, j;
    uint leafsize;
    real eta;

    dim = 1;
    points = 10;
    m = 10;
    cpos = 1;
    leafsize = 2 * m * m;
    eta = 1.0;

    /* Define kernelmatrix object */
    km = new_kernelmatrix(dim, points, m, cpos);
    for(i=0; i<points; i++) {
        for(j=0; j<dim; j++){
            km->x[i][j] = FIELD_RAND();	/* Random points in [-1,1]^2 */
        }
    }
    km->kernel = tps_kernel_1d;


    /* Construct block and clusterbasis */
    cg = creategeometry_kernelmatrix(km);
    idx = (uint *) allocmem(sizeof(uint) * points);
    for(i=0; i<points; i++)
        idx[i] = i;
    root = build_adaptive_cluster(cg, points, idx, leafsize);
    broot = build_strict_block(root, root, &eta, admissible_2_cluster);
    cb = build_from_cluster_clusterbasis(root);
    fill_clusterbasis_kernelmatrix(km, cb);

    /* Define blockkernelmatrix object */
    use_h2 = true;
    bkm_h2 = build_from_kernelmatrix_blockkernelmatrix(km, use_h2, broot, cb);
    use_h2 = false;
    bkm_a = build_from_kernelmatrix_blockkernelmatrix(km, use_h2, 0, 0);

    /* Check MVM */
    tvec = new_zero_avector(bkm_a->rows);
    random_avector(tvec);
    res_a = new_zero_avector(bkm_a->rows);
    res_h2 = new_zero_avector(bkm_a->rows);
    addeval_blockkernelmatrix_avector(1.0, bkm_h2, tvec, res_h2);
    addeval_blockkernelmatrix_avector(1.0, bkm_a, tvec, res_a);
    print_avector(res_a);
    print_avector(res_h2);




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
