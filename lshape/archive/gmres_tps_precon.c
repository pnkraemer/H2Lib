#include <stdio.h>

#include "basic.h"
#include "h2compression.h"
#include "h2matrix.h"
#include "matrixnorms.h"
#include "parameters.h"
#include "kernelmatrix.h"

#include "kernelfcts.h"
#include "addon_kernelmatrix.h"
#include "addon_krylovsolvers.h"

#include "../lshape/auxiliary.h"



real interpolant(real x, real y){
    return 0.5 * FIELD_RAND() + 0.5;
}

void
make_rhs(pavector rhsvec)
{
    uint    dim;
    uint    i;

    dim = rhsvec->dim;

    for(i = 0; i < dim - 3; i++)
        setentry_avector(rhsvec, i, 0.5 * FIELD_RAND() + 0.5);

    for(i=0; i<3; i++)
        setentry_avector(rhsvec, dim - 3 + i, 0.0);
}

int
main(int argc, char **argv)
{
    init_h2lib(&argc, &argv);

    pkernelmatrix       km;
    pclustergeometry    cg;
    pcluster            root;
    pblock              broot;
    pclusterbasis       cb;
    ph2matrix           Gh1, 
                        Gh2;
    pstopwatch          sw;
    ptruncmode          tm;
    pamatrix            pb, 
                        bigmat;
	pavector            testvec, 
                        res_h2, 
                        res_a; 
    pavector            preconval,
                        preconrow,
                        preconcol;
    psparsematrix       precon_sp;
    pamatrix            precon_am;
    pamatrix            res;
    pavector            rhs,
                        x0;

    uint                *idx;
    uint                points;
    uint                m, lsz;
    uint                dim;
    size_t              sz;
    real                eta;
    real                eps;
    real                t_setup;
    uint                i;
    uint                n;
    uint                iter;
    real                error;

    sw = new_stopwatch();
    points = atoi(argv[1]);     /* number of points*/
    n = atoi(argv[2]);
    m = atoi(argv[3]);          /* interpolation order */
    lsz = 2*m*m;                /* leafsize */
    eps = pow(10.0, -1.0 * m);  /* recompression tolerance */
    eta = 1.0;                  /* admissibility condition */
    assert(points<24000);       /* 24000 is all one can do with 8GB of RAM*/
    dim = 2;                    /* This script is 2D only*/

    (void) printf("\nCreating kernelmatrix object for %u points, interpolation order %u\n",	points, m);
    km = new_kernelmatrix(dim, points, m);
    km->kernel = tps_kernel_2d; /* Choose 2d-kernel*/
    construct_lattice_2d(km, 0);

    (void) printf("Creating clustergeometry, cluster and block tree\n");
    cg = creategeometry_kernelmatrix(km);
    idx = (uint *) allocmem(sizeof(uint) * (points));
    for(i=0; i<points; i++)
        idx[i] = i;
    root = build_adaptive_cluster(cg, points, idx, lsz);
    broot = build_strict_block(root, root, &eta, admissible_2_cluster);


    (void) printf("Creating and filling cluster basis\n");
    cb = build_from_cluster_clusterbasis(root);
    start_stopwatch(sw);
    fill_clusterbasis_kernelmatrix(km, cb);
    t_setup = stop_stopwatch(sw);
    sz = getsize_clusterbasis(cb);


    (void) printf("Creating, filling and recompressing H^2-matrix\n");
    Gh1 = build_from_block_h2matrix(broot, cb, cb);
    start_stopwatch(sw);
    fill_h2matrix_kernelmatrix(km, Gh1);
    tm = new_releucl_truncmode();
    Gh2 = compress_h2matrix_h2matrix(Gh1, false, false, tm, eps);
    t_setup = stop_stopwatch(sw);
    sz = getsize_h2matrix(Gh2);
    (void) printf("\t%.2f seconds\n"
		"\t%.1f MB\n"
		"\t%.1f KB/DoF\n",
		t_setup, sz / 1048576.0, sz / 1024.0 / points);


    (void) printf("Filling reference matrix\n");
    start_stopwatch(sw);
    bigmat = new_amatrix(points + 1 + dim, points + 1 + dim);
    assemble_big_kernelmatrix(km, bigmat);
    t_setup = stop_stopwatch(sw);
    sz = getsize_amatrix(bigmat);
    (void) printf("\t%.2f seconds\n"
        "\t%.1f MB\n"
        "\t%.1f KB/DoF\n",
        t_setup, sz / 1048576.0, sz / 1024.0 / points);
    pb = new_amatrix(points, 1 + dim);
    assemble_pblock(km, pb);


    (void) printf("Loading preconditioner\n");
    start_stopwatch(sw);

    preconval = new_avector(points*n);
    preconrow = new_avector(points*n);
    preconcol = new_avector(points*n);
    loadfromtxt_precon(preconval, preconrow, preconcol, points, n);
    precon_sp = make_precon_sparse(preconval, preconrow, preconcol, points, n);
    precon_am = make_precon_full(preconval, preconrow, preconcol, points, n);
    res = new_zero_amatrix(points + dim + 1, points + dim + 1);
    addmul_amatrix(1.0, 0, bigmat, 0, precon_am, res);
    t_setup = stop_stopwatch(sw);
    sz = getsize_sparsematrix(precon_sp);
    (void) printf("\t%.2f seconds\n"
        "\t%.1f MB\n"
        "\t%.1f KB/DoF\n",
        t_setup, sz / 1048576.0, sz / 1024.0 / points);


/*    
    print_pointset(km->x, points, dim);
    print_amatrix(bigmat);
    print_amatrix(precon_am);
*/


    (void) printf("Computing rel. MVM discrepancy\n");
    testvec = new_avector(points + dim + 1);
    res_h2 = new_zero_avector(points + dim + 1);
    res_a = new_zero_avector(points + dim + 1);
    random_real_avector(testvec);
    rhs = new_zero_avector(points + 1 + dim);
    addeval_sparsematrix_avector(1.0, precon_sp, testvec, rhs);
    addeval_amatrix_avector(1.0, bigmat, rhs, res_a);
    addeval_cond_kernelh2matrix_precon(1.0, Gh2, pb, precon_sp, testvec, res_h2);

 //   addeval_cond_kernelh2matrix(1.0, Gh2, pb, testvec, res_h2);
 //   addeval_amatrix_avector(1.0, bigmat, testvec, res_a);
    add_avector(-1.0, res_h2, res_a);
    real diff = norm2_avector(res_a) / norm2_avector(res_h2);
    printf("\t%.1e\n", diff);








    (void) printf("Computing GMRES\n");
    make_rhs(rhs);
    error = norm2_avector(rhs);
    x0 = new_zero_avector(points + dim + 1);
    iter = solve_gmres_h2precond_avector(Gh2, pb, precon_sp, rhs, x0, 1e-5, 10000, 20);
    (void) printf("\t%u iterations\n", iter);
    addeval_amatrix_avector(-1.0, res, x0, rhs);
    error = norm2_avector(rhs)/error;
    (void) printf("\t%.1e relative error\n\n", error);

/*
    error = norm2_avector(rhs);
    iter = solve_gmres_amatrix_avector(res, rhs, x0, 1e-7, 10000, 20);
    printf("Number of GMRES iterations: %u\n\n", iter);
    addeval_amatrix_avector(-1.0, bigmat, x0, rhs);
    error = norm2_avector(rhs);
    printf("ERROR of GMRES: %.1e\n\n", error);
*/

/*    
    cairo_surface_t *surface_h2mat = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 1920, 1080);
    cairo_t *cr_h2mat = cairo_create(surface_h2mat);
    draw_cairo_h2matrix(cr_h2mat, Gh2, 1, 0);
    cairo_surface_write_to_png (surface_h2mat, "./lshape/figures/h2m.png");
    cairo_destroy(cr_h2mat);
    cairo_surface_destroy(surface_h2mat);
*/

    del_avector(rhs);
    del_avector(x0);
    del_amatrix(res);
    del_avector(preconval);
    del_avector(preconrow);
    del_avector(preconcol);
    del_sparsematrix(precon_sp);
    del_amatrix(precon_am);
    del_amatrix(pb);
    del_avector(testvec);
    del_avector(res_h2);
    del_avector(res_a);
    del_amatrix(bigmat);
    del_kernelmatrix(km);
    del_clustergeometry(cg);
    del_cluster(root);
    del_block(broot);
    del_h2matrix(Gh1);
    del_h2matrix(Gh2);
    del_stopwatch(sw);
    del_truncmode(tm);
    freemem(idx);

    uninit_h2lib();
    return 0;
}







