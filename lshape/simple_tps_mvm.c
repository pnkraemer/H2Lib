/*
* Purpose of this script:
* provide basic checks for blockwise mvm of cond. pos. definite kernels
* with H2 matrices in terms of storage and relative errors
*
*
*/




#include <stdio.h>

#include "basic.h"
#include "h2compression.h"
#include "h2matrix.h"
#include "matrixnorms.h"
#include "parameters.h"
#include "kernelmatrix.h"

#include "kernelfcts.h"
#include "addon_kernelmatrix.h"


static void
print_pointset(real **pts, uint npts, uint dim);

int
main(int argc, char **argv)
{
    init_h2lib(&argc, &argv);
    pkernelmatrix km;
    pclustergeometry cg;
    pcluster root;
    pblock broot;
    pclusterbasis cb;
    ph2matrix Gh1, Gh2;
    pamatrix G;
    pstopwatch sw;
    ptruncmode tm;
    pamatrix pb, bigmat;
	pavector testvec, tv2, tv3;  
    uint points;
    uint m, lsz;
    uint *idx;
    uint dim;
    size_t sz;
    real eta;
    real eps;
    real t_setup, norm, error;
    uint i, j;


    sw = new_stopwatch();
    points = atoi(argv[1]);     /* number of points*/
    m = atoi(argv[2]);          /* interpolation order */
    lsz = 2*m*m;                /*  leafsize */
    eps = pow(10.0, -1.0 * m);  /* recompression tolerance */
    eta = 2.0;                  /* admissibility condition */
    assert(points<24000);       /* 24000 is all one can do with 8GB of RAM*/

    (void) printf("Creating kernelmatrix object for %u points, order %u\n",	points, m);
    dim = 1;
    km = new_kernelmatrix(dim, points, m);
    km->kernel = tps_kernel_1d;
    for(i=0; i<points; i++){
        for(j=0;j<dim;j++){
            km->x[i][j] = FIELD_RAND(); /* Random points in [-1,1]^2 */
        }
    }


    (void) printf("Creating clustergeometry object and cluster tree\n");
    cg = creategeometry_kernelmatrix(km);
    idx = (uint *) allocmem(sizeof(uint) * (points));
    for(i=0; i<points; i++)
        idx[i] = i;
    root = build_adaptive_cluster(cg, points, idx, lsz);
    broot = build_strict_block(root, root, &eta, admissible_2_cluster);



/*
    cairo_surface_t *surface_h2mat = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 1920, 1080);
    cairo_t *cr_h2mat = cairo_create(surface_h2mat);
    draw_cairo_block(cr_h2mat, broot,0);
    cairo_surface_write_to_png (surface_h2mat, "./lshape/figures/block.png");
    cairo_surface_destroy(surface_h2mat);
    cairo_destroy(cr_h2mat);
*/


    (void) printf("Creating and filling cluster basis\n");
    cb = build_from_cluster_clusterbasis(root);
    start_stopwatch(sw);
    fill_clusterbasis_kernelmatrix(km, cb);
    t_setup = stop_stopwatch(sw);
    sz = getsize_clusterbasis(cb);

    (void) printf("Creating and filling H^2-matrix\n");
    Gh1 = build_from_block_h2matrix(broot, cb, cb);
    sz = getsize_h2matrix(Gh1);
    (void) printf("    %.1f MB\n"
		"    %.1f KB/DoF\n",
		sz / 1048576.0, sz / 1024.0 / (points));

    start_stopwatch(sw);
    fill_h2matrix_kernelmatrix(km, Gh1);
    t_setup = stop_stopwatch(sw);
    
    (void) printf("Filling reference matrix\n");
    G = new_amatrix(points, points);
    start_stopwatch(sw);
    fillN_kernelmatrix(0, 0, km, G);
    t_setup = stop_stopwatch(sw);
    sz = getsize_amatrix(G);
    (void) printf("    %.2f seconds\n"
		"    %.1f MB\n"
		"    %.1f KB/DoF\n",
		t_setup, sz / 1048576.0, sz / 1024.0 / (points));
    norm = norm2_amatrix(G);

    (void) printf("Recompressing H^2-matrix, eps=%g\n",
		eps);
    start_stopwatch(sw);
    tm = new_releucl_truncmode();
    Gh2 = compress_h2matrix_h2matrix(Gh1, false, false, tm, eps);
    t_setup = stop_stopwatch(sw);
    sz = getsize_h2matrix(Gh2);
    (void) printf("    %.2f seconds\n"
		"    %.1f MB\n"
		"    %.1f KB/DoF\n",
		t_setup, sz / 1048576.0, sz / 1024.0 / points);


    (void) printf("Computing approximation error\n");
    error = norm2diff_amatrix_h2matrix(Gh2, G);
    (void) printf("    Relative spectral error %.3e\n",
		error/norm);

    (void) printf("Filling conditional part\n");
    bigmat = new_amatrix(points + 1 + dim, points + 1 + dim);
    assemble_big_kernelmatrix(km, bigmat);
    pb = assemble_pblock(km);




    testvec = new_avector(points + dim + 1);
    tv2 = new_avector(points + dim + 1);
    tv3 = new_avector(points + dim + 1);
    random_avector(testvec);
    copy_avector(testvec, tv2);
    copy_avector(testvec, tv3);

    addeval_cond_kernelh2matrix(1.0, Gh1, pb, testvec, tv2);
    addeval_amatrix_avector(1.0, bigmat, testvec, tv3);

    add_avector(-1.0, tv2, tv3);
    real diff = norm2_avector(tv3) / norm2_avector(tv2);
    printf("\nMVM discrepancy %.3e\n\n", diff);






    


    /* CHECK APPROXIMATION ERROR OF H2MATRIX*/
/*    pavector testvec = new_avector(points);
    random_avector(testvec);
    pavector residual_a = new_zero_avector(points);
    pavector residual_h2 = new_zero_avector(points);
    addeval_amatrix_avector(1.0, G, testvec, residual_a);
    addeval_h2matrix_avector(1.0, Gh1, testvec, residual_h2);
    add_avector(-1.0, residual_a, residual_h2);
    real diff = norm2_avector(residual_h2) / norm2_avector(residual_a);
    printf("\nMVM discrepancy %.3e\n\n", diff);
    del_avector(testvec);
    del_avector(residual_a);
    del_avector(residual_h2);

*/


    del_amatrix(pb);
    del_avector(testvec);
    del_avector(tv2);
    del_avector(tv3);
    del_amatrix(bigmat);


    del_kernelmatrix(km);
    del_clustergeometry(cg);
    del_cluster(root);
    del_block(broot);
    del_h2matrix(Gh1);
    del_h2matrix(Gh2);
    del_amatrix(G);
    del_stopwatch(sw);
    del_truncmode(tm);
    freemem(idx);
  

    uninit_h2lib();
    return 0;
}





static void
print_pointset(real **pts, uint npts, uint dim)
{
    printf("\n");
    for(uint i = 0; i < npts; i++){
        printf("(");
        printf("%.1e", pts[i][0]);
        for(uint j = 1; j < dim; j++){
            printf(", %.1e", pts[i][j]);
        }
        printf(")\n");
    }
    printf("\n");
}





