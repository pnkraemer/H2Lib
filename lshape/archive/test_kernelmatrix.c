
#include "kernelmatrix.h"

#include <stdio.h>

#include "basic.h"
#include "h2compression.h"
#include "h2matrix.h"
#include "matrixnorms.h"
#include "parameters.h"


static field
kernel_exp(const real *xx, const real *yy, void *data)
{
    real norm2;

    (void) data;

    norm2 = REAL_SQR(xx[0] - yy[0]) + REAL_SQR(xx[1] - yy[1]);

    return REAL_EXP(-norm2);
}

static field
kernel_matern15(const real *xx, const real *yy, void *data)
{
    real norm2;

    (void) data;

    norm2 = REAL_SQR(xx[0] - yy[0]) + REAL_SQR(xx[1] - yy[1]);

    return (1.0 + REAL_SQR(3) * norm2) * REAL_EXP(- REAL_SQR(3) * norm2);

}

static field
kernel_matern25(const real *xx, const real *yy, void *data)
{
    real norm2;

    (void) data;

    norm2 = REAL_SQR(xx[0] - yy[0]) + REAL_SQR(xx[1] - yy[1]);

    return (1 + REAL_SQR(5) * norm2 + 5.0 * norm2 * norm2/3.0  ) * REAL_EXP(- REAL_SQR(5) * norm2);

}

static field
kernel_tps(const real *xx, const real *yy, void *data)
{
    real norm2;

    (void) data;
    if(fabs(xx[0] * yy[0])>999999.0){
        return 0.0;
    }

    if(xx[0] == -1000.0){
        return 1.0;
    }
    if(xx[0] == 1000.0){
        return yy[0];
    }
    if(xx[0] == 2000.0){
        return yy[1];
    }

    if(yy[0] == -1000.0){
        return 1.0;
    }
    if(yy[0] == 1000.0){
        return xx[0];
    }
    if(yy[0] == 2000.0){
        return xx[1];
    }

    norm2 = REAL_SQR(xx[0] - yy[0]) + REAL_SQR(xx[1] - yy[1]);
    if(norm2>0.0){
        return norm2* REAL_LOG(norm2);
    }else{
        return 0.0;
    }
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


int
main(int argc, char **argv)
{
    pkernelmatrix km;
    pclustergeometry cg;
    pcluster root;
    pblock broot;
    pclusterbasis cb;
    ph2matrix Gh1, Gh2;
    pamatrix G;
    pstopwatch sw;
    ptruncmode tm;
    uint points;
    uint m, leafsize;
    uint *idx;
    uint dim;
    uint pblocksize;
    size_t sz;
    real eta;
    real eps;
    real t_setup, norm, error;
    uint i, j;

    init_h2lib(&argc, &argv);

    sw = new_stopwatch();
    points = atoi(argv[1]); /* number of points*/
    m = atoi(argv[2]);  /* interpolation order */
    leafsize = 2*m*m;
    eps = pow(10.0, -1.0 * m);  /* recompression tolerance */
    eta = 2.0;
    assert(points<24000); /* 24000 is all one can do with 8GB of RAM*/
    (void) printf("Creating kernelmatrix object for %u points, order %u\n",
		points, m);
    dim = 2;
    pblocksize = 1 + dim;
    km = new_kernelmatrix(dim, points + pblocksize, m);
    km->kernel = kernel_tps;
    for(i=0; i<points; i++) {
        for(j=0;j<dim;j++){
          km->x[i][j] = FIELD_RAND(); /* Random points in [-1,1]^2 */
        }
    }
    for(j=0;j<dim;j++){
            km->x[points][j] = -1000.0; /* becomes coordinate in cdm */
        for(i=points + 1; i<points + pblocksize; i++) {
                km->x[i][j] = (i - points) * 1000.0; /* becomes coordinate in cdm */
        }
    }


    (void) printf("Creating clustergeometry object and cluster tree\n");
    cg = creategeometry_kernelmatrix(km);
    idx = (uint *) allocmem(sizeof(uint) * (points + pblocksize));
    for(i=0; i<points + pblocksize; i++)
        idx[i] = i;
    root = build_adaptive_cluster(cg, points + pblocksize, idx, leafsize);
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
		sz / 1048576.0, sz / 1024.0 / (points + pblocksize));

    start_stopwatch(sw);
    fill_h2matrix_kernelmatrix(km, Gh1);
    t_setup = stop_stopwatch(sw);
    
    (void) printf("Filling reference matrix\n");
    G = new_amatrix(points + pblocksize, points + pblocksize);
    start_stopwatch(sw);
    fillN_kernelmatrix(0, 0, km, G);
    t_setup = stop_stopwatch(sw);
    sz = getsize_amatrix(G);
    (void) printf("    %.2f seconds\n"
		"    %.1f MB\n"
		"    %.1f KB/DoF\n",
		t_setup, sz / 1048576.0, sz / 1024.0 / (points + pblocksize));
    print_pointset(km->x, points + pblocksize, 2);
    print_amatrix(G);
    norm = norm2_amatrix(G);

    (void) printf("Computing approximation error\n");
    error = norm2diff_amatrix_h2matrix(Gh1, G);
    (void) printf("    Relative spectral error %.3e (norm %.1e)\n",
		error/norm, norm);

    
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

    /* CHECK APPROXIMATION ERROR OF H2MATRIX*/
    pavector testvec = new_avector(points + pblocksize);
    random_avector(testvec);
    pavector residual_a = new_zero_avector(points + pblocksize);
    pavector residual_h2 = new_zero_avector(points + pblocksize);
    addeval_amatrix_avector(1.0, G, testvec, residual_a);
    addeval_h2matrix_avector(1.0, Gh1, testvec, residual_h2);
    add_avector(-1.0, residual_a, residual_h2);
    real diff = norm2_avector(residual_h2) / norm2_avector(residual_a);
    printf("\nMVM discrepancy %.3e\n\n", diff);
    del_avector(testvec);
    del_avector(residual_a);
    del_avector(residual_h2);



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
