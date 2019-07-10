/*
* 
* In this file we attempt to cleverly construct conditionally positive definite kernel matrices
* and their h2-matrix approximation
*
*/
#include "kernelmatrix.h"

#include <stdio.h>

#include "basic.h"
#include "h2compression.h"
#include "h2matrix.h"
#include "matrixnorms.h"
#include "parameters.h"



static field
kernel_tps(const real *xx, const real *yy, void *data)
{
    real norm2;
    (void) data;


    if(xx[0] >=1.0){
    	if(xx[0] >=2.0){
    	    	if(yy[0] >= 1.0)
    	    		return 0.0;
    	    return yy[0];
    	}
    	if(yy[0] >= 1.0)
    	   	return 0.0;
    	return 1.0;
    }

    if(yy[0] >= 1.0){
    	if(yy[0] >= 2.0){
    		return xx[0];
	    }
	    return 1.0;
    }
//    printf("\n\nmade it\n\n");
    assert(xx[0]<1.0 && yy[0] < 1.0);
    norm2 = REAL_SQRT(REAL_SQR(xx[0] - yy[0]));
    if(norm2>0.0){
        return REAL_SQR(norm2)* REAL_LOG(norm2);
    }else{
        return 0.0;
    }
}

static field
kernel_exp(const real *xx, const real *yy, void *data)
{
    real norm2;
    (void) data;



    norm2 = REAL_SQRT(REAL_SQR(xx[0] - yy[0]));
    if(norm2>0.0){
        return REAL_EXP(-norm2);
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
    uint points;
    uint m, leafsize;
    uint *idx;
    uint dim;
    size_t sz;
    real eta;
    real eps;
    real t_setup, norm, error;
    uint i, j;
    uint msize;
    bool cpos;


    sw = new_stopwatch();
    points = atoi(argv[1]); /* number of points*/
    m = atoi(argv[2]);  /* interpolation order */
//    leafsize = m*m;
    leafsize = m;
    eps = pow(10.0, -1.0 * m);  /* recompression tolerance */
    eta = 2.0;
    cpos = false;
    assert(points<24000); /* 24000 is all one can do with 8GB of RAM*/
    (void) printf("Creating kernelmatrix object for %u (+ %u) points, order %u\n",
		points, (uint)(cpos) * (dim + 1), m);
    dim = 1;
    msize = points + (uint)(cpos) * (dim + 1);
    km = new_kernelmatrix(dim, msize, m);
    km->kernel = kernel_tps;
    for(j=0;j<dim;j++){
	    for(i=0; i<points; i++) {
            km->x[i][j] = (real)(i) * 1.0/points; /* Random points in [-1,1]^2 */
        }
	}
    printf("\nfine\n");

    (void) printf("Creating clustergeometry object and cluster tree\n");
    printf("%i", msize);
    cg = creategeometry_kernelmatrix(km);

//    print_pointset(cg->x, msize, dim);
    printf("%u \n\n\n", km->points);
    idx = (uint *) allocmem(sizeof(uint) * (msize));
    for(i=0; i<msize; i++)
        idx[i] = i;
    root = build_adaptive_cluster(cg, msize, idx, leafsize);

    printf("\nCluster info:\n");
    printf("\tnidcs: %i\n\tnsons: %i\n\tdim: %i\n\tndesc: %i\n\ttype: %i\n\n", root->size, root->sons, root->dim, root->desc, root->type);


    broot = build_strict_block(root, root, &eta, admissible_2_cluster);

    printf("\nBlock info:\n");
    printf("\tnrowsons: %i\n\tncolsons: %i\n\tdesc: %i\n\tadmis = %i\n\n", broot->rsons, broot->csons, broot->desc, broot->a);





    cairo_surface_t *surface_block = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 1600, 900);
    cairo_t *cr_block = cairo_create(surface_block);
    draw_cairo_block(cr_block, broot,0);
    cairo_surface_write_to_png (surface_block, "./lshape/figures/block.png");
    cairo_surface_destroy(surface_block);
    cairo_destroy(cr_block);


    /*GET INTO BLOCK */

    


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
		sz / 1048576.0, sz / 1024.0 / (msize));



    start_stopwatch(sw);
    fill_h2matrix_kernelmatrix(km, Gh1);
    t_setup = stop_stopwatch(sw);

/*
    cairo_surface_t *surface_h2mat = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 1600, 900);
    cairo_t *cr_h2mat = cairo_create(surface_h2mat);
    draw_cairo_h2matrix(cr_h2mat, Gh1, true, 0);
    cairo_surface_write_to_png (surface_h2mat, "./lshape/figures/h2mat.png");
    cairo_surface_destroy(surface_h2mat);
    cairo_destroy(cr_h2mat);
*/
    
    (void) printf("Filling reference matrix\n");
    G = new_amatrix(msize, msize);
    start_stopwatch(sw);
    fillN_kernelmatrix(0, 0, km, G);
    t_setup = stop_stopwatch(sw);
    sz = getsize_amatrix(G);
    (void) printf("    %.2f seconds\n"
		"    %.1f MB\n"
		"    %.1f KB/DoF\n",
		t_setup, sz / 1048576.0, sz / 1024.0 / (msize));
    norm = norm2_amatrix(G);
  
   // print_amatrix(G);

    (void) printf("Computing approximation error\n");
    error = norm2diff_amatrix_h2matrix(Gh1, G);
//    error = norm2_h2matrix(Gh1);
    (void) printf("    Relative spectral error %.3e (norm %.1e, error: %.1e)\n",
		error/norm, norm, error);




    
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
		t_setup, sz / 1048576.0, sz / 1024.0 / msize);

    (void) printf("Computing approximation error\n");
    error = norm2diff_amatrix_h2matrix(Gh2, G);
    (void) printf("    Relative spectral error %.3e\n",
		error/norm);

//    print_amatrix(G);

    /* CHECK APPROXIMATION ERROR OF H2MATRIX*/
    pavector testvec = new_avector(msize);
    random_avector(testvec);
    pavector residual_a = new_zero_avector(msize);
    pavector residual_h2 = new_zero_avector(msize);
    addeval_amatrix_avector(1.0, G, testvec, residual_a);
    addeval_h2matrix_avector(1.0, Gh1, testvec, residual_h2);
 //   print_avector(residual_h2);
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
