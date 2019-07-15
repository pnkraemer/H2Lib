/*
Run with either:

./lshape/lshape_gmres /home/kraemer/Programmieren/txts/uniform_lshape/ 21 23 15
./lshape/lshape_gmres /home/kraemer/Programmieren/txts/demlow/ 129 56 15

All rights reserved, Nicholas Krämer, 2019
*/

#include <stdio.h>
#include <string.h>

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
#include "tri2d.h" 



field
interpolant(field x, field y){
    return x + y;
}



/* 
real tpsKernel(real coord){
    if(coord > 0.0){
        return coord*coord * log(coord);
    }
    else{
        return 0;
    }
}
*/




void
make_rhs(pavector rhsvec, pclustergeometry cg)
{
    uint    dim;
    uint    i;
    real    pt;

    dim = rhsvec->dim;
    assert((rhsvec->dim) == (cg->nidx + 3));

    for(i = 0; i < dim - 3; i++)
    {
        pt = interpolant(cg->x[i][0], cg->x[i][1]);
        setentry_avector(rhsvec, i, pt);
    }

    for(i=0; i<3; i++)
        setentry_avector(rhsvec, dim - 3 + i, 0.0);

}


static void 
loadfromtxt_mesh(pkernelmatrix km, bool print_yes, char *string)
{
    uint points;
    char filename[50], filepath[250];
    FILE *meshfile;

    points = km->points;
    sprintf(filename, "mesh/mesh_N%u.txt", points);
    strcpy(filepath, string);
    strcat(filepath, filename);

    meshfile = fopen(filepath, "r");
    if(meshfile == NULL){
        printf("Error reading file\n");
        exit(0);
    }
    for(uint i = 0; i < points; i++){
            fscanf(meshfile, "%lf ", &(km->x[i][0]));      
            fscanf(meshfile, "%lf ", &(km->x[i][1]));
            if(print_yes==1)  
               (void) printf("(%.1f, %.1f)\n", km->x[i][0], km->x[i][1]);
    }
    fclose(meshfile);
}

/*
static void
make_random_lshape(pclustergeometry cg, bool printyes)
{
    uint    counter;
    uint    nidx;
    field   x1, x2;

    counter = 0;
    nidx = cg->nidx;

    while(counter < nidx)
    {
        x1 = FIELD_RAND();
        x2 = FIELD_RAND();

        if(x1 < 0 || x2 > 0)
        {
            cg->x[counter][0] = x1;
            cg->x[counter][1] = x2;
            counter++;
            if(printyes == 1)
                printf("(%.1f, %.1f)\n", x1, x2);
        }

    }

}



*/


/*


pamatrix new_kernelamatrix_tpssquare(clustergeometry *clGeom1, pkernelmatrix km, real shift) {

    uint numPts1 = clGeom1->nidx;
    uint numPts1kurz = clGeom1->nidx;
    uint dim1 = clGeom1->dim;
    uint numPts2 = km->points ;
    uint numPts2kurz = km->points;
    uint dim2 = km->dim;
    assert(dim1 == dim2);
    uint dim = dim1;

    pamatrix kernelMtrx = new_amatrix(numPts1 + dim + 1, numPts2 + dim + 1);
    real norm;
    for(uint i = 0; i < numPts1kurz; i++){
        real* X = allocreal(dim);
        for(uint d = 0; d < dim; d++){
            X[d] = (clGeom1->x)[i][d];
        }
        for(uint j = 0; j < numPts2kurz; j++){
            real* Y = allocreal(dim);
            for(uint d = 0; d < dim; d++){
                Y[d] = (km->x)[j][d];
            }
                        norm = 0;
            for(uint d = 0; d < dim; d++){
                norm += (X[d] - Y[d])*(X[d] - Y[d]);
            }
            setentry_amatrix(kernelMtrx, i, j, tpsKernel(sqrt(norm)));
            if(i==j){
                addentry_amatrix(kernelMtrx, i, j, shift);
            }
            freemem(Y);
        }
        setentry_amatrix(kernelMtrx, i, numPts2kurz, 1.0);
        setentry_amatrix(kernelMtrx, i, numPts2kurz + 1, clGeom1->x[i][0]);
        setentry_amatrix(kernelMtrx, i, numPts2kurz + 2, clGeom1->x[i][1]);

        freemem(X);
    }
    
    for(uint idx = 0; idx < numPts2kurz; idx++){
        setentry_amatrix(kernelMtrx, numPts1kurz, idx, 1.0);
        setentry_amatrix(kernelMtrx, numPts1kurz + 1, idx, km->x[idx][0]);
        setentry_amatrix(kernelMtrx, numPts1kurz + 2, idx, km->x[idx][1]);
    }
    return kernelMtrx;
}

*/

static real
rmse_onthefly(pavector sol, pkernelmatrix km, uint num_evalpts){

    pavector    rowofkm;
    uint        dim, points;
    real        rmse;
    uint        i, j;
    real        kij;
    real        pt[2];

    points = km->points;
    dim = km->dim;
    assert(dim == 2);       /* script not ready for higher dimensions yet */
    rowofkm = new_avector(points + 1 + dim);
    assert(rowofkm->dim == sol->dim);


    rmse = 0;
    for(i = 0; i < num_evalpts; i++){

        /* get random evaluation point*/
        do{
            pt[0] = FIELD_RAND();
            pt[1] = FIELD_RAND();
        }while(pt[0]>0 && pt[1]<0);   /* exclude bottom right corner */

        /* assemble row of kernelmatrix*/
        for(j = 0; j < points; j++){
            kij = km->kernel(pt, km->x[j], km->data);
            setentry_avector(rowofkm, j, kij);
        }
        setentry_avector(rowofkm, points, 1.0);
        setentry_avector(rowofkm, points + 1, pt[0]);
        setentry_avector(rowofkm, points + 2, pt[1]);

        kij = dotprod_avector(rowofkm, sol) - interpolant(pt[0], pt[1]);
        rmse += kij*kij;
    }

    del_avector(rowofkm);

    return REAL_SQRT(rmse / (1.0 * num_evalpts));
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
    ph2matrix           Gh1, Gh2;
    pstopwatch          sw;
    ptruncmode          tm;
    pamatrix            pb; 
/*  pavector            preconval,
                        preconrow,
                        preconcol;
*/  psparsematrix       precon_sp;
    pavector            rhs, x0;

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
    bool                print_yes;
    uint                evalpoints;
    real                currentRelError;
    char                filepath[250];
    real                gmres_tol;
    uint                gmres_kk, 
                        gmres_maxit;

    /* Parameters */
    sw = new_stopwatch();
    strcpy(filepath, argv[1]);    /* path to mesh and precon*/
    points = atoi(argv[2]);       /* number of points*/
    n = atoi(argv[3]);            /* number of neighbours */
    m = atoi(argv[4]);            /* interpolation order */
    lsz = 2*m*m;                  /* leafsize prop. to interpolation order */
    eps = pow(10.0, -1.0 * 2*m);  /* recompression tolerance proportional to interpolation order */
    eta = 2.0;                    /* generic admissibility condition */
//  assert(points<24000);         /* 24000 is all one can do with 8GB of RAM*/
    dim = 2;                      /* this script is 2D only*/
    evalpoints = 10000;           /* number of points for RMSE estimate*/
    gmres_tol = 1e-5;
    gmres_maxit = 500;
    gmres_kk = 20;



    (void) printf("\nCreating kernelmatrix object for %u points (%u neighbours), interpolation order %u\n", points, n, m);
    start_stopwatch(sw);
    km = new_kernelmatrix(dim, points, m);
    km->kernel = tps_kernel_2d; /* Choose 2d-kernel*/
    t_setup = stop_stopwatch(sw);
    t_setup = stop_stopwatch(sw);
    (void) printf("\t%.2f seconds\n", t_setup);

    (void) printf("Loading mesh from txt file\n");
    start_stopwatch(sw);
    print_yes = 0;
    loadfromtxt_mesh(km, print_yes, filepath);
    t_setup = stop_stopwatch(sw);
    (void) printf("\t%.2f seconds\n", t_setup);




    (void) printf("Creating clustergeometry, cluster and block tree\n");
    start_stopwatch(sw);
    cg = creategeometry_kernelmatrix(km);
    idx = (uint *) allocmem(sizeof(uint) * (points));
    for(i=0; i<points; i++)
        idx[i] = i;
    root = build_adaptive_cluster(cg, points, idx, lsz);
    broot = build_strict_block(root, root, &eta, admissible_2_cluster);
    t_setup = stop_stopwatch(sw);
    (void) printf("\t%.2f seconds\n", t_setup);


    (void) printf("Creating and filling cluster basis\n");
    start_stopwatch(sw);
    cb = build_from_cluster_clusterbasis(root);
    start_stopwatch(sw);
    fill_clusterbasis_kernelmatrix(km, cb);
    t_setup = stop_stopwatch(sw);
    sz = getsize_clusterbasis(cb);
    t_setup = stop_stopwatch(sw);
    (void) printf("\t%.2f seconds\n", t_setup);

    (void) printf("Loading preconditioner\n");
    start_stopwatch(sw);
//    preconval = new_avector(points*n);
//    preconrow = new_avector(points*n);
//    preconcol = new_avector(points*n);
    precon_sp = loadfromtxt_precon(points, n, filepath);
//    precon_sp = make_precon_sparse(preconval, preconrow, preconcol, points, n);
    t_setup = stop_stopwatch(sw);
    sz = getsize_sparsematrix(precon_sp);
    (void) printf("\t%.2f seconds\n"
        "\t%.1f MB\n"
        "\t%.1f KB/DoF\n",
        t_setup, sz / 1048576.0, sz / 1024.0 / points);



    (void) printf("Creating, filling and recompressing H^2-matrix (and polblock)\n");
    start_stopwatch(sw);
    Gh1 = build_from_block_h2matrix(broot, cb, cb);
    fill_h2matrix_kernelmatrix(km, Gh1);
    tm = new_releucl_truncmode();
    Gh2 = compress_h2matrix_h2matrix(Gh1, false, false, tm, eps);
	pb = new_amatrix(points, 1 + dim);
    assemble_pblock(km, pb);
    t_setup = stop_stopwatch(sw);
    sz = getsize_h2matrix(Gh2);
    (void) printf("\t%.2f seconds\n"
        "\t%.1f MB\n"
        "\t%.1f KB/DoF\n",
        t_setup, sz / 1048576.0, sz / 1024.0 / points);


/*  (void) printf("Filling reference matrix\n");
    start_stopwatch(sw);
    bigmat = new_amatrix(points + 1 + dim, points + 1 + dim);
    assemble_big_kernelmatrix(km, bigmat);
    t_setup = stop_stopwatch(sw);
    sz = getsize_amatrix(bigmat);
    (void) printf("\t%.2f seconds\n"
        "\t%.1f MB\n"
        "\t%.1f KB/DoF\n",
        t_setup, sz / 1048576.0, sz / 1024.0 / points);
*/






    (void) printf("Computing GMRES\n");
    start_stopwatch(sw);
    rhs = new_zero_avector(points + 1 + dim);
    make_rhs(rhs, cg);
    error = norm2_avector(rhs);
    x0 = new_zero_avector(points + dim + 1);
    iter = solve_gmres_h2precond_avector(Gh2, pb, precon_sp, rhs, x0, gmres_tol, gmres_maxit, gmres_kk);
    addeval_cond_kernelh2matrix_precon(-1.0, Gh2, pb, precon_sp, x0, rhs);
    error = norm2_avector(rhs)/error;
    pavector sol = new_zero_avector(points + dim + 1);
    addeval_sparsematrix_avector(1.0, precon_sp, x0, sol);
    t_setup = stop_stopwatch(sw);
    (void) printf("\t%.2f seconds\n", t_setup);
    (void) printf("\t%u iterations\n", iter);
    (void) printf("\t%.1e relative error\n", error);








    (void) printf("Approximating RMSE");
    start_stopwatch(sw);
    t_setup = stop_stopwatch(sw);
    (void) printf("\n\t%.2f seconds\n", t_setup);
    currentRelError = rmse_onthefly(sol, km, evalpoints);
    (void) printf("\t%.1e rmse\n\n", currentRelError);


/*  error = norm2_avector(rhs);
    x0 = new_zero_avector(points + dim + 1);
    iter = solve_gmres_amatrix_avector(res, rhs, x0, 1e-5, 10000, 20);
    printf("Number of GMRES iterations: %u\n\n", iter);
    addeval_amatrix_avector(-1.0, res, x0, rhs);
    error = norm2_avector(rhs);
    printf("ERROR of GMRES: %.1e\n\n", error);
*/

    
/*  cairo_surface_t *surface_h2mat = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 1920, 1080);
    cairo_t *cr_h2mat = cairo_create(surface_h2mat);
    draw_cairo_h2matrix(cr_h2mat, Gh2, 1, 0);
    cairo_surface_write_to_png (surface_h2mat, "./lshape/figures/h2m.png");
    cairo_destroy(cr_h2mat);
    cairo_surface_destroy(surface_h2mat);
*/

    del_avector(rhs);
    del_avector(x0);
    del_avector(sol);
//    del_amatrix(res);
//    del_avector(preconval);
//    del_avector(preconrow);
//    del_avector(preconcol);
    del_sparsematrix(precon_sp);
//    del_amatrix(precon_am);
    del_amatrix(pb);
//  del_avector(testvec);
//  del_avector(res_h2);
//  del_avector(res_a);
//    del_amatrix(bigmat);
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







