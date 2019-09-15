/*
 * All rights reserved, Nicholas Krämer, 2019
 * Use it as:
 * 
 * ./lshape/clean_lshape_gmres <path_to_mesh> <path_to_precon> <path_to_evalpts> <local_hmatrix_rank> <gmres_acc>
 * 
 * THIS SCRIPT IS 2D ONLY
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
#include "blockkernelmatrix.h"
#include "krylovsolvers.h"
#include "tri2d.h"



static void
print_avector_from_to(pcavector vec, uint from, uint to)
{
    uint i;
    assert(from<to);
    printf("\n( ");
    for(i = from; i < to; i++){
        printf("%e, ", getentry_avector(vec, i));
    }
    printf(")\n ");
}


static void
multi_sp(void *pdata, pavector r)
{
    psparsematrix  A = (psparsematrix) pdata;
    pavector rcopy = new_zero_avector(r->dim);
    assert(A->cols == rcopy->dim);
    addeval_sparsematrix_avector(1.0, A, r, rcopy);
    clear_avector(r);
    copy_avector(rcopy, r);
    del_avector(rcopy);
}


field
interpolant(field x, field y){
	return x;
}

void
make_rhs(pavector rhsvec, pclustergeometry cg)
{
	uint    dim;
	uint    i;
	real    pt;

	dim = rhsvec->dim;
	assert((rhsvec->dim) == (cg->nidx + 3));

	for(i = 0; i < dim - 3; i++){
	   pt = interpolant(cg->x[i][0], cg->x[i][1]);
	   setentry_avector(rhsvec, i, pt);
	}

	for(i=0; i<3; i++){
	   setentry_avector(rhsvec, dim - 3 + i, 0.0);
    }
}

static void
loadfromtxt_mesh(pkernelmatrix km, bool print_yes, char *filepath)
{
	uint points;
	FILE *meshfile;

	points = km->points;

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



HEADER_PREFIX psparsematrix
loadfromtxt_precon(uint N, char *filepath)
{
    real val;
    real colidx, rowidx;
    uint i;
    pavector preconval, preconrow, preconcol;
    psparsepattern sp;
    psparsematrix spm;
    char buffer[100];
    uint num_elements_precon;
    uint count;

    FILE *myFile;
    myFile = fopen(filepath, "r");

    if(myFile == NULL){
        printf("Error reading file\n");
        exit(0);
    }
    // Throw away the explanational line, i.e. the first one
    fgets(buffer, 100, myFile);

    // Count number of elements
    count = 0;
	while (fscanf(myFile, "%lf", &val) != EOF)
	{
		count++;
	}
	num_elements_precon = (uint)count / 3.0;
	fclose(myFile);


    myFile = fopen(filepath, "r");
    if(myFile == NULL){
        printf("Error reading file\n");
        exit(0);
    }
    fgets(buffer, 100, myFile);


    preconval = new_avector(num_elements_precon);
    preconrow = new_avector(num_elements_precon);
    preconcol = new_avector(num_elements_precon);


    sp = new_sparsepattern(N + 3, N + 3);
    i = 0;
    for(i = 0; i < num_elements_precon; i++){
    	fscanf(myFile, "%lf", &val);
    	fscanf(myFile, "%lf", &rowidx);
    	fscanf(myFile, "%lf", &colidx);

        setentry_avector(preconval, i, val);
        setentry_avector(preconrow, i, (uint)rowidx);
        setentry_avector(preconcol, i, (uint)colidx);
        addnz_sparsepattern(sp, rowidx, colidx);
    }

    for(uint i = 0; i < 3; i++){
        addnz_sparsepattern(sp, N + i, N + i);
    }
    fclose(myFile);

    spm = new_zero_sparsematrix(sp);
    for(i = 0; i < num_elements_precon; i++){
        val = getentry_avector(preconval, i);
        rowidx = getentry_avector(preconrow, i);
        colidx = getentry_avector(preconcol, i);
        setentry_sparsematrix(spm, rowidx, colidx, val);
    }
    for(uint i = 0; i < 3; i++){
        setentry_sparsematrix(spm, N+i, N+i, 0.0);
    }

    del_sparsepattern(sp);
    del_avector(preconval);
    del_avector(preconrow);
    del_avector(preconcol);
    return spm;
}


INLINE_PREFIX real /* Not happy with the name of the function */
rmse_onthefly(pavector sol, pkernelmatrix km, char *path_to_evalpts){

	pavector    rowofkm;
	uint        dim, points;
	real        rmse;
	uint        i, j;
	real        kij;
	real        pt[2];
	FILE 		*centers;
	uint 		count;
	bool 		print_yes;
	pkernelmatrix kmeval;
	real 		dum;
	uint num_evalpts;


	dim = 2;

	 //Reading number of evaluation points from file...
	centers = fopen(path_to_evalpts, "r");
	count = 0;
	while (fscanf(centers, "%lf", &dum) != EOF)
	{
		count++;
	}
	num_evalpts = (long int)count/dim;
	fclose(centers);


	kmeval = new_kernelmatrix(dim, num_evalpts, km->m, km->cpos);
	print_yes = 0;
	loadfromtxt_mesh(kmeval, print_yes, path_to_evalpts);


	points = km->points;
	rowofkm = new_avector(points + 1 + dim);
	assert(rowofkm->dim == sol->dim);


	rmse = 0;
	for(i = 0; i < num_evalpts; i++){

		// get evaluation point
		pt[0] = kmeval->x[i][0];
		pt[1] = kmeval->x[i][1];

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
	del_kernelmatrix(kmeval);

	return REAL_SQRT(rmse / (1.0 * num_evalpts));
}

/* Usage:
 *
 * ./lshape/clean_lshape_gmres <path_to_mesh> <path_to_precon> <path_to_evalpts_for_rmse> <local_hmatrix_rank> <gmres_acc>
 *
 */
int
main(int argc, char **argv)
{
	init_h2lib(&argc, &argv);

	pkernelmatrix       km;
	pclustergeometry    cg;
	pcluster            root;
	pblock              broot;
	pclusterbasis       cb;
	pstopwatch          sw;
	psparsematrix       precon_sp;
	pavector            rhs, x0;
	pavector            sol;
	pamatrix            fullmat;
	pavector            testvec, res_a, res_h2;
	pamatrix            precon_full, kp;
	pavector            sol2;
	pblockkernelmatrix	bkm, bkm_uncompressed, bkm_full;
	real                mvm_error;

	uint                *idx;
	uint                points;
	uint                m, lsz, cpos;
	uint                dim;
	size_t              sz;
	real                eta;
	real                t_setup;
	uint                i;
	uint                iter;
	real                error;
	bool                print_yes;
	real                currentRelError;
	char                path_to_mesh[250], path_to_precon[250], path_to_evalpts[250];
	real                gmres_tol;
	uint                gmres_kk, gmres_maxit;
	real                error_fullgmres;
    real                h2compression_acc;

    FILE 				*centers;
    real 				dummyvar;
    uint 				count;

	/* Parameters */
	sw = new_stopwatch();
	strcpy(path_to_mesh, argv[1]);    /* path to mesh*/
	strcpy(path_to_precon, argv[2]);    /* path to preconditioner*/
	strcpy(path_to_evalpts, argv[3]);    /* path to evaluation pointset*/
	m = atoi(argv[4]);            /* interpolation order */
	assert(m > 0);
	lsz = 2*m*m;                  /* leafsize prop. to interpolation order */
	eta = 2.0;                    /* generic admissibility condition */
//  assert(points<24000);         /* 24000 is all one can do with 8GB of RAM*/
	dim = 2;                      /* this script is 2D only*/
	gmres_tol = atof(argv[5]);    /* curiously: for ~1e-12, any gmres fails */

	// Some internal variables
	gmres_maxit = 5000;
	gmres_kk = 25;
	cpos = 1;
    h2compression_acc = 1e-3 * gmres_tol;


	(void) printf("\nReading number of evaluation points from file...");
	centers = fopen(path_to_mesh, "r");
	count = 0;
	while (fscanf(centers, "%lf", &dummyvar) != EOF)
	{
		count++;
	}
	points = (long int)count/dim;
	fclose(centers);






	(void) printf("\nCreating kernelmatrix object for %u points, interpolation order %u\n", points, m);
	start_stopwatch(sw);
	km = new_kernelmatrix(dim, points, m, cpos);
	km->kernel = tps_kernel_2d; /* Choose 2d-kernel*/
	t_setup = stop_stopwatch(sw);
	(void) printf("\t%.2f seconds\n", t_setup);


	(void) printf("Loading mesh from txt file\n");
	start_stopwatch(sw);
	print_yes = 0;
	loadfromtxt_mesh(km, print_yes, path_to_mesh);
	t_setup = stop_stopwatch(sw);
	(void) printf("\t%.2f seconds\n", t_setup);

	(void) printf("Creating clustergeometry, cluster and block tree\n");
	start_stopwatch(sw);
	cg = creategeometry_kernelmatrix(km);
	idx = (uint *) allocmem(sizeof(uint) * (points));
	for(i=0; i<points; i++)
		idx[i] = i;
	root = build_cluster(cg, points, idx, lsz, H2_ADAPTIVE);
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
	precon_sp = loadfromtxt_precon(points, path_to_precon);
	t_setup = stop_stopwatch(sw);
	sz = getsize_sparsematrix(precon_sp);
	(void) printf("\t%.2f seconds\n"
	"\t%.1f MB\n"
	"\t%.1f KB/DoF\n",
	t_setup, sz / 1048576.0, sz / 1024.0 / points);

	(void) printf("Creating, filling and compressing H^2-matrix (and polblock)\n");
	start_stopwatch(sw);
	bkm_uncompressed = build_from_kernelmatrix_h2_blockkernelmatrix(km, broot, cb);
    bkm = compress_blockkernelmatrix(bkm_uncompressed, h2compression_acc);
//    bkm_full = build_from_kernelmatrix_full_blockkernelmatrix(km);
	t_setup = stop_stopwatch(sw);
	sz = getsize_blockkernelmatrix(bkm);
	(void) printf("\t%.2f seconds\n"
	"\t%.1f MB (%.1f MB for full matrix)\n"
	"\t%.1f KB/DoF\n",
	t_setup, sz / 1048576.0, points*points*8.0/1048576.0, sz / 1024.0 / points);

	(void) printf("Preparing, solving and evaluating linear system\n");
	start_stopwatch(sw);
	rhs = new_zero_avector(points + 1 + dim);
	make_rhs(rhs, cg);
	error = norm2_avector(rhs);
	x0 = new_zero_avector(points + dim + 1);
	iter = solve_rpgmres_blockkernelmatrix_avector(bkm, multi_sp, precon_sp, rhs, x0, gmres_tol, gmres_maxit, gmres_kk);
    //print_avector_from_to(x0, x0->dim-3, x0->dim);
	sol = new_zero_avector(points + dim + 1);
	addeval_sparsematrix_avector(1.0, precon_sp, x0, sol);
	addeval_blockkernelmatrix_avector(-1.0, bkm, sol, rhs);
	error = norm2_avector(rhs)/error;
	t_setup = stop_stopwatch(sw);
	(void) printf("\t%.2f seconds\n", t_setup);
	(void) printf("\t%u GMRES iterations with the H2-matrix\n", iter);
	(void) printf("\t%.12e relative error", error);




	/* Check: if GMRES fails -> what was the issue?
	 * h2approx or bad preconditioning */
	sol2 = new_zero_avector(x0->dim);
	if((error > gmres_tol || error!=error )&& bkm->dof < 65000)  /* Dont want to assemble too large matrices */
//	if(1)
    {
		(void) printf(" ->TOO MUCH ERROR!");
		(void) printf("\n--------------------------------------------------");
		(void) printf("\nChecking reference matrix\n");
		start_stopwatch(sw);

		/* Assemble big matrix */
        bkm_full = build_from_kernelmatrix_full_blockkernelmatrix(km);
		fullmat = new_amatrix(bkm_full->dof, bkm_full->dof);
		convert_blockkernelmatrix_amatrix(bkm_full, fullmat);

		/* Check approximation quality of H2Matrix */
		testvec = new_avector(fullmat->rows);
		res_a = new_zero_avector(fullmat->rows);
		res_h2 = new_zero_avector(fullmat->rows);
		random_avector(testvec);
		addeval_amatrix_avector(1.0, fullmat, testvec, res_a);
		addeval_blockkernelmatrix_avector(1.0, bkm, testvec, res_h2);
		add_avector(-1.0, res_a, res_h2);
		mvm_error = norm2_avector(res_h2)/norm2_avector(res_a);


		/* Check GMRES functionality */
		precon_full = new_zero_amatrix(fullmat->rows, fullmat->rows);
		kp = new_zero_amatrix(fullmat->rows, fullmat->rows);
		add_sparsematrix_amatrix(1.0, 0, precon_sp, precon_full);
		addmul_amatrix(1.0, 0, fullmat, 0, precon_full, kp);
		clear_avector(x0);
		make_rhs(rhs, cg);
		error = norm2_avector(rhs);
		iter = solve_rpgmres_blockkernelmatrix_avector(bkm_full, multi_sp, precon_sp, rhs, x0, gmres_tol, gmres_maxit, gmres_kk);
//		iter = solve_rpgmres_amatrix_avector(fullmat, multi_sp, precon_sp, rhs, x0, gmres_tol, gmres_maxit, gmres_kk);
//		iter = solve_gmres_amatrix_avector(kp, rhs, x0, gmres_tol, gmres_maxit, gmres_kk);
    //    print_avector_from_to(x0, x0->dim-3, x0->dim);
        clear_avector(sol);
        addeval_sparsematrix_avector(1.0, precon_sp, x0, sol);
		addeval_amatrix_avector(-1.0, fullmat, sol, rhs);
		error_fullgmres = norm2_avector(rhs)/error;
		t_setup = stop_stopwatch(sw);
		(void) printf("\t%.2f seconds\n", t_setup);
		(void) printf("\t%.1e relative MVM error", mvm_error);
		if(mvm_error > gmres_tol){
			printf(" ->BAD!\n");
		} else{
			printf(" ->GOOD!\n");
		}
		(void) printf("\t%u GMRES iterations with the full matrix\n", iter);
		(void) printf("\t%.12e relative GMRES error", error_fullgmres);
		if(error_fullgmres > gmres_tol){
			printf(" ->BAD!");
		} else{
			printf("->GOOD!");
		}
		(void) printf("\n--------------------------------------------------");

		addeval_amatrix_avector(1.0, precon_full, x0, sol2);
		del_blockkernelmatrix(bkm_full);
		del_avector(testvec);
		del_avector(res_a);
		del_avector(res_h2);
		del_amatrix(fullmat);
		del_amatrix(precon_full);
		del_amatrix(kp);
	}

	(void) printf("\nApproximating RMSE");
	start_stopwatch(sw);
	currentRelError = rmse_onthefly(sol, km, path_to_evalpts);
	t_setup = stop_stopwatch(sw);
	(void) printf("\n\t%.2f seconds\n", t_setup);
	(void) printf("\t%.2e rmse\n", currentRelError);
	if(error > gmres_tol){
		currentRelError = rmse_onthefly(sol2, km, path_to_evalpts);
		(void) printf("\t(%.2e rmse with full matrix)\n\n", currentRelError);
	}





	if(bkm->h2kmat==false)
		del_clusterbasis(cb);

	del_blockkernelmatrix(bkm);
	del_blockkernelmatrix(bkm_uncompressed);
	del_avector(rhs);
	del_avector(x0);
	del_avector(sol);
	del_avector(sol2);
	del_sparsematrix(precon_sp);
	del_kernelmatrix(km);
	del_clustergeometry(cg);
	del_cluster(root);
	del_block(broot);
	del_stopwatch(sw);
	freemem(idx);

	uninit_h2lib();
	return 0;
}



/*  FOR POTENTIAL LATER USE:

    cairo_surface_t *surface_h2mat = cairo_image_surface_create (CAIRO_FORMAT_ARGB32, 1920, 1080);
	cairo_t *cr_h2mat = cairo_create(surface_h2mat);
	draw_cairo_h2matrix(cr_h2mat, Gh2, 1, 0);
	cairo_surface_write_to_png (surface_h2mat, "./lshape/figures/h2m.png");
	cairo_destroy(cr_h2mat);
	cairo_surface_destroy(surface_h2mat);
*/
