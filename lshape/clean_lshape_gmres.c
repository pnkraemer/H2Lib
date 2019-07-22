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
#include "blockkernelmatrix.h"
#include "krylovsolvers.h"
#include "tri2d.h"

static void
multi_sp(void *pdata, pavector r)
{
    psparsematrix  A = (psparsematrix) pdata;
    pavector rcopy = new_zero_avector(r->dim);
    assert(A->cols == rcopy->dim);
    addeval_sparsematrix_avector(1.0, A, r, rcopy);
    copy_avector(rcopy, r);
    del_avector(rcopy);
}


field
interpolant(field x, field y){
	return x + y;
}

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

HEADER_PREFIX psparsematrix
loadfromtxt_precon(uint N, uint n, char *filepath)
{
    real val;
    uint colidx, rowidx;
    uint i;
    char strVals[150], strRowIdx[150], strColIdx[150];
    char filepath_val[250], filepath_row[250], filepath_col[250];
    pavector preconval, preconrow, preconcol;
    psparsepattern sp;
    psparsematrix spm;

    preconval = new_avector(n*N);
    preconrow = new_avector(n*N);
    preconcol = new_avector(n*N);

    strcpy(filepath_val, filepath);
    strcpy(filepath_row, filepath);
    strcpy(filepath_col, filepath);
    sprintf(strVals, "precon/precon_val_N%d_n%d.txt", N, n);
    sprintf(strRowIdx, "precon/precon_row_N%d_n%d.txt", N, n);
    sprintf(strColIdx, "precon/precon_col_N%d_n%d.txt", N, n);
    strcat(filepath_val, strVals);
    strcat(filepath_row, strRowIdx);
    strcat(filepath_col, strColIdx);

    FILE *myFileVals, *myFileRowIdx, *myFileColIdx;
    myFileVals = fopen(filepath_val, "r");
    myFileRowIdx = fopen(filepath_row, "r");
    myFileColIdx = fopen(filepath_col, "r");
    if(myFileRowIdx == NULL || myFileColIdx == NULL || myFileVals == NULL){
        printf("Error reading file\n");
        exit(0);
    }

    sp = new_sparsepattern(N + 3, N + 3);
    for(i = 0; i < n*N; i++){
            fscanf(myFileVals, "%lf ", &val);
            fscanf(myFileRowIdx, "%u ", &rowidx);
            fscanf(myFileColIdx, "%u ", &colidx);
            setentry_avector(preconval, i, val);
            setentry_avector(preconrow, i, rowidx);
            setentry_avector(preconcol, i, colidx);
            addnz_sparsepattern(sp, rowidx, colidx);
    }
    for(uint i = 0; i < 3; i++){
        addnz_sparsepattern(sp, N + i, N + i);
    }
    fclose(myFileVals);
    fclose(myFileColIdx);
    fclose(myFileRowIdx);

    spm = new_zero_sparsematrix(sp);
    for(i = 0; i < n*N; i++){
        val = getentry_avector(preconval, i);
        rowidx = getentry_avector(preconrow, i);
        colidx = getentry_avector(preconcol, i);
        setentry_sparsematrix(spm, rowidx, colidx, val);
    }
    for(uint i = 0; i < 3; i++){
        setentry_sparsematrix(spm, N+i, N+i, 1.0);
    }

    del_sparsepattern(sp);
    del_avector(preconval);
    del_avector(preconrow);
    del_avector(preconcol);
    return spm;
}


INLINE_PREFIX real /* Not happy with the name of the function */
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
	uint                n;
	uint                iter;
	real                error;
	bool                print_yes;
	uint                evalpoints;
	real                currentRelError;
	char                filepath[250];
	real                gmres_tol;
	uint                gmres_kk, gmres_maxit;
	real                error_fullgmres;
    real                h2compression_acc;
	/* Parameters */
	sw = new_stopwatch();
	strcpy(filepath, argv[1]);    /* path to mesh and precon*/
	points = atoi(argv[2]);       /* number of points*/
	n = atoi(argv[3]);            /* number of neighbours */
	m = atoi(argv[4]);            /* interpolation order */
	assert(points>0 && n > 0 && m > 0);
	lsz = 2*m*m;                  /* leafsize prop. to interpolation order */
	eta = 2.0;                    /* generic admissibility condition */
//  assert(points<24000);         /* 24000 is all one can do with 8GB of RAM*/
	dim = 2;                      /* this script is 2D only*/
	evalpoints = 50000;           /* number of points for RMSE estimate*/
	gmres_tol = 1e-10;
	gmres_maxit = 5000;
	gmres_kk = 20;
	cpos = 1;
    h2compression_acc = 1e-2 * gmres_tol;

	(void) printf("\nCreating kernelmatrix object for %u points (%u neighbours), interpolation order %u\n", points, n, m);
	start_stopwatch(sw);
	km = new_kernelmatrix(dim, points, m, cpos);
	km->kernel = tps_kernel_2d; /* Choose 2d-kernel*/
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
	precon_sp = loadfromtxt_precon(points, n, filepath);
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
	sol = new_zero_avector(points + dim + 1);
	addeval_sparsematrix_avector(1.0, precon_sp, x0, sol);
	addeval_blockkernelmatrix_avector(-1.0, bkm, sol, rhs);
	error = norm2_avector(rhs)/error;
	t_setup = stop_stopwatch(sw);
	(void) printf("\t%.2f seconds\n", t_setup);
	(void) printf("\t%u GMRES iterations with the H2-matrix\n", iter);
	(void) printf("\t%.1e relative error", error);




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
//		iter = solve_rpgmres_amatrix_avector(fullmat, multi_sp, precon_sp, rhs, x0, gmres_tol, gmres_maxit, gmres_kk);
		iter = solve_gmres_amatrix_avector(kp, rhs, x0, gmres_tol, gmres_maxit, gmres_kk);
	// print_avector(x0);
		addeval_amatrix_avector(-1.0, kp, x0, rhs);
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
		(void) printf("\t%.1e relative GMRES error", error_fullgmres);
		if(error_fullgmres > gmres_tol){
			printf(" ->BAD!");
		} else{
			printf("->GOOD!");
		}
		(void) printf("\n--------------------------------------------------");

		addeval_amatrix_avector(1.0, precon_full, x0, sol2);
		del_avector(testvec);
		del_avector(res_a);
		del_avector(res_h2);
		del_amatrix(fullmat);
		del_amatrix(precon_full);
		del_amatrix(kp);
		del_blockkernelmatrix(bkm_full);
	}


	(void) printf("\nApproximating RMSE");
	start_stopwatch(sw);
	t_setup = stop_stopwatch(sw);
	(void) printf("\n\t%.2f seconds\n", t_setup);
	currentRelError = rmse_onthefly(sol, km, evalpoints);
	(void) printf("\t%.2e rmse\n", currentRelError);
	if(error > gmres_tol){
		currentRelError = rmse_onthefly(sol2, km, evalpoints);
		(void) printf("\t(%.2e rmse with full matrix)\n\n", currentRelError);
	}





	if(bkm->h2kmat==false)
		del_clusterbasis(cb);

	del_blockkernelmatrix(bkm_uncompressed);
	del_blockkernelmatrix(bkm);
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
