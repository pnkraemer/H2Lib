/*
* A few add-ons for the kernelmatrix.h module,
* mostly regarding conditionally positive definite kernels
*
* Author: Nicholas Krämer
*/

#include "addon_kernelmatrix.h"
#include "kernelmatrix.h"
#include "amatrix.h"

pamatrix 
assemble_pblock(pckernelmatrix km){

	pamatrix 	mtrx;
	uint 		points, dim;
	uint 		i, j;

	points = km->points;
	dim = km->dim;
	mtrx = new_amatrix(points, 1 + dim);

	for(i=0;i<points;i++){
		setentry_amatrix(mtrx, i, 0, 1.0);
		for(j=0; j<dim; j++){
			setentry_amatrix(mtrx, i, 1 + j, km->x[i][j]);
		}
	}
	return mtrx;
}




pamatrix
assemble_pblock_trans(pckernelmatrix km){

	pamatrix 	mtrx;
	uint 		points, dim;
	uint 		i, j;

	points = km->points;
	dim = km->dim;
	mtrx = new_amatrix(1 + dim, points);

	for(i=0;i<points;i++){
		setentry_amatrix(mtrx, 0, i, 1.0);
		for(j=0; j<dim; j++){
			setentry_amatrix(mtrx, 1 + j, i, km->x[i][j]);
		}
	}
	return mtrx;
}

void 
addeval_cond_kernelh2matrix(field alpha, ph2matrix h2km, pamatrix pb, pavector src, pavector trg){

    pavector 	top_subvec_src, bottom_subvec_src, top_subvec_trg, bottom_subvec_trg;
    uint 		pts, pbsize;

    pts = pb->rows;
    pbsize = pb->cols;
    top_subvec_src = new_sub_avector(src, pts, 0);
    bottom_subvec_src = new_sub_avector(src, pbsize, pts);
    top_subvec_trg = new_sub_avector(trg, pts, 0);
    bottom_subvec_trg = new_sub_avector(trg, pbsize, pts);

    addeval_h2matrix_avector(alpha, h2km, top_subvec_src, top_subvec_trg);
    addeval_amatrix_avector(alpha, pb, bottom_subvec_src, top_subvec_trg);
    addevaltrans_amatrix_avector(alpha, pb, top_subvec_src, bottom_subvec_trg);

    del_avector(top_subvec_src);
    del_avector(bottom_subvec_src);
    del_avector(top_subvec_trg);
    del_avector(bottom_subvec_trg);
}


void
assemble_big_kernelmatrix(pkernelmatrix km, pamatrix mat){

	uint points, dim;
	pamatrix topleft, bottomleft, topright, bottomright;
	pamatrix pb;


	points = km->points;
	dim = km->dim;

	topleft = new_sub_amatrix(mat, points, 0, points, 0);
	topright = new_sub_amatrix(mat, points,  0, 1 + dim, points);
	bottomleft = new_sub_amatrix(mat, 1 + dim, points, points, 0);
	bottomright = new_sub_amatrix(mat, 1 + dim, points, 1 + dim, points);


    fillN_kernelmatrix(0, 0, km, topleft);
    pb = assemble_pblock(km);
    copy_amatrix(0, pb, topright);
    copy_amatrix(1, pb, bottomleft);
    clear_amatrix(bottomright);


    del_amatrix(topleft);
    del_amatrix(bottomleft);
    del_amatrix(topright);
    del_amatrix(bottomright);
    del_amatrix(pb);

}

