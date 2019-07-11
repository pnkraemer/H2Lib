/*
* A few add-ons for the kernelmatrix.h module,
* mostly regarding conditionally positive definite kernels
*
* Author: Nicholas Krämer
*/

#include "addon_kernelmatrix.h"
#include "kernelmatrix.h"
#include "amatrix.h"

void 
assemble_pblock(pckernelmatrix km, pamatrix mtrx){

	uint 		points, dim;
	uint 		i, j;

	points = km->points;
	dim = km->dim;
	assert(mtrx->rows == points && mtrx->cols == (1+dim));

	for(i=0;i<points;i++){
		setentry_amatrix(mtrx, i, 0, 1.0);
		for(j=0; j<dim; j++){
			setentry_amatrix(mtrx, i, 1 + j, km->x[i][j]);
		}
	}
}




void
assemble_pblock_trans(pckernelmatrix km, pamatrix mtrx){

	uint 		points, dim;
	uint 		i, j;

	points = km->points;
	dim = km->dim;
	assert(mtrx->cols == points && mtrx->rows == (1+dim));

	for(i=0;i<points;i++){
		setentry_amatrix(mtrx, 0, i, 1.0);
		for(j=0; j<dim; j++){
			setentry_amatrix(mtrx, 1 + j, i, km->x[i][j]);
		}
	}
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
    pb = new_amatrix(points, 1 + dim);
    assemble_pblock(km, pb);
    copy_amatrix(0, pb, topright);
    copy_amatrix(1, pb, bottomleft);
    clear_amatrix(bottomright);


    del_amatrix(topleft);
    del_amatrix(bottomleft);
    del_amatrix(topright);
    del_amatrix(bottomright);
    del_amatrix(pb);

}

