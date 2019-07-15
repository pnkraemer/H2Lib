/*
NAME: 'kernelmat.h'

PURPOSE: Auxiliary functions to assemble kernel matrices

AUTHOR: kraemer(at)ins.uni-bonn.de
*/

#ifndef KERNELMAT_H_
#define KERNELMAT_H_
typedef struct {
	uint leafSize;
	real admPar;
	real truncAcc;
	real luAcc;
}HPar;



real expKernel(real coord){
	return exp(-0.5 * coord);
}


real maternKernel(real coord){
	return (1+ sqrt(3) * coord  ) * exp(- sqrt(3) * coord);
}

real maternKernel_15(real coord){
	return (1+ sqrt(3) * coord  ) * exp(- sqrt(3) * coord);
}

real tpsKernel(real coord){
    if(coord > 0.0){
        return coord * log(coord);
    }
    else{
        return 0;
    }
}


real maternKernel_25(real coord){
	return (1 + sqrt(5) * coord + 5.0 * coord * coord/3.0  ) * exp(- sqrt(5) * coord);
}


pamatrix new_kernelamatrix(clustergeometry *clGeom1, clustergeometry *clGeom2, real shift) {

	uint numPts1 = clGeom1->nidx;
	uint dim1 = clGeom1->dim;
	uint numPts2 = clGeom2->nidx;
	uint dim2 = clGeom2->dim;
	assert(dim1 == dim2);
	uint dim = dim1;

	pamatrix kernelMtrx = new_amatrix(numPts1, numPts2);
	real norm;
	for(uint i = 0; i < numPts1; i++){
		real* X = allocreal(dim);
		for(uint d = 0; d < dim; d++){
			X[d] = (clGeom1->x)[i][d];
		}
		for(uint j = 0; j < numPts2; j++){
			real* Y = allocreal(dim);
			for(uint d = 0; d < dim; d++){
				Y[d] = (clGeom2->x)[j][d];
			}
			norm = 0;
			for(uint d = 0; d < dim; d++){
				norm += (X[d] - Y[d])*(X[d] - Y[d]);
			}
			kernelMtrx->a[i + (kernelMtrx->ld) * j]  = maternKernel(sqrt(norm));
			freemem(Y);
			if(i==j){
				kernelMtrx->a[i + (kernelMtrx->ld) * j] += shift;
			}
		}
		freemem(X);
	}
	return kernelMtrx;
}








pamatrix new_kernelamatrix_tpssphere(clustergeometry *clGeom1, clustergeometry *clGeom2, real shift) {

	uint numPts1 = clGeom1->nidx;
	uint numPts1kurz = clGeom1->nidx - 4;
	uint dim1 = clGeom1->dim;
	uint numPts2 = clGeom2->nidx ;
	uint numPts2kurz = clGeom2->nidx - 4;
	uint dim2 = clGeom2->dim;
	assert(dim1 == dim2);
	uint dim = dim1;

	pamatrix kernelMtrx = new_amatrix(numPts1, numPts2);
	real norm;
	for(uint i = 0; i < numPts1kurz; i++){
		real* X = allocreal(dim);
		for(uint d = 0; d < dim; d++){
			X[d] = (clGeom1->x)[i][d];
		}
		for(uint j = 0; j < numPts2kurz; j++){
			real* Y = allocreal(dim);
			for(uint d = 0; d < dim; d++){
				Y[d] = (clGeom2->x)[j][d];
			}
			norm = 0;
			for(uint d = 0; d < dim; d++){
				norm += (X[d] * Y[d]);
			}
            setentry_amatrix(kernelMtrx, i, j, tpsKernel(1.0 - norm));
			if(i==j){
                addentry_amatrix(kernelMtrx, i, j, shift);
            }
			freemem(Y);
		}
        setentry_amatrix(kernelMtrx, i, numPts2kurz, 1.0);
        setentry_amatrix(kernelMtrx, i, numPts2kurz + 1, clGeom1->x[i][0]);
        setentry_amatrix(kernelMtrx, i, numPts2kurz + 2, clGeom1->x[i][1]);
        setentry_amatrix(kernelMtrx, i, numPts2kurz + 3, clGeom1->x[i][2]);

		freemem(X);
	}
	
	for(uint idx = 0; idx < numPts2kurz; idx++){
        setentry_amatrix(kernelMtrx, numPts1kurz, idx, 1.0);
        setentry_amatrix(kernelMtrx, numPts1kurz + 1, idx, clGeom2->x[idx][0]);
        setentry_amatrix(kernelMtrx, numPts1kurz + 2, idx, clGeom2->x[idx][1]);
        setentry_amatrix(kernelMtrx, numPts1kurz + 3, idx, clGeom2->x[idx][2]);
    }
	return kernelMtrx;
}










pamatrix new_kernelamatrix_tpssquare(clustergeometry *clGeom1, clustergeometry *clGeom2, real shift) {

	uint numPts1 = clGeom1->nidx;
	uint numPts1kurz = clGeom1->nidx - 3;
	uint dim1 = clGeom1->dim;
	uint numPts2 = clGeom2->nidx ;
	uint numPts2kurz = clGeom2->nidx - 3;
	uint dim2 = clGeom2->dim;
	assert(dim1 == dim2);
	uint dim = dim1;

	pamatrix kernelMtrx = new_amatrix(numPts1, numPts2);
	real norm;
	for(uint i = 0; i < numPts1kurz; i++){
		real* X = allocreal(dim);
		for(uint d = 0; d < dim; d++){
			X[d] = (clGeom1->x)[i][d];
		}
		for(uint j = 0; j < numPts2kurz; j++){
			real* Y = allocreal(dim);
			for(uint d = 0; d < dim; d++){
				Y[d] = (clGeom2->x)[j][d];
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
        setentry_amatrix(kernelMtrx, numPts1kurz + 1, idx, clGeom2->x[idx][0]);
        setentry_amatrix(kernelMtrx, numPts1kurz + 2, idx, clGeom2->x[idx][1]);
    }
	return kernelMtrx;
}






pamatrix new_kernelamatrix_maternsquare(clustergeometry *clGeom1, clustergeometry *clGeom2, real shift) {

	uint numPts1 = clGeom1->nidx;
	uint dim1 = clGeom1->dim;
	uint numPts2 = clGeom2->nidx ;
	uint dim2 = clGeom2->dim;
	assert(dim1 == dim2);
	uint dim = dim1;

	pamatrix kernelMtrx = new_amatrix(numPts1, numPts2);
	real norm;
	for(uint i = 0; i < numPts1; i++){
		real* X = allocreal(dim);
		for(uint d = 0; d < dim; d++){
			X[d] = (clGeom1->x)[i][d];
		}
		for(uint j = 0; j < numPts2; j++){
			real* Y = allocreal(dim);
			for(uint d = 0; d < dim; d++){
				Y[d] = (clGeom2->x)[j][d];
			}
						norm = 0;
			for(uint d = 0; d < dim; d++){
				norm += (X[d] - Y[d])*(X[d] - Y[d]);
			}
            setentry_amatrix(kernelMtrx, i, j, maternKernel(sqrt(norm)));
			if(i==j){
                addentry_amatrix(kernelMtrx, i, j, shift);
            }
			freemem(Y);
		}
		freemem(X);
	}
	return kernelMtrx;
}






pamatrix new_kernelamatrix_matern25(clustergeometry *clGeom1, clustergeometry *clGeom2, real shift) {

	uint numPts1 = clGeom1->nidx;
	uint dim1 = clGeom1->dim;
	uint numPts2 = clGeom2->nidx;
	uint dim2 = clGeom2->dim;
	assert(dim1 == dim2);
	uint dim = dim1;

	pamatrix kernelMtrx = new_amatrix(numPts1, numPts2);
	real norm;
	for(uint i = 0; i < numPts1; i++){
		real* X = allocreal(dim);
		for(uint d = 0; d < dim; d++){
			X[d] = (clGeom1->x)[i][d];
		}
		for(uint j = 0; j < numPts2; j++){
			real* Y = allocreal(dim);
			for(uint d = 0; d < dim; d++){
				Y[d] = (clGeom2->x)[j][d];
			}
			norm = 0;
			for(uint d = 0; d < dim; d++){
				norm += (X[d] - Y[d])*(X[d] - Y[d]);
			}
			kernelMtrx->a[i + (kernelMtrx->ld) * j] = maternKernel_25(sqrt(norm));
			freemem(Y);
			if(i==j){
				kernelMtrx->a[i + (kernelMtrx->ld) * j] += shift;
			}
		}
		freemem(X);
	}
	return kernelMtrx;
}


pamatrix new_kernelamatrix_exp(clustergeometry *clGeom1, clustergeometry *clGeom2, real shift) {

	uint numPts1 = clGeom1->nidx;
	uint dim1 = clGeom1->dim;
	uint numPts2 = clGeom2->nidx;
	uint dim2 = clGeom2->dim;
	assert(dim1 == dim2);
	uint dim = dim1;

	pamatrix kernelMtrx = new_amatrix(numPts1, numPts2);
	real norm;
	for(uint i = 0; i < numPts1; i++){
		real* X = allocreal(dim);
		for(uint d = 0; d < dim; d++){
			X[d] = (clGeom1->x)[i][d];
		}
		for(uint j = 0; j < numPts2; j++){
			real* Y = allocreal(dim);
			for(uint d = 0; d < dim; d++){
				Y[d] = (clGeom2->x)[j][d];
			}
			norm = 0;
			for(uint d = 0; d < dim; d++){
				norm += (X[d] - Y[d])*(X[d] - Y[d]);
			}
			kernelMtrx->a[i + (kernelMtrx->ld) * j] = expKernel(sqrt(norm));
			freemem(Y);
			if(i==j){
				kernelMtrx->a[i + (kernelMtrx->ld) * j] += shift;
			}
		}
		freemem(X);
	}
	return kernelMtrx;
}



phmatrix makeHMat(pamatrix kernelMtrx, pblock blockTree, HPar *hPar) {
	real truncAcc = hPar->truncAcc;

	ptruncmode truncMode = new_truncmode();
	ph2matrix h2KernelMtrx = compress_amatrix_h2matrix(kernelMtrx, blockTree, truncMode, truncAcc);
	phmatrix hKernelMtrx = convert_h2matrix_hmatrix(h2KernelMtrx);

	del_h2matrix(h2KernelMtrx);
	del_truncmode(truncMode);

	return hKernelMtrx;
}

ph2matrix makeH2Mat(pamatrix kernelMtrx, pblock blockTree, HPar *hPar) {
	real truncAcc = hPar->truncAcc;

	ptruncmode truncMode = new_abseucl_truncmode();
	ph2matrix h2KernelMtrx = compress_amatrix_h2matrix(kernelMtrx, blockTree, truncMode, truncAcc);

	del_truncmode(truncMode);

	return h2KernelMtrx;
}

ph2matrix makeH2Matfrob(pamatrix kernelMtrx, pblock blockTree, HPar *hPar) {
	real truncAcc = hPar->truncAcc;

	ptruncmode truncMode = new_relfrob_truncmode();
	ph2matrix h2KernelMtrx = compress_amatrix_h2matrix(kernelMtrx, blockTree, truncMode, truncAcc);

	del_truncmode(truncMode);

	return h2KernelMtrx;
}

#endif
