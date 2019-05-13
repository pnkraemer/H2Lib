
#include "helpfcts.h"

typedef struct  {
	uint maxIt;
	real acc;
	real shift;
}RileyPar;


typedef struct {
	uint numit;
	real lastRelError;
	pavector relError;
	pavector solution;
}RileyOutput;





RileyOutput* rileyAlgo(RileyPar *rileyPar, ph2matrix h2KernelMtrx, ph2matrix L_h2KeSh, ph2matrix R_h2KeSh, pavector rhsVec, pavector startVec) {

	startVec = clone_avector(rhsVec);
	lrsolve_h2matrix_avector(L_h2KeSh, R_h2KeSh, startVec);

	uint numPts = startVec->dim;
	RileyOutput *rilOut = malloc(sizeof(RileyOutput));
	rilOut->solution = new_zero_avector(numPts);
	rilOut->relError = new_zero_avector(rileyPar->maxIt);

	pavector rhsVecCopyForErrorAn = clone_avector(rhsVec);
	real normOfRhs = norm2_avector(startVec);

	/* SET PARAMETERS */
	pavector currIt = clone_avector(rhsVec);
	int counter = 0;
	real currentRelError = 100.0;
	pavector relError = new_zero_avector(rileyPar->maxIt);
	/* ITERATION */
	while(currentRelError >= rileyPar->acc && counter < rileyPar->maxIt){

		
		/* SOLVE SYSTEM*/
		pavector solVec = clone_avector(currIt);
		lrsolve_h2matrix_avector(L_h2KeSh, R_h2KeSh, solVec);
		if(isnan(norm2_avector(solVec))==1){
			printf("NAN!!! %i\n", counter);
			printf("currentRelError: %.1e\n", currentRelError);
			printf("lastRelError: %.1e\n", relError->v[counter-2]);
		}

		/* RILEY ITERATION */	
		pavector startVecCopy = clone_avector(startVec);
		add_avector(rileyPar->shift, solVec, startVecCopy);
		copy_avector(startVecCopy, currIt);

		/* COMPUTE ERROR */
		pavector matVecProd = new_zero_avector(numPts);
		addeval_h2matrix_avector(1.0, h2KernelMtrx, currIt, matVecProd);
		add_avector(-1.0, rhsVecCopyForErrorAn, matVecProd);
		currentRelError = norm2_avector(matVecProd) / normOfRhs;
		if(isinf(currentRelError)==1){
			printf("INFINITY!!! (at counter=%i)\n", counter);
			printf("lastRelError: %.1e\n", relError->v[counter-1]);
		}
		relError->v[counter] = currentRelError;
		counter += 1;
		del_avector(startVecCopy);
		del_avector(matVecProd);
		del_avector(solVec);
	}

	rilOut->numit = counter;
	rilOut->lastRelError = currentRelError;
	copy_avector(currIt, rilOut->solution);
	copy_avector(relError, rilOut->relError);


	del_avector(startVec);
	del_avector(relError);
	del_avector(rhsVecCopyForErrorAn);
	del_avector(currIt);

	return rilOut;
}







RileyOutput* rileyAlgo2(RileyPar *rileyPar, ph2matrix h2KernelMtrx, ph2matrix L_h2KeSh, ph2matrix R_h2KeSh, pavector rhsVec, pavector startVec) {

	startVec = clone_avector(rhsVec);
	lrsolve_h2matrix_avector(L_h2KeSh, R_h2KeSh, startVec);

	uint numPts = startVec->dim;
	RileyOutput *rilOut = malloc(sizeof(RileyOutput));
	rilOut->solution = new_zero_avector(numPts);
	rilOut->relError = new_zero_avector(rileyPar->maxIt);

	pavector rhsVecCopyForErrorAn = clone_avector(rhsVec);
	real normOfRhs = norm2_avector(startVec);

	/* SET PARAMETERS */
	pavector currIt = clone_avector(rhsVec);
	int counter = 0;
	real currentRelError = 100.0;
	pavector relError = new_zero_avector(rileyPar->maxIt);
	/* ITERATION */
	while(counter < rileyPar->maxIt){

		
		/* SOLVE SYSTEM*/
		pavector solVec = clone_avector(currIt);
		lrsolve_h2matrix_avector(L_h2KeSh, R_h2KeSh, solVec);
		if(isnan(norm2_avector(solVec))==1){
			printf("NAN!!! %i\n", counter);
			printf("currentRelError: %.1e\n", currentRelError);
			printf("lastRelError: %.1e\n", relError->v[counter-2]);
		}

		/* RILEY ITERATION */	
		pavector startVecCopy = clone_avector(startVec);
		add_avector(rileyPar->shift, solVec, startVecCopy);
		copy_avector(startVecCopy, currIt);

		/* COMPUTE ERROR */
		pavector matVecProd = new_zero_avector(numPts);
		addeval_h2matrix_avector(1.0, h2KernelMtrx, currIt, matVecProd);
		add_avector(-1.0, rhsVecCopyForErrorAn, matVecProd);
		currentRelError = norm2_avector(matVecProd) / normOfRhs;
		if(isinf(currentRelError)==1){
			printf("INFINITY!!! (at counter=%i)\n", counter);
			printf("lastRelError: %.1e\n", relError->v[counter-1]);
		}
		relError->v[counter] = currentRelError;
		counter += 1;
		del_avector(startVecCopy);
		del_avector(matVecProd);
		del_avector(solVec);
	}

	rilOut->numit = counter;
	rilOut->lastRelError = currentRelError;
	copy_avector(currIt, rilOut->solution);
	copy_avector(relError, rilOut->relError);


	del_avector(startVec);
	del_avector(relError);
	del_avector(rhsVecCopyForErrorAn);
	del_avector(currIt);

	return rilOut;
}






RileyOutput* rileyAlgo3(RileyPar *rileyPar, pcamatrix h2KernelMtrx, ph2matrix L_h2KeSh, ph2matrix R_h2KeSh, pavector rhsVec, pavector startVec) {

	startVec = clone_avector(rhsVec);
	lrsolve_h2matrix_avector(L_h2KeSh, R_h2KeSh, startVec);

	uint numPts = startVec->dim;
	RileyOutput *rilOut = malloc(sizeof(RileyOutput));
	rilOut->solution = new_zero_avector(numPts);
	rilOut->relError = new_zero_avector(rileyPar->maxIt);

	pavector rhsVecCopyForErrorAn = clone_avector(rhsVec);
	real normOfRhs = norm2_avector(startVec);

	/* SET PARAMETERS */
	pavector currIt = clone_avector(rhsVec);
	int counter = 0;
	real currentRelError = 100.0;
	pavector relError = new_zero_avector(rileyPar->maxIt);
	/* ITERATION */
	while(currentRelError >= rileyPar->acc && counter < rileyPar->maxIt){

		
		/* SOLVE SYSTEM*/
		pavector solVec = clone_avector(currIt);
		lrsolve_h2matrix_avector(L_h2KeSh, R_h2KeSh, solVec);
		if(isnan(norm2_avector(solVec))==1){
			printf("NAN!!! %i\n", counter);
			printf("currentRelError: %.1e\n", currentRelError);
			printf("lastRelError: %.1e\n", relError->v[counter-2]);
		}

		/* RILEY ITERATION */	
		pavector startVecCopy = clone_avector(startVec);
		add_avector(rileyPar->shift, solVec, startVecCopy);
		copy_avector(startVecCopy, currIt);

		/* COMPUTE ERROR */
		pavector matVecProd = new_zero_avector(numPts);
		addeval_amatrix_avector(1.0, h2KernelMtrx, currIt, matVecProd);
		add_avector(-1.0, rhsVecCopyForErrorAn, matVecProd);
		currentRelError = norm2_avector(matVecProd) / normOfRhs;
		if(isinf(currentRelError)==1){
			printf("INFINITY!!! (at counter=%i)\n", counter);
			printf("lastRelError: %.1e\n", relError->v[counter-1]);
		}
		relError->v[counter] = currentRelError;
		counter += 1;
		del_avector(startVecCopy);
		del_avector(matVecProd);
		del_avector(solVec);
	}

	rilOut->numit = counter;
	rilOut->lastRelError = currentRelError;
	copy_avector(currIt, rilOut->solution);
	copy_avector(relError, rilOut->relError);


	del_avector(startVec);
	del_avector(relError);
	del_avector(rhsVecCopyForErrorAn);
	del_avector(currIt);

	return rilOut;
}



RileyOutput* rileyAlgo4(RileyPar *rileyPar, pcavector trueSol, ph2matrix L_h2KeSh, ph2matrix R_h2KeSh, pavector rhsVec, pavector startVec) {

	startVec = clone_avector(rhsVec);
	lrsolve_h2matrix_avector(L_h2KeSh, R_h2KeSh, startVec);

	uint numPts = startVec->dim;
	RileyOutput *rilOut = malloc(sizeof(RileyOutput));
	rilOut->solution = new_zero_avector(numPts);
	rilOut->relError = new_zero_avector(rileyPar->maxIt);

	pavector rhsVecCopyForErrorAn = clone_avector(rhsVec);
	real normOfRhs = norm2_avector(trueSol);

	/* SET PARAMETERS */
	pavector currIt = clone_avector(rhsVec);
	int counter = 0;
	real currentRelError = 100.0;
	pavector relError = new_zero_avector(rileyPar->maxIt);
	/* ITERATION */
	while(currentRelError >= rileyPar->acc && counter < rileyPar->maxIt){

		
		/* SOLVE SYSTEM*/
		pavector solVec = clone_avector(currIt);
		lrsolve_h2matrix_avector(L_h2KeSh, R_h2KeSh, solVec);
		if(isnan(norm2_avector(solVec))==1){
			printf("NAN!!! %i\n", counter);
			printf("currentRelError: %.1e\n", currentRelError);
			printf("lastRelError: %.1e\n", relError->v[counter-2]);
		}

		/* RILEY ITERATION */	
		pavector startVecCopy = clone_avector(startVec);
		add_avector(rileyPar->shift, solVec, startVecCopy);
		copy_avector(startVecCopy, currIt);

		/* COMPUTE ERROR */
		pavector matVecProd = clone_avector(currIt);
	//	addeval_amatrix_avector(1.0, h2KernelMtrx, currIt, matVecProd);
		add_avector(-1.0, trueSol, matVecProd);
		currentRelError = norm2_avector(matVecProd) / normOfRhs;
		if(isinf(currentRelError)==1){
			printf("INFINITY!!! (at counter=%i)\n", counter);
			printf("lastRelError: %.1e\n", relError->v[counter-1]);
		}
		relError->v[counter] = currentRelError;
		counter += 1;
		del_avector(startVecCopy);
		del_avector(matVecProd);
		del_avector(solVec);
	}

	rilOut->numit = counter;
	rilOut->lastRelError = currentRelError;
	copy_avector(currIt, rilOut->solution);
	copy_avector(relError, rilOut->relError);


	del_avector(startVec);
	del_avector(relError);
	del_avector(rhsVecCopyForErrorAn);
	del_avector(currIt);

	return rilOut;
}








RileyOutput* rileyAlgo5(RileyPar *rileyPar, pcavector trueSol, ph2matrix L_h2KeSh, ph2matrix R_h2KeSh, pavector rhsVec, pavector startVec) {

	startVec = clone_avector(rhsVec);
	lrsolve_h2matrix_avector(L_h2KeSh, R_h2KeSh, startVec);

	uint numPts = startVec->dim;
	RileyOutput *rilOut = malloc(sizeof(RileyOutput));
	rilOut->solution = new_zero_avector(numPts);
	rilOut->relError = new_zero_avector(rileyPar->maxIt);

	pavector rhsVecCopyForErrorAn = clone_avector(rhsVec);
	real normOfRhs = norm2_avector(trueSol);

	/* SET PARAMETERS */
	pavector currIt = clone_avector(rhsVec);
	int counter = 0;
	real currentRelError = 100.0;
	pavector relError = new_zero_avector(rileyPar->maxIt);
	/* ITERATION */
	while(counter < rileyPar->maxIt){

		
		/* SOLVE SYSTEM*/
		pavector solVec = clone_avector(currIt);
		lrsolve_h2matrix_avector(L_h2KeSh, R_h2KeSh, solVec);
		if(isnan(norm2_avector(solVec))==1){
			printf("NAN!!! %i\n", counter);
			printf("currentRelError: %.1e\n", currentRelError);
			printf("lastRelError: %.1e\n", relError->v[counter-2]);
		}

		/* RILEY ITERATION */	
		pavector startVecCopy = clone_avector(startVec);
		add_avector(rileyPar->shift, solVec, startVecCopy);
		copy_avector(startVecCopy, currIt);

		/* COMPUTE ERROR */
		pavector matVecProd = clone_avector(currIt);
	//	addeval_amatrix_avector(1.0, h2KernelMtrx, currIt, matVecProd);
		add_avector(-1.0, trueSol, matVecProd);
		currentRelError = norm2_avector(matVecProd) / normOfRhs;
		if(isinf(currentRelError)==1){
			printf("INFINITY!!! (at counter=%i)\n", counter);
			printf("lastRelError: %.1e\n", relError->v[counter-1]);
		}
		relError->v[counter] = currentRelError;
		counter += 1;
		del_avector(startVecCopy);
		del_avector(matVecProd);
		del_avector(solVec);
	}

	rilOut->numit = counter;
	rilOut->lastRelError = currentRelError;
	copy_avector(currIt, rilOut->solution);
	copy_avector(relError, rilOut->relError);


	del_avector(startVec);
	del_avector(relError);
	del_avector(rhsVecCopyForErrorAn);
	del_avector(currIt);

	return rilOut;
}







