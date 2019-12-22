/*
 * author: pedro.borges.melo@gmail.com
 * date: October/2019
 */

#include <stdio.h>
#include <stdlib.h>

#include "LinearSystemSolver.h"

LinearSystemSolver::LinearSystemSolver(int n_, int* ia_, int* ja_, double* a_, bool symmetricReal) {
	// set data
	this->n = n_;
	this->ia = ia_;
	this->ja = ja_;
	this->a = a_;
	// symmetric indefinite matrix
	if (symmetricReal == true) {
		this->mtype = LinearSystemSolver::TYPE_SYS_REAL_SYMMETRIC_INDEFINITE;
	} else {
		this->mtype = LinearSystemSolver::TYPE_SYS_REAL_NON_SYMMETRIC;
	}
	// Number of right hand sides
	this->nrhs = 1;
	// Setup Pardiso control parameters
	MKL_INT i;
	for ( i = 0; i < 64; i++ ) {
		this->iparm[i] = 0;
	}
	// No solver default
	this->iparm[0] = 1;
	// Fill-in reordering from METIS 
	this->iparm[1] = 2;
	// No iterative-direct algorithm
	this->iparm[3] = 0;
	// No user fill-in reducing permutation 
	this->iparm[4] = 0;
	// Write solution into x 
	this->iparm[5] = 0;
	// Not in use 
	this->iparm[6] = 0;
	// Max numbers of iterative refinement steps 
	this->iparm[7] = 2;
	// Not in use
	this->iparm[8] = 0;
	// Perturb the pivot elements with 1E-13 
	this->iparm[9] = 13;
	// Use nonsymmetric permutation and scaling MPS 
	this->iparm[10] = 1;
	// Not in use 
	this->iparm[11] = 0;
	// Maximum weighted matching algorithm is switched-off (default for symmetric)
	// Try iparm[12] = 1 in case of inappropriate accuracy
	this->iparm[12] = 0;
	// Output: Number of perturbed pivots 
	this->iparm[13] = 0;
	// Not in use 
	this->iparm[14] = 0;
	// Not in use 
	this->iparm[15] = 0;
	// Not in use 
	this->iparm[16] = 0;
	// Output: Number of nonzeros in the factor LU 
	this->iparm[17] = -1;
	// Output: Mflops for LU factorization 
	this->iparm[18] = -1;
	// Output: Numbers of CG Iterations 
	this->iparm[19] = 0;
	// Maximum number of numerical factorizations. 
	this->maxfct = 1;
	// Which factorization to use. 
	this->mnum = 1;
	// Print statistical information in file 
	this->msglvl = 0;
	// Initialize error flag 
	this->error = 0;
	// Initialize the internal solver memory pointer. This is only
	// necessary for the FIRST call of the PARDISO solver.
	for ( i = 0; i < 64; i++ ) {
		this->pt[i] = 0;
	}
	// Reordering and Symbolic Factorization. This step also allocates
	// all memory that is necessary for the factorization.
	this->phase = 11;
	PARDISO (this->pt, &(this->maxfct), &(this->mnum), 
		&(this->mtype), &(this->phase),
		&n, a, ia, ja, &(this->idum), 
		&(this->nrhs), this->iparm, 
		&(this->msglvl), &(this->ddum), 
		&(this->ddum), &(this->error));
	if ( this->error != 0 ) {
		printf ("\nERROR during symbolic factorization: %d", this->error);
		exit(-1);
	}
	// Numerical factorization
	this->phase = 22;
	PARDISO (this->pt, &(this->maxfct), &(this->mnum), 
		&(this->mtype), &(this->phase),
		&n, a, ia, ja, &(this->idum), 
		&(this->nrhs), (this->iparm), &(this->msglvl),
		&(this->ddum), &(this->ddum), &(this->error));
	// Check for error
	if ( this->error != 0 ) {
		printf ("\nERROR during numerical factorization: %d", this->error);
		exit(-1);
	}
}

LinearSystemSolver::~LinearSystemSolver() {
	// Release internal memory
	this->phase = -1;
	PARDISO(this->pt, &(this->maxfct), &(this->mnum), 
		&(this->mtype), &(this->phase),
		&(this->n), &(this->ddum), 
		this->ia, this->ja, &(this->idum), &(this->nrhs),
		this->iparm, &(this->msglvl), &(this->ddum), 
		&(this->ddum), &(this->error));
}

int LinearSystemSolver::solveSystemWithPardiso(double* b, double* x) {
	// Back substitution and iterative refinement
	this->phase = 33;
	// Max numbers of iterative refinement steps
	this->iparm[7] = 2;
	PARDISO (this->pt, &(this->maxfct), 
		&(this->mnum), &(this->mtype), 
		&(this->phase), &(this->n),
		this->a, this->ia, this->ja, &(this->idum), &(this->nrhs), 
		this->iparm, &(this->msglvl), 
		b, x, &(this->error));
	// Check for error
	if ( this->error != 0 ) {
		printf ("\nERROR during solution: %d", this->error);
		exit(-1);
	}
	// return
	return 0;
}
