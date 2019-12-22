/*
 * author: pedro.borges.melo@gmail.com
 * date: October/2019
 */

#ifndef LINEAR_SYSTEM_SOLVER_H_
#define LINEAR_SYSTEM_SOLVER_H_

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

class LinearSystemSolver {

	private:

		// constants
		const int TYPE_SYS_REAL_SYMMETRIC_INDEFINITE = -2;
		const int TYPE_SYS_REAL_NON_SYMMETRIC = 11;
		// data
		int n;
		int* ia;
		int* ja;
		double* a;
		// pardiso data
		void *pt[64];
		MKL_INT maxfct;
		MKL_INT mnum;
		MKL_INT mtype;
		MKL_INT phase;
		MKL_INT idum;
		MKL_INT nrhs;
		MKL_INT iparm[64];
		MKL_INT msglvl;
		MKL_INT error;
		double ddum;

	public:
	
		LinearSystemSolver(int n_, int* ia_, int* ja_, double* a_, bool symmetricReal);

		~LinearSystemSolver();

		int solveSystemWithPardiso(double* b, double* x);
};

#endif
