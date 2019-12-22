/*
 * author: pedro.borges.melo@gmail.com
 * date: October/2019
 */

#ifndef LINEAR_SYSTEM_H_
#define LINEAR_SYSTEM_H_

#include <vector>
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

using namespace std;

class LinearSystem {

	private:

		int totalEqLinearSystem;
		int totalNzEqLinearSystem;
		MKL_INT* iaPardiso;
		MKL_INT* jaPardiso;
		double* valsPardiso;
		double* rhsPardiso;
		double* solPardiso;

	public:
	
		// NOTE: vector< pair<int, int> >* is assumed to be sorted lexicographically
		LinearSystem(int totalEqLinearSystem_, int totalNzEqLinearSystem_, vector< pair<int, int> > *nzPositionsLinSys );

		void createSparsityStructureWithPairs( vector< pair<int, int> > *nzPositionsLinSys_ );

		~LinearSystem();

		int getTotalEqLinearSystem();
	
		int getTotalNzEqLinearSystem();
	
		MKL_INT* getIaPardiso();
		
		MKL_INT* getJaPardiso();
		
		double* getValsPardiso();
		
		double* getRhsPardiso();

		double* getSolPardiso();
};

#endif
