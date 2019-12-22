/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#ifndef LINEAR_CBK_H_
#define LINEAR_CBK_H_

#include <map>

#include "Cbk.h"

using namespace std;

class LinearCbk : public Cbk {
		
	private:
		
		int totNz;
		int* nzRowsOneBasedLinearCost;
		double* nzValsLinearCost;
		double freeTerm;
		double scalingLinearCost;
		map<int, int> *nzRowToValPos;
		
	public:
		
		LinearCbk(int totalVariables_, int totNz_,
			int* nzRowsOneBasedLinearCost_,
			double* nzValsLinearCost_,
			double freeTerm_,
			double scalingLinearCost_);
		
		~LinearCbk();
		
		double getValueFunc(double* x);
		
		void getFullGradFunc(double* x, double* grad);
		
		double getPartialFunc(int jOneBased, double* x);
};

#endif
