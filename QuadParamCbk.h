/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#ifndef QUAD_PARAM_CBK_H_
#define QUAD_PARAM_CBK_H_

#include <vector>

#include "MySparseMatrix.h"
#include "ParamCbk.h"

class QuadParamCbk : public ParamCbk {
	
	private:
		
		// dimensions
		int totalDecisionVars;
		int totalParameters;
		// terms only on the decisions
		int totNzSymQuadPart;
		int* nzRowsOneBasedSymQuadPart;
		int* nzColsOneBasedSymQuadPart;
		double* nzValsSymQuadPart;
		double* linearCostX;
		double freeTerm;
		// linear parametric perturation
		double* linearCostP;
		// sparsity patterns grad
		vector<int> *nzPosOneBasedGradX;
		vector<int> *nzPosOneBasedGradP;
		// sparsity pattern hess x
		vector<int> *nzRowsOneBasedHessX;
		vector<int> *nzColsOneBasedHessX;
		// sparsity pattern partial p grad x
		vector<int> *nzRowsOneBasedPartialPGradX;
	
	public:
		
		QuadParamCbk(
			int totalDecisionVars_,
			int totalParameters_,
			int totNzSymQuadPart_,
			int* nzRowsOneBasedSymQuadPart_,
			int* nzColsOneBasedSymQuadPart_,
			double* nzValsSymQuadPart_,
			double* linearCostX_,
			double freeTerm_,
			double* linearCostP_);

		~QuadParamCbk();

		double getValueFunc(double* x, double* p);
		
		void addScaledGradXFunc(double* x, double* p, int* dict, double* values, double scaling);
		
		void addScaledGradPFunc(double* x, double* p, int* dict, double* values, double scaling);

		void addScaledHessXFunc(double* x, double* p, int* dict, double* values, double scaling);
		
		void addScaledPartialPGradXFunc(int jOneBased, double* x, double* p, int* dict, double* values, double scaling);
		
		double getPartialPFunc(int jOneBased, double* x, double*p);

		vector<int>* getNzPosOneBasedGradX();
		
		int getTotalNzGradX();
		
		vector<int>* getNzPosOneBasedGradP();
		
		int getTotalNzGradP();
		
		vector<int>* getNzRowsOneBasedHessX();
		
		vector<int>* getNzColsOneBasedHessX();
		
		int getTotalNzHessX();

		vector<int>* getNzRowsOneBasedPartialPGradX();
		
		int getTotalNzRowsOneBasedPartialPGradX();
};

#endif
