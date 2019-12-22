/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#ifndef MTX_PARAM_CBK_H_
#define MTX_PARAM_CBK_H_

#include <map>
#include <vector>

#include "MySparseMatrix.h"
#include "Cbk.h"

using namespace std;

class MtxParamCbk {
	
	private:
		
		// non parametric part
		MySparseMatrix* fixedPart;
		// parametric part
		vector<int> *nzRowsOneBasedParamPart;
		vector<int> *nzColsOneBasedParamPart;
		vector<Cbk*> *nzValsParamPart;
		// final sparsity pattern
		vector<int> *finalNzRows;
		vector<int> *finalNzCols;
		// inner filling instructions 
		int* fillInstructFixedPart;
		int* fillInstructParamPart;
		// filling instructions for the whole mtx param cbk
		int* fillWithJumpsLinSys;
		int* fillWithJumpsFixedPartLinSys;
		int* fillSequentially;

	public:
		
		MtxParamCbk(int totNzFixedPart_, int* nzRowsFixedPart_, int* nzColsFixedPart_, double* nzValsFixedPart_);
		
		~MtxParamCbk();
		
		MySparseMatrix* getFixedPart();
		
		int getTotalNz();

		int getTotalNzParamPart();
		
		void addTimesVec(double* x, double* p, double* result);
		
		void writeNzRowsAndColsSequentiallyInZeroBasedForm(int* iRow, int* jCol);
		
		void addWithJumpsNzValsLinSys(double* p, double* values);
		
		void writeSequentiallyNzVals(double* p, double* values);
		
		void calcFinalSparsityPatternAndCalcFillingInstructions();
		
		vector<int>* getFinalNzRows();
		
		vector<int>* getFinalNzCols();
		
		void addScaledPartialPTimesVec(int jOneBased, double* x, double* p, double* result, double scaling);
		
		void addScaledTransposePartialPTimesVec(int jOneBased, double* mul, double* p, double* result, double scaling);
		
		void fillFillingInstructionsForMtxCbkWithJumpsLinSys(int totalDecisionVars, map< pair<int, int> , int>* dictOfGlobalPositions);
};

#endif
