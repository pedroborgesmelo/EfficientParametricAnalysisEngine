/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#ifndef PARAM_CBK_H_
#define PARAM_CBK_H_

#include <map>
#include <vector>
#include <functional>

using namespace std;

class ParamCbk {
	
	private:
		
		// instructions to fill grad x
		int* fillWithJumpsGradX;
		int* fillSequentiallyGradX;
		int* fillWithJumpsGradXGradXtranspose;
		// instructions to fill grad p
		int* fillWithJumpsGradP;
		int* fillSequentiallyGradP;
		// instructions to fill hess x
		int* fillWithJumpsHessX;
		int* fillSequentiallyHessX;
		int* fillWithJumpsLinSysHessX;
		// instructions to fill partial p hess x
		int* fillWithJumpsPartialPGradX;
		int* fillSequentiallyPartialPGradX;
		// buffers
		double* bufferGradX;
		
	public:
		
		ParamCbk();
		
		~ParamCbk();
		
		void initInnerSparsityRelatedStructures();

		virtual double getValueFunc(double* x, double* p) = 0;
		
		virtual void addScaledGradXFunc(double* x, double* p, int* dict, double* values, double scaling) = 0;
		
		virtual void addScaledGradPFunc(double* x, double* p, int* dict, double* values, double scaling) = 0;

		virtual void addScaledHessXFunc(double* x, double* p, int* dict, double* values, double scaling) = 0;
		
		virtual void addScaledPartialPGradXFunc(int jOneBased, double* x, double* p, int* dict, double* values, double scaling) = 0;
		
		virtual double getPartialPFunc(int jOneBased, double* x, double*p) = 0;
		
		double getValue(double* x, double* p);
		
		void addSequentiallyGradX(double* x, double* p, double* gradX);

		void addWithJumpsGradX(double* x, double* p, double* gradX);

		void addScaledWithJumpsGradX(double* x, double* p, double* gradX, double scaling);
		
		void addScaledWithJumpsLinSysGradXGradXTranspose(double* x, double* p, double* sortOfHessX, double scaling);
		
		void addSequentiallyGradP(double* x, double* p, double* gradP);
		
		void addWithJumpsGradP(double* x, double* p, double* gradP);
		
		void addSequentiallyScaledHessX(double* x, double* p, double* hessX, double scaling);
		
		void addWithJumpsScaledHessX(double* x, double* p, double* hessX, double scaling);
		
		void addWithJumpsLinSysScaledHessX(double* x, double* p, double* hessX, double scaling);
		
		void addScaledWithJumpsPartialPGradX(int jOneBased, double* x, double* p, double* partialPGradX, double scaling);
		
		double getPartialP(int jOneBased, double* x, double* p);
		
		virtual vector<int>* getNzPosOneBasedGradX() = 0;
		
		virtual int getTotalNzGradX() = 0;
		
		virtual vector<int>* getNzPosOneBasedGradP() = 0;
		
		virtual int getTotalNzGradP() = 0;
		
		virtual vector<int>* getNzRowsOneBasedHessX() = 0;
		
		virtual vector<int>* getNzColsOneBasedHessX() = 0;
		
		virtual int getTotalNzHessX() = 0;

		virtual vector<int>* getNzRowsOneBasedPartialPGradX() = 0;
		
		virtual int getTotalNzRowsOneBasedPartialPGradX() = 0;
		
		double computeInnerProductVecWithGradX(double* x, double* p, double* vec);
		
		void fillFillingInstructionsForHessWithJumps(map< pair<int, int> , int>* dictOfGlobalPositions);
		
		void fillFillingInstructionsForHessWithJumpsLinSys(map< pair<int, int> , int>* dictOfGlobalPositions);
		
		void fillFillingInstructionsForGradXGradXTransposeWithJumpsLinSys(map< pair<int, int> , int>* dictOfGlobalPositions);

		double* getBufferGradX();
};

#endif
