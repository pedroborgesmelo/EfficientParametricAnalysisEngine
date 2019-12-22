/*
 * author: pedro.borges.melo@gmail.com
 * date: October/2019
 */

#ifndef MY_VECTOR_H_
#define MY_VECTOR_H_

#include <random>

using namespace std;

class MyVector {

	public:

		static int getNegativePositionVec(double* vec, int size);

		static void initWithRandomValues(double* vec, int size, default_random_engine *gen, double min, double max);

		static void sumScaledVectorIntoAnother(double* from, double* to, int size, double scaling);

		static void copyScaledVecFromTo(double* from, double* to, int size, double scaling);

		static double getNormOfDiffVec(double* vec1, double* vec2, int size);

		static double normVec(double* vec, int size);
		
		static void initVecWithVal(int* vec, int size, int val);

		static void initVecWithVal(double* vec, int size, double val);

		static void copyVecFromTo(double* from, double* to, int size);

		static void shiftVecToMinVal(double* vec, int size, double minVal);

		static double computeInnerProduct(double* vec1, double* vec2, int size);

		static double computeSquaredNorm(double* vec1, int size);
		
		static void mulVecByScalar(double* vec, int size, double scalar);

		static void printVec(char* vecName, double* vec, int size);

		// computes vec1'mtx*vec2 with sparse information
		// NOTE: rows and cols are one-based
		static double prodVecMtxVec(double* vec1, double* mtx, double* vec2, int* rows, int* cols, int totNz);

		static double prodVecMtxVecFull(double* vec1, double* mtx, double* vec2, int size);

		// computes vec1'mtx*(e_j) with sparse information
		// NOTE: rows and cols are one-based
		static double prodVecColMtx(double* vec1, double* mtx, int* rows, int* cols, int totNz, int colOneBased);
};

#endif
