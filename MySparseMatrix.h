/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#ifndef SPARSE_MATRIX_H_
#define SPARSE_MATRIX_H_

class MySparseMatrix {

	private:

		int totalNz;
		int* nzRowsOneBased;
		int* nzColsOneBased;
		double* nzVals;

	public:
		
		MySparseMatrix(int totalNz_, int* nzRowsOneBased_, int* nzColsOneBased_, double* nzVals_);
		
		~MySparseMatrix();
		
		int getTotalNz();
		
		int* getNzRowsOneBased();
		
		int* getNzColsOneBased();
		
		double* getNzVals();

		void addTimesVec(double* vec, double* result);

		void writeNzRowsAndColsSequentiallyInZeroBasedForm(int* iRow, int* jCol);
		
		void writeSequentiallyNzVals(double* values);

		void writeGeneralNzVals(int* dict, double* values);
		
		void addSequentiallyNzVals(double* values);

		void addGeneralNzVals(int* dict, double* values);

		double vecMtxVec(double* vec1, double* vec2);
};

#endif
