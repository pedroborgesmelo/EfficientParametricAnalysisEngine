/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#include "MySparseMatrix.h"

MySparseMatrix::MySparseMatrix(int totalNz_, int* nzRowsOneBased_, int* nzColsOneBased_, double* nzVals_) {
	this->totalNz = totalNz_;
	this->nzRowsOneBased = nzRowsOneBased_;
	this->nzColsOneBased = nzColsOneBased_;
	this->nzVals = nzVals_;
}

MySparseMatrix::~MySparseMatrix() {
	// NOTHING
}

int MySparseMatrix::getTotalNz() {
	return this->totalNz;
}

int* MySparseMatrix::getNzRowsOneBased() {
	return this->nzRowsOneBased;
}

int* MySparseMatrix::getNzColsOneBased() {
	return this->nzColsOneBased;
}

double* MySparseMatrix::getNzVals() {
	return this->nzVals;
}

void MySparseMatrix::addTimesVec(double* vec, double* result) {
	// iterate
	for (int i=0; i<(this->totalNz); ++i) {
		// get data
		int r = this->nzRowsOneBased[i];
		int c = this->nzColsOneBased[i];
		double v = this->nzVals[i];
		// add
		result[r-1] += (v) * (vec[c-1]);
	}
}

void MySparseMatrix::writeNzRowsAndColsSequentiallyInZeroBasedForm(int* iRow, int* jCol) {
	// iterate
	for (int i=0; i<(this->totalNz); ++i) {
		iRow[i] = this->nzRowsOneBased[i] - 1;
		jCol[i] = this->nzColsOneBased[i] - 1;
	}
}

void MySparseMatrix::writeSequentiallyNzVals(double* values) {
	// iterate
	for (int i=0; i<(this->totalNz); ++i) {
		values[i] = this->nzVals[i];
	}
}

void MySparseMatrix::writeGeneralNzVals(int* dict, double* values) {
	// iterate
	for (int i=0; i<(this->totalNz); ++i) {
		values[ dict[i] ] = this->nzVals[i];
	}
}

void MySparseMatrix::addSequentiallyNzVals(double* values) {
	// iterate
	for (int i=0; i<(this->totalNz); ++i) {
		values[i] += this->nzVals[i];
	}
}

void MySparseMatrix::addGeneralNzVals(int* dict, double* values) {
	// iterate
	for (int i=0; i<(this->totalNz); ++i) {
		values[ dict[i] ] += this->nzVals[i];
	}
}

double MySparseMatrix::vecMtxVec(double* vec1, double* vec2) {
	// init
	double total = 0.0;
	// iterate
	for (int i=0; i<(this->totalNz); ++i) {
		int rowOneBased = nzRowsOneBased[i];
		int colOneBased = nzColsOneBased[i];
		double value = nzVals[i];
		total += value * ( vec1[ rowOneBased-1 ] ) * ( vec2[ colOneBased-1 ] );
	}
	// return
	return total;
}
