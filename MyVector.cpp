/*
 * author: pedro.borges.melo@gmail.com
 * date: October/2019
 */

#include <math.h>
#include <iostream>
#include <random>

#include "MyVector.h"

using namespace std;

int MyVector::getNegativePositionVec(double* vec, int size) {
	// search
	for (int i=0; i<size; ++i) {
		if (vec[i] < 0.0) {
			return i;
		}
	}
	// return
	return -1;
}

void MyVector::initWithRandomValues(double* vec, int size, default_random_engine *gen, double min, double max) {
	// gen
	uniform_real_distribution<double> dist(min, max);
	// sample
	for (int i=0;i<size; ++i) {
		vec[i] = dist(*gen);
	}
}

void MyVector::sumScaledVectorIntoAnother(double* from, double* to, int size, double scaling) {
	// sum
	for (int i=0; i<size; ++i) {
		to[i] += (scaling * from[i]);
	}
}

void MyVector::copyScaledVecFromTo(double* from, double* to, int size, double scaling) {
	// copy
	for (int i=0; i<size; ++i) {
		to[i] = (scaling * from[i]);
	}
}

double MyVector::getNormOfDiffVec(double* vec1, double* vec2, int size) {
	// init
	double total = 0.0;
	// sum
	for (int i=0; i<size; ++i) {
		total += (vec1[i]-vec2[i])*(vec1[i]-vec2[i]);
	}
	// return
	return sqrt(total);
}

double MyVector::normVec(double* vec, int size) {
	// init
	double total = 0.0;
	// norm
	for (int i=0; i<size; ++i) {
		total += (vec[i])*(vec[i]);
	}
	// return
	return sqrt(total);
}

void MyVector::initVecWithVal(int* vec, int size, int val) {
	// init
	for (int i=0; i<size; ++i) {
		vec[i] = val;
	}
}

void MyVector::initVecWithVal(double* vec, int size, double val) {
	// init
	for (int i=0; i<size; ++i) {
		vec[i] = val;
	}
}

void MyVector::copyVecFromTo(double* from, double* to, int size) {
	// copy
	for (int i=0; i<size; ++i) {
		to[i] = from[i];
	}
}

void MyVector::shiftVecToMinVal(double* vec, int size, double minVal) {
	// shift
	for (int i=0; i<size; ++i) {
		if (vec[i] < minVal) {
			vec[i] = minVal;
		}
	}
}

double MyVector::computeInnerProduct(double* vec1, double* vec2, int size) {
	// init
	double total = 0.0;
	// calc
	for(int i=0; i<size; ++i) {
		total += (vec1[i] * vec2[i]);
	}
	// return
	return total;
}

void MyVector::mulVecByScalar(double* vec, int size, double scalar) {
	// multiply
	for(int i=0; i<size; ++i) {
		vec[i] *= scalar;
	}
}

void MyVector::printVec(char* vecName, double* vec, int size) {
	// print
	cout << vecName << endl;
	for(int i=0; i<size; ++i) {
		cout << vec[i] << endl; 
	}
	cout << endl;
}

double MyVector::prodVecMtxVec(double* vec1, double* mtx, double* vec2, int* rows, int* cols, int totNz) {
	// init
	double total = 0.0;
	// calc
	for (int i=0; i<totNz; ++i) {
		// get data
		int r = rows[i] - 1;
		int c = cols[i] - 1;
		double v = mtx[i];
		// add
		total += (v * vec1[r] * vec2[c]);
	}
	// return
	return total;
}

double MyVector::prodVecColMtx(double* vec1, double* mtx, int* rows, int* cols, int totNz, int colOneBased) {
	// init
	double total = 0.0;
	// calc
	for (int i=0; i<totNz; ++i) {
		// get data
		int r = rows[i] - 1;
		int c = cols[i];
		double v = mtx[i];
		// add
		if(c == colOneBased) {
			total += (v * vec1[r]);
		}
	}
	// return
	return total;
}

double MyVector::prodVecMtxVecFull(double* vec1, double* mtx, double* vec2, int size) {
	// init
	double total = 0.0;
	// calc
	int i = 0;
	for (int r=0; r<size; ++r) {
		for (int c=r; c<size; ++c) {
			// get data
			double v = mtx[i];
			// add
			total += (v * vec1[r] * vec2[c]);
			i += 1;
		}
	}
	// add lower dial part
	total *= 2.0;
	// return
	return total;
}

double MyVector::computeSquaredNorm(double* vec1, int size) {
	// return
	return MyVector::computeInnerProduct(vec1, vec1, size);
}
