/*
 * author: pedro.borges.melo@gmail.com
 * date: October/2019
 */

#include "QuadParamCbk.h"
#include "MyVector.h"

QuadParamCbk::QuadParamCbk(int totalDecisionVars_,
		int totalParameters_,
		int totNzSymQuadPart_,
		int* nzRowsOneBasedSymQuadPart_,
		int* nzColsOneBasedSymQuadPart_,
		double* nzValsSymQuadPart_,
		double* linearCostX_,
		double freeTerm_,
		double* linearCostP_) : ParamCbk() {
	// set normal data
	this->totalDecisionVars = totalDecisionVars_;
	this->totalParameters = totalParameters_;
	this->totNzSymQuadPart = totNzSymQuadPart_;
	this->nzRowsOneBasedSymQuadPart = nzRowsOneBasedSymQuadPart_,
	this->nzColsOneBasedSymQuadPart = nzColsOneBasedSymQuadPart_;
	this->nzValsSymQuadPart = nzValsSymQuadPart_;
	this->linearCostX = linearCostX_;
	this->freeTerm = freeTerm_;
	this->linearCostP = linearCostP_;
	// alloc sparsity vectors
	this->nzPosOneBasedGradX = new vector<int>();
	this->nzPosOneBasedGradP = new vector<int>();
	this->nzRowsOneBasedHessX = new vector<int>();
	this->nzColsOneBasedHessX = new vector<int>();
	this->nzRowsOneBasedPartialPGradX = new vector<int>();
	// fill grad x sparsity
	for ( int i=0; i<totalDecisionVars_; ++i ) {
		this->nzPosOneBasedGradX->push_back( i+1 );
	}
	// fill grad p sparsity
	for ( int i=0; i<totalParameters_; ++i ) {
		if ( this->linearCostP[i] != 0) {
			this->nzPosOneBasedGradP->push_back( i+1 );
		}
	}
	// fill hess x
	for ( int i=0; i<totNzSymQuadPart_; ++i ) {
		this->nzRowsOneBasedHessX->push_back( nzRowsOneBasedSymQuadPart_[i] );
		this->nzColsOneBasedHessX->push_back( nzColsOneBasedSymQuadPart_[i] );
	}
	// final step
	initInnerSparsityRelatedStructures();
}

QuadParamCbk::~QuadParamCbk() {
	delete this->nzPosOneBasedGradX;
	delete this->nzPosOneBasedGradP;
	delete this->nzRowsOneBasedHessX;
	delete this->nzColsOneBasedHessX;
	delete this->nzRowsOneBasedPartialPGradX;
}

double QuadParamCbk::getValueFunc(double* x, double* p) {
	// init
	double total = this->freeTerm;
	// linear term x
	total += MyVector::computeInnerProduct(x, this->linearCostX, this->totalDecisionVars);
	// quadratic part
	for ( int i=0; i<(this->totNzSymQuadPart); ++i) {
		// get data
		int rowOneBased = this->nzRowsOneBasedSymQuadPart[i];
		int colOneBased = this->nzColsOneBasedSymQuadPart[i];
		double val = this->nzValsSymQuadPart[i];
		// add
		if (rowOneBased == colOneBased) {
			total += val * x[ rowOneBased-1] * x[ colOneBased-1 ];
		} else {
			total += 2.0 * val * x[ rowOneBased-1] * x[ colOneBased-1 ];
		}
	}
	// linear term x
	total += MyVector::computeInnerProduct(p, this->linearCostP, this->totalParameters);
	// return
	return total;
}

void QuadParamCbk::addScaledGradXFunc(double* x, double* p, int* dict, double* values, double scaling) {
	// add linear cost
	for ( int i=0; i<(this->totalDecisionVars) ; ++i ) {
		values[ dict[i] ] += scaling * (this->linearCostX[i]);
	}
	// add contrib from quadratic term
	for ( int i=0; i<(this->totNzSymQuadPart); ++i) {
		// get data
		int rowOneBased = this->nzRowsOneBasedSymQuadPart[i];
		int colOneBased = this->nzColsOneBasedSymQuadPart[i];
		double val = this->nzValsSymQuadPart[i];
		// add
		if (rowOneBased == colOneBased) {
			values[ dict[ rowOneBased-1 ] ] += scaling * 2.0 * val * x[ colOneBased-1 ];
		} else {
			values[ dict[ rowOneBased-1 ] ] += scaling * 2.0 * val * x[ colOneBased-1 ];
			values[ dict[ colOneBased-1 ] ] += scaling * 2.0 * val * x[ rowOneBased-1 ];
		}
	}
}

void QuadParamCbk::addScaledGradPFunc(double* x, double* p, int* dict, double* values, double scaling) {
	int totNzGradP = this->nzPosOneBasedGradP->size();
	for ( int i=0; i<totNzGradP; ++i ) {
		// get data
		int jOneBased = (*nzPosOneBasedGradP)[i];
		double partial = getPartialPFunc(jOneBased, x, p);
		// add
		values[ dict[i] ] += scaling * partial;
	}
}

void QuadParamCbk::addScaledHessXFunc(double* x, double* p, int* dict, double* values, double scaling) {
	// add contrib from quadratic term
	for ( int i=0; i<(this->totNzSymQuadPart); ++i) {
		// get data
		int rowOneBased = this->nzRowsOneBasedSymQuadPart[i];
		int colOneBased = this->nzColsOneBasedSymQuadPart[i];
		double val = this->nzValsSymQuadPart[i];
		// add
		values[ dict[i] ] += scaling * 2.0 * val;
	}
}

void QuadParamCbk::addScaledPartialPGradXFunc(int jOneBased, double* x, double* p, int* dict, double* values, double scaling) {
	// NOTHING
}

double QuadParamCbk::getPartialPFunc(int jOneBased, double* x, double*p) {
	if ( linearCostP != NULL) {
		return this->linearCostP[jOneBased-1];
	} else {
		return 0.0;
	}
}

vector<int>* QuadParamCbk::getNzPosOneBasedGradX() {
	return this->nzPosOneBasedGradX;
}

int QuadParamCbk::getTotalNzGradX() {
	return this->nzPosOneBasedGradX->size();
}

vector<int>* QuadParamCbk::getNzPosOneBasedGradP() {
	return this->nzPosOneBasedGradP;
}

int QuadParamCbk::getTotalNzGradP() {
	return this->nzPosOneBasedGradP->size();
}

vector<int>* QuadParamCbk::getNzRowsOneBasedHessX() {
	return this->nzRowsOneBasedHessX;
}

vector<int>* QuadParamCbk::getNzColsOneBasedHessX() {
	return this->nzColsOneBasedHessX;
}

int QuadParamCbk::getTotalNzHessX() {
	return this->totNzSymQuadPart;
}

vector<int>* QuadParamCbk::getNzRowsOneBasedPartialPGradX() {
	return this->nzRowsOneBasedPartialPGradX;
}

int QuadParamCbk::getTotalNzRowsOneBasedPartialPGradX() {
	return this->nzRowsOneBasedPartialPGradX->size();
}
