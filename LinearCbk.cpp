/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#include <stdlib.h>
#include <map>

#include "LinearCbk.h"
#include "MyVector.h"

LinearCbk::LinearCbk(int totalVariables_,
		int totNz_, 
		int* nzRowsOneBasedLinearCost_,
		double* nzValsLinearCost_,
		double freeTerm_,
		double scalingLinearCost_) : Cbk(totalVariables_) {
	// set normal data
	this->totNz = totNz_;
	this->nzRowsOneBasedLinearCost = nzRowsOneBasedLinearCost_;
	this->nzValsLinearCost = nzValsLinearCost_;
	this->freeTerm = freeTerm_;
	this->scalingLinearCost = scalingLinearCost_;
	// create map
	this->nzRowToValPos = new map<int, int>();
	for ( int i=0; i<totNz; ++i ) {
		int rowOneBased = this->nzRowsOneBasedLinearCost[i];
		(*nzRowToValPos)[ rowOneBased ] = i;
	}
}

LinearCbk::~LinearCbk() {
	delete this->nzRowToValPos;
}

double LinearCbk::getValueFunc(double* x) {
	// init
	double total = this->freeTerm;
	// linear cost
	for ( int i=0; i<totNz; ++i ) {
		// get data
		int rowOneBased = this->nzRowsOneBasedLinearCost[i];
		double val = this->nzValsLinearCost[i];
		// add
		total += scalingLinearCost * val * ( x[ rowOneBased-1 ] );
	}
	// return
	return total;
}

void LinearCbk::getFullGradFunc(double* x, double* grad) {
	// init with zeros
	MyVector::initVecWithVal(grad, this->getTotalVariables(), 0.0);
	// linear cost
	for ( int i=0; i<totNz; ++i ) {
		// get data
		int rowOneBased = this->nzRowsOneBasedLinearCost[i];
		double val = this->nzValsLinearCost[i];
		// add
		grad[ rowOneBased-1 ] = scalingLinearCost * val;
	}
}

double LinearCbk::getPartialFunc(int jOneBased, double* x) {
	map<int, int>::iterator it = this->nzRowToValPos->find(jOneBased);
	if ( it != this->nzRowToValPos->end() ) {
		return scalingLinearCost * (this->nzValsLinearCost[ it->second ]);
	} else {
		return 0.0;
	}
}
