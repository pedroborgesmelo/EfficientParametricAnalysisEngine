/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#include "Cbk.h"
#include "MyVector.h"

Cbk::Cbk(int totalVariables_) {
	// set normal data
	this->totalVariables = totalVariables_;
	// value buffer
	this->valueLastX = (double*) malloc(sizeof(double) * totalVariables_);
	this->valueBkp = 1e29;
	// grad buffer
	this->gradLastX = (double*) malloc(sizeof(double) * totalVariables_);
	this->gradBkp = (double*) malloc(sizeof(double) * totalVariables_);
	// init buffers
	for (int i=0; i<totalVariables_; ++i) {
		this->valueLastX[i] = 1e29;
		this->gradLastX[i] = 1e29;
		this->gradBkp[i] = 1e29;
	}
}

Cbk::~Cbk() {
	free(this->valueLastX);
	free(this->gradLastX);
	free(this->gradBkp);
}

int Cbk::getTotalVariables() {
	return this->totalVariables;
}

double Cbk::getValue(double* x) {
	// check
	if (MyVector::getNormOfDiffVec(x, this->valueLastX, this->totalVariables) > BUFFER_TOL) {
		MyVector::copyVecFromTo(x, this->valueLastX, this->totalVariables);
		this->valueBkp = this->getValueFunc(x);
	}
	// return
	return this->valueBkp;
}

double Cbk::getPartial(int jOneBased, double* x) {
	// check value buffer
	this->getValue(x);
	// check grad buffer
	if (MyVector::getNormOfDiffVec(x, this->gradLastX, this->totalVariables) > BUFFER_TOL) {
		MyVector::copyVecFromTo(x, this->gradLastX, this->totalVariables);
		this->getFullGradFunc(x, this->gradBkp);
	}
	// return
	return this->gradBkp[jOneBased-1];
}
