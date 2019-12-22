/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#include "MyVector.h"

#include "ParamCbk.h"

ParamCbk::ParamCbk() {
	// NOTHING
}

void ParamCbk::initInnerSparsityRelatedStructures() {
	// filling patterns for grad x
	int totNzGradX = this->getTotalNzGradX();
	this->fillWithJumpsGradX = (int*) malloc(sizeof(int) * totNzGradX);
	this->fillSequentiallyGradX = (int*) malloc(sizeof(int) * totNzGradX);
	this->fillWithJumpsGradXGradXtranspose = (int*) malloc(sizeof(int) * totNzGradX * (totNzGradX + 1) * 0.5);
	// filling patterns for grad p
	int totNzGradP = this->getTotalNzGradP();
	this->fillWithJumpsGradP = (int*) malloc(sizeof(int) * totNzGradP);
	this->fillSequentiallyGradP = (int*) malloc(sizeof(int) * totNzGradP);
	// instructions to fill hess x
	int totNzHessX = this->getTotalNzHessX();
	this->fillWithJumpsHessX = (int*) malloc(sizeof(int) * totNzHessX);
	this->fillSequentiallyHessX = (int*) malloc(sizeof(int) * totNzHessX);
	this->fillWithJumpsLinSysHessX = (int*) malloc(sizeof(int) * totNzHessX);
	// instructions to fill partial p hess x
	int totNzPartialPGradX = this->getTotalNzRowsOneBasedPartialPGradX();
	this->fillWithJumpsPartialPGradX = (int*) malloc(sizeof(int) * totNzPartialPGradX);
	this->fillSequentiallyPartialPGradX = (int*) malloc(sizeof(int) * totNzPartialPGradX);
	// grad x filling
	vector<int>* nzPosOneBasedGradX = this->getNzPosOneBasedGradX();
	for (int i=0; i<totNzGradX; ++i) {
		// get data
		int posOneBased = (*nzPosOneBasedGradX)[i];
		// fill
		this->fillWithJumpsGradX[i] = posOneBased-1;
		this->fillSequentiallyGradX[i] = i;
	}
	// grad p filling
	vector<int>* nzPosOneBasedGradP = this->getNzPosOneBasedGradP();
	for (int i=0; i<totNzGradP; ++i) {
		// get data
		int posOneBased = (*nzPosOneBasedGradP)[i];
		// fill
		this->fillWithJumpsGradP[i] = posOneBased-1;
		this->fillSequentiallyGradP[i] = i;
	}
	// hess x filling
	for ( int i=0; i<totNzHessX; ++i ) {
		this->fillWithJumpsHessX[i] = i;
		this->fillSequentiallyHessX[i] = i;
		this->fillWithJumpsLinSysHessX[i] = i;
	}
	// instructions to fill partial p hess x
	vector<int>* nzRowsOneBasedPartialPGradX = this->getNzRowsOneBasedPartialPGradX();
	for ( int i=0; i<totNzPartialPGradX; ++i ) {
		// get data
		int posOneBased = (*nzRowsOneBasedPartialPGradX)[i];
		// fill
		this->fillWithJumpsPartialPGradX[i] = posOneBased - 1;
		this->fillSequentiallyPartialPGradX[i] = i;
	}
	// alloc buffers
	this->bufferGradX = (double*) malloc(sizeof(double) * totNzGradX);
}

ParamCbk::~ParamCbk() {
	free(this->fillWithJumpsGradX);
	free(this->fillSequentiallyGradX);
	free(this->fillWithJumpsGradXGradXtranspose);
	free(this->fillWithJumpsGradP);
	free(this->fillSequentiallyGradP);
	free(this->fillWithJumpsHessX);
	free(this->fillSequentiallyHessX);
	free(this->fillWithJumpsLinSysHessX);
	free(this->fillWithJumpsPartialPGradX);
	free(this->fillSequentiallyPartialPGradX);
	free(this->bufferGradX);
}

double ParamCbk::getValue(double* x, double* p) {
	return this->getValueFunc(x, p);
}

void ParamCbk::addSequentiallyGradX(double* x, double* p, double* gradX) {
	this->addScaledGradXFunc(x, p, this->fillSequentiallyGradX, gradX, 1.0);
}

void ParamCbk::addWithJumpsGradX(double* x, double* p, double* gradX) {
	this->addScaledGradXFunc(x, p, this->fillWithJumpsGradX, gradX, 1.0);
}

void ParamCbk::addScaledWithJumpsGradX(double* x, double* p, double* gradX, double scaling) {
	this->addScaledGradXFunc(x, p, this->fillWithJumpsGradX, gradX, scaling);
}

void ParamCbk::addScaledWithJumpsLinSysGradXGradXTranspose(double* x, double* p, double* sortOfHessX, double scaling) {
	// dimensions
	int totNzGradX = this->getTotalNzGradX();
	vector<int>* nzPosOneBasedGradX = this->getNzPosOneBasedGradX();
	// calc grad x
	MyVector::initVecWithVal(this->bufferGradX, totNzGradX, 0.0);
	this->addSequentiallyGradX(x, p, this->bufferGradX);
	// fill with jumps
	int currentPos = 0;
	for (int k=0; k<totNzGradX; ++k) {
		for (int l=k; l<totNzGradX; ++l) {
			int posOneBased1 = (*nzPosOneBasedGradX)[k];
			int posOneBased2 = (*nzPosOneBasedGradX)[l];
			int pos = this->fillWithJumpsGradXGradXtranspose[currentPos];
			sortOfHessX[pos] += scaling * (this->bufferGradX[ posOneBased1-1 ]) * (this->bufferGradX[ posOneBased2-1 ]);
			currentPos += 1;
		}
	}
}

double* ParamCbk::getBufferGradX() {
	return this->bufferGradX;
}

void ParamCbk::addSequentiallyGradP(double* x, double* p, double* gradP) {
	this->addScaledGradPFunc(x, p, this->fillSequentiallyGradP, gradP, 1.0);
}

void ParamCbk::addWithJumpsGradP(double* x, double* p, double* gradP) {
	this->addScaledGradPFunc(x, p, this->fillWithJumpsGradP, gradP, 1.0);
}

void ParamCbk::addSequentiallyScaledHessX(double* x, double* p, double* hessX, double scaling) {
	this->addScaledHessXFunc(x, p, this->fillSequentiallyHessX, hessX, scaling);
}

void ParamCbk::addWithJumpsScaledHessX(double* x, double* p, double* hessX, double scaling) {
	this->addScaledHessXFunc(x, p, this->fillWithJumpsHessX, hessX, scaling);
}

void ParamCbk::addWithJumpsLinSysScaledHessX(double* x, double* p, double* hessX, double scaling) {
	this->addScaledHessXFunc(x, p, this->fillWithJumpsLinSysHessX, hessX, scaling);
}

void ParamCbk::addScaledWithJumpsPartialPGradX(int jOneBased, double* x, double* p, double* partialPGradX, double scaling) {
	this->addScaledPartialPGradXFunc(jOneBased, x, p, this->fillWithJumpsPartialPGradX, partialPGradX, scaling);
}

double ParamCbk::getPartialP(int jOneBased, double* x, double* p) {
	return this->getPartialPFunc(jOneBased, x, p);
}

double ParamCbk::computeInnerProductVecWithGradX(double* x, double* p, double* vec) {
	// get sparsity
	vector<int>* nzPosOneBasedGradX = this->getNzPosOneBasedGradX();
	int totNzGradX = this->getTotalNzGradX();
	// store gradient
	MyVector::initVecWithVal(this->bufferGradX, totNzGradX, 0.0);
	this->addSequentiallyGradX(x, p, this->bufferGradX);
	// inner product
	double total = 0.0;
	for (int i=0; i<totNzGradX; ++i) {
		int rowOneBased = (*nzPosOneBasedGradX)[i];
		total += (this->bufferGradX[ rowOneBased-1 ] ) * ( vec[ rowOneBased-1 ]);
	}
	// return
	return total;
}

void ParamCbk::fillFillingInstructionsForHessWithJumps(map< pair<int, int> , int>* dictOfGlobalPositions) {
	vector<int>* nzRowsOneBasedHessX = this->getNzRowsOneBasedHessX();
	vector<int>* nzColsOneBasedHessX = this->getNzColsOneBasedHessX();
	int totNzHessX = nzRowsOneBasedHessX->size();
	for (int i=0; i<totNzHessX; ++i) {
		// get data
		int row = (*nzRowsOneBasedHessX)[i];
		int col = (*nzColsOneBasedHessX)[i];
		// pair
		pair<int, int> pIter = make_pair(row, col);
		// fill instructions
		this->fillWithJumpsHessX[i] = (*dictOfGlobalPositions)[pIter];
	}
}

void ParamCbk::fillFillingInstructionsForHessWithJumpsLinSys(map< pair<int, int> , int>* dictOfGlobalPositions) {
	vector<int>* nzRowsOneBasedHessX = this->getNzRowsOneBasedHessX();
	vector<int>* nzColsOneBasedHessX = this->getNzColsOneBasedHessX();
	int totNzHessX = nzRowsOneBasedHessX->size();
	for (int i=0; i<totNzHessX; ++i) {
		// get data
		int row = (*nzRowsOneBasedHessX)[i];
		int col = (*nzColsOneBasedHessX)[i];
		// pair
		pair<int, int> pIter = make_pair(row, col);
		// fill instructions
		this->fillWithJumpsLinSysHessX[i] = (*dictOfGlobalPositions)[pIter];
	}
}

void ParamCbk::fillFillingInstructionsForGradXGradXTransposeWithJumpsLinSys(map< pair<int, int> , int>* dictOfGlobalPositions) {
	vector<int>* nzPosOneBasedGradX = this->getNzPosOneBasedGradX();
	int totNzGradX = this->getTotalNzGradX();
	int currentPos = 0;
	for (int k=0; k<totNzGradX; ++k) {
		for (int l=k; l<totNzGradX; ++l) {
			// get data
			int row = (*nzPosOneBasedGradX)[k];
			int col = (*nzPosOneBasedGradX)[l];
			// make pair
			pair<int, int> pIter = make_pair(row, col);
			// fill instruction
			fillWithJumpsGradXGradXtranspose[currentPos] = (*dictOfGlobalPositions)[pIter];
			currentPos += 1;
		}
	}
}
