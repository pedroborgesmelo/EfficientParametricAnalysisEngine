/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#include <set>
#include <algorithm>

#include "MtxParamCbk.h"

MtxParamCbk::MtxParamCbk(int totNzFixedPart_, int* nzRowsFixedPart_, int* nzColsFixedPart_, double* nzValsFixedPart_) {
	// fixed part
	this->fixedPart = new MySparseMatrix(totNzFixedPart_, 
		nzRowsFixedPart_, nzColsFixedPart_, nzValsFixedPart_);
	// oaram part
	this->nzRowsOneBasedParamPart = NULL;
	this->nzColsOneBasedParamPart = NULL;
	this->nzValsParamPart = NULL;
	// dimensions
	int totNzFixedPart = this->fixedPart->getTotalNz();
	int totNzParamPart = this->getTotalNzParamPart();
	// inner filling instructions
	this->fillInstructFixedPart = (int*) malloc(sizeof(int) * totNzFixedPart);
	this->fillInstructParamPart = (int*) malloc(sizeof(int) * totNzParamPart);
	// sparsity pattern
	calcFinalSparsityPatternAndCalcFillingInstructions();
	// filling instructions for the whole mtx param cbk
	int totNzMtxCbk = this->finalNzRows->size();
	this->fillWithJumpsLinSys = (int*) malloc(sizeof(int) * totNzMtxCbk);
	this->fillWithJumpsFixedPartLinSys = (int*) malloc(sizeof(int) * totNzFixedPart);
	this->fillSequentially = (int*) malloc(sizeof(int) * totNzMtxCbk);
	// fill default values
	for ( int i=0; i<totNzMtxCbk; ++i ) {
		this->fillWithJumpsLinSys[i] = i;
		this->fillSequentially[i] = i;
	}
	// fill default for fixed part
	for ( int i=0; i<totNzFixedPart; ++i ) {
		this->fillWithJumpsFixedPartLinSys[i] = i;
	}
}

MtxParamCbk::~MtxParamCbk() {
	// fixed part
	delete fixedPart;
	// final sparsity pattern
	delete this->finalNzRows;
	delete this->finalNzCols;
	// inner filling instructions
	free(this->fillInstructFixedPart);
	free(this->fillInstructParamPart);
	// outter filling instructions
	free(this->fillWithJumpsLinSys);
	free(this->fillWithJumpsFixedPartLinSys);
	free(this->fillSequentially);
}

MySparseMatrix* MtxParamCbk::getFixedPart() {
	return this->fixedPart;
}

int MtxParamCbk::getTotalNz() {
	return this->finalNzRows->size();
}

int MtxParamCbk::getTotalNzParamPart() {
	if (this->nzRowsOneBasedParamPart == NULL) return 0; 
	return this->nzRowsOneBasedParamPart->size();
}

void MtxParamCbk::addTimesVec(double* x, double* p, double* result) {
	// add fixed part
	this->fixedPart->addTimesVec(x, result);
	// add parametric part
	if (this->nzRowsOneBasedParamPart != NULL && this->nzColsOneBasedParamPart != NULL && this->nzValsParamPart != NULL) {
		int totNzParamPart = this->getTotalNzParamPart();
		for ( int i=0; i<totNzParamPart; ++i ) {
			// get data
			int rowOneBased = (*nzRowsOneBasedParamPart)[i];
			int colOneBased = (*nzColsOneBasedParamPart)[i];
			double value = (*nzValsParamPart)[i]->getValue(p);
			// add result
			result[ rowOneBased - 1 ] += value * ( x[ colOneBased - 1] );
		}
	}
}

void MtxParamCbk::writeNzRowsAndColsSequentiallyInZeroBasedForm(int* iRow, int* jCol) {
	int totNzAll = finalNzRows->size();
	for (int i=0; i<totNzAll; ++i) {
		iRow[i] = (*finalNzRows)[i] - 1;
		jCol[i] = (*finalNzCols)[i] - 1;
	}
}

void MtxParamCbk::addWithJumpsNzValsLinSys(double* p, double* values) {
	// fixed part
	this->fixedPart->addGeneralNzVals( this->fillWithJumpsFixedPartLinSys, values );
	// add parametric part
	if (this->nzRowsOneBasedParamPart != NULL && this->nzColsOneBasedParamPart != NULL && this->nzValsParamPart != NULL) {
		int totNzParamPart = this->getTotalNzParamPart();
		for (int i=0; i<totNzParamPart; ++i) {
			// get data
			int rowOneBased = (*(this->nzRowsOneBasedParamPart))[i];
			int colOneBased = (*(this->nzColsOneBasedParamPart))[i];
			int pos = this->fillInstructParamPart[i];
			int posToWrite = this->fillWithJumpsLinSys[pos];
			values[posToWrite] += (*nzValsParamPart)[i]->getValue(p);
		}
	}
}

void MtxParamCbk::writeSequentiallyNzVals(double* p, double* values) {
	// fixed part
	this->fixedPart->writeGeneralNzVals(this->fillInstructFixedPart, values);
	// add parametric part
	if (this->nzRowsOneBasedParamPart != NULL && this->nzColsOneBasedParamPart != NULL && this->nzValsParamPart != NULL) {
		int totNzParamPart = this->getTotalNzParamPart();
		for (int i=0; i<totNzParamPart; ++i) {
			// get data
			int rowOneBased = (*(this->nzRowsOneBasedParamPart))[i];
			int colOneBased = (*(this->nzColsOneBasedParamPart))[i];
			int pos = this->fillInstructParamPart[i];
			int posToWrite = this->fillSequentially[pos];
			values[posToWrite] = (*nzValsParamPart)[i]->getValue(p);
		}
	}
}

void MtxParamCbk::calcFinalSparsityPatternAndCalcFillingInstructions() {
	// set (without repetitions)
	set< pair<int, int> > setTmp;
	// fixed part data
	int totNzFixedPart = this->fixedPart->getTotalNz();
	int* nzRowsOneBased = this->fixedPart->getNzRowsOneBased();
	int* nzColsOneBased = this->fixedPart->getNzColsOneBased();
	// add fixed part pattern
	for (int i=0; i<totNzFixedPart; ++i) {
		setTmp.insert( make_pair(nzRowsOneBased[i], nzColsOneBased[i]));
	}
	// parametric part
	if (this->nzRowsOneBasedParamPart != NULL && this->nzColsOneBasedParamPart != NULL && this->nzValsParamPart != NULL) {
		int totNzParamPart = this->getTotalNzParamPart();
		for (int i=0; i<totNzParamPart; ++i) {
			// get data
			int rowOneBased = (*(this->nzRowsOneBasedParamPart))[i];
			int colOneBased = (*(this->nzColsOneBasedParamPart))[i];
			// insert
			setTmp.insert( make_pair(rowOneBased, colOneBased) );
		}
	}
	// put set in vector and sort
	vector< pair<int, int> > vecToSort;
	for (set< pair<int, int> >::iterator it = setTmp.begin(); it != setTmp.end(); ++it) {
		vecToSort.push_back( *it );
	}
	// sort lexicographically
	sort( vecToSort.begin(), vecToSort.end() );
	// put vector in map to calc indexes efficiently
	map< pair<int, int>, int > mapToStoreIndices;
	for (int pos = 0; pos< vecToSort.size() ; ++pos) {
		mapToStoreIndices[ vecToSort[pos] ] = pos; 
	}
	// alloc
	this->finalNzRows = new vector<int>();
	this->finalNzCols = new vector<int>();
	// get data from set without repetitions
	for ( vector< pair<int, int>>::iterator it = vecToSort.begin(); it != vecToSort.end(); ++it ) {
		this->finalNzRows->push_back( it->first );
		this->finalNzCols->push_back( it->second );
	}
	// filling instructions for fixed part
	int* nzRowsFixedPart = this->fixedPart->getNzRowsOneBased();
	int* nzColsFixedPart = this->fixedPart->getNzColsOneBased();
	for ( int i=0; i<totNzFixedPart; ++i ) {
		int rowOneBased = nzRowsFixedPart[i];
		int colOneBased = nzColsFixedPart[i];
		pair<int, int> pIter = make_pair(rowOneBased, colOneBased);
		this->fillInstructFixedPart[i] = mapToStoreIndices[pIter];
	}
	// filling instructions for param part
	int totNzParamPart = this->getTotalNzParamPart();
	for ( int i=0; i<totNzParamPart; ++i ) {
		int rowOneBased = (*nzRowsOneBasedParamPart)[i];
		int colOneBased = (*nzColsOneBasedParamPart)[i];
		pair<int, int> pIter = make_pair(rowOneBased, colOneBased);
		this->fillInstructParamPart[i] = mapToStoreIndices[pIter];
	}
}

vector<int>* MtxParamCbk::getFinalNzRows() {
	return this->finalNzRows;
}

vector<int>* MtxParamCbk::getFinalNzCols() {
	return this->finalNzCols;
}

void MtxParamCbk::addScaledPartialPTimesVec(int jOneBased, double* x, double* p, double* result, double scaling) {
	int totNzParamPart = this->getTotalNzParamPart();
	for ( int i=0; i<totNzParamPart; ++i ) {
		// get data
		int rowOneBased = (*nzRowsOneBasedParamPart)[i];
		int colOneBased = (*nzColsOneBasedParamPart)[i];
		double partial = (*nzValsParamPart)[i]->getPartial(jOneBased, p);
		// add result
		result[ rowOneBased-1 ] += scaling * partial * ( x[ colOneBased-1 ] );
	}
}

void MtxParamCbk::addScaledTransposePartialPTimesVec(int jOneBased, double* mul, double* p, double* result, double scaling) {
	int totNzParamPart = this->getTotalNzParamPart();
	for ( int i=0; i<totNzParamPart; ++i ) {
		// get data
		int rowOneBased = (*nzRowsOneBasedParamPart)[i];
		int colOneBased = (*nzColsOneBasedParamPart)[i];
		double partial = (*nzValsParamPart)[i]->getPartial(jOneBased, p);
		// add result
		result[ colOneBased-1 ] += scaling * partial * ( mul[ rowOneBased-1 ] );
	}
}

void MtxParamCbk::fillFillingInstructionsForMtxCbkWithJumpsLinSys(int totalDecisionVars, map< pair<int, int> , int>* dictOfGlobalPositions) {
	// NOTE: for the linear system it is the transpose that enters
	// for the whole mtx cbk
	int totNzAll = this->finalNzRows->size();
	for ( int i=0; i<totNzAll; ++i ) {
		int rowOneBased = (*finalNzCols)[i];
		int colOneBased = totalDecisionVars + (*finalNzRows)[i];
		pair<int, int> pIter = make_pair( rowOneBased, colOneBased );
		this->fillWithJumpsLinSys[i] = (*dictOfGlobalPositions)[pIter];
	}
	// for the fixed part only
	int totNzFixedPart = this->fixedPart->getTotalNz();
	for ( int i=0; i<totNzFixedPart; ++i ) {
		int posToStoreLinSys = this->fillInstructFixedPart[i];
		this->fillWithJumpsFixedPartLinSys[i] = this->fillWithJumpsLinSys[posToStoreLinSys];	
	}
}
