/*
 * author: pedro.borges.melo@gmail.com
 * date: October/2019
 */

#include <stdlib.h>
#include <vector>

#include "LinearSystem.h"

LinearSystem::LinearSystem(int totalEqLinearSystem_, int totalNzEqLinearSystem_, vector< pair<int, int> > *nzPositionsLinSys) {
	// set dimensions
	this->totalEqLinearSystem = totalEqLinearSystem_;
	this->totalNzEqLinearSystem = totalNzEqLinearSystem_;
	// alloc
	this->iaPardiso = (MKL_INT*) malloc(sizeof(MKL_INT) * (totalEqLinearSystem_ + 1));
	this->jaPardiso = (MKL_INT*) malloc(sizeof(MKL_INT) * totalNzEqLinearSystem_);
	this->valsPardiso = (double*) malloc(sizeof(double) * totalNzEqLinearSystem_);
	this->rhsPardiso = (double*) malloc(sizeof(double) * totalEqLinearSystem_);
	this->solPardiso = (double*) malloc(sizeof(double) * totalEqLinearSystem_);
	// fill data
	this->createSparsityStructureWithPairs( nzPositionsLinSys );
}

void LinearSystem::createSparsityStructureWithPairs( vector< pair<int, int> > *nzPositionsLinSys_ ) {
	// index counters
	int globalIndexCol = 0;
	// init lines
	this->iaPardiso[0] = 1;
	this->iaPardiso[totalEqLinearSystem] = totalEqLinearSystem + 1;
	// iterate
	for ( vector< pair<int, int>>::iterator it = nzPositionsLinSys_->begin(); it != nzPositionsLinSys_->end(); ++it ) {
		// get data
		int row = it->first;
		int col = it->second;
		// set sparsity
		this->jaPardiso[globalIndexCol] = col;
		globalIndexCol += 1;
		this->iaPardiso[row] = globalIndexCol + 1;
	}
}

LinearSystem::~LinearSystem() {
	free(this->iaPardiso);
	free(this->jaPardiso);
	free(this->valsPardiso);
	free(this->rhsPardiso);
	free(this->solPardiso);
}

int LinearSystem::getTotalEqLinearSystem() {
	return this->totalEqLinearSystem;
}

int LinearSystem::getTotalNzEqLinearSystem() {
	return this->totalNzEqLinearSystem;
}

MKL_INT* LinearSystem::getIaPardiso() {
	return this->iaPardiso;
}

MKL_INT* LinearSystem::getJaPardiso() {
	return this->jaPardiso;
}

double* LinearSystem::getValsPardiso() {
	return this->valsPardiso;
}

double* LinearSystem::getRhsPardiso() {
	return this->rhsPardiso;
}

double* LinearSystem::getSolPardiso() {
	return this->solPardiso;
}
