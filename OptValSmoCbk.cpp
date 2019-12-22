/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#include <iostream>

#include "OptValSmoCbk.h"

using namespace std;

OptValSmoCbk::OptValSmoCbk(ParamConvexProb* paramConvexProb_,
		double epsilon_, double mu_) : Cbk( paramConvexProb_->getTotalParameters() ) {
	// set normal data
	this->paramConvexProb = paramConvexProb_;
	this->epsilon = epsilon_;
	this->mu = mu_;
}

OptValSmoCbk::~OptValSmoCbk() {
	// NOTHING
}

double OptValSmoCbk::getValueFunc(double* x) {
	// update buffer
	this->paramConvexProb->updateRegBuffer(
		this->epsilon, this->mu, x);
	// get smoothed opt val
	double total = this->paramConvexProb->getLastRegOptVal();
	// return
	return total;
}

double OptValSmoCbk::getSqpApproxValueFunc() {
	// return
	return this->paramConvexProb->getSqpOptVal();
}

void OptValSmoCbk::getFullGradFunc(double* x, double* grad) {
	// get pointer to objective
	ParamCbk* paramObjFunc = this->paramConvexProb->getParamObjFunction();
	// create linear system
	this->paramConvexProb->buildSysForDerOfSolutionMappings();
	// compute partials
	int totalParameters = this->paramConvexProb->getTotalParameters();
	for (int jOneBased=1; jOneBased <= totalParameters; ++jOneBased) {
		// compute derivative of solution mappings with respect to the j-th parameter
		this->paramConvexProb->calcPartialDerSolMap(jOneBased);
		// get reg primal sol
		double* lastRegPrimalSol = this->paramConvexProb->getLastRegPrimalSol();
		// get der reg sol
		double* derSol = this->paramConvexProb->getLastDerPrimalSol();
		// init
		double total = paramObjFunc->getPartialP( jOneBased, lastRegPrimalSol, x);
		// grad contrib
		total += paramObjFunc->computeInnerProductVecWithGradX(lastRegPrimalSol, x, derSol); 
		// set grad
		grad[ jOneBased-1 ] = total;
	}
	// clean up linear system
	this->paramConvexProb->destroySysForDerOfSolutionMappings();
}
