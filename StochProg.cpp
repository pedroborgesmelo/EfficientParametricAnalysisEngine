/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "StochProg.h"
#include "MyVector.h"
#include "OptValSmoCbk.h"

using namespace std;

StochProg::StochProg(int totalDecisionVars_, double* linearCost_, int totalEqConstr_,
		double* rhsValues_, int totalNzEqConstr_, int* nzRowsEqConstr_,
		int* nzColsEqConstr_, double* nzValsEqConstr_, 
		int totalScenarios_, double* probScenariosVec_, vector<Cbk*> *recourseFunctions_,
		char* pathToWrtObjIterates_) {
	this->totalDecisionVars = totalDecisionVars_;
	this->linearCost = linearCost_;
	this->totalEqConstr = totalEqConstr_;
	this->rhsValues = rhsValues_;
	this->eqConstrMtx = new MySparseMatrix(totalNzEqConstr_, nzRowsEqConstr_, nzColsEqConstr_, nzValsEqConstr_);
	this->totalScenarios = totalScenarios_;
	this->probScenariosVec = probScenariosVec_;
	this->recourseFunctions = recourseFunctions_;
	this->totalIterations = 0;
	this->totalSecondsSolvingSubproblems = 0;
	this->pathToWrtObjIterates = pathToWrtObjIterates_;
	this->lastObjMasterProblem = 1e29;
	this->totalTimesVariationWasBelowTarget = 0;
}

StochProg::~StochProg() {
	delete (this->eqConstrMtx);
}

void StochProg::addOneToTotalIterations() {
	this->totalIterations += 1;
}

int StochProg::getTotalIterations() {
	return this->totalIterations;
}

int StochProg::getTotalDecisionVars() {
	return this->totalDecisionVars;
}

double* StochProg::getLinearCost() {
	return this->linearCost;
}

int StochProg::getTotalEqConstr() {
	return this->totalEqConstr;
}

char* StochProg::getPathToWrtObjIterates() {
	return this->pathToWrtObjIterates;
}

double StochProg::getStartTime() {
	return this->startTime;
}

void StochProg::addTimeSolvingSubproblems(double toAdd) {
	this->totalSecondsSolvingSubproblems += toAdd;
}

double* StochProg::getRhsValues() {
	return this->rhsValues;
}

MySparseMatrix* StochProg::getEqConstrMtx() {
	return this->eqConstrMtx;
}

double StochProg::getProbabilitiesScenario(int s) {
	return this->probScenariosVec[s];
}

vector<Cbk*>* StochProg::getRecourseFunctions() {
	return this->recourseFunctions;
}

double StochProg::getLastObjMasterProblem() {
	return this->lastObjMasterProblem;
}

void StochProg::setLastObjMasterProblem(double newVal) {
	this->lastObjMasterProblem = newVal;
}

int StochProg::getTotalTimesVariationWasBelowTarget() {
	return this->totalTimesVariationWasBelowTarget;
}

void StochProg::addOneToTotalTimesVariationWasBelowTarget() {
	this->totalTimesVariationWasBelowTarget += 1;
}

void StochProg::setTotalTimesVariationWasBelowTarget(double newVal) {
	this->totalTimesVariationWasBelowTarget = newVal;
}

double StochProg::getBelowThisVarItIsConsideredNoise() {
	return this->belowThisVarItIsConsideredNoise;
}

int StochProg::getMaxTimesVariationWasBelowTarget() {
	return this->maxTimesVariationWasBelowTarget;
}

Bool StochProg::intermediate_callback_master_prob(
		Index alg_mod, Index iter_count, Number obj_value, Number inf_pr, 
		Number inf_du, Number mu, Number d_norm, Number regularization_size, 
		Number alpha_du, Number alpha_pr, Index ls_trials, UserDataPtr user_data_all) {
	// cast
	StochProg* stochProg = (StochProg*) user_data_all;
	// count number of times below target
	double lastObj = stochProg->getLastObjMasterProblem();
	double var = abs(obj_value - lastObj);
	double relativeVar = var / (abs(lastObj) + 1);
	if (var < stochProg->getBelowThisVarItIsConsideredNoise()) {
		stochProg->addOneToTotalTimesVariationWasBelowTarget();
	} else {
		stochProg->setTotalTimesVariationWasBelowTarget( 0 );
	}
	// update last val
	stochProg->setLastObjMasterProblem( obj_value );
	// decide if should stop
	if (stochProg->getTotalTimesVariationWasBelowTarget() > stochProg->getMaxTimesVariationWasBelowTarget()) {
		return false;
	}
	// return
	return TRUE;
}

Bool StochProg::eval_f_master_prob(Index n, Number* x, Bool new_x, Number* obj_value, UserDataPtr user_data_all) {
	// get data
	StochProg* stochProg = (StochProg*) user_data_all;
	// total types
	double total = 0.0;
	double sqpTotal = 0.0;
	// first stage cost
	total += MyVector::computeInnerProduct(stochProg->getLinearCost(), x, n);
	// copy first stage cost
	sqpTotal += total;
	// second stage cost
	vector<Cbk*>* recourseFunctionsVec = stochProg->getRecourseFunctions();
	int totalScenarios = recourseFunctionsVec->size();
	for (int s=0; s<totalScenarios; ++s) {
		// cast
		OptValSmoCbk* castedCbk = (OptValSmoCbk*) (*recourseFunctionsVec)[s];
		// begin time
		double startedSolvingSubproblem = clock();
		// solve subproblem
		double prob = stochProg->getProbabilitiesScenario(s);
		// get sqp
		double valRecFunc = castedCbk->getValue(x);
		double valSqpRecFunc = castedCbk->getSqpApproxValueFunc();
		// update totals
		total += (prob) * (valRecFunc);
		sqpTotal += (prob) * (valSqpRecFunc);
		// log total time
		double time = (clock() - startedSolvingSubproblem) / CLOCKS_PER_SEC;
		stochProg->addTimeSolvingSubproblems(time);
	}
	// set final value
	*obj_value = total;
	// log
	stochProg->addOneToTotalIterations();
	FILE* fp = fopen(stochProg->getPathToWrtObjIterates(), "a");
	double totTime = (clock() - stochProg->getStartTime()) / CLOCKS_PER_SEC;
	fprintf(fp, "UPPOPTGENERALSQP, %.6f, UPPTIME, %.6f\n", sqpTotal, totTime);
	fclose(fp);
	// return
	return TRUE;
}

Bool StochProg::eval_grad_f_master_prob(Index n, Number* x, Bool new_x, Number* grad_f, UserDataPtr user_data_all) {
	StochProg* stochProg = (StochProg*) user_data_all;
	// first stage grad
	MyVector::copyVecFromTo(stochProg->getLinearCost(), grad_f, n);
	// second stage grad
	vector<Cbk*>* recourseFunctionsVec = stochProg->getRecourseFunctions();
	int totalScenarios = recourseFunctionsVec->size();
	for (int s=0; s<totalScenarios; ++s) {
		double prob = stochProg->getProbabilitiesScenario(s);
		for (int iOneBased=1; iOneBased<=n; ++iOneBased) {
			// begin time
			double startedSolvingSubproblem = clock();
			// calc
			double partialDer = (*recourseFunctionsVec)[s]->getPartial(iOneBased, x);
			grad_f[iOneBased-1] += (prob) * (partialDer);
			// log total time
			double time = (clock() - startedSolvingSubproblem) / CLOCKS_PER_SEC;
			stochProg->addTimeSolvingSubproblems(time);
		}
	}
	// return
	return TRUE;
}

Bool StochProg::eval_g_master_prob(Index n, Number* x, Bool new_x, Index m, Number* g, UserDataPtr user_data_all) {
	// get data
	StochProg* stochProg = (StochProg*) user_data_all;
	MySparseMatrix* stochProgMtx = stochProg->getEqConstrMtx();
	// copy - rhs to output
	MyVector::copyScaledVecFromTo(stochProg->getRhsValues(), g, stochProg->getTotalEqConstr(), -1.0);
	// add mtx * x to the output
	stochProgMtx->addTimesVec(x, g);
	// return
	return TRUE;
}

Bool StochProg::eval_jac_g_master_prob(Index n, Number *x, Bool new_x, Index m, Index nele_jac, 
		Index *iRow, Index *jCol, Number *values, UserDataPtr user_data_all) {
	// get data
	StochProg* stochProg = (StochProg*) user_data_all;
	MySparseMatrix* stochProgMtx = stochProg->getEqConstrMtx();
	// sparsity pattern
	if (values == NULL) {
		// get data
		int totNz = stochProgMtx->getTotalNz();
		int* nzRows = stochProgMtx->getNzRowsOneBased();
		int* nzCols = stochProgMtx->getNzColsOneBased();
		// copy
		for (int i=0; i < totNz; ++i) {
			iRow[i] = nzRows[i] - 1;
			jCol[i] = nzCols[i] - 1;
		}
	// values
	} else {
		// get data
		int totNz = stochProgMtx->getTotalNz();
		double* vals = stochProgMtx->getNzVals();
		// copy
		for (int i=0; i < totNz; ++i) { 
			values[i] = vals[i];
		}
	}
	// return
	return TRUE;
}

Bool StochProg::eval_h_master_prob(Index n, Number *x, Bool new_x, Number obj_factor, Index m, 
		Number *lambda, Bool new_lambda, Index nele_hess, 
		Index *iRow, Index *jCol, Number *values, UserDataPtr user_data_all) {
	return TRUE;
}

void StochProg::solveStochProg(double* objValOut, double* initiaPrimalIterate, double* initialDualIterateEqConstr,
		double* optMutXvarLowerConstrOut, double* optMutXvarUpperConstrOut) {
	// log
	FILE* fp = fopen(this->getPathToWrtObjIterates(), "a");
	fprintf(fp, "UPPOPTGENERAL, starting\n");
	fclose(fp);
	// initial clock
	this->startTime = clock();
	// number of variables 
	Index n = this->totalDecisionVars;
	 // number of constraints 
	Index m = this->totalEqConstr;
	// lower bounds on g 
	Number* g_L = NULL;
	// upper bounds on g 
	Number* g_U = NULL;
	// IpoptProblem 
	IpoptProblem nlp = NULL;
	// Solve return code 
	enum ApplicationReturnStatus status;
	// generic counter
	Index i;
	// Number of nonzeros in the Jacobian of the constraints 
	Index nele_jac = this->eqConstrMtx->getTotalNz();
	// Number of nonzeros in the Hessian of the Lagrangian (lower or upper triangual part only) 
	Index nele_hess = 0;
	// indexing style for matrices 
	Index index_style = 0;
	// set the number of variables and allocate space for the bounds 
	Number* x_L = (Number*)malloc(sizeof(Number)*n);
	Number* x_U = (Number*)malloc(sizeof(Number)*n);
	// set the values for the variable bounds 
	for (i=0; i<n; i++) {
		x_L[i] = 0.0;
		x_U[i] = 2e19;
	}
	// set the number of constraints and allocate space for the bounds 
	if (m > 0) {
		g_L = (Number*)malloc(sizeof(Number)*m);
		g_U = (Number*)malloc(sizeof(Number)*m);
		// set the values of the constraint bounds 
		for (int i=0; i<m; ++i) {
			g_L[i] = 0.0;
			g_U[i] = 0.0;
		}
	}
	// create the IpoptProblem
	nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
		index_style, &(StochProg::eval_f_master_prob), 
		&(StochProg::eval_g_master_prob), 
		&(StochProg::eval_grad_f_master_prob), 
		&(StochProg::eval_jac_g_master_prob),
		&(StochProg::eval_h_master_prob));
	// We can free the memory now - the values for the bounds have been
	// copied internally in CreateIpoptProblem
	free(x_L);
	free(x_U);
	if (m > 0) {
		free(g_L);
		free(g_U);
	}
	// Set some options.  Note the following ones are only examples,
	// they might not be suitable for your problem.
	AddIpoptIntOption(nlp, "print_level", 5);
	AddIpoptIntOption(nlp, "max_iter", 200);
	AddIpoptStrOption(nlp, "warm_start_init_point", "yes");
	AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
	AddIpoptNumOption(nlp, "tol", 1e-9);
	// intermediate callback
	SetIntermediateCallback(nlp, &(StochProg::intermediate_callback_master_prob));
	// solve the problem
	status = IpoptSolve(nlp, initiaPrimalIterate, NULL, objValOut, 
		initialDualIterateEqConstr, optMutXvarLowerConstrOut, 
		optMutXvarUpperConstrOut, this);
	// free
	FreeIpoptProblem(nlp);
	// log
	fp = fopen( this->getPathToWrtObjIterates() , "a");
	fprintf(fp, "UPPOPTGENERAL, TIME SUBPROBLEMS = %.6f\n", totalSecondsSolvingSubproblems);
	fprintf(fp, "UPPOPTGENERAL, leaving\n");
	fclose(fp);
}
