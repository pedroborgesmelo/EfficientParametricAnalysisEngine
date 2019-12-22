/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#include <set>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <ilcp/cp.h>
#include <ilcplex/ilocplex.h>

#include "MyVector.h"
#include "ParamConvexProb.h"

using namespace std;

ParamConvexProb::ParamConvexProb(int totalDecisionVars_, 
		int totalParameters_, int totalEqualityConstraints_,
		ParamCbk* paramObjFunction_, MtxParamCbk* mtxParamCbk_, 
		vector<Cbk*> *paramRhsVec_, vector<ParamCbk*> *paramIneqConstrVec_) {
	// opt problem data
	this->totalDecisionVars = totalDecisionVars_;
	this->totalParameters = totalParameters_;
	this->totalEqualityConstraints = totalEqualityConstraints_;
	this->paramObjFunction = paramObjFunction_;
	this->mtxParamCbk = mtxParamCbk_;
	this->paramRhsVec = paramRhsVec_;
	this->paramIneqConstrVec = paramIneqConstrVec_;
	// sparsity patterns for ipopt
	this->totalNzJacobianAllConstraints = this->getTotalNzForJacobianOfAllConstraints();
	this->nzRowNzColToPosValMap = new map< pair<int,int> , int>(); 
	this->positionsToAddTikhonovTerms = (int*) malloc(sizeof(int) * totalDecisionVars);
	this->fillNzFoHessianOfLagrangian();
	// sqp buffer
	this->sqpOptVal = 1e29;
	this->sqpBufferHessX = (double*) malloc(sizeof(double) * (this->nzRowNzColToPosValMap->size()));
	this->sqpNzRowZeroBasedEqConstr = (int*) malloc( sizeof(int) * (totalNzJacobianAllConstraints) );
	this->sqpNzColZeroBasedEqConstr = (int*) malloc( sizeof(int) * (totalNzJacobianAllConstraints) );
	this->sqpNzValEqConstr = (double*) malloc( sizeof(double) * (totalNzJacobianAllConstraints) );
	// alloc buffer
	int totalInequalities = paramIneqConstrVec_->size();
	this->lastEpsilon = -1.0;
	this->lastMu = -1.0;
	this->lastParameterVecRegBuffer = (double*) malloc(sizeof(double) * totalParameters_);
	this->lastRegOptVal = 1e29;
	this->lastRegPrimalSol = (double*) malloc(sizeof(double) * totalDecisionVars_);
	this->lastRegDualSolEqConstraints = (double*) malloc(sizeof(double) * totalEqualityConstraints_);
	this->lastRegLowerIneqMul = (double*) malloc(sizeof(double) * totalDecisionVars_);
	this->lastRegUpperIneqMul = (double*) malloc(sizeof(double) * totalDecisionVars_);
	this->lastRegIneqConstraintValues = (double*) malloc(sizeof(double) * totalInequalities);
	// init buffer
	for (int i=0; i<totalDecisionVars_; ++i) {
		this->lastRegPrimalSol[i] = -1.0;
		this->lastRegLowerIneqMul[i] = -1.0;
		this->lastRegUpperIneqMul[i] = -1.0;
	}
	for (int i=0; i<totalEqualityConstraints_; ++i) {
		this->lastRegDualSolEqConstraints[i] = -1.0;
	}
	for (int i=0; i<totalParameters_; ++i) {
		this->lastParameterVecRegBuffer[i] = -1.0;
	}
	// calc non trivial sparsity pattern for linear system
	this->nzPositionsLinSys = new vector< pair<int, int> >();
	this->calcNonTrivialSparsityPatternForLinSys();
	// lin sys
	int totalEqLinearSystem = totalDecisionVars + totalEqualityConstraints;
	int totalNzEqLinearSystem = this->nzPositionsLinSys->size();
	this->lSysDerSolMap = new LinearSystem(totalEqLinearSystem, 
		totalNzEqLinearSystem, this->nzPositionsLinSys);
	// set filling structions
	this->setFillingStructionsForParamCbks();
}

ParamConvexProb::~ParamConvexProb() {
	// sparsity patterns for ipopt
	delete nzRowNzColToPosValMap;
	free(this->positionsToAddTikhonovTerms);
	// sqp buffer
	free(this->sqpBufferHessX);
	free(this->sqpNzRowZeroBasedEqConstr);
	free(this->sqpNzColZeroBasedEqConstr);
	free(this->sqpNzValEqConstr);
	// buffer
	free(this->lastParameterVecRegBuffer);
	free(this->lastRegPrimalSol);
	free(this->lastRegDualSolEqConstraints);
	free(this->lastRegLowerIneqMul);
	free(this->lastRegUpperIneqMul);
	free(this->lastRegIneqConstraintValues);
	// lin sys structure
	delete this->nzPositionsLinSys;
	delete this->lSysDerSolMap;
}

double ParamConvexProb::getSqpOptVal() {
	return this->sqpOptVal;
}

int ParamConvexProb::getTotalDecisionVars() {
	return this->totalDecisionVars;
}

int ParamConvexProb::getTotalParameters() {
	return this->totalParameters;
}

int ParamConvexProb::getTotalInequalityConstraints() {
	return this->paramIneqConstrVec->size();
}

int ParamConvexProb::getTotalEqualityConstraints() {
	return this->totalEqualityConstraints;
}

ParamCbk* ParamConvexProb::getParamObjFunction() {
	return this->paramObjFunction;
}

MtxParamCbk* ParamConvexProb::getMtxParamCbk() {
	return this->mtxParamCbk;
}

int* ParamConvexProb::getPositionsToAddTikhonovTerms() {
	return this->positionsToAddTikhonovTerms;
}

vector<Cbk*>* ParamConvexProb::getParamRhsVec() {
	return this->paramRhsVec;
}

ParamCbk* ParamConvexProb::getInequalityConstraint(int i) {
	return (*(this->paramIneqConstrVec))[i];
}

map< pair<int,int> , int>* ParamConvexProb::getNzRowNzColToPosValMap() {
	return this->nzRowNzColToPosValMap;
}

double ParamConvexProb::getLastEpsilon() {
	return this->lastEpsilon;
}

double ParamConvexProb::getLastMu() {
	return this->lastMu;
}

double ParamConvexProb::getLastRegOptVal() {
	return this->lastRegOptVal;
}

double* ParamConvexProb::getLastRegPrimalSol() {
	return this->lastRegPrimalSol;
}

double* ParamConvexProb::getLastRegDualSolEqConstraints() {
	return this->lastRegDualSolEqConstraints;
}

double* ParamConvexProb::getLastParameterVecRegBuffer() {
	return this->lastParameterVecRegBuffer;
}

double* ParamConvexProb::getLastDerPrimalSol() {
	return this->lSysDerSolMap->getSolPardiso();
}

void ParamConvexProb::calcWarmStartForIpoptWithCplex() {
	// get common data
	MtxParamCbk* mtxParamCbk = this->getMtxParamCbk();
	double* lastRegPrimalSol = getLastRegPrimalSol();
	double* lastRegParam = getLastParameterVecRegBuffer();
	// get dimensions of quadratic problem
	int totDecVarIter = this->getTotalDecisionVars();
	int totEqConstrIter = this->getTotalEqualityConstraints();
	int totalInequalityConstraints = this->getTotalInequalityConstraints();
	int totalNzMatrixA = mtxParamCbk->getTotalNz();
	// build model
	IloEnv env;
	IloModel model(env);
	// create var
	IloNumVarArray x(env);
	for(int i=0; i<totDecVarIter; ++i) {
		x.add(IloNumVar(env));
	}
	// create objective for cplex solver
	IloExpr objExpr(env);
	// grad contrib of objective
	ParamCbk* objFunc = this->getParamObjFunction();
	vector<int>* nzPosOneBasedGradX = objFunc->getNzPosOneBasedGradX();
	double* bufferObjGradX = objFunc->getBufferGradX();
	MyVector::initVecWithVal(bufferObjGradX, nzPosOneBasedGradX->size(), 0.0);
	objFunc->addSequentiallyGradX( getLastRegPrimalSol(), getLastParameterVecRegBuffer(), bufferObjGradX );
	for (int i=0; i< nzPosOneBasedGradX->size(); ++i) {
		int rowOneBased = (*nzPosOneBasedGradX)[i];
		objExpr += (bufferObjGradX[i]) * (x[ rowOneBased-1 ] - lastRegPrimalSol[ rowOneBased-1 ]);
	}
	// get hess data for contrib of objective
	vector<int>* nzRowsOneBasedHessX = objFunc->getNzRowsOneBasedHessX();
	vector<int>* nzColsOneBasedHessX = objFunc->getNzColsOneBasedHessX();
	objFunc->addSequentiallyScaledHessX( getLastRegPrimalSol(), getLastParameterVecRegBuffer(), sqpBufferHessX, 1.0 );
	for (int i=0; i<nzRowsOneBasedHessX->size(); ++i) {
		// get data
		int row = (*nzRowsOneBasedHessX)[i];
		int col = (*nzColsOneBasedHessX)[i];
		double val = sqpBufferHessX[i];
		// add (this way because only the upper diagonal comes)
		if ( row == col) {
			objExpr += 0.5 * val * (x[row-1] - lastRegPrimalSol[row-1]) * (x[col-1] - lastRegPrimalSol[col-1]);
		} else {
			objExpr += 0.5 * 2.0 * val * (x[row-1] - lastRegPrimalSol[row-1]) * (x[col-1] - lastRegPrimalSol[col-1]);
		}
	}
	// set cplex objective
	IloObjective obj = IloMinimize(env, objExpr);
	model.add(obj);
	// set non-negativity constraints for variables of second stage problems
	IloRangeArray nonNegConstraints(env);
	for (int i=0; i<totDecVarIter; ++i) {
		nonNegConstraints.add(x[i] >= 0);
	}
	model.add(nonNegConstraints);
	// get data for equality constraints
	mtxParamCbk->writeNzRowsAndColsSequentiallyInZeroBasedForm( 
		this->sqpNzRowZeroBasedEqConstr, this->sqpNzColZeroBasedEqConstr );
	mtxParamCbk->writeSequentiallyNzVals(
		this->getLastParameterVecRegBuffer(), this->sqpNzValEqConstr );
	// create equality constraints
	int globalIndexScenarioMatrixA = 0;
	IloRangeArray eqConstraints(env);
	for(int i=0; i<totEqConstrIter; ++i) {
		// compute expression
		IloExpr linExpr(env);
		int counter = 0;
		while(counter < totalNzMatrixA && sqpNzRowZeroBasedEqConstr[globalIndexScenarioMatrixA] == i) {
			int rowIter = sqpNzRowZeroBasedEqConstr[globalIndexScenarioMatrixA];
			int colIter = sqpNzColZeroBasedEqConstr[globalIndexScenarioMatrixA];
			double valIter = sqpNzValEqConstr[globalIndexScenarioMatrixA];
			linExpr += valIter * x[colIter];
			globalIndexScenarioMatrixA++;
			counter++;
		}
		// get rhs val
		double rhsVal = (*paramRhsVec)[i]->getValue( getLastParameterVecRegBuffer() );
		// add constraint
		eqConstraints.add(linExpr == rhsVal);
	}
	// add eq constraints
	model.add(eqConstraints);
	// inequality constraints
	IloRangeArray ineqConstraints(env);
	int totalInequalities = this->getTotalInequalityConstraints();
	for ( int i=0; i<totalInequalities; ++i ) {
		// get inequality
		ParamCbk* ineqParamCbk = getInequalityConstraint(i);
		// create linear expression
		IloExpr expr(env);
		// get value
		double valueIneq = ineqParamCbk->getValue( getLastRegPrimalSol(), getLastParameterVecRegBuffer() );
		expr += valueIneq;
		// gradient part
		vector<int>* ineqNzPosOneBasedGradX = ineqParamCbk->getNzPosOneBasedGradX();
		double* bufferIneqGradX = ineqParamCbk->getBufferGradX();
		MyVector::initVecWithVal(bufferIneqGradX, ineqNzPosOneBasedGradX->size(), 0.0);
		ineqParamCbk->addSequentiallyGradX( getLastRegPrimalSol(), getLastParameterVecRegBuffer(), bufferIneqGradX );
		for ( int i=0; i<ineqNzPosOneBasedGradX->size(); ++i ) {
			int rowOneBased = (*ineqNzPosOneBasedGradX)[i];
			expr += (bufferIneqGradX[i]) * (x[ rowOneBased-1 ] - lastRegPrimalSol[ rowOneBased-1 ]);
		}
		// quadratic part
		vector<int>* ineqNzRowsOneBasedHessX = ineqParamCbk->getNzRowsOneBasedHessX();
		vector<int>* ineqNzColsOneBasedHessX = ineqParamCbk->getNzColsOneBasedHessX();
		ineqParamCbk->addSequentiallyScaledHessX( getLastRegPrimalSol(), getLastParameterVecRegBuffer(), sqpBufferHessX, 1.0 );
		for (int i=0; i<ineqNzRowsOneBasedHessX->size(); ++i) {
			// get data
			int row = (*ineqNzRowsOneBasedHessX)[i];
			int col = (*ineqNzColsOneBasedHessX)[i];
			double val = sqpBufferHessX[i];
			// add (this way because only the upper diagonal comes)
			if ( row == col) {
				expr += 0.5 * val * (x[row-1] - lastRegPrimalSol[row-1]) * (x[col-1] - lastRegPrimalSol[col-1]);
			} else {
				expr += 0.5 * 2.0 * val * (x[row-1] - lastRegPrimalSol[row-1]) * (x[col-1] - lastRegPrimalSol[col-1]);
			}
		}
		// add constraint
		ineqConstraints.add(expr <= 0.0);
	}
	// add ineq constraints
	model.add(ineqConstraints);
	// solve scenario lp
	IloCplex cplex(model);
	// set parameters for the quadratic case
	if ( this->lastMu > 0.0 || this->nzRowNzColToPosValMap->size() > 0) {
		cplex.setParam(IloCplex::RootAlg, IloCplex::Param::Barrier::Algorithm);
		cplex.setParam(IloCplex::Param::Barrier::ConvergeTol, 0.0000000001);
		cplex.setParam(IloCplex::BarCrossAlg, -1);
	}
	// normal parameters
	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::Param::Tune::TimeLimit, 100.0);
	// solve
	IloBool isFeasible = cplex.solve();
	// get primal solution
	for(int i=0; i<totDecVarIter; ++i) {
		this->lastRegPrimalSol[i] = cplex.getValue(x[i]);
	}
	// get dual solution eq constraints
	for(int i=0; i<totEqConstrIter; ++i) {
		this->lastRegDualSolEqConstraints[i] = cplex.getDual(eqConstraints[i]);
		this->lastRegDualSolEqConstraints[i] *= -1.0;
	}
	// get dual solution non neg constraints
	for(int i=0; i<totDecVarIter; ++i) {
		this->lastRegLowerIneqMul[i] = cplex.getDual(nonNegConstraints[i]);
		this->lastRegUpperIneqMul[i] = 1e29;
	}
	// get opt value
	this->sqpOptVal = this->paramObjFunction->getValue(
		getLastRegPrimalSol(), getLastParameterVecRegBuffer() );
	// finish
	env.end();
}

Bool ParamConvexProb::intermediate_callback_param_convex_prob(
		Index alg_mod, Index iter_count, Number obj_value, Number inf_pr, 
		Number inf_du, Number mu, Number d_norm, Number regularization_size, 
		Number alpha_du, Number alpha_pr, Index ls_trials, UserDataPtr user_data_all) {
	// get data
	ParamConvexProb* paramConvexProb = (ParamConvexProb*) user_data_all;
	
	// return
	return TRUE;
}

Bool ParamConvexProb::eval_f_param_convex_prob(Index n, Number* x, Bool new_x, Number* obj_value, UserDataPtr user_data_all) {
	// get data
	ParamConvexProb* paramConvexProb = (ParamConvexProb*) user_data_all;
	ParamCbk* paramObjFunction = paramConvexProb->getParamObjFunction();
	double lastEpsilon = paramConvexProb->getLastEpsilon();
	double lastMu = paramConvexProb->getLastMu();
	double* lastParamVec = paramConvexProb->getLastParameterVecRegBuffer();
	// init
	double total = 0.0;
	// normal value
	total += paramObjFunction->getValue(x, lastParamVec);
	// tikhonov term
	double normSquared = MyVector::computeSquaredNorm(x, n);
	total += 0.5 * lastEpsilon * lastMu * normSquared;
	// set
	*obj_value = total;
	// return
	return TRUE;
}

Bool ParamConvexProb::eval_grad_f_param_convex_prob(Index n, Number* x, Bool new_x, Number* grad_f, UserDataPtr user_data_all) {
	// get data
	ParamConvexProb* paramConvexProb = (ParamConvexProb*) user_data_all;
	ParamCbk* paramObjFunction = paramConvexProb->getParamObjFunction();
	double lastEpsilon = paramConvexProb->getLastEpsilon();
	double lastMu = paramConvexProb->getLastMu();
	double* lastParamVec = paramConvexProb->getLastParameterVecRegBuffer();
	// init
	MyVector::initVecWithVal(grad_f, n, 0.0);
	// normal grad
	paramObjFunction->addWithJumpsGradX(x, lastParamVec, grad_f);
	// tikhonov part
	for (int i=0; i<n; ++i) {
		grad_f[i] += lastEpsilon * lastMu * (x[i]);
	}
	// return
	return TRUE;
}

Bool ParamConvexProb::eval_g_param_convex_prob(Index n, Number* x, Bool new_x, Index m, Number* g, UserDataPtr user_data_all) {
	// get data
	ParamConvexProb* paramConvexProb = (ParamConvexProb*) user_data_all;
	ParamCbk* paramObjFunction = paramConvexProb->getParamObjFunction();
	int totalEqualityConstraints = paramConvexProb->getTotalEqualityConstraints();
	double lastEpsilon = paramConvexProb->getLastEpsilon();
	double lastMu = paramConvexProb->getLastMu();
	double* lastParamVec = paramConvexProb->getLastParameterVecRegBuffer();
	// init with - rhs
	vector<Cbk*>* paramRhsVec = paramConvexProb->getParamRhsVec();
	for (int i=0; i<totalEqualityConstraints; ++i) {
		Cbk* rhsCbkIter = (*paramRhsVec)[i];
		double value = rhsCbkIter->getValue( lastParamVec );
		g[i] = - value;
	}
	// linear equality constraints
	MtxParamCbk* mtxParamCbk = paramConvexProb->getMtxParamCbk();
	mtxParamCbk->addTimesVec(x, lastParamVec, g);
	// position
	int posToStore = totalEqualityConstraints;
	// inequality constraints
	int totalInequalityConstraints = paramConvexProb->getTotalInequalityConstraints();
	for (int i=0; i<totalInequalityConstraints; ++i) {
		ParamCbk* paramCbk = paramConvexProb->getInequalityConstraint(i);
		double value = paramCbk->getValue(x, lastParamVec);
		g[posToStore] = value;
		posToStore += 1;
	}
	// return
	return TRUE;
}

Bool ParamConvexProb::eval_jac_g_param_convex_prob(Index n, Number *x, Bool new_x, Index m, Index nele_jac, 
		Index *iRow, Index *jCol, Number *values, UserDataPtr user_data_all) {
	// get data
	ParamConvexProb* paramConvexProb = (ParamConvexProb*) user_data_all;
	int totalEqualityConstraints = paramConvexProb->getTotalEqualityConstraints();
	// sparsity pattern
	if (values == NULL) {
		// equality constraints
		paramConvexProb->getMtxParamCbk()->writeNzRowsAndColsSequentiallyInZeroBasedForm(iRow, jCol);
		// init
		int currentPos = paramConvexProb->getMtxParamCbk()->getTotalNz();
		// iterate inequalities
		int totalInequalityConstraints = paramConvexProb->getTotalInequalityConstraints();
		for (int j=0; j<totalInequalityConstraints; ++j) {
			ParamCbk* paramCbkIneq = paramConvexProb->getInequalityConstraint(j);
			vector<int>* nzPosOneBasedGradX = paramCbkIneq->getNzPosOneBasedGradX();
			int totNzIneq = nzPosOneBasedGradX->size();
			for (int i=0; i<totNzIneq; ++i) {
				iRow[currentPos] = totalEqualityConstraints + j;
				jCol[currentPos] = (*nzPosOneBasedGradX)[i] - 1;
				currentPos += 1;
			}
		}
	// values
	} else {
		// equality constraints
		paramConvexProb->getMtxParamCbk()->writeSequentiallyNzVals(
			paramConvexProb->getLastParameterVecRegBuffer(), values);
		// init
		int currentPos = paramConvexProb->getMtxParamCbk()->getTotalNz();
		// iterate inequalities
		int totalInequalityConstraints = paramConvexProb->getTotalInequalityConstraints();
		for (int j=0; j<totalInequalityConstraints; ++j) {
			ParamCbk* paramCbkIneq = paramConvexProb->getInequalityConstraint(j);
			int totNzIneq = paramCbkIneq->getTotalNzGradX();
			double* valuesIneq = &(values[currentPos]);
			MyVector::initVecWithVal(valuesIneq, totNzIneq, 0.0);
			paramCbkIneq->addSequentiallyGradX(x, paramConvexProb->getLastParameterVecRegBuffer(), valuesIneq);
			currentPos += totNzIneq;
		}
	}
	// return
	return TRUE;
}

Bool ParamConvexProb::eval_h_param_convex_prob(Index n, Number *x, Bool new_x, Number obj_factor, Index m, 
		Number *lambda, Bool new_lambda, Index nele_hess, 
		Index *iRow, Index *jCol, Number *values, UserDataPtr user_data_all) {
	// get data
	ParamConvexProb* paramConvexProb = (ParamConvexProb*) user_data_all;
	double lastEpsilon = paramConvexProb->getLastEpsilon();
	double lastMu = paramConvexProb->getLastMu();
	// sparsity pattern
	if (values == NULL) {
		map< pair<int,int> , int> *nzRowNzColToPosValMap = paramConvexProb->getNzRowNzColToPosValMap();
		for(map< pair<int,int> ,int>::iterator it = nzRowNzColToPosValMap->begin(); it != nzRowNzColToPosValMap->end(); ++it) {
			// get data
			pair<int, int> pIter = it->first;
			int posZeroBased = it->second;
			// split pair
			int rowOneBased = get<0>(pIter);
			int colOneBased = get<1>(pIter);
			// assign
			iRow[posZeroBased] = rowOneBased - 1;
			jCol[posZeroBased] = colOneBased - 1;
		}
	// values
	} else {
		// dict to fill values
		int totalEqualityConstraints = paramConvexProb->getTotalEqualityConstraints();
		map< pair<int,int> , int> *nzRowNzColToPosValMap = paramConvexProb->getNzRowNzColToPosValMap();
		// init with zeros
		MyVector::initVecWithVal(values, nele_hess, 0.0);
		// fill obj data
		ParamCbk* paramObjFunction = paramConvexProb->getParamObjFunction();
		paramObjFunction->addWithJumpsScaledHessX(x,  paramConvexProb->getLastParameterVecRegBuffer(), 
			values, obj_factor);
		// iterate inequalities
		int totalInequalityConstraints = paramConvexProb->getTotalInequalityConstraints();
		for (int i=0; i<totalInequalityConstraints; ++i) {
			ParamCbk* paramCbkIneq = paramConvexProb->getInequalityConstraint(i);
			paramCbkIneq->addWithJumpsScaledHessX(x, paramConvexProb->getLastParameterVecRegBuffer(), 
				values, lambda[totalEqualityConstraints + i]);
		}
		// tikhonov part
		int* positionsToAddTikhonovTerms = paramConvexProb->getPositionsToAddTikhonovTerms();
		for (int i=0; i<n; ++i) {
			int posToAdd = positionsToAddTikhonovTerms[i];
			values[posToAdd] += obj_factor * lastEpsilon * lastMu;
		}
	}
	// return
	return TRUE;	
}

int ParamConvexProb::getTotalNzForJacobianOfAllConstraints() {
	// init
	int totalNzJacobian = 0;
	// equality constraints
	totalNzJacobian += this->mtxParamCbk->getTotalNz();
	// inequality constraints
	int totalInequalityConstraints = this->paramIneqConstrVec->size();
	for (int i=0; i<totalInequalityConstraints; ++i) {
		ParamCbk* paramCbk = (*(this->paramIneqConstrVec))[i];
		totalNzJacobian += paramCbk->getTotalNzGradX();
	}
	// return
	return totalNzJacobian;
}

void ParamConvexProb::fillNzFoHessianOfLagrangian() {
	// init
	int posToStore = 0;
	// for the tikhonov perturbation or the inequality constraints
	int totalDecisionVars = this->getTotalDecisionVars();
	for (int i=0; i<totalDecisionVars; ++i) {
		(*nzRowNzColToPosValMap)[make_pair(i+1, i+1)] = posToStore;
		posToStore += 1;
	} 
	// objective function
	ParamCbk* paramObjFunction = this->getParamObjFunction();
	vector<int>* objNzRowsOneBasedHessX = paramObjFunction->getNzRowsOneBasedHessX();
	vector<int>* objNzColsOneBasedHessX = paramObjFunction->getNzColsOneBasedHessX();
	int totalNzHessObj = objNzRowsOneBasedHessX->size();
	for (int i=0; i<totalNzHessObj; ++i) {
		// get data
		int row = (*objNzRowsOneBasedHessX)[i];
		int col = (*objNzColsOneBasedHessX)[i];
		// check if value is not present
		if (nzRowNzColToPosValMap->find(make_pair(row, col)) == nzRowNzColToPosValMap->end()) {
			(*nzRowNzColToPosValMap)[make_pair(row, col)] = posToStore;
			posToStore += 1;
		}
	}
	// inequality constraints
	int totalInequalityConstraints = this->paramIneqConstrVec->size();
	for (int j=0; j<totalInequalityConstraints; ++j) {
		// get ineq data
		ParamCbk* paramCbkIneq = (*(this->paramIneqConstrVec))[j];
		vector<int>* ineqNzRowsOneBasedHessX = paramCbkIneq->getNzRowsOneBasedHessX();
		vector<int>* ineqNzColsOneBasedHessX = paramCbkIneq->getNzColsOneBasedHessX();
		int totalNzHessIneq = ineqNzRowsOneBasedHessX->size();
		// iterate
		for (int i=0; i<totalNzHessIneq; ++i) {
			// get data
			int row = (*ineqNzRowsOneBasedHessX)[i];
			int col = (*ineqNzColsOneBasedHessX)[i];
			// check if value is not present
			if (nzRowNzColToPosValMap->find(make_pair(row, col)) == nzRowNzColToPosValMap->end()) {
				(*nzRowNzColToPosValMap)[make_pair(row, col)] = posToStore;
				posToStore += 1;
			}
		}
	}
	// fill positions to add tikhonov terms
	for (int i=0; i<totalDecisionVars; ++i) {
		pair<int, int> pIter = make_pair(i+1, i+1);
		int posToAdd = (*nzRowNzColToPosValMap)[pIter];
		positionsToAddTikhonovTerms[i] = posToAdd;
	}
}

void ParamConvexProb::calcNonTrivialSparsityPatternForLinSys() {
	// init undordered set
	set< pair<int, int> > setTmp;
	// part of the pattern repeats
	for (map< pair<int, int>, int>::iterator it = nzRowNzColToPosValMap->begin(); it != nzRowNzColToPosValMap->end(); ++it) {
		// put element into vector
		setTmp.insert(it->first);
	}
	// diagonal of multiplier part should be there even if zero
	for (int i=0; i<totalEqualityConstraints; ++i) {
		// idx multiplier
		int idxMulOneBased = totalDecisionVars + i + 1;
		// put element into vector
		setTmp.insert(make_pair(idxMulOneBased, idxMulOneBased));
	}
	// iterate inequalities
	int totalInequalityConstraints = this->paramIneqConstrVec->size();
	for (int j=0; j<totalInequalityConstraints; ++j) {
		// get ineq data
		ParamCbk* paramCbkIneq = (*(this->paramIneqConstrVec))[j];
		vector<int> *nzPosOneBasedGradX = paramCbkIneq->getNzPosOneBasedGradX();
		int totNzGradX = nzPosOneBasedGradX->size();
		// NOTE: need to add the upper part of the non zero positions of gradx * (gradx')
		for (int k=0; k<totNzGradX; ++k) {
			for (int l=k; l<totNzGradX; ++l) {
				// get data
				int rowOneBased = (*nzPosOneBasedGradX)[k];
				int colOneBased = (*nzPosOneBasedGradX)[l];
				// put element into vector
				pair<int, int> pIter = make_pair(rowOneBased, colOneBased);
				setTmp.insert(pIter);
			}
		}
	}
	// get equality constraints data (what enters is the transpose)
	int totNzEqualityConstraints = this->mtxParamCbk->getTotalNz();
	vector<int>* nzRowsEqConstr = this->mtxParamCbk->getFinalNzRows();
	vector<int>* nzColsEqConstr = this->mtxParamCbk->getFinalNzCols();
	// equality constraints
	for (int i=0; i<totNzEqualityConstraints; ++i) {
		// get data ( *the transpose* )
		int rowOneBased = (*nzColsEqConstr)[i];
		int colOneBased = (*nzRowsEqConstr)[i] + totalDecisionVars;
		// put element into vector
		pair<int, int> pIter = make_pair(rowOneBased, colOneBased);
		setTmp.insert(pIter);
	}
	// put set to vector
	for ( set< pair<int, int> >::iterator it = setTmp.begin() ; it != setTmp.end() ; ++it) {
		this->nzPositionsLinSys->push_back( *it );
	}
	// sort lexicographically
	sort( this->nzPositionsLinSys->begin(), this->nzPositionsLinSys->end() );
}

void ParamConvexProb::setFillingStructionsForParamCbks() {
	// get obj
	ParamCbk* objFunc = this->getParamObjFunction();
	MtxParamCbk* mtxParamCbk = this->getMtxParamCbk();
	// using the global reference set the sparsity patterns of the parametric callbacks
	// NOTE: these fillings intructions is for Ipopt
	objFunc->fillFillingInstructionsForHessWithJumps(nzRowNzColToPosValMap);
	int totalInequalities = this->getTotalInequalityConstraints();
	for (int i=0; i<totalInequalities; ++i) {
		ParamCbk* ineqParamCbk = getInequalityConstraint(i);
		ineqParamCbk->fillFillingInstructionsForHessWithJumps(nzRowNzColToPosValMap);
	}
	// create map for lin sys filling positions
	// NOTE: assume nzPositionsLinSys is lexicographically sorted
	map < pair<int, int>, int> dictOfPositionsForLinSys;
	int currentPosition = 0;
	for ( vector< pair<int, int>>::iterator it = nzPositionsLinSys->begin(); it != nzPositionsLinSys->end(); ++it) {
		dictOfPositionsForLinSys[*it] = currentPosition;
		currentPosition += 1;
	}
	// using the global reference set the sparsity patterns of the parametric callbacks
	// NOTE: these fillings intructions are the the linear system
	mtxParamCbk->fillFillingInstructionsForMtxCbkWithJumpsLinSys(this->totalDecisionVars, &dictOfPositionsForLinSys);
	objFunc->fillFillingInstructionsForHessWithJumpsLinSys(&dictOfPositionsForLinSys);
	for (int i=0; i<totalInequalities; ++i) {
		ParamCbk* ineqParamCbk = getInequalityConstraint(i);
		ineqParamCbk->fillFillingInstructionsForHessWithJumpsLinSys(&dictOfPositionsForLinSys);
		ineqParamCbk->fillFillingInstructionsForGradXGradXTransposeWithJumpsLinSys(&dictOfPositionsForLinSys);
	}
}

int ParamConvexProb::getTotalNzFoHessianOfLagrangian() {
	return this->nzRowNzColToPosValMap->size();
}

void ParamConvexProb::recalcRegBuffer() {
	// warm start
	if ( this->WARM_START_IPOPT_STEP == true ) {
		calcWarmStartForIpoptWithCplex();
	}
	// number of variables 
	Index n = this->totalDecisionVars;
	 // number of constraints 
	Index m = this->totalEqualityConstraints + this->getTotalInequalityConstraints();
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
	Index nele_jac = this->totalNzJacobianAllConstraints;
	// Number of nonzeros in the Hessian of the Lagrangian (lower or upper triangual part only) 
	Index nele_hess = this->getTotalNzFoHessianOfLagrangian();
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
	int mToAlloc = m;
	if (m == 0) mToAlloc = 1;
	g_L = (Number*)malloc(sizeof(Number)*mToAlloc);
	g_U = (Number*)malloc(sizeof(Number)*mToAlloc);
	// set equality constraints
	for ( int i=0; i<(this->totalEqualityConstraints); ++i ) {
		g_L[i] = 0.0;
		g_U[i] = 0.0;
	}
	// set bounds for inequality constranints
	for ( int i=0; i<(this->getTotalInequalityConstraints()); ++i ) {
		g_L[this->totalEqualityConstraints + i] = -2e19;
		g_U[this->totalEqualityConstraints + i] = 0.0;
	}
	// create the IpoptProblem
	nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess,
		index_style, &(ParamConvexProb::eval_f_param_convex_prob), 
		&(ParamConvexProb::eval_g_param_convex_prob), 
		&(ParamConvexProb::eval_grad_f_param_convex_prob), 
		&(ParamConvexProb::eval_jac_g_param_convex_prob),
		&(ParamConvexProb::eval_h_param_convex_prob));
	// We can free the memory now - the values for the bounds have been
	// copied internally in CreateIpoptProblem
	free(x_L);
	free(x_U);
	free(g_L);
	free(g_U);
	// Set some options.  Note the following ones are only examples,
	// they might not be suitable for your problem.
	AddIpoptIntOption(nlp, "print_level", 0);
	AddIpoptIntOption(nlp, "max_iter", 200);
	AddIpoptNumOption(nlp, "mu_target", this->lastEpsilon);
	AddIpoptStrOption(nlp, "warm_start_init_point", "yes");
	AddIpoptStrOption(nlp, "hessian_approximation", "exact");
	AddIpoptNumOption(nlp, "tol", 1e-9);
	// intermediate callback
	SetIntermediateCallback(nlp, &(ParamConvexProb::intermediate_callback_param_convex_prob));
	// solve the problem
	status = IpoptSolve(nlp, this->lastRegPrimalSol, NULL, &(this->lastRegOptVal), 
		this->lastRegDualSolEqConstraints, this->lastRegLowerIneqMul, 
		this->lastRegUpperIneqMul, this);
	// free
	FreeIpoptProblem(nlp);
	// fill reg buffer for values of inequalities
	int totalInequalityConstraints = this->getTotalInequalityConstraints();
	for (int i=0; i<totalInequalityConstraints; ++i) {
		ParamCbk* paramCbk = this->getInequalityConstraint(i);
		lastRegIneqConstraintValues[i] = paramCbk->getValue(
			this->lastRegPrimalSol, this->getLastParameterVecRegBuffer());
	}
	// truncate primal solutions if needed
	for (int i=0; i<totalDecisionVars; ++i) {
		if ( lastRegPrimalSol[i] <= SAFE_GUARD_INTERIOR_COND_TRUNCATION) {
			lastRegPrimalSol[i] = SAFE_GUARD_INTERIOR_COND_TRUNCATION;
		}
	}
	// truncate inequality values if needed
	for (int i=0; i<totalInequalityConstraints; ++i) {
		if ( lastRegIneqConstraintValues[i] >= -SAFE_GUARD_INTERIOR_COND_TRUNCATION ) {
			lastRegIneqConstraintValues[i] = -SAFE_GUARD_INTERIOR_COND_TRUNCATION;
		}
	}
	// change sign of multiplier of equality constraints
//	for (int i=0; i<(this->totalEqualityConstraints); ++i ) {
//		this->lastRegDualSolEqConstraints[i] *= -1.0;
//	}
}

void ParamConvexProb::updateRegBuffer(double epsilon, double mu, double* p) {
	// check if should update buffer
	double distParameters = MyVector::getNormOfDiffVec(p, this->lastParameterVecRegBuffer, this->totalParameters);
	if (abs(epsilon - this->lastEpsilon) > BUFFER_TOL || abs(mu - this->lastMu) > BUFFER_TOL || distParameters > BUFFER_TOL) {
		// update targets for the buffer
		this->lastEpsilon = epsilon;
		this->lastMu = mu;
		MyVector::copyVecFromTo(p, this->lastParameterVecRegBuffer, this->totalParameters);
		// calc again
		this->recalcRegBuffer();
	}	
}

void ParamConvexProb::buildSysForDerOfSolutionMappings() {
	// get data
	int n_ = lSysDerSolMap->getTotalEqLinearSystem();
	int totNz = lSysDerSolMap->getTotalNzEqLinearSystem();
	int* ia_ = lSysDerSolMap->getIaPardiso();
	int* ja_ = lSysDerSolMap->getJaPardiso();
	double* a_ = lSysDerSolMap->getValsPardiso();
	bool symmetricReal_ = true;
	// init
	MyVector::initVecWithVal(a_, totNz, 0.0);
	// fill equality constraints
	MtxParamCbk* mtxParamCbk = this->getMtxParamCbk();
	mtxParamCbk->addWithJumpsNzValsLinSys(
		getLastParameterVecRegBuffer(), a_);
	// objective part
	ParamCbk* objFunc = this->getParamObjFunction();
	objFunc->addWithJumpsLinSysScaledHessX(
		getLastRegPrimalSol(), getLastParameterVecRegBuffer(), a_, 1.0);
	// iterate inequalities (hessian part)
	int totalInequalityConstraints = this->paramIneqConstrVec->size();
	for (int j=0; j<totalInequalityConstraints; ++j) {
		// get ineq data
		ParamCbk* paramCbkIneq = (*(this->paramIneqConstrVec))[j];
		// calc scaling
		double scaling = - lastEpsilon / lastRegIneqConstraintValues[j];
		// add contrib
		paramCbkIneq->addWithJumpsLinSysScaledHessX(
			getLastRegPrimalSol(), getLastParameterVecRegBuffer(), a_, scaling);
	}
	// iterate inequalities ( gradx * (gradx') part)
	for (int j=0; j<totalInequalityConstraints; ++j) {
		// get ineq data
		ParamCbk* paramCbkIneq = (*(this->paramIneqConstrVec))[j];
		// calc scaling
		double gSquared = ( lastRegIneqConstraintValues[j] ) * ( lastRegIneqConstraintValues[j] );
		double scaling = lastEpsilon / gSquared;
		// add contrib
		paramCbkIneq->addScaledWithJumpsLinSysGradXGradXTranspose(
			getLastRegPrimalSol(), getLastParameterVecRegBuffer(), a_, scaling);
	}
	// add eps * tikhonov and eps * X^(-2)
	double* lastRegPrimalSol = this->getLastRegPrimalSol();
	int totalDecisionVars = this->getTotalDecisionVars();
	for ( int i=0; i<totalDecisionVars; ++i ) {
		int beginLineI = ia_[i] - 1;
		a_[ beginLineI ] += (this->lastEpsilon) * pow(lastRegPrimalSol[i], -2.0);
		a_[ beginLineI ] += (this->lastEpsilon) * (this->lastMu);
	}
	// create new solver instance (accounts for the factorization)
	this->lSysSolver = new LinearSystemSolver(n_, ia_, ja_, a_, symmetricReal_);
}

void ParamConvexProb::destroySysForDerOfSolutionMappings() {
	// delete
	if (this->lSysSolver != NULL) {
		delete this->lSysSolver;
	}
	// set
	this->lSysSolver = NULL;
}

void ParamConvexProb::calcRhsForDerSolutionMappingsCalculation(int jOneBased) {
	// get data
	int sizeLinSys = this->lSysDerSolMap->getTotalEqLinearSystem();
	double* rhs = this->lSysDerSolMap->getRhsPardiso();
	// init
	MyVector::initVecWithVal(rhs, sizeLinSys, 0.0);
	// theta_j, part 1
	ParamCbk* objFunc = this->getParamObjFunction();
	objFunc->addScaledWithJumpsPartialPGradX( jOneBased, getLastRegPrimalSol(), getLastParameterVecRegBuffer(), rhs, -1.0 );
	// theta_j, part 2
	MtxParamCbk* mtxParamCbk = this->getMtxParamCbk();
	mtxParamCbk->addScaledTransposePartialPTimesVec(jOneBased,
		this->getLastRegDualSolEqConstraints(), getLastParameterVecRegBuffer(), rhs, 1.0);
	// eps * varphi_j, part 1
	int totalInequalityConstraints = this->paramIneqConstrVec->size();
	for (int j=0; j<totalInequalityConstraints; ++j) {
		// get ineq data
		ParamCbk* paramCbkIneq = (*(this->paramIneqConstrVec))[j];
		// calc scaling
		double scaling = lastEpsilon / lastRegIneqConstraintValues[j];
		// add partial p grad x
		paramCbkIneq->addScaledWithJumpsPartialPGradX(jOneBased,
			getLastRegPrimalSol(), getLastParameterVecRegBuffer(), rhs, scaling);
	}
	// eps * varphi_j, part 2
	for (int j=0; j<totalInequalityConstraints; ++j) {
		// get ineq data
		ParamCbk* paramCbkIneq = (*(this->paramIneqConstrVec))[j];
		// calc scaling to the gradient
		double gSquared = ( lastRegIneqConstraintValues[j] ) * ( lastRegIneqConstraintValues[j] );
		double partialP = paramCbkIneq->getPartialP( jOneBased, getLastRegPrimalSol(), getLastParameterVecRegBuffer() );
		double scaling = - ( partialP * lastEpsilon ) / gSquared;
		// add scaled gradient
		paramCbkIneq->addScaledWithJumpsGradX(
			getLastRegPrimalSol(), getLastParameterVecRegBuffer(), rhs, scaling);
	}
	// space for the second part
	int totEqCosntraints = this->getTotalEqualityConstraints();
	if ( totEqCosntraints > 0 ) {
		double* shiftedRhs = &(rhs[ this->totalDecisionVars ]);
		// phi_j, part 1
		for (int i=0; i<totEqCosntraints; ++i) {
			Cbk* rhsIter = (*paramRhsVec)[i];
			shiftedRhs[i] = rhsIter->getPartial(
				jOneBased, getLastParameterVecRegBuffer());
		}
		// phi_j, part 2
		mtxParamCbk->addScaledPartialPTimesVec(jOneBased, 
			getLastRegPrimalSol(), getLastParameterVecRegBuffer(), shiftedRhs, -1.0);
	}
}

void ParamConvexProb::calcPartialDerSolMap(int jOneBased) {
	// set rhs
	this->calcRhsForDerSolutionMappingsCalculation(jOneBased);
	// solve linear system
	// NOTE: assume factorization is already computed
	this->lSysSolver->solveSystemWithPardiso( lSysDerSolMap->getRhsPardiso(), lSysDerSolMap->getSolPardiso() );
}
