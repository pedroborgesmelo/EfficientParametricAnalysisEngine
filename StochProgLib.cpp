/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#include <iostream>

#include "StochProgLib.h"
#include "StochProg.h"
#include "ParamConvexProb.h"
#include "Cbk.h"
#include "LinearCbk.h"
#include "OptValSmoCbk.h"
#include "ParamCbk.h"
#include "QuadParamCbk.h"
#include "MtxParamCbk.h"

int solveNonDiffProbUpperSmoothingWithQuadObjQuadConstraints(
		int maxIterPerSubProblem,
		double epsilon,
		double tikhonovMu,
		int totalDecisionVarsFirstStage,
		int totalScenarios,
		int* totalDecisionVarsPerScenario,
		int totalEqConsFirstStage,
		double* rhsFirstStage,
		int* totalEqConstrPerScenario,
		double* costFirstStage,
		double* probVecScenarios,
		double* costPerScenario,
		double* rhsPerScenario,
		int totalNonZeroFirstStageW,
		int* nzRowsFirstStageW,
		int* nzColsFirstStageW,
		double* nzValsFirstStageW,
		int* totalNonZeroPerScenarioW,
		int* nzRowsPerScenarioW,
		int* nzColsPerScenarioW,
		double* nzValsPerScenarioW,
		int* totalNonZerPerScenarioT,
		int* nzRowsPerScenarioT,
		int* nzColsPerScenarioT,
		double* nzValsPerScenarioT,
		int* totalNzSymQuadPerturbationObjPerScenario,
		int* nzRowsSymQuadPerturbationObjPerScenario,
		int* nzColsSymQuadPerturbationObjPerScenario,
		double* nzValsSymQuadPerturbationObjPerScenario,
		int* totalIneqPerScenario,
		double* linearPartIneqPerScenario,
		double* freeTermIneqPerScenario,
		int* totNzSymQuadPartIneqPerScenario,
		int* nzRowsSymQuadPartIneqPerScenario,
		int* nzColsSymQuadPartIneqPerScenario,
		double* nzValsSymQuadPartIneqPerScenario,
		double* linearCostParamPerIneqPerScenario,
		double* optSolOut,
		double* optValOut,
		double* optMulEqConsOut,
		double* optMutXvarLowerConstrOut,
		double* optMutXvarUpperConstrOut,
		char* pathToWrtObjIterates) {

	// to free at the end
	vector<Cbk*> cbkToFree;
	vector<LinearCbk*> linearCbkToFree;
	vector<MtxParamCbk*> mtxParamCbkVec;
	vector<ParamCbk*> paramCbkToFree;
	vector<QuadParamCbk*> quadParamCbkToFree;
	vector<vector<ParamCbk*>*> vectorsParamCbkToFree;
	vector<vector<Cbk*>*> vectorsCbkToFree;
	
	// init param convex prob vec
	vector<ParamConvexProb*> paramConvexProbVec;
	
	// create parametric problems
	int globalCounterMatrixTperScenario = 0;
	int globalIndexRhsPerScenario = 0;
	int globalIndexFixedPart = 0;
	int globalIndexQuadraticPerturbation = 0;
	int globalIndexCostPerScenario = 0;
	int globalIndexIneqQuadScen = 0;
	int globalIndexFreeTermIneq = 0;
	int globalIndexLinearPartIneq = 0;
	int globalIndexLinearCostParamIneq = 0;
	for (int s=0; s<totalScenarios; ++s) {
		
		// objective function
		int totNzSymQuadPart_ = totalNzSymQuadPerturbationObjPerScenario[s];
		int* nzRowsOneBasedSymQuadPart_ = &(nzRowsSymQuadPerturbationObjPerScenario[globalIndexQuadraticPerturbation]);
		int* nzColsOneBasedSymQuadPart_ = &(nzColsSymQuadPerturbationObjPerScenario[globalIndexQuadraticPerturbation]);
		double* nzValsSymQuadPart_ = &(nzValsSymQuadPerturbationObjPerScenario[globalIndexQuadraticPerturbation]);
		double* linearCostX_ = &(costPerScenario[globalIndexCostPerScenario]);
		QuadParamCbk* paramObjFunction = new QuadParamCbk(
			totalDecisionVarsPerScenario[s],
			0,
			totNzSymQuadPart_,
			nzRowsOneBasedSymQuadPart_,
			nzColsOneBasedSymQuadPart_,
			nzValsSymQuadPart_,
			linearCostX_,
			0.0,
			NULL);
		quadParamCbkToFree.push_back(paramObjFunction);
		
		// parametric equality constraints
		int totNzFixedPart_ = totalNonZeroPerScenarioW[s];
		int* nzRowsFixedPart_ = &(nzRowsPerScenarioW[globalIndexFixedPart]);
		int* nzColsFixedPart_ = &(nzColsPerScenarioW[globalIndexFixedPart]);
		double* nzValsFixedPart_= &(nzValsPerScenarioW[globalIndexFixedPart]);
		MtxParamCbk* mtxParamCbk = new MtxParamCbk(
			totNzFixedPart_, nzRowsFixedPart_, 
			nzColsFixedPart_, nzValsFixedPart_);
		
		// param rhs vec
		int totEqConstrScen = totalEqConstrPerScenario[s];
		vector<Cbk*> *paramRhsVec = new vector<Cbk*>();
		vectorsCbkToFree.push_back(paramRhsVec);
		int counterMtxTscen = 0;
		int maxMtxT = totalNonZerPerScenarioT[s];

		for (int i=0; i<totEqConstrScen; ++i) {

			// calc total nz on the line
			
			int totNzLine_ = 0;
			while (counterMtxTscen < maxMtxT && nzRowsPerScenarioT[counterMtxTscen] == i+1) {
				counterMtxTscen += 1;
				totNzLine_ += 1;
			}

			// get data

			int totalVariables_ = totalDecisionVarsFirstStage;
			int* nzRowsOneBasedLinearCost_ = &(nzColsPerScenarioT[globalCounterMatrixTperScenario]);
			double* nzValsLinearCost_ = &(nzValsPerScenarioT[globalCounterMatrixTperScenario]);
			double freeTerm_ = rhsPerScenario[globalIndexRhsPerScenario];
			
			// create cbk
			
			LinearCbk* rhsIter = new LinearCbk(
				totalVariables_, totNzLine_,
				nzRowsOneBasedLinearCost_,
				nzValsLinearCost_,
				freeTerm_, -1.0);

			// add
			
			paramRhsVec->push_back( rhsIter );
			linearCbkToFree.push_back( rhsIter );
			
			// controls
			
			globalIndexRhsPerScenario += 1;
			globalCounterMatrixTperScenario += totNzLine_;
		}

		// param ineq constraints
		int totalInequalitiesScen = totalIneqPerScenario[s];
		vector<ParamCbk*> *paramIneqConstrVec = new vector<ParamCbk*>();
		vectorsParamCbkToFree.push_back(paramIneqConstrVec);

		for ( int i=0; i<totalInequalitiesScen ; ++i ) {
			
			// get data
			
			int totNzSymQuadPart_ = totNzSymQuadPartIneqPerScenario[globalIndexFreeTermIneq];
			int* nzRowsOneBasedSymQuadPart_ = &(nzRowsSymQuadPartIneqPerScenario[globalIndexIneqQuadScen]);
			int* nzColsOneBasedSymQuadPart_ = &(nzColsSymQuadPartIneqPerScenario[globalIndexIneqQuadScen]);
			double* nzValsSymQuadPart_ = &(nzValsSymQuadPartIneqPerScenario[globalIndexIneqQuadScen]);
			double* linearCostX_ = &(linearPartIneqPerScenario[globalIndexLinearPartIneq]);
			double freeTerm_ = freeTermIneqPerScenario[globalIndexFreeTermIneq];
			double* linearCostP_ = &(linearCostParamPerIneqPerScenario[globalIndexLinearCostParamIneq]);

			// quad param cbk
			
			QuadParamCbk* ineqQuadParamCbk = new QuadParamCbk(
				totalDecisionVarsPerScenario[s],
				totalDecisionVarsFirstStage,
				totNzSymQuadPart_,
				nzRowsOneBasedSymQuadPart_,
				nzColsOneBasedSymQuadPart_,
				nzValsSymQuadPart_,
				linearCostX_,
				freeTerm_,
				linearCostP_);

			// append to free or to use
			
			quadParamCbkToFree.push_back(ineqQuadParamCbk);
			paramIneqConstrVec->push_back(ineqQuadParamCbk);

			// controls

			globalIndexFreeTermIneq += 1;
			globalIndexLinearPartIneq += totalDecisionVarsPerScenario[s];
			globalIndexLinearCostParamIneq += totalDecisionVarsFirstStage;
			globalIndexIneqQuadScen += totNzSymQuadPart_;
		}

		// alloc
		ParamConvexProb* probIter = new ParamConvexProb(
			totalDecisionVarsPerScenario[s],
			totalDecisionVarsFirstStage,
			totalEqConstrPerScenario[s],
			paramObjFunction,
			mtxParamCbk,
			paramRhsVec,
			paramIneqConstrVec);
		
		// add
		paramConvexProbVec.push_back(probIter);
		
		// controls
		globalIndexQuadraticPerturbation += totNzSymQuadPart_;
		globalIndexCostPerScenario += totalDecisionVarsPerScenario[s];
		globalIndexFixedPart += totNzFixedPart_;
		mtxParamCbkVec.push_back(mtxParamCbk);
	}
	
	// init callback vector
	vector<Cbk*> recourseFunctions;
	
	// create callback functions
	for (int s=0; s<totalScenarios; ++s) {
		
		// create recourse function
		
		OptValSmoCbk* recourseIter = new OptValSmoCbk(
			paramConvexProbVec[s], epsilon, tikhonovMu);
		
		// add recourse function
		
		recourseFunctions.push_back( recourseIter );
	}
	
	// init stoch prog
	StochProg stochProg(totalDecisionVarsFirstStage, 
		costFirstStage, totalEqConsFirstStage,
		rhsFirstStage, totalNonZeroFirstStageW, nzRowsFirstStageW, 
		nzColsFirstStageW, nzValsFirstStageW,
		totalScenarios, probVecScenarios, &recourseFunctions,
		pathToWrtObjIterates);
	
	// solve
	stochProg.solveStochProg(optValOut, optSolOut, optMulEqConsOut, 
		optMutXvarLowerConstrOut, optMutXvarUpperConstrOut);
	
	// free param convex problems
	int totalProblems = paramConvexProbVec.size();
	for (int i=0; i<totalProblems; ++i) {
		delete paramConvexProbVec[i];
	}
	
	// free cbks
	int totalCbks = recourseFunctions.size();
	for (int i=0; i<totalCbks; ++i) {
		delete recourseFunctions[i];
	}
	
	// free mtx param cbk
	int totalMtxCbk = mtxParamCbkVec.size();
	for (int i=0; i<totalMtxCbk; ++i) {
		delete mtxParamCbkVec[i];
	}
	
	// free param cbk to free
	for (int i=0; i<paramCbkToFree.size(); ++i) {
		delete paramCbkToFree[i];
	}
	
	// free quad param cbk
	for (int i=0; i<quadParamCbkToFree.size(); ++i) {
		delete quadParamCbkToFree[i];
	}
	
	// cbk to free
	for (int i=0; i<cbkToFree.size(); ++i) {
		delete cbkToFree[i];
	}
	
	// linear cbk to free
	for (int i=0; i<linearCbkToFree.size(); ++i) {
		delete linearCbkToFree[i];
	}
	
	// free vectors to free
	for (int i=0; i<vectorsParamCbkToFree.size(); ++i) {
		delete vectorsParamCbkToFree[i];
	}
	
	// free vectors of cbks
	for (int i=0; i<vectorsCbkToFree.size(); ++i) {
		delete vectorsCbkToFree[i];
	}
	
	// return
	return stochProg.getTotalIterations();
}
