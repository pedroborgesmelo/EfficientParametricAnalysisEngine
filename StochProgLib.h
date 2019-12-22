/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#ifndef STOCH_PROG_LIB_H_
#define STOCH_PROG_LIB_H_

#ifdef __cplusplus
extern "C" {
#endif

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
	char* pathToWrtObjIterates);

#ifdef __cplusplus
}
#endif

#endif
