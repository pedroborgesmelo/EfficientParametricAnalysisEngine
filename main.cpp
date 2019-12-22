/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#include <vector>
#include <algorithm>
#include <iostream>
#include <stdlib.h>

#include "StochProgLib.h"
#include "LinearSystem.h"
#include "LinearSystemSolver.h"

using namespace std;

// test linear system
void test00() {
		
	vector< pair<int, int> > nzPositionsLinSys;
	nzPositionsLinSys.push_back( make_pair(1, 1) );
	nzPositionsLinSys.push_back( make_pair(1, 3) );
	nzPositionsLinSys.push_back( make_pair(2, 2) );
	nzPositionsLinSys.push_back( make_pair(3, 3) );
	sort( nzPositionsLinSys.begin(), nzPositionsLinSys.end() );

	int totalEqLinearSystem_ = 3;
	int totalNzEqLinearSystem_ = nzPositionsLinSys.size();
	LinearSystem lSys(totalEqLinearSystem_, totalNzEqLinearSystem_, &nzPositionsLinSys );
	
	int n_ = lSys.getTotalEqLinearSystem();
	int* ia_ = lSys.getIaPardiso();
	int* ja_ = lSys.getJaPardiso();
	double* a_ = lSys.getValsPardiso();
	bool symmetricReal = true;

	double* sol = lSys.getSolPardiso();
	double* rhs = lSys.getRhsPardiso();
	
	rhs[0] = 1.0;
	rhs[1] = 1.0;
	rhs[2] = 2.0;

	a_[0] = 1.0;
	a_[1] = 1.0;
	a_[2] = 1.0;

	LinearSystemSolver lSysSolver(n_, ia_, ja_, a_, symmetricReal);

	lSysSolver.solveSystemWithPardiso( rhs, sol );	
	
	cout << "solution linear system = " << endl;
	for (int i=0; i<totalEqLinearSystem_; ++i) {
		cout << "sol " << i << " = " << sol[i] << endl;
	}
}

// simple test
// non-reg optimal sol = 1.0
void test01() {
	
	int maxIterPerSubProblem = 40;
	double epsilon = 0.1;
	double tikhonovMu = 0.1;
	int totalDecisionVarsFirstStage = 1;
	int totalScenarios = 1;
	int totalDecisionVarsPerScenario[] = {2};
	int totalEqConsFirstStage = 1;
	double rhsFirstStage[] = {1.0};
	int totalEqConstrPerScenario[] = {1};
	double costFirstStage[] = {1.0};
	double probVecScenarios[] = {1.0};
	double costPerScenario[] = {-1.0, -1.0};
	double rhsPerScenario[] = {2.0};
	int totalNonZeroFirstStageW = 1;
	int nzRowsFirstStageW[] = {1};
	int nzColsFirstStageW[] = {1};
	double nzValsFirstStageW[] = {1.0};
	int totalNonZeroPerScenarioW[] = {1};
	int nzRowsPerScenarioW[] = {1};
	int nzColsPerScenarioW[] = {1};
	double nzValsPerScenarioW[] = {1.0};
	int totalNonZerPerScenarioT[] = {1};
	int nzRowsPerScenarioT[] = {1};
	int nzColsPerScenarioT[] = {1};
	double nzValsPerScenarioT[] = {1};
	int totalNzSymQuadPerturbationObjPerScenario[] = {2};
	int nzRowsSymQuadPerturbationObjPerScenario[] = {1, 2};
	int nzColsSymQuadPerturbationObjPerScenario[] = {1, 2};
	double nzValsSymQuadPerturbationObjPerScenario[] = {1.0, 1.0};
	int totalIneqPerScenario[] = {1};
	double linearPartIneqPerScenario[] = {0.0, -1.0};
	double freeTermIneqPerScenario[] = {1.0};
	int totNzSymQuadPartIneqPerScenario[] = {2};
	int nzRowsSymQuadPartIneqPerScenario[] = {1, 2};
	int nzColsSymQuadPartIneqPerScenario[] = {1, 2};
	double nzValsSymQuadPartIneqPerScenario[] = {0.0, 0.0};
	double linearCostParamPerIneqPerScenario[] = {0.0};
	double optSolOut[] = {1.0};
	double optValOut = 1e29;
	double optMulEqConsOut[] = {100.0};
	double optMutXvarLowerConstrOut[] = {100.0};
	double optMutXvarUpperConstrOut[] = {100.0};
	char* pathToWrtObjIterates = "out.txt";
	solveNonDiffProbUpperSmoothingWithQuadObjQuadConstraints(
		maxIterPerSubProblem,
		epsilon,
		tikhonovMu,
		totalDecisionVarsFirstStage,
		totalScenarios,
		totalDecisionVarsPerScenario,
		totalEqConsFirstStage,
		rhsFirstStage,
		totalEqConstrPerScenario,
		costFirstStage,
		probVecScenarios,
		costPerScenario,
		rhsPerScenario,
		totalNonZeroFirstStageW,
		nzRowsFirstStageW,
		nzColsFirstStageW,
		nzValsFirstStageW,
		totalNonZeroPerScenarioW,
		nzRowsPerScenarioW,
		nzColsPerScenarioW,
		nzValsPerScenarioW,
		totalNonZerPerScenarioT,
		nzRowsPerScenarioT,
		nzColsPerScenarioT,
		nzValsPerScenarioT,
		totalNzSymQuadPerturbationObjPerScenario,
		nzRowsSymQuadPerturbationObjPerScenario,
		nzColsSymQuadPerturbationObjPerScenario,
		nzValsSymQuadPerturbationObjPerScenario,
		totalIneqPerScenario,
		linearPartIneqPerScenario,
		freeTermIneqPerScenario,
		totNzSymQuadPartIneqPerScenario,
		nzRowsSymQuadPartIneqPerScenario,
		nzColsSymQuadPartIneqPerScenario,
		nzValsSymQuadPartIneqPerScenario,
		linearCostParamPerIneqPerScenario,
		optSolOut,
		&optValOut,
		optMulEqConsOut,
		optMutXvarLowerConstrOut,
		optMutXvarUpperConstrOut,
		pathToWrtObjIterates);

	cout << endl;
	cout << "test01" << endl;
	cout << "opt val = " << optValOut << endl;
	cout << "opt sol = " << optSolOut[0] << endl;
	cout << endl;
}

// test, recourse = |x| using inequality constraints
// non-reg opt val = 0
void test02() {
	
	int maxIterPerSubProblem = 40;
	double epsilon = 0.1;
	double tikhonovMu = 0.1;
	int totalDecisionVarsFirstStage = 1;
	int totalScenarios = 1;
	int totalDecisionVarsPerScenario[] = {1};
	int totalEqConsFirstStage = 0;
	double* rhsFirstStage = NULL;
	int totalEqConstrPerScenario[] = {0};
	double costFirstStage[] = {0.0};
	double probVecScenarios[] = {1.0};
	double costPerScenario[] = {1.0};
	double* rhsPerScenario = NULL;
	int totalNonZeroFirstStageW = 0;
	int* nzRowsFirstStageW = NULL;
	int* nzColsFirstStageW = NULL;
	double* nzValsFirstStageW = NULL;
	int totalNonZeroPerScenarioW[] = {0};
	int* nzRowsPerScenarioW = NULL;
	int* nzColsPerScenarioW = NULL;
	double* nzValsPerScenarioW = NULL;
	int totalNonZerPerScenarioT[] = {0};
	int* nzRowsPerScenarioT = NULL;
	int* nzColsPerScenarioT = NULL;
	double* nzValsPerScenarioT = NULL;
	int totalNzSymQuadPerturbationObjPerScenario[] = {0};
	int* nzRowsSymQuadPerturbationObjPerScenario = NULL;
	int* nzColsSymQuadPerturbationObjPerScenario = NULL;
	double* nzValsSymQuadPerturbationObjPerScenario = NULL;
	int totalIneqPerScenario[] = {2};
	double linearPartIneqPerScenario[] = {-1.0, -1.0};
	double freeTermIneqPerScenario[] = {0.0, 0.0};
	int totNzSymQuadPartIneqPerScenario[] = {0};
	int* nzRowsSymQuadPartIneqPerScenario = NULL;
	int* nzColsSymQuadPartIneqPerScenario = NULL;
	double* nzValsSymQuadPartIneqPerScenario = NULL;
	double linearCostParamPerIneqPerScenario[] = {1.0, -1.0};
	double optSolOut[] = {100.0};
	double optValOut = 1e29;
	double optMulEqConsOut[] = {100.0};
	double optMutXvarLowerConstrOut[] = {100.0};
	double optMutXvarUpperConstrOut[] = {100.0};
	char* pathToWrtObjIterates = "out.txt";
	
	solveNonDiffProbUpperSmoothingWithQuadObjQuadConstraints(
		maxIterPerSubProblem,
		epsilon,
		tikhonovMu,
		totalDecisionVarsFirstStage,
		totalScenarios,
		totalDecisionVarsPerScenario,
		totalEqConsFirstStage,
		rhsFirstStage,
		totalEqConstrPerScenario,
		costFirstStage,
		probVecScenarios,
		costPerScenario,
		rhsPerScenario,
		totalNonZeroFirstStageW,
		nzRowsFirstStageW,
		nzColsFirstStageW,
		nzValsFirstStageW,
		totalNonZeroPerScenarioW,
		nzRowsPerScenarioW,
		nzColsPerScenarioW,
		nzValsPerScenarioW,
		totalNonZerPerScenarioT,
		nzRowsPerScenarioT,
		nzColsPerScenarioT,
		nzValsPerScenarioT,
		totalNzSymQuadPerturbationObjPerScenario,
		nzRowsSymQuadPerturbationObjPerScenario,
		nzColsSymQuadPerturbationObjPerScenario,
		nzValsSymQuadPerturbationObjPerScenario,
		totalIneqPerScenario,
		linearPartIneqPerScenario,
		freeTermIneqPerScenario,
		totNzSymQuadPartIneqPerScenario,
		nzRowsSymQuadPartIneqPerScenario,
		nzColsSymQuadPartIneqPerScenario,
		nzValsSymQuadPartIneqPerScenario,
		linearCostParamPerIneqPerScenario,
		optSolOut,
		&optValOut,
		optMulEqConsOut,
		optMutXvarLowerConstrOut,
		optMutXvarUpperConstrOut,
		pathToWrtObjIterates);
	
	cout << endl;
	cout << "test02, opt val = " << optValOut << endl;
	cout << endl;
}

// simple test, recourse = |x| using only equality constraints
// non-reg optimal sol = 1.0
void test03() {
	
	int maxIterPerSubProblem = 40;
	double epsilon = 0.1;
	double tikhonovMu = 0.1;
	int totalDecisionVarsFirstStage = 1;
	int totalScenarios = 1;
	int totalDecisionVarsPerScenario[] = {3};
	int totalEqConsFirstStage = 0;
	double* rhsFirstStage = NULL;
	int totalEqConstrPerScenario[] = {2};
	double costFirstStage[] = {1.0};
	double probVecScenarios[] = {1.0};
	double costPerScenario[] = {1.0, 0.0, 0.0};
	double rhsPerScenario[] = {0.0, 0.0};
	int totalNonZeroFirstStageW = 0;
	int* nzRowsFirstStageW = NULL;
	int* nzColsFirstStageW = NULL;
	double* nzValsFirstStageW = NULL;
	int totalNonZeroPerScenarioW[] = {4};
	int nzRowsPerScenarioW[] = {1, 1, 2, 2};
	int nzColsPerScenarioW[] = {1, 2, 1, 3};
	double nzValsPerScenarioW[] = {1.0, -1.0, 1.0, -1.0};
	int totalNonZerPerScenarioT[] = {2};
	int nzRowsPerScenarioT[] = {1, 2};
	int nzColsPerScenarioT[] = {1, 1};
	double nzValsPerScenarioT[] = {1.0, 1.0};
	int totalNzSymQuadPerturbationObjPerScenario[] = {0};
	int* nzRowsSymQuadPerturbationObjPerScenario = NULL;
	int* nzColsSymQuadPerturbationObjPerScenario = NULL;
	double* nzValsSymQuadPerturbationObjPerScenario = NULL;
	int totalIneqPerScenario[] = {0};
	double* linearPartIneqPerScenario = NULL;
	double* freeTermIneqPerScenario = NULL;
	int* totNzSymQuadPartIneqPerScenario = NULL;
	int* nzRowsSymQuadPartIneqPerScenario = NULL;
	int* nzColsSymQuadPartIneqPerScenario = NULL;
	double* nzValsSymQuadPartIneqPerScenario = NULL;
	double* linearCostParamPerIneqPerScenario = NULL;
	double optSolOut[] = {1.0};
	double optValOut = 1e29;
	double optMulEqConsOut[] = {100.0};
	double optMutXvarLowerConstrOut[] = {100.0};
	double optMutXvarUpperConstrOut[] = {100.0};
	char* pathToWrtObjIterates = "out.txt";
	solveNonDiffProbUpperSmoothingWithQuadObjQuadConstraints(
		maxIterPerSubProblem,
		epsilon,
		tikhonovMu,
		totalDecisionVarsFirstStage,
		totalScenarios,
		totalDecisionVarsPerScenario,
		totalEqConsFirstStage,
		rhsFirstStage,
		totalEqConstrPerScenario,
		costFirstStage,
		probVecScenarios,
		costPerScenario,
		rhsPerScenario,
		totalNonZeroFirstStageW,
		nzRowsFirstStageW,
		nzColsFirstStageW,
		nzValsFirstStageW,
		totalNonZeroPerScenarioW,
		nzRowsPerScenarioW,
		nzColsPerScenarioW,
		nzValsPerScenarioW,
		totalNonZerPerScenarioT,
		nzRowsPerScenarioT,
		nzColsPerScenarioT,
		nzValsPerScenarioT,
		totalNzSymQuadPerturbationObjPerScenario,
		nzRowsSymQuadPerturbationObjPerScenario,
		nzColsSymQuadPerturbationObjPerScenario,
		nzValsSymQuadPerturbationObjPerScenario,
		totalIneqPerScenario,
		linearPartIneqPerScenario,
		freeTermIneqPerScenario,
		totNzSymQuadPartIneqPerScenario,
		nzRowsSymQuadPartIneqPerScenario,
		nzColsSymQuadPartIneqPerScenario,
		nzValsSymQuadPartIneqPerScenario,
		linearCostParamPerIneqPerScenario,
		optSolOut,
		&optValOut,
		optMulEqConsOut,
		optMutXvarLowerConstrOut,
		optMutXvarUpperConstrOut,
		pathToWrtObjIterates);

	cout << endl;
	cout << "test01" << endl;
	cout << "opt val = " << optValOut << endl;
	cout << "opt sol = " << optSolOut[0] << endl;
	cout << endl;
}


int main(int argc, char* argv[]) {
	
	// run test
	// test00();
	// test01();
	// test02();
	test03();
	
	// return
	return 0;
}
