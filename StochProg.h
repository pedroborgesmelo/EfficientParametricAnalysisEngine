/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#ifndef STOCH_PROG_H_
#define STOCH_PROG_H_

#include <vector>
#include "IpStdCInterface.h"

#include "Cbk.h"
#include "MySparseMatrix.h"

using namespace std;

class StochProg {

	private:
		
		// first stage objective
		int totalDecisionVars;
		double *linearCost;
		// linear equality constraints
		int totalEqConstr;
		double *rhsValues;
		MySparseMatrix *eqConstrMtx;
		// recourse functions in the objective
		int totalScenarios;
		double* probScenariosVec;
		vector<Cbk*> *recourseFunctions;
		// intermediate callback control
		double lastObjMasterProblem;
		int totalTimesVariationWasBelowTarget;
		const double belowThisVarItIsConsideredNoise = 1e-2;
		const int maxTimesVariationWasBelowTarget = 5;
		// outputs
		int totalIterations;
		double startTime;
		double totalSecondsSolvingSubproblems;
		char* pathToWrtObjIterates;
		
	public:
		
		StochProg(int totalDecisionVars_, double *linearCost_, int totalEqConstr_,
			double *rhsValues_, int totalNzEqConstr_, int *nzRowsEqConstr_, int *nzColsEqConstr_, double *nzValsEqConstr_,
			int totalScenarios_, double* probScenariosVec_, vector<Cbk*> *recourseFunctions_,
			char* pathToWrtObjIterates_);
		
		~StochProg();

		void addOneToTotalIterations();

		int getTotalIterations();

		int getTotalDecisionVars();

		double* getLinearCost();

		int getTotalEqConstr();

		char* getPathToWrtObjIterates();

		double getStartTime();

		double* getRhsValues();

		MySparseMatrix* getEqConstrMtx();

		double getProbabilitiesScenario(int s);

		vector<Cbk*>* getRecourseFunctions();

		void addTimeSolvingSubproblems(double toAdd);

		double getLastObjMasterProblem();

		void setLastObjMasterProblem(double newVal);
		
		int getTotalTimesVariationWasBelowTarget();

		void addOneToTotalTimesVariationWasBelowTarget();

		void setTotalTimesVariationWasBelowTarget(double newVal);
		
		double getBelowThisVarItIsConsideredNoise();
		
		int getMaxTimesVariationWasBelowTarget();
		
		static Bool intermediate_callback_master_prob(
			Index alg_mod, Index iter_count, Number obj_value, Number inf_pr, 
			Number inf_du, Number mu, Number d_norm, Number regularization_size, 
			Number alpha_du, Number alpha_pr, Index ls_trials, UserDataPtr user_data_all);
		
		static Bool eval_f_master_prob(Index n, Number* x, Bool new_x, Number* obj_value, UserDataPtr user_data_all);
		
		static Bool eval_grad_f_master_prob(Index n, Number* x, Bool new_x, Number* grad_f, UserDataPtr user_data_all);
		
		static Bool eval_g_master_prob(Index n, Number* x, Bool new_x, Index m, Number* g, UserDataPtr user_data_all);
		
		static Bool eval_jac_g_master_prob(Index n, Number *x, Bool new_x, Index m, Index nele_jac, 
			Index *iRow, Index *jCol, Number *values, UserDataPtr user_data_all);
		
		static Bool eval_h_master_prob(Index n, Number *x, Bool new_x, Number obj_factor, Index m, 
			Number *lambda, Bool new_lambda, Index nele_hess, 
			Index *iRow, Index *jCol, Number *values, UserDataPtr user_data_all);
		
		void solveStochProg(double* objValOut, double* initiaPrimalIterate, double* initialDualIterateEqConstr,
			double* optMutXvarLowerConstrOut, double* optMutXvarUpperConstrOut);
};

#endif
