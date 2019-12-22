/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#ifndef PARAM_CONVEX_PROB_H_
#define PARAM_CONVEX_PROB_H_

#include <map>
#include <vector>
#include "IpStdCInterface.h"

#include "Cbk.h"
#include "ParamCbk.h"
#include "QuadParamCbk.h"
#include "MtxParamCbk.h"
#include "LinearSystem.h"
#include "LinearSystemSolver.h"

using namespace std;

class ParamConvexProb {

	private:
		
		// dimensions
		int totalDecisionVars;
		int totalParameters;
		int totalEqualityConstraints;
		// objective data
		ParamCbk* paramObjFunction;
		// linear equality constraints
		MtxParamCbk *mtxParamCbk;
		// rhs parametric vector
		vector<Cbk*> *paramRhsVec;
		// inequality constraints
		vector<ParamCbk*> *paramIneqConstrVec;
		// resulting sparsity patterns
		int totalNzJacobianAllConstraints;
		map< pair<int,int> , int> *nzRowNzColToPosValMap;
		int* positionsToAddTikhonovTerms;
		// buffer
		const double SAFE_GUARD_INTERIOR_COND_TRUNCATION = 1e-9;
		const double BUFFER_TOL = 1e-15;
		double lastEpsilon;
		double lastMu;
		double* lastParameterVecRegBuffer;
		double lastRegOptVal;
		double* lastRegPrimalSol;
		double* lastRegDualSolEqConstraints;
		double* lastRegLowerIneqMul;
		double* lastRegUpperIneqMul;
		double* lastRegIneqConstraintValues;
		// buffer to find starting point with SQP model
		const bool WARM_START_IPOPT_STEP = true;
		double sqpOptVal;
		double* sqpBufferHessX;
		int* sqpNzRowZeroBasedEqConstr;
		int* sqpNzColZeroBasedEqConstr;
		double* sqpNzValEqConstr;
		// data for the lin sys of the derivatives of solution mappings
		vector< pair<int, int> > *nzPositionsLinSys;
		LinearSystem *lSysDerSolMap;
		LinearSystemSolver *lSysSolver;
		
	public:
		
		ParamConvexProb(int totalDecisionVars_, int totalParameters_, int totalEqualityConstraints_,
			ParamCbk* paramObjFunction_, MtxParamCbk* mtxParamCbk_, 
			vector<Cbk*> *paramRhsVec_, vector<ParamCbk*> *paramIneqConstrVec_);

		ParamConvexProb();
		
		~ParamConvexProb();

		double getSqpOptVal();
		
		int getTotalDecisionVars();
		
		int getTotalParameters();

		int getTotalEqualityConstraints();
		
		int getTotalInequalityConstraints();

		int* getPositionsToAddTikhonovTerms();

		ParamCbk* getParamObjFunction();

		MtxParamCbk* getMtxParamCbk();

		vector<Cbk*>* getParamRhsVec();

		ParamCbk* getInequalityConstraint(int i);

		map< pair<int,int> , int>* getNzRowNzColToPosValMap();

		double getLastEpsilon();

		double getLastMu();

		double getLastRegOptVal();

		double* getLastRegPrimalSol();

		double* getLastRegDualSolEqConstraints();

		double* getLastParameterVecRegBuffer();

		double* getLastDerPrimalSol();

		void calcWarmStartForIpoptWithCplex();

		static Bool intermediate_callback_param_convex_prob(
			Index alg_mod, Index iter_count, Number obj_value, Number inf_pr, 
			Number inf_du, Number mu, Number d_norm, Number regularization_size, 
			Number alpha_du, Number alpha_pr, Index ls_trials, UserDataPtr user_data_all);
		
		static Bool eval_f_param_convex_prob(Index n, Number* x, Bool new_x, Number* obj_value, UserDataPtr user_data_all);
		
		static Bool eval_grad_f_param_convex_prob(Index n, Number* x, Bool new_x, Number* grad_f, UserDataPtr user_data_all);
		
		static Bool eval_g_param_convex_prob(Index n, Number* x, Bool new_x, Index m, Number* g, UserDataPtr user_data_all);
		
		static Bool eval_jac_g_param_convex_prob(Index n, Number *x, Bool new_x, Index m, Index nele_jac, 
			Index *iRow, Index *jCol, Number *values, UserDataPtr user_data_all);
		
		static Bool eval_h_param_convex_prob(Index n, Number *x, Bool new_x, Number obj_factor, Index m, 
			Number *lambda, Bool new_lambda, Index nele_hess, 
			Index *iRow, Index *jCol, Number *values, UserDataPtr user_data_all);
		
		int getTotalNzForJacobianOfAllConstraints();
		
		void fillNzFoHessianOfLagrangian();
		
		void calcNonTrivialSparsityPatternForLinSys();

		void setFillingStructionsForParamCbks();
		
		int getTotalNzFoHessianOfLagrangian();

		void recalcRegBuffer();
		
		void updateRegBuffer(double epsilon, double mu, double* p);
		
		void buildSysForDerOfSolutionMappings();
		
		void destroySysForDerOfSolutionMappings();

		void calcRhsForDerSolutionMappingsCalculation(int jOneBased);
		
		void calcPartialDerSolMap(int jOneBased);
};

#endif
