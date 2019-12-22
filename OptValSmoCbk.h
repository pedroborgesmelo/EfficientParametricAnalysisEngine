/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#ifndef OPT_VAL_SMO_CBK_H_
#define OPT_VAL_SMO_CBK_H_

#include <map>

#include "Cbk.h"
#include "ParamConvexProb.h"

using namespace std;

class OptValSmoCbk : public Cbk {
		
	private:
		
		ParamConvexProb* paramConvexProb;
		double epsilon;
		double mu;
		
	public:
		
		OptValSmoCbk(ParamConvexProb* paramConvexProb_,
			double epsilon_, double mu_);
		
		~OptValSmoCbk();
		
		double getValueFunc(double* x);
		
		double getSqpApproxValueFunc();
		
		void getFullGradFunc(double* x, double* grad);
};

#endif
