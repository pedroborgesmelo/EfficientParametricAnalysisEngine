/*
 * author: pedro.borges.melo@gmail.com
 * date: December/2019
 */

#ifndef CBK_H_
#define CBK_H_

#include <functional>

using namespace std;

class Cbk {
		
	private:
		
		// target to update buffer
		const double BUFFER_TOL = 1e-10;
		// data
		int totalVariables;
		// value buffer
		double* valueLastX;
		double valueBkp;
		// grad buffer
		double* gradLastX;
		double* gradBkp;
		
	public:
		
		Cbk(int totalVariables_);
		
		~Cbk();
		
		virtual double getValueFunc(double* x) = 0;
		
		virtual void getFullGradFunc(double* x, double* grad) = 0;

		int getTotalVariables();

		double getValue(double* x);

		double getPartial(int jOneBased, double* x);

		void getSparseHessVals(double* x, double* vals);
};

#endif
