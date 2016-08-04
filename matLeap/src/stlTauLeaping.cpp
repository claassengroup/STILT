/*
 * stlTauLeaping.cpp
 *
 *  Created on: Jun 29, 2016
 *      Author: fjustin
 */


#include "mex.h"
#include "tauLeaping.h"
#include <iostream>
#include <iomanip>
#include "model_def.h"

#include <map>

#if defined (_WIN32)
    #include <windows.h>
#elif defined (__linux__)
    #include <unistd.h>
#endif

#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
#else
    extern bool utIsInterruptPending();
#endif

using namespace std;

//#define dalloc(N) (double*)mxMalloc(sizeof(double) * N)

inline bool mxIsScalar_(const mxArray* a) {
	return mxGetNumberOfElements(a)==size_t(1);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
  const mxArray *prhs[]) {


	/* extract model parameters */

	// X0: current state of the system SPNR * PNr
	// Theta: matrix of parameter values for the simulations THNR * PNr
	// dt: timestep
	// NIntervals: number of intervals for intermediate output
	// options struct

	// validate model parameters
	if (nrhs < 4)
		mexErrMsgTxt("Not enough arguments supplied");

	if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]))
		mexErrMsgTxt("Incorrect type of input arguments!");

	if (nrhs > 4 && !mxIsStruct(prhs[4]))
		mexErrMsgTxt("Incorrect type of input arguments!");

	if (!mxIsScalar_(prhs[2])  || !mxIsScalar_(prhs[3]))
		mexErrMsgTxt("Incorrect dimensions of input arguments!");

	// validate options struct
	enum DATA_TYPES {
		SCALAR_DOUBLE, SCALAR_INT, LOGICAL
	};

	map<string, DATA_TYPES> optionsFormat;
	optionsFormat["DELTA_PE"] = SCALAR_DOUBLE;
	optionsFormat["N_STIFF"] = SCALAR_INT;
	optionsFormat["NEWTON_EPS"] = SCALAR_DOUBLE;
	optionsFormat["MAX_ITER"] = SCALAR_INT;
	optionsFormat["SSA_CONDITION_NUMBER"] = SCALAR_INT;
	optionsFormat["CRITICAL_NUMBER"] = SCALAR_INT;
	optionsFormat["USE_IMPLICIT"] = LOGICAL;
	optionsFormat["N_SSA_1"] = SCALAR_INT;
	optionsFormat["N_SSA_2"] = SCALAR_INT;
	optionsFormat["EPSILON"] = SCALAR_DOUBLE;
	optionsFormat["N_C"] = SCALAR_INT;
	optionsFormat["SSA_ONLY"] = LOGICAL;


	RUN_OPTIONS run_options;
	if (nrhs > 4) {
		const mxArray* p_options = prhs[4];

		for (auto it = optionsFormat.begin(); it != optionsFormat.end(); ++it) {
//			cout << "searching " << it->first.c_str() << endl;

			const mxArray *tok = mxGetField(p_options, 0, it->first.c_str());
			if ( tok != NULL) {
				// validate
				stringstream ss;
				ss << "Incorrect type of options." << it->first;
//				cout << "Found " << it->first << " = " <<  << endl;


				switch (optionsFormat[it->first] ) {
				case SCALAR_DOUBLE:
					if (!mxIsDouble(tok))
						mexErrMsgTxt(ss.str().c_str());
					if (!mxIsScalar_(tok))
						mexErrMsgTxt(ss.str().c_str());

					run_options.set(it->first, *mxGetPr(tok));
				break;
				case SCALAR_INT:
				{
					if (!mxIsDouble(tok))
						mexErrMsgTxt(ss.str().c_str());
					if (!mxIsScalar_(tok))
						mexErrMsgTxt(ss.str().c_str());
					int arg = int(*mxGetPr(tok));
					run_options.set(it->first, arg);
				}
				break;

				case LOGICAL:
					if (!mxIsLogical(tok))
						mexErrMsgTxt(ss.str().c_str());
					run_options.set(it->first, *mxGetLogicals(tok) );
				break;
				}
			}

		}
	}

//	run_options.print();

	const mxArray *p_X0 = prhs[0], *p_Theta = prhs[1];

	int PNr;
	double *X0, *Theta;

	if (mxGetN(p_X0) == 1 && mxGetN(p_Theta) > 1) {
		// single intial condition, multiple values for theta
		double *X0_ = mxGetPr(p_X0);
		PNr = mxGetN(p_Theta);
		X0 = (double*)mxMalloc(sizeof(double) * SPNR * PNr);
		for (int k = 0; k < PNr; ++k) {
			copy(X0_, X0_+SPNR, X0+k*SPNR);
		}
		Theta = mxGetPr(p_Theta);
	} else if (mxGetN(p_Theta) == 1 && mxGetN(p_X0) > 1)  {
		// multiple initial conditions, single value for theta
		double *Theta_ = mxGetPr(p_Theta);
		PNr = mxGetN(p_X0);
		Theta = (double*)mxMalloc(sizeof(double) * THNR * PNr);
		for (int k = 0; k < PNr; ++k) {
			copy(Theta_, Theta_+SPNR, Theta+k*THNR);
		}
		X0 = mxGetPr(p_X0);
	} else {
		if (!(mxGetN(p_X0)==mxGetN(p_Theta)))
			mexErrMsgTxt("X0 and Theta must have same number of cols!");
		PNr = mxGetN(p_X0); // same as for Theta
		X0 = mxGetPr(p_X0); 
		Theta = mxGetPr(p_Theta);
	}


	if (mxGetM(p_X0) != SPNR)
		mexErrMsgTxt("X0 must be species number x number of replicates");

	if (mxGetM(p_Theta) != THNR)
		mexErrMsgTxt("Theta must be parameters number x number of replicates");


	double dt = *mxGetPr(prhs[2]);
	int NIntervals = static_cast<int>(*(mxGetPr(prhs[3])));
//	int SPNR = mxGetM(p_X0), THNR = mxGetM(p_Theta);


	double *X_out = static_cast<double*>(mxMalloc(sizeof(double)*PNr * SPNR * (NIntervals+1))); // state of the system
	double *R_out = static_cast<double*>(mxMalloc(sizeof(double)*PNr * THNR * (NIntervals+1))); // state of the system
	double *G_out = static_cast<double*>(mxMalloc(sizeof(double)*PNr * THNR * (NIntervals+1))); // state of the system

	fill_n(R_out, PNr * THNR * (NIntervals+1), 0.0);
	fill_n(G_out, PNr * THNR * (NIntervals+1), 0.0);

	if (nlhs==0)
		return;
	std::vector<OUTPUT_HANDLE> fileHandles;

	CALL_STATS *call_stats = (CALL_STATS*)mxMalloc(sizeof(CALL_STATS) * PNr);

	for (int k=0; k<PNr; k++)
	{
		try {
			if (utIsInterruptPending()) {
				mexErrMsgTxt("Execution interrupted");
				return;
			}
			CALL_STATS local_call_stats = tauLeaping(&X0[k*SPNR], &Theta[k*THNR], NIntervals, dt, run_options,
				&X_out[k*SPNR*(NIntervals+1)], &R_out[k*THNR*(NIntervals+1)], &G_out[k*THNR*(NIntervals+1)],
				false, fileHandles);
			call_stats[k].ssaCalls = local_call_stats.ssaCalls;
			call_stats[k].expCalls = local_call_stats.expCalls;
			call_stats[k].impCalls = local_call_stats.impCalls;

		} catch (const std::exception &e) {
			mexErrMsgTxt(e.what());
		}
	}

	/* output simulated states X */

	if (nlhs > 0)
	{
//		plhs[0] = mxCreateDoubleMatrix( SPNR, PNr, mxREAL);
		int dims[] = {SPNR, NIntervals+1, PNr};
		plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
		mxFree(mxGetPr(plhs[0]));
		mxSetPr(plhs[0], X_out);
	} else
	{
		mxFree(X_out);
		mxFree(R_out);
		mxFree(G_out);
		mxFree(call_stats);
		return;
	}

	/* output time */
	if (nlhs > 1)
	{
		plhs[1] = mxCreateDoubleMatrix(1, NIntervals+1, mxREAL);
		double *t_out = mxGetPr(plhs[1]);
		t_out[0] = 0.0;

		double delta_t = dt/(double)NIntervals;

		for (int k = 1; k <= NIntervals; ++k) {
			t_out[k] = k*delta_t;
		}

	}
	else {
		mxFree(R_out);
		mxFree(G_out);
		mxFree(call_stats);
		return;
	}

	/* output r */
	if (nlhs > 2)
	{
		int dims[] = {THNR, NIntervals+1, PNr};
		plhs[2] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
		mxFree(mxGetPr(plhs[2]));
		mxSetPr(plhs[2], R_out);
	}
	else {
		mxFree(R_out);
		mxFree(G_out);
		mxFree(call_stats);
		return;
	}

	/* output G */
	if (nlhs > 3)
	{
		int dims[] = {THNR, NIntervals+1, PNr};
		plhs[3] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
		mxFree(mxGetPr(plhs[3]));
		mxSetPr(plhs[3], G_out);
	}
	else {
		mxFree(G_out);
		mxFree(call_stats);
		return;
	}

	/* output call stats */
	if (nlhs > 4) {
		// output the call statistics
		char a[] = "ssaCalls", b[] = "expCalls", c[] = "impCalls";
		const char *fieldNames[] = {a, b, c};

		mxArray *s_callStats = mxCreateStructMatrix(1, PNr, 3, &fieldNames[0]);
		for (int k = 0; k < PNr; ++k) {
			mxSetField(s_callStats, k, "ssaCalls", mxCreateDoubleScalar(call_stats[k].ssaCalls));
			mxSetField(s_callStats, k, "expCalls", mxCreateDoubleScalar(call_stats[k].expCalls));
			mxSetField(s_callStats, k, "impCalls", mxCreateDoubleScalar(call_stats[k].impCalls));
		}

		plhs[4] = s_callStats;

	} else {
		mxFree(call_stats);
		return;
	}


}


