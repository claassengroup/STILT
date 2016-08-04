/*
 * performTauLeaping.cpp
 *
 *  Created on: Jun 30, 2016
 *      Author: fjustin
 */


#include "tauLeaping.h"
#include <random>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <float.h>
#include <sstream>
#include <array>

#include <ctime>
#include "Eigen/Core"
#include "Eigen/LU"

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/generator_iterator.hpp>

using namespace std;

#ifdef MEX
#define dalloc(N) (double*)mxMalloc(sizeof(double) * N)
#define FREE mxFree
#include "mex.h"
#else
#define dalloc(N) (double*)malloc(sizeof(double) * N)
#define FREE free
#endif

#define S(i,j) StoichMatrix[THNR*i+j]
#define Sed(i,j) StoichMatrix_ed[THNR*i+j]

typedef boost::minstd_rand base_generator_type;

static array<vector<int>, THNR> reactionMap;

CALL_STATS tauLeaping (const double *X0, const double *Theta, int NIntervals, const double dt,
		RUN_OPTIONS run_options, double *X, double *R, double *G,
		bool writeOutput, vector<OUTPUT_HANDLE> &fileHandles) {
	// perform the tau leaping, outputting results at each of the intervals

	// copy the left time point to the output
	copy(X0, X0+SPNR, X);
	fill_n(R, (NIntervals+1)*THNR, 0.0);
	fill_n(G, (NIntervals+1)*THNR, 0.0);

	CALL_STATS call_stats;

	double ddt = dt/static_cast<double>(NIntervals);
	for (int k = 0; k < NIntervals; ++k) {
		// use previous endpoint as initial condition for this interval
		CALL_STATS tmp_stats = tauLeapingInterval(&X[SPNR*k], &R[THNR*k], &G[THNR*k],
				Theta, ddt, run_options, &X[SPNR*(k+1)], &R[THNR*(k+1)], &G[THNR*(k+1)]);
		call_stats.expCalls += tmp_stats.expCalls;
		call_stats.impCalls += tmp_stats.impCalls;
		call_stats.ssaCalls += tmp_stats.ssaCalls;

		if (!writeOutput)
			continue;

		// write to file
		for (int i = 0; i < SPNR; ++i) {
			fileHandles[k].X << setw(10) << X[SPNR*k+i] ;
		}
		fileHandles[k].X << endl;

		for (int j = 0; j < THNR; ++j) {
			fileHandles[k].R << setw(10) << R[THNR*k+j];
			fileHandles[k].G << setw(10) << G[THNR*k+j];
		}
		fileHandles[k].R << setw(10) << endl;
		fileHandles[k].G << setw(10) << endl;
	}

	return call_stats;

}


CALL_STATS tauLeapingInterval (const double *X0, const double *R0, const double *G0,
		const double *Theta, double dt, RUN_OPTIONS run_options, double *X, double *R, double *G)
{
	// extract run options
	int n_c = run_options.N_C;
	double eps = run_options.EPSILON;
	bool SSA_ONLY = run_options.SSA_ONLY;
	int ssaSteps1 = SSA_ONLY ? std::numeric_limits<int>::max() : run_options.N_SSA_1;
	int ssaSteps2 = run_options.N_SSA_2;
	double delta_pe = run_options.DELTA_PE;
	int N_stiff = run_options.N_STIFF;
	int SSA_CONDITION_NUMBER = run_options.SSA_CONDITION_NUMBER;
	int CRITICAL_NUMBER = run_options.CRITICAL_NUMBER;
	bool USE_IMPLICIT = run_options.USE_IMPLICIT;

	// initialize RNGs
	static unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
//	static default_random_engine generator (seed);
//	static uniform_real_distribution<double> uniformDist(0.0, 1.0);
	static base_generator_type generator(seed);
	static boost::uniform_real<> uni_dist(0,1);
	static boost::variate_generator<base_generator_type&, boost::uniform_real<> > uniformDist(generator, uni_dist);

	static bool computeReactions = true;
	static int StoichMatrix[SPNR * THNR];
	static int StoichMatrix_ed[SPNR * THNR];

	if (computeReactions) {
		// get stoichiometric matrix
		computeStoichiometricMatrix(StoichMatrix, StoichMatrix_ed);

		for (int j = 0; j < THNR; ++j) {
			for (int i = 0; i < SPNR; ++i) {
				if (S(i,j) != 0)
					reactionMap[j].push_back(i);
			}
		}

		computeReactions = false;
	}

//	for (int j = 0; j < THNR; ++j) {
//		cout << "reaction " << j;
//		for (auto it=reactionMap[j].begin(); it != reactionMap[j].end(); it++)
//		{
//			cout << " " << *it;
//		}
//		cout << endl;
//	}

	double t=0.0;

	copy(X0, X0+SPNR, X);
	copy(R0, R0+THNR, R);
	copy(G0, G0+THNR, G);

	double X_tmp[SPNR];
	double a[THNR], a_tmp[THNR];
	double *a_ptr[] = {a}, *a_tmp_ptr[] = {a_tmp};

	// statistics of simulation
	int imLeap=0, exLeap=0, SSA_call=0;

	bool critical[THNR], reactantSpecies[SPNR], pe[THNR];
	fill_n(reactantSpecies, SPNR, false);

	/* properties of the stoichiometric matrix */
	// reactant species
	for (int i = 0; i < SPNR; i++) {
		for (int j = 0; j < THNR; ++j) {
			if (Sed(i,j) != 0)
			{
				reactantSpecies[i] = true;
				break;
			}
		}
	}

	// reversible reaction pairs
	list<pair<int,int> > reversiblePairs;

	for (int j = 0; j < THNR; ++j) {
		for (int j2 = j+1; j2 < THNR; ++j2) {
			bool isPair=true;
			for (int i = 0; i < SPNR; ++i) {
				isPair = isPair && (S(i,j)==-S(i,j2));
			}
			if (isPair)
				reversiblePairs.push_back(pair<int,int>(j,j2));
		}
	}
	bool SSA_flag = SSA_ONLY;

	// number of reaction firings
	int k[THNR];

	/* specific to the current state below */

	// initialize variables
	double thisG[THNR];
	double Xnew[SPNR];


	while (t < dt)
	{
		// current propensities
		computePropensities(X, Theta, 0, a_ptr);
		// compute total propensity
		double a0 = 0.0;
		for (int j = 0; j < THNR; ++j) {
			a0 += a[j];
		}

		// if total propensity is zero then exit the simulation
		if (a0 < DBL_EPSILON) {
			t = dt;
			break;
		}

		// Step 1: determine critical reactions
		for (int j = 0; j < THNR; ++j) {
			critical[j] = false;
		}
		for (int i = 0; i < SPNR; ++i) {
			for (int j = 0; j < THNR; ++j) {
				if (a[j] == 0)
					continue;

				int nu_ij = S(i,j);
				if ( nu_ij < 0) {
					int L = (int)floor(X[i] / abs(double(nu_ij)));
					critical[j] = critical[j] || (L < n_c);
				}
			}
		}

		// determine reactions in partial equilibrium
		for (int j = 0; j < THNR; ++j) {
			pe[j] = false;
		}

		for (list<pair<int,int> >::iterator it = reversiblePairs.begin(); it != reversiblePairs.end(); it++) {
			int j1=it->first, j2=it->second;
			if (abs(a[j1] - a[j2]) <= delta_pe * min(a[j1],a[j2]))
			{
				pe[j1] = true;
				pe[j2] = true;
			}
		}


		double mu_ex[SPNR], sigma2_ex[SPNR],
			mu_im[SPNR], sigma2_im[SPNR];

		// compute explicit mu and sigma
		fill_n(mu_ex, SPNR, 0.0);
		fill_n(sigma2_ex, SPNR, 0.0);

		for (int j = 0; j < THNR; ++j) {
			if (critical[j])
				continue;
			for (auto it = reactionMap[j].begin(); it != reactionMap[j].end(); it++) {
				int i = *it;
				mu_ex[i] += static_cast<double>(S(i,j))*a[j];
				sigma2_ex[i] += pow(static_cast<double>(S(i,j)),2)*a[j];
			}
		}

//		for (int i = 0; i < SPNR; ++i) {
//			mu_ex[i] = 0;
//			sigma2_ex[i] = 0;
//			for (int j = 0; j < THNR; ++j) {
//				if (critical[j])
//					continue;
//				mu_ex[i] += static_cast<double>(S(i,j))*a[j];
//				sigma2_ex[i] += pow(static_cast<double>(S(i,j)),2)*a[j];
//			}
//		}




		// compute reaction orders g_i
		double g[SPNR];
		computeG(X, StoichMatrix, StoichMatrix_ed, g);

		// step 3: compute the candidate leaps
		// explicit tau leap: Eq 8
		double tau_ex = INFINITY;
		for (int i = 0; i < SPNR; i++) {
			if (!reactantSpecies[i])
				continue;
			tau_ex = min(tau_ex, max(eps*X[i]/g[i],1.0)/abs(mu_ex[i]));
			tau_ex = min(tau_ex, pow(max(eps*X[i]/g[i],1.0),2)/sigma2_ex[i]);
		}

		// implicit tau leap: Eq 14
		double tau_im = INFINITY;

		fill_n(mu_im, SPNR, 0.0);
		fill_n(sigma2_im, SPNR, 0.0);

		for (int j = 0; j < THNR; ++j) {
			if (critical[j] || pe[j])
				continue;

			for (auto it = reactionMap[j].begin(); it != reactionMap[j].end(); it++) {
				int i = *it;
				mu_im[i] += static_cast<double>(S(i,j))*a[j];
				sigma2_im[i] += pow(static_cast<double>(S(i,j)),2)*a[j];
			}
		}

//		for (int i = 0; i < SPNR; ++i) {
//			mu_im[i] = 0;
//			sigma2_im[i] = 0;
//			for (int j = 0; j < THNR; ++j) {
//				if (critical[j] || pe[j])
//					continue;
//
//				mu_im[i] += static_cast<double>(S(i,j))*a[j];
//				sigma2_im[i] += pow(static_cast<double>(S(i,j)),2)*a[j];
//			}
//		}

		for (int i = 0; i < SPNR; i++) {
			if (!reactantSpecies[i])
				continue;
			tau_im = min(tau_im, max(eps*X[i]/g[i],1.0)/abs(mu_im[i]));
			tau_im = min(tau_im, pow(max(eps*X[i]/g[i],1.0),2)/sigma2_im[i]);
		}

		bool stiff = tau_im > N_stiff*tau_ex;

		double tau_1 = stiff ? tau_im : tau_ex;

		bool invalidLeap = true;
		// loop until the proposed leap is valid
		while (invalidLeap) {
			// reset temp variables
			fill_n(thisG, THNR, 0.0);
			fill_n(k, THNR, 0);

			for (int i = 0; i < SPNR; i++)
				X_tmp[i] = X[i];

			// total propensity of critical reactions only
			double ac0 = 0.0;
			for (int j = 0; j < THNR; ++j) {
				if (!critical[j])
					continue;
				ac0 += a[j];
			}

			// step 4: check if need to use SSA now
			if (SSA_ONLY || tau_1 < SSA_CONDITION_NUMBER * 1.0/a0 || ac0 > a0/CRITICAL_NUMBER)
			{
				int SSA_steps = (SSA_flag) ? ssaSteps1 : ssaSteps2;
				SSA_flag = true;
				// directly update R and G!
#ifdef DEBUG
//				cout << "SSA call @t = " << t << endl;
#endif
				performSSA(X, Theta, StoichMatrix, dt, SSA_steps, t, R, G);
#ifdef DEBUG
//				cout << "new time " << t << endl;
//				for (int j=0; j<SPNR; j++)
//					cout << setw(10) << R[j];
//				cout << endl;
#endif
//				for (int i = 0; i < SPNR; ++i) {
//					assert(X[i]>=0);
//				};


				SSA_call++;
				break;
			} else {
				SSA_flag = false;
			}

			// step 5: sum of critical reaction propensities

			double tau_2;
			if (ac0 != 0) {
				boost::random::exponential_distribution<double> expDist(ac0);
				boost::variate_generator< base_generator_type&, boost::random::exponential_distribution<double> > rvt(generator, expDist);
				tau_2 = rvt();
			} else {
				tau_2 = INFINITY;
			}
			double tau;

			// step 6: choose tau and fire reactions
			if (tau_2 > tau_1)
			{
				tau = tau_1;
				if (t+tau > dt)
				{
					tau=dt-t;
				}

				if (!stiff || !USE_IMPLICIT)
				{
					// not stiff: regular tau leaping
					// no critical reactions firing
					for (int j = 0; j < THNR; ++j) {
						if (critical[j])
							continue;

						if (a[j]==0) {
							k[j] = 0;
						} else {
							boost::random::poisson_distribution<int> distribution(a[j]*tau);
							boost::variate_generator<base_generator_type&, boost::random::poisson_distribution<int> > rvt(generator, distribution);
							k[j] = rvt();
						}
						for (auto it = reactionMap[j].begin(); it != reactionMap[j].end(); it++)
							X_tmp[*it] += k[j]*S(*it,j);
					}
					exLeap++;
					SSA_flag = true;

					// update states for propensity computation
					computePropensities(X_tmp, Theta, 0, a_tmp_ptr);
					for (int j=0; j<THNR; j++)
						thisG[j] = (a[j]+a_tmp[j])/(2.0*Theta[j])*tau;
				} else {
					performStiffTau(X, Theta, StoichMatrix, critical, tau, k, thisG, run_options);
					SSA_flag = false;
					imLeap++;
				}

			} else { // tau_2 <= tau_1
				tau = tau_2;
				if (t+tau > dt)
				{
					// next critical reaction is after the end of the simulation
					// don't fire it
					tau=dt-t;
				}
				else {
					// fire a critical reaction
					// but only if we actually reached the time for it to fire, not before then

					double r = uniformDist(); // (generator)
					double cumsum=0.0;
					int j;
					for (j = 0; j < THNR; ++j) {
						if (!critical[j])
							continue;
						cumsum += a[j];
						if ((cumsum/ac0) >= r)
							break;
					}
					k[j]++;
					// don't update state until the end of the interval
				}

				if (!USE_IMPLICIT || !stiff || (stiff && tau_2 <= tau_ex))
				{
					// use explicit tau leaping
					// tau-leap the non-critical reactions

					for (int j = 0; j < THNR; ++j) {
						if (critical[j])
							continue;

						if (a[j] == 0) {
							k[j] = 0;
						} else {
							boost::random::poisson_distribution<int> distribution(a[j]*tau);
							boost::variate_generator<base_generator_type&, boost::random::poisson_distribution<int> > rvt(generator, distribution);
							k[j] = rvt();
						}
						for (auto it = reactionMap[j].begin(); it != reactionMap[j].end(); it++)
							X_tmp[*it] += k[j]*S(*it,j);
					}
					exLeap++;
					SSA_flag = true;
					computePropensities(X_tmp, Theta, 0, a_tmp_ptr);
					for (int j=0; j<THNR; j++)
						thisG[j] = (a[j]+a_tmp[j])/(2.0*Theta[j])*tau;
				} else {
					performStiffTau(X, Theta, StoichMatrix, critical, tau, k, thisG, run_options);
					SSA_flag = false;
					imLeap++;
				}
			}

			// update states
			invalidLeap = false;
			copy(X, X+SPNR, Xnew);
			for (int j = 0; j < THNR; ++j) {
				for (auto it = reactionMap[j].begin(); it != reactionMap[j].end(); it++)
					Xnew[*it] += k[j]*S(*it,j);
			}
//			invalidLeap = any_of(Xnew, Xnew+SPNR, [](int i){return i<0;});
			bool invalidLeap = false;
			for (int i = 0; i < SPNR; ++i) {
				if (Xnew[i]<0) {
					invalidLeap = true;
					break;
				}
			}

			if (!invalidLeap)
			{
				for (int i=0; i< SPNR; i++)
					X[i] = Xnew[i];
				for (int j = 0; j < THNR; ++j) {
					R[j] += (double)k[j];
					G[j] += thisG[j];
#ifdef DEBUG
//				cout << setw(10) << R[j] ;
#endif
				}
#ifdef DEBUG
//				cout << endl;
#endif
				t += tau;
			}
			else
			{
				tau_1 /= 2;
			}
		}



	}

	CALL_STATS call_stats;
	call_stats.expCalls = exLeap;
	call_stats.impCalls = imLeap;
	call_stats.ssaCalls = SSA_call;

	return call_stats;
}

void performSSA(double *X,const double *Theta,const int *StoichMatrix,const double dt, const int SSA_steps, double &t, double *R, double *G) {
	static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//	static default_random_engine generator (seed);
//	static uniform_real_distribution<double> uniformDist(0.0, 1.0);
	static boost::uniform_real<> uni_dist(0,1);
	static base_generator_type generator(seed);
	static boost::variate_generator<base_generator_type&, boost::uniform_real<> > uniformDist(generator, uni_dist);


	double a[THNR];
	double *a_ptr[] = {a};
    int lastReaction = -1;
    
	for (int step=0; step<SSA_steps && t < dt; step++) {
		computePropensities(X, Theta, lastReaction+1, a_ptr);
		// compute a0
		double a0 = 0.0;
		for (int j=0; j<THNR; j++)
			a0 += a[j];

		if (a0 < DBL_EPSILON) {
			// propensities are zero
			// no reactions fire
			t = dt;
			break;
		}

		double r = uniformDist(); // uniformDist(generator);
		double cumsum = 0.0;
		int j;
		for (j=0; j<THNR; j++)
		{
			// check for rounding error?
			if (a[j] < DBL_EPSILON)
				a[j] = 0.0;

			cumsum += a[j];
			if (cumsum/a0 > r)
				break;
		}

		boost::random::exponential_distribution<double> expDist(a0);
		boost::variate_generator< base_generator_type&, boost::random::exponential_distribution<double> > rvt(generator, expDist);
		double tau = rvt();
		if (t+tau <= dt)
		{
			// reaction fires
			for (auto it=reactionMap[j].begin(); it != reactionMap[j].end(); it++) {
				X[*it] += S(*it,j);
			}

//			for (int i = 0; i < SPNR; ++i) {
//				X[i] += S(i,j);
//			}
			R[j]++;
		} else {
			// time runs out
			tau = dt-t;
			t = dt;
			for (int j = 0; j < THNR; j++) {
				G[j]+=a[j]*tau/Theta[j];
			}
			break;
		}

		for (int j = 0; j < THNR; j++) {
			G[j]+=a[j]*tau/Theta[j];
		}
		t += tau;
        lastReaction = j;
	}
}
void performStiffTau(const double *X0, const double *Theta, const int *StoichMatrix, const bool *critical, const double tau,
		int *R, double *G, const RUN_OPTIONS run_options)
{
	static unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//	static default_random_engine generator (seed);
	static base_generator_type generator(seed);
	static boost::uniform_real<> uni_dist(0,1);
	static boost::random::variate_generator< base_generator_type&, boost::uniform_real<> > uniformDist(generator, uni_dist);
//	static uniform_real_distribution<double> uniformDist(0.0, 1.0);

	int MAX_ITER = run_options.MAX_ITER;
	double NEWTON_EPS = run_options.NEWTON_EPS;

	// have to perform implicit tau leaping using e.g. trapezoidal rule
	// requires solving for X(t+tau) using Newton's method

	// x'(n+1) = x'(n) - J_F(x'(n))^-1 * F(x'(n))
	// J_F(x'(n)) = I - tau/2 S*da/dx'
	// compute da/dx' from the propensity functions ahead of time and link it in

//	double *a = dalloc(THNR);
	double a[THNR];
	double *a_ptr[] = {a};
	double X[SPNR];
	copy(X0,X0+SPNR,X);

	computePropensities(X, Theta, 0, a_ptr);

	double a0[THNR];
	copy(a,a+THNR,a0);

	// get Poisson jumps
	int P[THNR];
	for (int j = 0; j < THNR; ++j) {
		if (critical[j])
		{
			P[j] = 0;
		}
		else {
			if (a[j]==0) {
				P[j]=0;
			} else {
				boost::random::poisson_distribution<int> distribution(a[j]*tau);
				boost::variate_generator< base_generator_type&, boost::random::poisson_distribution<int> > rvt(generator, distribution);
				P[j] = rvt();
			}
		}
	}

	// first guess for X(t+tau) = X(t) + S*P
	for (int i = 0; i < SPNR; ++i) {
		for (int j=0; j < THNR; j++) {
			X[i] += P[j]*S(i,j);
		}
	}

	int iter=0;
	double deltaX = 1E6;
	double k1[THNR], k2[THNR], F[SPNR];
	double Sk1[SPNR], Sk2[SPNR];
//	double invJF[SPNR];

	// compute k1
	for (int j = 0; j < THNR; ++j) {
		k1[j] = static_cast<double>(P[j]) - tau/2.0 * a[j];
	}
	// compute S*k1
	MVprod(SPNR, THNR, StoichMatrix, k1, Sk1);

	computePropensities(X, Theta, 0, a_ptr);


	// solve for next value of X
	while (iter++ < MAX_ITER && deltaX > NEWTON_EPS) {
		// compute k2
		for (int i = 0; i < THNR; ++i) {
			k2[i] = tau/2.0*a[i];
		}
		// compute S*k2
		MVprod(SPNR, THNR, StoichMatrix, k2, Sk2);
		// compute F
		for (int i = 0; i < SPNR; ++i) {
			F[i] = -(X[i] - X0[i] - Sk1[i] - Sk2[i]); // evaluated at X'_n
		}


		double dX[SPNR];
#ifdef INVJ
		double invJF[SPNR];
		double  invJ[SPNR*SPNR];
		fill_n(invJF, SPNR, 0.0);
		fill_n(invJ, SPNR*SPNR, 0.0);
		inverseJ(X, Theta, &tau, invJ);
		MVprod(SPNR, SPNR, invJ, F, invJF);
		copy(invJF, invJF+SPNR, dX);
#else
		// F'
		double dFdX[SPNR*SPNR];
		fill_n(dFdX, SPNR*SPNR, 0.0);
		J(X, Theta, &tau, dFdX);
		Eigen::Map<SPMatrix> map(dFdX);
		SPVector b(F);
		SPMatrix A = map;
		SPVector y = A.lu().solve(b);
		for (int i = 0; i < SPNR; ++i) {
			dX[i] = y[i];
		}
#endif

		// piv LU solver
//		Eigen::FullPivLU<SPMatrix> lu(A);
//		SPVector y = lu.solve(b);

		// compute invJ * F at X'_n
		// comput next iteration of Xnew
		deltaX = 0.0;
		for (int i = 0; i < SPNR; ++i) {
			X[i] += dX[i];
			deltaX += pow(dX[i], 2.0);
		}
		// new propensity at X_n+1
		computePropensities(X, Theta, 0, a_ptr);
		deltaX = sqrt(deltaX);
	}

	// get k2 from the final state
	for (int i = 0; i < THNR; ++i) {
		k2[i] = tau/2.0*a[i];
	}

	// approximate integral of G using trapezoidal approximation
	for (int j = 0; j < THNR; ++j) {
		G[j] += (a0[j] + a[j])*tau/2.0/Theta[j];
	}

	// update R
	for (int j = 0; j < THNR; j++) {
		R[j] += static_cast<int>(round(k1[j] + k2[j]));
	}

}

void MVprod(const int N,const  int M, const double *A, const double *B, double *Z) {
	for (int i = 0; i < N; ++i) {
		Z[i] = 0.0;
		for (int j = 0; j < M; ++j) {
			Z[i] += A[i*M+j]*B[j];
		}
	}
}

//void MVprod(const int N,const  int M,const  int *A,const  int *B, int *Z) {
//	for (int i = 0; i < N; ++i) {
//		Z[i] = 0.0;
//		for (int j = 0; j < M; ++j) {
//			Z[i] += A[i*M+j]*B[j];
//		}
//	}
//}

void MVprod(const int N,const  int M, const int *A, const double *B, double *Z) {
	for (int i = 0; i < N; ++i) {
		Z[i] = 0.0;
		for (int j = 0; j < M; ++j) {
			Z[i] += static_cast<double>(A[i*M+j])*B[j];
		}
	}
}

void computeG(const double *X, int *StoichMatrix, int *StoichMatrix_ed, double *g) {
	// determine highest order reaction of each species
	static bool computeOrders = true;
	static int HOR[SPNR], HEO[SPNR], RO[THNR];

	if (computeOrders) {

		fill_n(RO, THNR, 0);
		fill_n(HOR, SPNR, 0);
		fill_n(HEO, SPNR, 0);

		// compute reaction orders
		for (int j = 0; j < THNR; ++j) {
			for (int i = 0; i < SPNR; ++i) {
				RO[j] += Sed(i,j) != 0 ? 1:0;
			}
		}

		for (int i = 0; i < SPNR; ++i) {
			HOR[i] = 0;
			for (int j = 0; j < THNR; ++j) {
				HOR[i] = max(HOR[i], (Sed(i,j) != 0) ? RO[j] : 0);
			}
		}

		// highest order educt reaction
		for (int i = 0; i < SPNR; ++i) {
			for (int j = 0; j < THNR; ++j) {
				HEO[i] = max(HEO[i], abs(Sed(i,j)));
			}
		}
		computeOrders = false;
	}

	for (int i = 0; i < SPNR; ++i) {
		if (HOR[i] == 1)
			g[i] = 1;
		else if (HOR[i] == 2 && HEO[i] < 2)
			g[i] = 2;
		else if (HOR[i] == 2 && HEO[i] == 2)
			g[i] = 2 + (1.0/(X[i]-1.0));
		else if (HOR[i] == 3 && HEO[i] < 2)
			g[i] = 3;
		else if (HOR[i] == 3 && HEO[i] == 2)
			g[i] = 1.5*(2.0+1.0/(X[i]-1.0));
		else if (HOR[i] == 3 && HEO[i] == 3)
			g[i] = 3.0 + 1.0/(X[i]-1.0) + 2.0/(X[i]-2.0);
	}
}


/*
void inverseJ(const double *X, const double *Theta, const double *tau, double *invJ) {
	// compute the inverse of J for the Newton iteration to compute the update for the implicit tau-leaping scheme

	// get J from Matlab computation
	J(X, Theta, tau, J_mem);
	int k=0;

//	alglib::real_2d_array m_J;
//	m_J.setlength(SPNR, SPNR);
//	m_J.setcontent(SPNR, SPNR, J_mem);

	// get inverse
//	alglib::ae_int_t info;
//	alglib::matinvreport rep;


//	alglib::rmatrixinverse(m_J, info, rep);

	// use Eigen
	// create matrix with J
	Eigen::Map<Eigen::Matrix2d> J_map(J_mem);

	// solve inv(F')(x_n+1 - x_n) = -F(x)





	// m_J contains the inverse
//	int k=0;
	k=0;
	for (int i = 0; i < SPNR; ++i) {
		for (int i2 = 0; i2 < SPNR; ++i2)  {
			invJ[k] = m_J[i][i2];
			k++;
		}
	}
//	copy(m_J[0], m_J[0]+SPNR*SPNR*m_J.getstride(), invJ);

}
*/
