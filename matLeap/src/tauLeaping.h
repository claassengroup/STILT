/*
 * stlTauLeaping.h
 *
 *  Created on: Jun 30, 2016
 *      Author: fjustin
 */

#ifndef TAULEAPING_H_
#define TAULEAPING_H_

#include "math.h"
#include <algorithm>
#include <cmath>
#include <utility>
#include <array>
#include <vector>
#include <list>
#include <chrono>
#include <iostream>
#include "model_def.h"
#include <fstream>
#include "Eigen/Core"
#include "Eigen/LU"

using namespace std;

struct CALL_STATS {
	int ssaCalls;
	int expCalls;
	int impCalls;

	CALL_STATS() {
		ssaCalls = 0;
		expCalls = 0;
		impCalls = 0;
	}
};

class OUTPUT_HANDLE {
public:
	std::ofstream X;
	std::ofstream R;
	std::ofstream G;

	OUTPUT_HANDLE() {

	}
	~OUTPUT_HANDLE() {
		X.close();
		R.close();
		G.close();
	}
};

struct RUN_OPTIONS {

	double DELTA_PE;
	int N_STIFF;
	double NEWTON_EPS;
	int MAX_ITER;
	int SSA_CONDITION_NUMBER;
	int CRITICAL_NUMBER;
	bool USE_IMPLICIT;
	bool SSA_ONLY;
	int N_SSA_1;
	int N_SSA_2;
	double EPSILON;
	int N_C;

	RUN_OPTIONS() {
		DELTA_PE = 0.05;
		N_STIFF  = 100;
		NEWTON_EPS = 0.1;
		MAX_ITER  = 30;
		SSA_CONDITION_NUMBER = 10;
		CRITICAL_NUMBER = 10;
		USE_IMPLICIT = true;
		SSA_ONLY = false;
		N_SSA_1 = 100;
		N_SSA_2 = 10;
		EPSILON = 0.03;
		N_C = 10;
	}

	template <typename T>
	void set(std::string field, T value) {
		if (field.compare("DELTA_PE") == 0) {
			DELTA_PE = value;
		} else if(field.compare("N_STIFF") == 0) {
			N_STIFF = value;
		} else if(field.compare("NEWTON_EPS") == 0) {
			NEWTON_EPS = value;
		} else if(field.compare("MAX_ITER") == 0) {
			MAX_ITER = value;
		} else if(field.compare("SSA_CONDITION == 0_NUMBER") == 0) {
			SSA_CONDITION_NUMBER = value;
		} else if(field.compare("CRITICAL_NUMBER") == 0) {
			CRITICAL_NUMBER = value;
		} else if(field.compare("USE_IMPLICIT") == 0) {
			USE_IMPLICIT = value;
		} else if(field.compare("N_SSA_1") == 0) {
			N_SSA_1 = value;
		} else if(field.compare("N_SSA_2") == 0) {
			N_SSA_2 = value;
		} else if(field.compare("EPSILON") == 0) {
			EPSILON = value;
		} else if(field.compare("N_C") == 0) {
			N_C = value;
		} else if(field.compare("SSA_ONLY") == 0) {
			SSA_ONLY = value;
		} else {
//			cout << "error!" << endl;
		}
	}

	void print() {
		cout << DELTA_PE << endl << N_STIFF << endl << NEWTON_EPS << endl
				<< MAX_ITER << endl << SSA_CONDITION_NUMBER << endl << CRITICAL_NUMBER << endl
				<< USE_IMPLICIT << endl << N_SSA_1 << endl << N_SSA_2 << endl << EPSILON  << endl
				<< N_C << endl << SSA_ONLY << endl;
	}
};

void computePropensities(const double *X, const double *theta, const int lastReacIdx, double **propensities);
void computeStoichiometricMatrix(int *S, int *Sed);

void performSSA(double *X, const double *Theta, const int *StoichMatrix, const double dt, const int SSA_steps, double &t, double *R, double *G);
void performStiffTau(const double *X, const double *Theta, const int *StoichMatrix, const bool *critical, const double tau, int *R, double *G, const RUN_OPTIONS run_options);
void J(const double *in1, const double *in2, const double *in3 , double *p_J);
void inverseJ(const double *X, const double *Theta, const double *tau, double *J);

void MVprod(const int N, const int M,const   double *A, const double *B, double *Z);
void MVprod(const int N, const int M,const  int *A,const  int *B, int *Z);
void MVprod(const int N, const int M,const  int *A,const  double *B, double *Z);

void computeG(const double *X, int *StoichMatrix, int *StoichMatrix_ed, double *g);
CALL_STATS tauLeapingInterval (const double *X0, const double *R0, const double *G0,
		const double *Theta, const double dt, RUN_OPTIONS run_options, double *X, double *R, double *G);
CALL_STATS tauLeaping (const double *X0, const double *Theta, int NIntervals, const double dt,
		RUN_OPTIONS run_options, double *X, double *R, double *G, bool writeOutput, std::vector<OUTPUT_HANDLE> &output_handles);

typedef Eigen::Matrix<double, SPNR, SPNR, Eigen::RowMajor> SPMatrix;
typedef Eigen::Matrix<double, SPNR, 1> SPVector;

#endif /* TAULEAPING_H_ */
