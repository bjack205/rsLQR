#pragma once

#include "lqr_data.h"

typedef struct {
  int nhorizon;
  double* x0;
  LQRData** lqrdata;
} LQRProblem;

/**
 * @brief Initialize the problem with an initial state and the LQR data
 * 
 * @param lqrproblem 
 * @param x0          Initial state vector. The data is copied into the problem.
 * @param lqrdata     A vector of LQR data. Each element is copied into the problem. 
 * @return int 
 */
int ndlqr_InitializeLQRProblem(LQRProblem* lqrproblem, double* x0, LQRData** lqrdata);
LQRProblem* ndlqr_NewLQRProblem(int nstates, int ninputs, int nhorizon);
int ndlqr_FreeLQRProblem(LQRProblem* lqrprob);
