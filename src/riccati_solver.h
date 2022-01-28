#pragma once

#include "matrix.h"
#include "lqr_problem.h"

typedef struct {
  LQRProblem* prob;
  int nhorizon;
  int nstates;
  int ninputs;
  int nvars;
  double* data;
  Matrix* K;
  Matrix* d;
  Matrix* P;
  Matrix* p;
  Matrix* X;
  Matrix* U;
  Matrix* Y;  // lagrange multiplier
  Matrix* Qx;
  Matrix* Qu;
  Matrix* Qxx;
  Matrix* Qux;
  Matrix* Quu;
  double t_solve_ms;
  double t_backward_pass_ms;
  double t_forward_pass_ms;
} RiccatiSolver;

RiccatiSolver* ndlqr_NewRiccatiSolver(LQRProblem* lqrprob);
int ndlqr_FreeRiccatiSolver(RiccatiSolver* solver);

int ndlqr_PrintRiccatiSummary(RiccatiSolver* solver);

Matrix ndlqr_GetRiccatiSolution(RiccatiSolver* solver);

int ndlqr_CopyRiccatiSolution(RiccatiSolver* solver, double* soln);

int ndlqr_GetRiccatiSolveTimes(RiccatiSolver* solver, double* t_solve, double* t_bp, double* t_fp);
