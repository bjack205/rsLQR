#pragma once

#include "binary_tree.h"
#include "linalg.h"
#include "lqr_problem.h"
#include "nddata.h"
#include "cholesky_factors.h"

typedef struct {
  double t_total_ms;
  double t_leaves_ms;
  double t_products_ms;
  double t_cholesky_ms;
  double t_cholsolve_ms;
  double t_shur_ms;
  int num_threads;
} NdLqrProfile;

NdLqrProfile ndlqr_NewNdLqrProfile();
void ndlqr_ResetProfile(NdLqrProfile* prof);
void ndlqr_CopyProfile(NdLqrProfile *dest, NdLqrProfile *src);
void ndlqr_PrintProfile(NdLqrProfile* profile);
void ndlqr_CompareProfile(NdLqrProfile* base, NdLqrProfile* prof);

typedef struct {
  int nstates;
  int ninputs;
  int nhorizon;
  int depth;
  int nvars;  // number of decision variables (size of the linear system)
  OrderedBinaryTree tree;
  Matrix* diagonals;  // (nhorizon,2) array of diagonal blocks (Q,R)
  NdData* data;       // original matrix data
  NdData* fact;       // factorization
  NdData* soln;       // solution vector (also the initial RHS)
  NdLqrCholeskyFactors* cholfacts;
  double solve_time_ms;
  double linalg_time_ms;
  NdLqrProfile profile;
  int num_threads;
} NdLqrSolver;

NdLqrSolver* ndlqr_NewNdLqrSolver(int nstates, int ninputs, int nhorizon);

int ndlqr_FreeNdLqrSolver(NdLqrSolver* solver);

int ndlqr_InitializeWithLQRProblem(const LQRProblem* lqrprob,
                                   NdLqrSolver* solver);

void ndlqr_ResetSolver(NdLqrSolver* solver);

void ndlqr_PrintSolveSummary(NdLqrSolver* solver);

int ndlqr_GetNumVars(NdLqrSolver* solver);

int ndlqr_SetNumThreads(NdLqrSolver* solver, int num_threads);

int ndlqr_GetNumThreads(NdLqrSolver* solver);

int ndlqr_PrintSolveProfile(NdLqrSolver* solver);

NdLqrProfile ndlqr_GetProfile(NdLqrSolver* solver);
