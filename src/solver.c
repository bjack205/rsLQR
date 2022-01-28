#include "solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "omp.h"
#include "linalg_utils.h"
#include "utils.h"

NdLqrProfile ndlqr_NewNdLqrProfile() {
  NdLqrProfile prof = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1};
  return prof;
}

void ndlqr_ResetProfile(NdLqrProfile *prof) {
  prof->t_total_ms = 0.0;
  prof->t_leaves_ms = 0.0;
  prof->t_products_ms = 0.0;
  prof->t_cholesky_ms = 0.0;
  prof->t_cholsolve_ms = 0.0;
  prof->t_shur_ms = 0.0;
}

void ndlqr_CopyProfile(NdLqrProfile *dest, NdLqrProfile *src) {
  dest->num_threads = src->num_threads;
  dest->t_total_ms = src->t_total_ms;
  dest->t_leaves_ms = src->t_leaves_ms;
  dest->t_products_ms = src->t_products_ms;
  dest->t_cholesky_ms = src->t_cholesky_ms;
  dest->t_cholsolve_ms = src->t_cholsolve_ms;
  dest->t_shur_ms = src->t_shur_ms;
}

void ndlqr_PrintProfile(NdLqrProfile* profile) {
  printf("Solved with %d threads\n", profile->num_threads);
  printf("Solve Total:    %.3f ms\n", profile->t_total_ms);
  printf("Solve Leaves:   %.3f ms\n", profile->t_leaves_ms);
  printf("Solve Products: %.3f ms\n", profile->t_products_ms);
  printf("Solve Cholesky: %.3f ms\n", profile->t_cholesky_ms);
  printf("Solve Solve:    %.3f ms\n", profile->t_cholsolve_ms);
  printf("Solve Shur:     %.3f ms\n", profile->t_shur_ms);
}

void PrintComp(double base, double new) {
  printf("%.3f / %.3f (%.2f speedup)\n", base, new, base / new);
}

void ndlqr_CompareProfile(NdLqrProfile* base, NdLqrProfile* prof) {
  printf("Num Threads:     %d / %d\n", base->num_threads, prof->num_threads);
  printf("Solve Total:     "); PrintComp(base->t_total_ms, prof->t_total_ms);
  printf("Solve Leaves:    "); PrintComp(base->t_leaves_ms, prof->t_leaves_ms);
  printf("Solve Products:  "); PrintComp(base->t_products_ms, prof->t_products_ms);
  printf("Solve Cholesky:  "); PrintComp(base->t_cholesky_ms, prof->t_cholesky_ms);
  printf("Solve CholSolve: "); PrintComp(base->t_cholsolve_ms, prof->t_cholsolve_ms);
  printf("Solve Shur Comp: "); PrintComp(base->t_shur_ms, prof->t_shur_ms);
}

NdLqrSolver *ndlqr_NewNdLqrSolver(int nstates, int ninputs, int nhorizon) {
  OrderedBinaryTree tree = ndlqr_BuildTree(nhorizon);
  NdLqrSolver *solver = (NdLqrSolver *)malloc(sizeof(NdLqrSolver));
  int nvars = (2 * nstates + ninputs) * nhorizon - ninputs;

  int diag_size = (nstates * nstates + ninputs * ninputs) * nhorizon;
  double *diag_data = (double *)malloc(diag_size * sizeof(double));
  Matrix *diagonals = (Matrix *)malloc(2 * nhorizon * sizeof(Matrix));
  for (int k = 0; k < nhorizon; ++k)
  {
    int blocksize = nstates * nstates + ninputs * ninputs;
    diagonals[2 * k].rows = nstates;
    diagonals[2 * k].cols = nstates;
    diagonals[2 * k].data = diag_data + k * blocksize;
    diagonals[2 * k + 1].rows = ninputs;
    diagonals[2 * k + 1].cols = ninputs;
    diagonals[2 * k + 1].data = diag_data + k * blocksize + nstates * nstates;
  }
  NdLqrCholeskyFactors *cholfacts = ndlqr_NewCholeskyFactors(tree.depth, nhorizon);

  solver->nstates = nstates;
  solver->ninputs = ninputs;
  solver->nhorizon = nhorizon;
  solver->depth = tree.depth;
  solver->nvars = nvars;
  solver->tree = tree;
  solver->diagonals = diagonals;
  solver->data = ndlqr_NewNdData(nstates, ninputs, nhorizon - 1, nstates);
  solver->fact = ndlqr_NewNdData(nstates, ninputs, nhorizon - 1, nstates);
  solver->soln = ndlqr_NewNdData(nstates, ninputs, nhorizon - 1, 1);
  solver->cholfacts = cholfacts;
  solver->solve_time_ms = 0.0;
  solver->linalg_time_ms = 0.0;
  solver->profile = ndlqr_NewNdLqrProfile();
  solver->num_threads = omp_get_num_procs() / 2;
  return solver;
}

void ndlqr_ResetSolver(NdLqrSolver *solver) {
  ndlqr_ResetNdData(solver->data);
  ndlqr_ResetNdData(solver->fact);
  ndlqr_ResetNdData(solver->soln);
  ndlqr_ResetProfile(&solver->profile);
  for (int i = 0; i < 2 * solver->nhorizon; ++i)
  {
    MatrixSetConst(&solver->diagonals[i], 0.0);
  }
}

int ndlqr_FreeNdLqrSolver(NdLqrSolver *solver) {
  if (!solver)
    return -1;
  ndlqr_FreeTree(&(solver->tree));
  ndlqr_FreeNdData(solver->data);
  ndlqr_FreeNdData(solver->fact);
  ndlqr_FreeNdData(solver->soln);
  ndlqr_FreeCholeskyFactors(solver->cholfacts);
  free(solver->diagonals[0].data);
  free(solver->diagonals);
  free(solver);
  return 0;
}

int ndlqr_InitializeWithLQRProblem(const LQRProblem *lqrprob,
                                   NdLqrSolver *solver) {
  int nstates = solver->nstates;
  int ninputs = solver->ninputs;

  // Create a minux identity matrix for copying into the original matrix
  Matrix minus_identity = NewMatrix(nstates, nstates);
  MatrixSetConst(&minus_identity, 0);
  for (int i = 0; i < nstates; ++i)
  {
    MatrixSetElement(&minus_identity, i, i, -1);
  }

  // Loop over the knot points, copying the LQR data into the matrix data
  // and populating the right-hand-side vector
  NdFactor *Cfactor;
  NdFactor *zfactor;
  ndlqr_GetNdFactor(solver->soln, 0, 0, &zfactor);
  memcpy(zfactor->lambda.data, lqrprob->x0, nstates * sizeof(double));
  int k;
  for (k = 0; k < solver->nhorizon - 1; ++k)
  {
    // Copy data into C factors and rhs vector from LQR data
    int level = ndlqr_GetIndexLevel(&(solver->tree), k);
    ndlqr_GetNdFactor(solver->data, k, level, &Cfactor);
    ndlqr_GetNdFactor(solver->soln, k, 0, &zfactor);
    Matrix A = {nstates, nstates, lqrprob->lqrdata[k]->A};
    Matrix B = {nstates, ninputs, lqrprob->lqrdata[k]->B};
    MatrixCopyTranspose(&Cfactor->state, &A);
    MatrixCopyTranspose(&Cfactor->input, &B);
    // memcpy(Cfactor->state.data, lqrprob->lqrdata[k]->A, nstates * nstates *
    // sizeof(double)); memcpy(Cfactor->input.data, lqrprob->lqrdata[k]->B,
    // nstates * ninputs * sizeof(double));
    memcpy(zfactor->state.data, lqrprob->lqrdata[k]->q,
           nstates * sizeof(double));
    memcpy(zfactor->input.data, lqrprob->lqrdata[k]->r,
           ninputs * sizeof(double));

    // Copy Q and R into diagonals
    Matrix Q = solver->diagonals[2 * k];
    Matrix R = solver->diagonals[2 * k + 1];
    MatrixSetConst(&Q, 0);
    MatrixSetConst(&R, 0);
    for (int i = 0; i < nstates; ++i)
    {
      MatrixSetElement(&Q, i, i, lqrprob->lqrdata[k]->Q[i]);
    }
    for (int i = 0; i < ninputs; ++i)
    {
      MatrixSetElement(&R, i, i, lqrprob->lqrdata[k]->R[i]);
    }

    // Next time step
    ndlqr_GetNdFactor(solver->data, k + 1, level, &Cfactor);
    ndlqr_GetNdFactor(solver->soln, k + 1, 0, &zfactor);
    memcpy(Cfactor->state.data, minus_identity.data,
           nstates * nstates * sizeof(double));
    MatrixSetConst(&Cfactor->input, 0.0);
    memcpy(zfactor->lambda.data, lqrprob->lqrdata[k]->d,
           nstates * sizeof(double));
  }

  // Terminal step
  memcpy(zfactor->state.data, lqrprob->lqrdata[k]->q, nstates * sizeof(double));
  Matrix Q = solver->diagonals[2 * k];
  MatrixSetConst(&Q, 0);
  for (int i = 0; i < nstates; ++i)
  {
    MatrixSetElement(&Q, i, i, lqrprob->lqrdata[k]->Q[i]);
  }

  // Negate the entire rhs vector
  for (int i = 0; i < solver->nvars; ++i)
  {
    solver->soln->data[i] *= -1;
  }

  FreeMatrix(&minus_identity);
  return 0;
}

void ndlqr_PrintSolveSummary(NdLqrSolver *solver) {
  printf("Solve time:  %f ms\n", solver->solve_time_ms);
  if (kMatrixLinearAlgebraTimingEnabled)
  {
    printf("LinAlg time: %f ms (%.1f%% of total)\n", solver->linalg_time_ms,
           100.0 * solver->linalg_time_ms / solver->solve_time_ms);
  }
  printf("Solved with %d threads.\n", solver->num_threads);
}

int ndlqr_GetNumVars(NdLqrSolver *solver) { return solver->nvars; }

int ndlqr_SetNumThreads(NdLqrSolver *solver, int num_threads) {
  if (!solver)
    return -1;
  solver->num_threads = num_threads;
  return 0;
}

int ndlqr_GetNumThreads(NdLqrSolver *solver) {
  if (!solver)
    return -1;
  return solver->num_threads;
}

int ndlqr_PrintSolveProfile(NdLqrSolver* solver) { 
  if (!solver) return -1;
  ndlqr_PrintProfile(&solver->profile); 
  return 0;
}

NdLqrProfile ndlqr_GetProfile(NdLqrSolver* solver) { return solver->profile; }
