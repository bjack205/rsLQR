#include "riccati_solver.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

RiccatiSolver* ndlqr_NewRiccatiSolver(LQRProblem* lqrprob) {
  int nhorizon = lqrprob->nhorizon;
  int nstates = lqrprob->lqrdata[0]->nstates;
  int ninputs = lqrprob->lqrdata[0]->ninputs;
  int nvars = (2 * nstates + ninputs) * nhorizon - ninputs;

  int dim_K = ninputs * nstates;
  int dim_d = ninputs;
  int dim_P = nstates * nstates;
  int dim_p = nstates;
  int dim_X = nstates;
  int dim_U = ninputs;
  int dim_Y = nstates;
  int dim_Q = nstates * nstates + ninputs * ninputs + ninputs * nstates + nstates + ninputs; 
  int len_Q = 2; 

  int num_K = dim_K * (nhorizon - 1);
  int num_d = dim_d * (nhorizon - 1);
  int num_P = dim_P * nhorizon;
  int num_p = dim_p * nhorizon;
  int num_X = dim_X * nhorizon;
  int num_U = dim_U * (nhorizon - 1);
  int num_Y = dim_Y * nhorizon;
  int num_Q = dim_Q * len_Q;
  int total_size = num_K + num_d + num_P + num_p + num_X + num_U + num_Y + num_Q;
  double* data = (double*) malloc(total_size * sizeof(double));
  if (!data) return NULL;
  memset(data, 0, total_size * sizeof(double));

  RiccatiSolver* solver = (RiccatiSolver*) malloc(sizeof(RiccatiSolver));
  if (!solver) {
    free(data);
  }

  Matrix* K = (Matrix*) malloc((nhorizon - 1) * sizeof(Matrix));
  Matrix* d = (Matrix*) malloc((nhorizon - 1) * sizeof(Matrix));
  Matrix* P = (Matrix*) malloc(nhorizon * sizeof(Matrix));
  Matrix* p = (Matrix*) malloc(nhorizon * sizeof(Matrix));
  Matrix* X = (Matrix*) malloc(nhorizon * sizeof(Matrix));
  Matrix* U = (Matrix*) malloc((nhorizon - 1) * sizeof(Matrix));
  Matrix* Y = (Matrix*) malloc(nhorizon * sizeof(Matrix));
  int offset = 0;
  for (int k = 0; k < nhorizon; ++k) {
    P[k].rows = nstates;
    P[k].cols = nstates;
    p[k].rows = nstates;
    p[k].cols = 1;
    X[k].rows = nstates;
    X[k].cols = 1;
    Y[k].rows = nstates;
    Y[k].cols = 1;

    P[k].data = data + offset; offset += nstates * nstates;
    p[k].data = data + offset; offset += nstates;

    if (k < nhorizon - 1) {
      K[k].rows = ninputs;
      K[k].cols = nstates;
      d[k].rows = ninputs;
      d[k].cols = 1;
      U[k].rows = ninputs;
      U[k].cols = 1;

      K[k].data = data + offset; offset += ninputs * nstates;
      d[k].data = data + offset; offset += ninputs;
    }
  }

  // Place the solution vector [y0, x0, u0, y1, x1, u1, ..., yn, xn] in contiguous memory
  for (int k = 0; k < nhorizon; ++k) {
    Y[k].data = data + offset; offset += nstates;
    X[k].data = data + offset; offset += nstates;
    if (k < nhorizon - 1) {
      U[k].data = data + offset; offset += ninputs; 
    }
  }

  // Initialize the temporary Q matrices
  offset = total_size - num_Q;
  Matrix* Q = (Matrix*) malloc(5 * len_Q * sizeof(Matrix));
  Matrix* Qx = Q + 0 * len_Q;
  Matrix* Qu = Q + 1 * len_Q;
  Matrix* Qxx = Q + 2 * len_Q;
  Matrix* Qux = Q + 3 * len_Q;
  Matrix* Quu = Q + 4 * len_Q;
  for (int i = 0; i < len_Q; ++i) {
    Qx[i].rows = nstates;
    Qx[i].cols = 1;
    Qu[i].rows = ninputs;
    Qu[i].cols = 1;
    Qxx[i].rows = nstates;
    Qxx[i].cols = nstates;
    Qux[i].rows = ninputs;
    Qux[i].cols = nstates;
    Quu[i].rows = ninputs;
    Quu[i].cols = ninputs;

    Qx[i].data = data + offset;
    Qu[i].data = data + offset + nstates;
    Qxx[i].data = data + offset + nstates + ninputs;
    Qux[i].data = data + offset + nstates + ninputs + nstates * nstates;
    Quu[i].data = data + offset + nstates + ninputs + nstates * nstates + ninputs * nstates;
    offset += dim_Q;
  }

  // Initialize the solver
  solver->prob = lqrprob;
  solver->nhorizon = nhorizon;
  solver->nstates = nstates;
  solver->ninputs = ninputs;
  solver->nvars = nvars;
  solver->data = data;
  solver->K = K;
  solver->d = d;
  solver->P = P;
  solver->p = p;
  solver->X = X;
  solver->U = U;
  solver->Y = Y;
  solver->t_solve_ms = 0.0;

  solver->Qx = Qx;
  solver->Qu = Qu;
  solver->Qxx = Qxx;
  solver->Qux = Qux;
  solver->Quu = Quu;
  return solver;
}

int ndlqr_FreeRiccatiSolver(RiccatiSolver* solver) {
  if (!solver) return -1;
  free(solver->data);
  free(solver->K);
  free(solver->d);
  free(solver->P);
  free(solver->p);
  free(solver->X);
  free(solver->U);
  free(solver->Y);
  free(solver->Qx);
  free(solver);
  solver = NULL;
  return 0;
}

int ndlqr_PrintRiccatiSummary(RiccatiSolver* solver) {
  if (!solver) return -1;
  printf("NDLQR Riccati Solve Summary\n");
  double t_solve = solver->t_solve_ms;
  double t_bp = solver->t_backward_pass_ms;
  double t_fp = solver->t_forward_pass_ms;
  printf("  Solve time:    %.2f ms\n", t_solve);
  printf("  Backward Pass: %.2f ms (%.1f %% of total)\n", t_bp, t_bp / t_solve * 100.0);
  printf("  Foward Pass:   %.2f ms (%.1f %% of total)\n", t_fp, t_fp / t_solve * 100.0);
  return 0;
}

Matrix ndlqr_GetRiccatiSolution(RiccatiSolver* solver) {
  if (!solver) {
    Matrix nullmat = {0, 0, NULL};
    return nullmat;
  }
  // This works because of the way the memory is laid out
  // The solution vector matches that of ndlqr:
  //   [y0, x0, u0, y1, x1, y1, .., yn, xn]
  Matrix soln = {solver->nvars, 1, solver->Y->data};
  return soln;
}

int ndlqr_CopyRiccatiSolution(RiccatiSolver* solver, double* soln) {
  if (!solver) return -1;
  memcpy(soln, solver->Y->data, solver->nvars * sizeof(double));
  return solver->nvars;
}

int ndlqr_GetRiccatiSolveTimes(RiccatiSolver* solver, double* t_solve, double* t_bp, double* t_fp) {
  if (!solver) return -1;
  *t_solve = solver->t_solve_ms;
  *t_bp = solver->t_backward_pass_ms;
  *t_fp = solver->t_forward_pass_ms;
  printf("Solve time: %f\n", *t_solve);
  return 0;
}
