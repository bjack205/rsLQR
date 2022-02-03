#include "test/minunit.h"
#include "solver.h"
#include "test/test_problem.h"

mu_test_init

int NewSolver() {
  NdLqrSolver* solver = ndlqr_NewNdLqrSolver(6, 3, 8);
  int level = ndlqr_GetIndexLevel(&solver->tree, 0);
  mu_assert(level == 0);
  level = ndlqr_GetIndexLevel(&solver->tree, 3);
  mu_assert(level == 2);
  ndlqr_FreeNdLqrSolver(solver);
  return 1;
}

int InitializeWithLQRProblem() {
  LQRProblem* lqrprob = ndlqr_ReadTestLQRProblem();
  NdLqrSolver* solver = ndlqr_GenTestSolver();
  mu_assert(solver->nvars == 117);
  int nstates = solver->nstates;
  int ninputs = solver->ninputs; 

  // Create an identity matrix
  Matrix minus_identity = NewMatrix(nstates, nstates);
  MatrixSetConst(&minus_identity, 0);
  for (int i = 0; i < nstates; ++i) {
    MatrixSetElement(&minus_identity, i, i, -1);
  }

  NdData* nddata = solver->data;
  NdData* rhs = solver->soln;
  NdFactor* C;
  NdFactor* z;
  Matrix A = {nstates, nstates, NULL};
  Matrix B = {nstates, ninputs, NULL};
  Matrix d = {nstates, 1, NULL};
  Matrix q = {nstates, 1, NULL};
  Matrix r = {ninputs, 1, NULL};
  Matrix x0 = {nstates, 1, lqrprob->x0};
  MatrixScaleByConst(&x0, -1);
  Matrix Cx, Cu, yd, yx, yu;
  Matrix At = NewMatrix(nstates, nstates);
  Matrix Bt = NewMatrix(ninputs, nstates);

  ndlqr_GetNdFactor(nddata, 0, 0, &C);
  ndlqr_GetNdFactor(rhs, 0, 0, &z);
  Cx = ndlqr_GetStateFactor(C);
  Cu = ndlqr_GetInputFactor(C);
  yd = ndlqr_GetLambdaFactor(z);
  yx = ndlqr_GetStateFactor(z);
  yu = ndlqr_GetInputFactor(z);

  A.data = lqrprob->lqrdata[0]->A;
  B.data = lqrprob->lqrdata[0]->B;
  q.data = lqrprob->lqrdata[0]->q;
  r.data = lqrprob->lqrdata[0]->r;
  MatrixCopyTranspose(&At, &A);
  MatrixCopyTranspose(&Bt, &B);
  MatrixScaleByConst(&q, -1);
  MatrixScaleByConst(&r, -1);

  mu_assert(MatrixNormedDifference(&At, &Cx) < 1e-6);
  mu_assert(MatrixNormedDifference(&Bt, &Cu) < 1e-6);
  mu_assert(MatrixNormedDifference(&yd, &x0) < 1e-6);
  mu_assert(MatrixNormedDifference(&yx, &q) < 1e-6);
  mu_assert(MatrixNormedDifference(&yu, &r) < 1e-6);

  // Check the -I
  ndlqr_GetNdFactor(nddata, 1, 0, &C);
  Cx = ndlqr_GetStateFactor(C);
  mu_assert(MatrixNormedDifference(&minus_identity, &Cx) < 1e-6);

  // Check the next step
  ndlqr_GetNdFactor(nddata, 1, 1, &C);
  ndlqr_GetNdFactor(rhs, 1, 0, &z);
  Cx = ndlqr_GetStateFactor(C);
  Cu = ndlqr_GetInputFactor(C);
  yd = ndlqr_GetLambdaFactor(z);
  yx = ndlqr_GetStateFactor(z);
  yu = ndlqr_GetInputFactor(z);

  A.data = lqrprob->lqrdata[1]->A;
  B.data = lqrprob->lqrdata[1]->B;
  d.data = lqrprob->lqrdata[0]->d;
  q.data = lqrprob->lqrdata[1]->q;
  r.data = lqrprob->lqrdata[1]->r;
  MatrixCopyTranspose(&At, &A);
  MatrixCopyTranspose(&Bt, &B);
  MatrixScaleByConst(&d, -1);
  MatrixScaleByConst(&q, -1);
  MatrixScaleByConst(&r, -1);

  mu_assert(MatrixNormedDifference(&At, &Cx) < 1e-6);
  mu_assert(MatrixNormedDifference(&Bt, &Cu) < 1e-6);
  mu_assert(MatrixNormedDifference(&yd, &d) < 1e-6);
  mu_assert(MatrixNormedDifference(&yx, &q) < 1e-6);
  mu_assert(MatrixNormedDifference(&yu, &r) < 1e-6);

  // Check the next step, make sure it goes back to level 0
  ndlqr_GetNdFactor(nddata, 2, 0, &C);
  ndlqr_GetNdFactor(rhs, 2, 0, &z);
  Cx = ndlqr_GetStateFactor(C);
  Cu = ndlqr_GetInputFactor(C);
  yd = ndlqr_GetLambdaFactor(z);
  yx = ndlqr_GetStateFactor(z);
  yu = ndlqr_GetInputFactor(z);

  A.data = lqrprob->lqrdata[2]->A;
  B.data = lqrprob->lqrdata[2]->B;
  d.data = lqrprob->lqrdata[1]->d;
  q.data = lqrprob->lqrdata[2]->q;
  r.data = lqrprob->lqrdata[2]->r;
  MatrixCopyTranspose(&At, &A);
  MatrixCopyTranspose(&Bt, &B);
  MatrixScaleByConst(&d, -1);
  MatrixScaleByConst(&q, -1);
  MatrixScaleByConst(&r, -1);

  mu_assert(MatrixNormedDifference(&At, &Cx) < 1e-6);
  mu_assert(MatrixNormedDifference(&Bt, &Cu) < 1e-6);
  mu_assert(MatrixNormedDifference(&yd, &d) < 1e-6);
  mu_assert(MatrixNormedDifference(&yx, &q) < 1e-6);
  mu_assert(MatrixNormedDifference(&yu, &r) < 1e-6);

  // Check the terminal knot point
  ndlqr_GetNdFactor(nddata, 7, 0, &C);
  ndlqr_GetNdFactor(rhs, 7, 0, &z);
  Cx = ndlqr_GetStateFactor(C);
  yx = ndlqr_GetStateFactor(z);
  q.data = lqrprob->lqrdata[7]->q;
  MatrixScaleByConst(&q, -1);
  mu_assert(MatrixNormedDifference(&minus_identity, &Cx) < 1e-6);
  mu_assert(MatrixNormedDifference(&yx, &q) < 1e-6);

  int k = solver->nhorizon - 1;
  for (int i = 0; i < nstates; ++i) {
    mu_assert(fabs(*MatrixGetElement(&solver->diagonals[2 * k], i, i) - lqrprob->lqrdata[k]->Q[i]) < 1e-6);
  }

  FreeMatrix(&minus_identity);
  FreeMatrix(&At);
  FreeMatrix(&Bt);
  ndlqr_FreeNdLqrSolver(solver);
  ndlqr_FreeLQRProblem(lqrprob);
  return 1;
}

void AllTests() {
  mu_run_test(NewSolver);
  mu_run_test(InitializeWithLQRProblem);
}
 
mu_test_main
