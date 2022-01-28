#include <time.h>

#include "ndlqr.h"

#include "solver.h"
#include "test/minunit.h"
#include "nested_dissection.h"
#include "test/test_problem.h"
#include "linalg.h"

int SolveLeaves() {
  LQRProblem* lqrprob = ndlqr_ReadTestLQRProblem();
  int nstates = lqrprob->lqrdata[0]->nstates;
  int ninputs = lqrprob->lqrdata[0]->ninputs;
  Matrix A = {nstates, nstates, NULL};
  Matrix B = {nstates, ninputs, NULL};
  Matrix Anull = NewMatrix(nstates, nstates);
  MatrixSetConst(&Anull, 0.0);
  NdLqrSolver* solver = ndlqr_GenTestSolver();
  ndlqr_SolveLeaves(solver);

  NdFactor* C;
  NdFactor* F;
  NdFactor* z;

  // Check first time step
  int k = 0;
  A.data = lqrprob->lqrdata[0]->A;
  B.data = lqrprob->lqrdata[0]->B;
  Matrix At = NewMatrix(nstates, nstates);
  Matrix Bt = NewMatrix(ninputs, nstates);
  ndlqr_GetNdFactor(solver->data, k, 0, &C);
  ndlqr_GetNdFactor(solver->fact, k, 0, &F);
  ndlqr_GetNdFactor(solver->soln, k, 0, &z);

  Matrix Fy = ndlqr_GetLambdaFactor(F);
  Matrix Fx = ndlqr_GetStateFactor(F);
  Matrix Fu = ndlqr_GetInputFactor(F);
  MatrixScaleByConst(&A, -1);
  MatrixScaleByConst(&B, 1 / lqrprob->lqrdata[0]->R[0]);  // same as left-multiplying by inv(R)
  MatrixCopyTranspose(&At, &A);
  MatrixCopyTranspose(&Bt, &B);
  mu_assert(MatrixNormedDifference(&Fy, &At) < 1e-6);
  mu_assert(MatrixNormedDifference(&Fx, &Anull) < 1e-6);
  mu_assert(MatrixNormedDifference(&Fu, &Bt) < 1e-6);

  double zdata0[15] = {-1.0, -2.2, 1.6, -1.6, 4.2, -1.0, 1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 100.0, -0.0, -100.0};
  Matrix zans = {2 * nstates + ninputs, 1, zdata0};
  Matrix z0 = {2 * nstates + ninputs, 1, z->lambda.data};
  mu_assert(MatrixNormedDifference(&z0, &zans) < 1e-6);

  // Check the next time step
  k = 1;
  ndlqr_GetNdFactor(solver->fact, k, 1, &F);
  ndlqr_GetNdFactor(solver->soln, k, 0, &z);

  MatrixScaleByConst(&At, -1);
  mu_assert(MatrixNormedDifference(&F->state, &At) < 1e-6);
  mu_assert(MatrixNormedDifference(&F->input, &Bt) < 1e-6);

  ndlqr_GetNdFactor(solver->fact, k, 0, &F);
  Matrix* Q = &solver->diagonals[0];
  for (int i = 0; i < nstates; ++i) {
    double x = *MatrixGetElement(Q, i, i);
    MatrixSetElement(Q, i, i, -1 / (x * x));
  }
  mu_assert(MatrixNormedDifference(&F->state, Q) < 1e-6); // Check the Q \ -I  term

  double zdata1[15] = {-1.5, -1.5, -1.5, -1.5, -1.5, -1.5, 4.0, 2.4, 0.8, -0.8, -2.4, -4.0, 200.0, -0.0, -200.0};
  zans.data = zdata1;
  z0.data = z->lambda.data;
  mu_assert(MatrixNormedDifference(&z0, &zans) < 1e-6);

  // Check last index
  k = solver->nhorizon - 2;
  ndlqr_GetNdFactor(solver->fact, k, 0, &F);
  ndlqr_GetNdFactor(solver->soln, k, 0, &z);
  mu_assert(MatrixNormedDifference(&F->state, &At) < 1e-6);
  mu_assert(MatrixNormedDifference(&F->input, &Bt) < 1e-6);

  ndlqr_GetNdFactor(solver->fact, k, 1, &F);
  Q = &solver->diagonals[2 * k];
  for (int i = 0; i < nstates; ++i) {
    double x = *MatrixGetElement(Q, i, i);
    MatrixSetElement(Q, i, i, -1 / (x * x));
  }
  mu_assert(MatrixNormedDifference(&F->state, Q) < 1e-6); // Check the Q \ -I  term

  ndlqr_GetNdFactor(solver->fact, k+1, 0, &F);
  Q = &solver->diagonals[2 * (k + 1)];
  for (int i = 0; i < nstates; ++i) {
    double x = *MatrixGetElement(Q, i, i);
    MatrixSetElement(Q, i, i, -1 / (x * x));
  }
  // PrintMatrix(Q);
  // PrintMatrix(&F->state);
  mu_assert(MatrixNormedDifference(&F->state, Q) < 1e-6); // Check the Q \ -I  term

  // Check right-hand side vector
  Matrix b_ans = ReadMatrixJSONFile(SAMPLEPROBFILE, "b");
  Matrix b = {b_ans.rows, 1, solver->soln->data};
  mu_assert(MatrixNormedDifference(&b, &b_ans) < 1e-6);
  FreeMatrix(&b_ans);

  ndlqr_FreeNdLqrSolver(solver);
  ndlqr_FreeLQRProblem(lqrprob);
  FreeMatrix(&Anull);
  FreeMatrix(&At);
  FreeMatrix(&Bt);
  return 1;
}

int FactorInnerProduct() {
  NdLqrSolver* solver = ndlqr_GenTestSolver();
  ndlqr_SolveLeaves(solver);
  ndlqr_FactorInnerProduct(solver->data, solver->fact, 0, 0, 0);
  NdFactor* F; 
  ndlqr_GetNdFactor(solver->fact, 1, 0, &F);
  Matrix S = ndlqr_GetLambdaFactor(F);
  PrintMatrix(&S);
  double Sdata[36] = {1.0025, 0.0, 0.0, 0.05, 0.0, 0.0, 
                      0.0, 1.0025, 0.0, 0.0, 0.05, 0.0, 
                      0.0, 0.0, 1.0025, 0.0, 0.0, 0.05, 
                      0.05, 0.0, 0.0, 2.0, 0.0, 0.0, 
                      0.0, 0.05, 0.0, 0.0, 2.0, 0.0, 
                      0.0, 0.0, 0.05, 0.0, 0.0, 2.0};
  Matrix Sans = {6, 6, Sdata};
  mu_assert(MatrixNormedDifference(&S, &Sans) < 1e-6);
  ndlqr_FreeNdLqrSolver(solver);
  return 1;
}

int ShurCompliment() {
  NdLqrSolver* solver = ndlqr_GenTestSolver();
  ndlqr_SolveLeaves(solver);
  ndlqr_FactorInnerProduct(solver->data, solver->fact, 0, 0, 0);

  int index = 0;
  int level = 0;
  NdFactor* F;
  ndlqr_GetNdFactor(solver->fact, index + 1, level, &F);
  Matrix Sbar = F->lambda;
  MatrixCholeskyFactorize(&Sbar);

  int upper_level = level + 1;
  ndlqr_FactorInnerProduct(solver->data, solver->fact, index, level, upper_level);
  NdFactor* G;
  ndlqr_GetNdFactor(solver->fact, index + 1, upper_level, &G);
  Matrix f = G->lambda;
  double fdata[36] = {  // transposed since the data is colwise
   -1.0,   0.0,   0.0,  -0.1,   0.0,   0.0,
    0.0,  -1.0,   0.0,   0.0,  -0.1,   0.0,
    0.0,   0.0,  -1.0,   0.0,   0.0,  -0.1,
    0.0,   0.0,   0.0,  -1.0,   0.0,   0.0,
    0.0,   0.0,   0.0,   0.0,  -1.0,   0.0,
    0.0,   0.0,   0.0,   0.0,   0.0,  -1.0,
  };
  Matrix fans = {6, 6, fdata};
  mu_assert(MatrixNormedDifference(&fans, &f) < 1e-6);

  MatrixCholeskySolve(&Sbar, &f);
  double fdata2[36] = {  // transposed since the data is colwise
    -0.996255,   0.0,       0.0,       -0.0250936,  0.0,        0.0,
    0.0,       -0.996255,   0.0,        0.0,       -0.0250936,  0.0,
    0.0,       0.0,       -0.996255,   0.0,        0.0,       -0.0250936,
    0.0249688,  0.0,        0.0,       -0.500624,   0.0,        0.0,
    0.0,       0.0249688,  0.0,        0.0,       -0.500624,   0.0,
    0.0,        0.0,        0.0249688,  0.0,        0.0,       -0.500624,
  };
  fans.data = fdata2;
  mu_assert(MatrixNormedDifference(&fans, &f) < 1e-6);

  ndlqr_ComputeShurCompliment(solver, index, level, upper_level);
  NdFactor* e;
  ndlqr_GetNdFactor(solver->fact, 0, 1, &e);

  char name[5] = "E01y\0";
  for (int i = 0; i < 2; ++i) {
    ndlqr_GetNdFactor(solver->fact, i, upper_level, &e);

    name[1] =  i + '0';
    name[3] = 'y';
    Matrix F_lambda = ReadMatrixJSONFile(SAMPLEPROBFILE, name);
    name[3] = 'x';
    Matrix F_state = ReadMatrixJSONFile(SAMPLEPROBFILE, name);
    name[3] = 'u';
    Matrix F_input = ReadMatrixJSONFile(SAMPLEPROBFILE, name);
    mu_assert(MatrixNormedDifference(&e->lambda, &F_lambda) < 1e-6);
    mu_assert(MatrixNormedDifference(&e->state, &F_state) < 1e-6);
    mu_assert(MatrixNormedDifference(&e->input, &F_input) < 1e-6);
    FreeMatrix(&F_lambda);
    FreeMatrix(&F_state);
    FreeMatrix(&F_input);
  }

  // Next level 
  {
    ++upper_level;
    ndlqr_FactorInnerProduct(solver->data, solver->fact, index, level, upper_level);
    NdFactor* G;
    ndlqr_GetNdFactor(solver->fact, index + 1, upper_level, &G);
    Matrix f = G->lambda;

    MatrixCholeskySolve(&Sbar, &f);
    ndlqr_ComputeShurCompliment(solver, index, level, upper_level);
  }
  for (int i = 0; i < 2; ++i) {
    ndlqr_GetNdFactor(solver->fact, i, upper_level, &e);

    name[2] = upper_level + '0';
    name[1] =  i + '0';
    name[3] = 'y';
    Matrix F_lambda = ReadMatrixJSONFile(SAMPLEPROBFILE, name);
    name[3] = 'x';
    Matrix F_state = ReadMatrixJSONFile(SAMPLEPROBFILE, name);
    name[3] = 'u';
    Matrix F_input = ReadMatrixJSONFile(SAMPLEPROBFILE, name);
    mu_assert(MatrixNormedDifference(&e->lambda, &F_lambda) < 1e-6);
    mu_assert(MatrixNormedDifference(&e->state, &F_state) < 1e-6);
    mu_assert(MatrixNormedDifference(&e->input, &F_input) < 1e-6);
    FreeMatrix(&F_lambda);
    FreeMatrix(&F_state);
    FreeMatrix(&F_input);
  }

  ndlqr_FreeNdLqrSolver(solver);
  return 1;
}

int RunSolve() {
  NdLqrSolver* solver = ndlqr_GenTestSolver();

  clock_t start = clock();
  ndlqr_Solve(solver);
  ndlqr_SetNumThreads(solver, 1);
  clock_t diff = clock() - start;
  double msec = diff * 1000.0 / (double)CLOCKS_PER_SEC;
  printf("Solve time: %f milliseconds\n", msec);

  char name[5] = "F01y\0";
  int upper_level = solver->depth - 1;
  NdFactor* e;
  NdFactor* z;
  for (int k = 0; k < solver->nhorizon; ++k) {
    ndlqr_GetNdFactor(solver->fact, k, upper_level, &e);
    ndlqr_GetNdFactor(solver->soln, k, 0, &z);

    name[2] = upper_level + '0';
    name[1] =  k + '0';
    name[3] = 'y';
    Matrix F_lambda = ReadMatrixJSONFile(SAMPLEPROBFILE, name);
    name[3] = 'x';
    Matrix F_state = ReadMatrixJSONFile(SAMPLEPROBFILE, name);
    name[3] = 'u';
    Matrix F_input = ReadMatrixJSONFile(SAMPLEPROBFILE, name);
    // mu_assert(MatrixNormedDifference(&e->lambda, &F_lambda) < 1e-6);
    // mu_assert(MatrixNormedDifference(&e->state, &F_state) < 1e-6);
    // mu_assert(MatrixNormedDifference(&e->input, &F_input) < 1e-6);
    FreeMatrix(&F_lambda);
    FreeMatrix(&F_state);
    FreeMatrix(&F_input);
  }

  Matrix x_ans = ReadMatrixJSONFile(SAMPLEPROBFILE, "soln");
  Matrix x = {x_ans.rows, 1, solver->soln->data};
  double err = MatrixNormedDifference(&x, &x_ans);
  printf("Accuracy of final solution: %e\n", err);
  mu_assert(MatrixNormedDifference(&x, &x_ans) < 1e-6);
  FreeMatrix(&x_ans);
  ndlqr_PrintSolveSummary(solver);

  ndlqr_FreeNdLqrSolver(solver);
  return 1;
}

int SolveTwice() {
  LQRProblem* lqrprob = ndlqr_ReadTestLQRProblem();
  int nstates = lqrprob->lqrdata[0]->nstates;
  int ninputs = lqrprob->lqrdata[0]->ninputs;
  int nhorizon = lqrprob->nhorizon;


  NdLqrSolver* solver = ndlqr_NewNdLqrSolver(nstates, ninputs, nhorizon);
  ndlqr_InitializeWithLQRProblem(lqrprob, solver);
  ndlqr_Solve(solver);

  Matrix x_ans = ReadMatrixJSONFile(SAMPLEPROBFILE, "soln");
  Matrix x = {x_ans.rows, 1, solver->soln->data};
  double err = MatrixNormedDifference(&x, &x_ans);
  printf("Accuracy of final solution: %e\n", err);
  mu_assert(MatrixNormedDifference(&x, &x_ans) < 1e-6);

  ndlqr_ResetNdData(solver->fact);
  ndlqr_InitializeWithLQRProblem(lqrprob, solver);
  ndlqr_Solve(solver);

  err = MatrixNormedDifference(&x, &x_ans);
  printf("Accuracy of final solution (2nd solve): %e\n", err);
  mu_assert(MatrixNormedDifference(&x, &x_ans) < 1e-6);
  FreeMatrix(&x_ans);

  ndlqr_FreeLQRProblem(lqrprob);
  ndlqr_FreeNdLqrSolver(solver);
  return 1;
}


void AllTests() {
  // mu_run_test(SolveLeaves);
  mu_run_test(RunSolve);
  // mu_run_test(SolveTwice);
  // mu_run_test(FactorInnerProduct);
  // mu_run_test(ShurCompliment);
}

mu_test_main
