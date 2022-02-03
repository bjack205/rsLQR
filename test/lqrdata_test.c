#include <stdlib.h>
#include <string.h>

#include "lqr/lqr_problem.h"
#include "test/minunit.h"
#include "test/test_problem.h"

#ifndef LQRDATAFILE
#define LQRDATAFILE "../lqr_data.json"
#endif

mu_test_init

    int
    CheckLQRData(LQRData* lqrdata) {
  mu_assert(lqrdata->nstates == 6);
  mu_assert(lqrdata->ninputs == 3);
  double h = 0.1;

  Matrix Q = ndlqr_GetQ(lqrdata);
  for (int i = 0; i < 6; ++i) {
    mu_assert(*MatrixGetElement(&Q, i, 0) == 1.0);
  }

  Matrix R = ndlqr_GetR(lqrdata);
  for (int i = 0; i < 3; ++i) {
    mu_assert(*MatrixGetElement(&R, i, 0) == 0.01);
  }

  Matrix A = ndlqr_GetA(lqrdata);
  Matrix B = ndlqr_GetB(lqrdata);
  for (int i = 0; i < 6; ++i) {
    mu_assert(*MatrixGetElement(&A, i, i) == 1.0);
  }
  for (int i = 0; i < 3; ++i) {
    mu_assert(fabs(*MatrixGetElement(&A, i, i + 3) - h) < 1e-8);
    mu_assert(fabs(*MatrixGetElement(&B, i + 3, i) - h) < 1e-8);
  }
  return 1;
}

int ReadLQRDataFileTest() {
  const char* filename = LQRDATAFILE;
  LQRData* lqrdata = ndlqr_ReadLQRDataJSONFile(filename);
  mu_assert(CheckLQRData(lqrdata) == 1);
  ndlqr_FreeLQRData(lqrdata);
  return 1;
}

int ReadTestDataFile() {
  LQRData* lqrdata = ndlqr_ReadTestLQRData();
  mu_assert(lqrdata->nstates == 6);
  mu_assert(lqrdata->ninputs == 3);
  ndlqr_FreeLQRData(lqrdata);
  return 1;
}

int ReadProblemFile() {
  LQRProblem* lqrprob = ndlqr_ReadTestLQRProblem();
  mu_assert(lqrprob->nhorizon == 8);

  // Check initial state
  mu_assert(lqrprob->x0[0] == +1);
  mu_assert(lqrprob->x0[1] == -1);
  mu_assert(lqrprob->x0[2] == +2);
  mu_assert(lqrprob->x0[3] == -2);
  mu_assert(lqrprob->x0[4] == +3);
  mu_assert(lqrprob->x0[5] == -3);

  // Check each knot point
  for (int k = 0; k < 7; ++k) {
    LQRData* lqrdata = lqrprob->lqrdata[0];
    mu_assert(CheckLQRData(lqrdata) == 1);
  }
  LQRData* lqrdata = lqrprob->lqrdata[7];

  // Check Terminal knot point
  mu_assert(lqrdata->nstates == 6);
  mu_assert(lqrdata->ninputs == 3);

  Matrix Q = ndlqr_GetQ(lqrdata);
  for (int i = 0; i < 6; ++i) {
    mu_assert(*MatrixGetElement(&Q, i, 0) == 10.0);
  }

  Matrix R = ndlqr_GetR(lqrdata);
  for (int i = 0; i < 3; ++i) {
    mu_assert(*MatrixGetElement(&R, i, 0) == 0.0);
  }

  Matrix d = ndlqr_Getd(lqrdata);
  for (int i = 0; i < 6; ++i) {
    mu_assert(fabs(*MatrixGetElement(&d, i, 0)) < 1e-8);
  }

  ndlqr_FreeLQRProblem(lqrprob);
  return 1;
}

int NewLQRProblem() {
  LQRProblem* prob = ndlqr_NewLQRProblem(6, 3, 8);
  mu_assert(prob->nhorizon == 8);
  for (int k = 0; k < 8; ++k) {
    mu_assert(prob->lqrdata[k] != NULL);
  }
  ndlqr_FreeLQRProblem(prob);
  return 1;
}

int InitializeLQRProblem() {
  int nhorizon = 8;
  int nstates = 6;
  int ninputs = 3;
  LQRProblem* prob = ndlqr_NewLQRProblem(nstates, ninputs, nhorizon);
  Matrix Q = NewMatrix(nstates, 1);
  Matrix R = NewMatrix(ninputs, 1);
  Matrix q = NewMatrix(nstates, 1);
  Matrix r = NewMatrix(ninputs, 1);
  Matrix A = NewMatrix(nstates, nstates);
  Matrix B = NewMatrix(nstates, ninputs);
  Matrix d = NewMatrix(nstates, 1);
  LQRData** alldata = (LQRData**)malloc(nhorizon * sizeof(LQRData*));
  for (int k = 0; k < nhorizon; ++k) {
    LQRData* lqrdata = ndlqr_NewLQRData(nstates, ninputs);
    MatrixSetConst(&Q, 1.2);
    MatrixSetConst(&R, 1.4);
    MatrixSetConst(&q, 0.2 + 3 * k);
    MatrixSetConst(&r, 0.4 + 2 * k);
    double c = -0.5 + k;
    MatrixSetConst(&A, 10.1);
    MatrixSetConst(&B, 11.1);
    MatrixSetConst(&d, -20.2 + k);
    ndlqr_InitializeLQRData(lqrdata, Q.data, R.data, q.data, r.data, c, A.data, B.data,
                            d.data);
    alldata[k] = lqrdata;
  }

  double* x0 = (double*)malloc(nstates * sizeof(double));
  for (int i = 0; i < nstates; ++i) {
    x0[i] = i * i - 0.2 * i;
  }
  printf("nstates = %d\n", alldata[0]->nstates);
  ndlqr_InitializeLQRProblem(prob, x0, alldata);

  for (int k = 0; k < nhorizon; ++k) {
    MatrixSetConst(&Q, 1.2);
    MatrixSetConst(&R, 1.4);
    MatrixSetConst(&q, 0.2 + 3 * k);
    MatrixSetConst(&r, 0.4 + 2 * k);
    double c = -0.5 + k;
    MatrixSetConst(&A, 10.1);
    MatrixSetConst(&B, 11.1);
    MatrixSetConst(&d, -20.2 + k);
    Matrix A2 = ndlqr_GetA(prob->lqrdata[k]);
    Matrix B2 = ndlqr_GetB(prob->lqrdata[k]);
    Matrix d2 = ndlqr_Getd(prob->lqrdata[k]);
    Matrix Q2 = ndlqr_GetQ(prob->lqrdata[k]);
    Matrix R2 = ndlqr_GetR(prob->lqrdata[k]);
    Matrix q2 = ndlqr_Getq(prob->lqrdata[k]);
    Matrix r2 = ndlqr_Getr(prob->lqrdata[k]);
    double c2 = *(prob->lqrdata[k]->c);
    mu_assert(MatrixNormedDifference(&A, &A2) < 1e-6);
    mu_assert(MatrixNormedDifference(&B, &B2) < 1e-6);
    mu_assert(MatrixNormedDifference(&d, &d2) < 1e-6);
    mu_assert(MatrixNormedDifference(&Q, &Q2) < 1e-6);
    mu_assert(MatrixNormedDifference(&R, &R2) < 1e-6);
    mu_assert(MatrixNormedDifference(&q, &q2) < 1e-6);
    mu_assert(MatrixNormedDifference(&r, &r2) < 1e-6);
    mu_assert(fabs(c - c2) < 1e-8);
  }

  for (int k = 0; k < nhorizon; ++k) {
    ndlqr_FreeLQRData(alldata[k]);
  }
  free(alldata);
  free(x0);
  ndlqr_FreeLQRProblem(prob);
  FreeMatrix(&Q);
  FreeMatrix(&R);
  FreeMatrix(&q);
  FreeMatrix(&r);
  FreeMatrix(&A);
  FreeMatrix(&B);
  FreeMatrix(&d);
  return 1;
}

int ReadLongProb() {
  LQRProblem* prob = ndlqr_ReadLongTestLQRProblem();
  mu_assert(prob->nhorizon == 256);
  ndlqr_FreeLQRProblem(prob);
  return 1;
}

void AllTests() {
  mu_run_test(NewLQRProblem);
  mu_run_test(ReadLQRDataFileTest);
  mu_run_test(ReadTestDataFile);
  mu_run_test(ReadProblemFile);
  mu_run_test(InitializeLQRProblem);
  mu_run_test(ReadLongProb);
}

mu_test_main
