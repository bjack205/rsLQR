#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "cholesky_factors.h"
#include "linalg.h"
#include "eigen_c/eigen_c.h"
#include "matrix.h"
#include "test/minunit.h"

mu_test_init

    Matrix* A;
Matrix* B;
Matrix* C;

void InitMatrices(int N, int n) {
  A = (Matrix*)malloc(N * sizeof(Matrix));
  B = (Matrix*)malloc(N * sizeof(Matrix));
  C = (Matrix*)malloc(N * sizeof(Matrix));
  for (int i = 0; i < N; ++i) {
    A[i] = NewMatrix(n, n);
    B[i] = NewMatrix(n, n);
    C[i] = NewMatrix(n, n);
    for (int j = 0; j < n * n; ++j) {
      B[i].data[j] = cos(0.01 * j);
      A[i].data[j] = (2.1 * j - 3.2 * j * j) / 100.0;
    }
    MatrixSetConst(&C[i], 0.0);
  }
}

void MakePSD(int N) {
  for (int i = 0; i < N; ++i) {
    MatrixCopyTranspose(B + i, A + i);
    MatrixMultiply(A + i, B + i, C + i, 0, 0, 1.0, 0.0);
    MatrixCopy(A + i, C + i);
    MatrixSetConst(C + i, 0.0);
    for (int j = 0; j < A->rows; ++j) {
      double* djj = MatrixGetElement(A + i, j, j);
      *djj += 100;
    }
  }
}

void FreeMatrices(int N) {
  for (int i = 0; i < N; ++i) {
    FreeMatrix(&A[i]);
    FreeMatrix(&B[i]);
    FreeMatrix(&C[i]);
  }

  free(A);
  free(B);
  free(C);
}

int CholeskyFree() {
  int n = 8;
  int N = 1;
  InitMatrices(N, n);
  MakePSD(N);
  CholeskyInfo cholinfo;
  MatrixCholeskyFactorizeWithInfo(A, &cholinfo);
  eigen_MatrixMultiply(n, n, n, A->data, A->data, A->data, 1, 0, 1.0, 0.0);
  MatrixScaleByConst(A, 1e-6);
  MatrixCholeskyFactorizeWithInfo(A, &cholinfo);
  FreeFactorization(&cholinfo);
  FreeMatrices(N);
  return 1;
}

int CholeskyFactors() {
  int n = 8;
  int N = 2;
  InitMatrices(N, n);
  MakePSD(N);
  bool usingeigen = MatrixGetLinearAlgebraLibrary() == libEigen;

  NdLqrCholeskyFactors* cholfacts = ndlqr_NewCholeskyFactors(3, 8);
  CholeskyInfo* cholinfo = NULL;
  ndlqr_GetQFactorizon(cholfacts, 0, &cholinfo);
  mu_assert(cholinfo == cholfacts->cholinfo);
  MatrixCholeskyFactorizeWithInfo(A, cholinfo);
  MatrixCholeskySolveWithInfo(A, B, cholinfo);
  if (usingeigen) {
    mu_assert(cholinfo->is_freed == false);
    mu_assert(cholinfo->lib == 'E');
  }

  ndlqr_GetRFactorizon(cholfacts, 0, &cholinfo);
  if (usingeigen) {
    mu_assert(cholinfo->is_freed == true);
  }
  MatrixCholeskyFactorizeWithInfo(A + 1, cholinfo);
  if (usingeigen) {
    mu_assert(cholinfo->is_freed == false);
  }

  int res = ndlqr_FreeCholeskyFactors(cholfacts);
  mu_assert(cholinfo->is_freed == true);
  FreeMatrices(N);
  mu_assert(res == 0);
  return 1;
}

void AllTests() {
  mu_run_test(CholeskyFree);
  mu_run_test(CholeskyFactors);
}

mu_test_main
