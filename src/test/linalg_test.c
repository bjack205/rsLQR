#include "matrix.h"
#include "test/minunit.h"
#include "linalg.h"

int MatMul() {
  // Matrix-matrix
  Matrix A = NewMatrix(3,4);
  Matrix B = NewMatrix(4,5);
  Matrix C = NewMatrix(3,5);
  Matrix Cans = NewMatrix(3,5);
  Matrix Bans = NewMatrix(4,5);
  MatrixSetConst(&A, 4);
  MatrixSetConst(&B, 3);
  MatrixSetConst(&C, 2);
  MatrixSetConst(&Cans, 50);
  MatrixSetConst(&Bans, 0);
  MatrixMultiply(&A, &B, &C, 0, 0, 1.0, 1.0);
  mu_assert(MatrixNormedDifference(&C, &Cans) < 1e-6);

  MatrixMultiply(&A, &C, &B, 1, 0, 1.0, -200);
  mu_assert(MatrixNormedDifference(&B, &Bans) < 1e-6);

  // Matrix-vector
  double xdata[4] = {1,2,3,4};
  double bdata[3] = {40, 40, 40};
  Matrix x = {4, 1, xdata};
  Matrix bans = {3, 1, bdata};
  Matrix b = NewMatrix(3, 1);
  MatrixMultiply(&A, &x, &b, 0, 0, 1.0, 0.0);
  mu_assert(MatrixNormedDifference(&b, &bans) < 1e-6);

  FreeMatrix(&A);
  FreeMatrix(&B);
  FreeMatrix(&C);
  FreeMatrix(&Cans);
  FreeMatrix(&Bans);
  FreeMatrix(&b);

  // printf("Using IntelMKL library? %d\n", kUseIntelMKL);
  return 1;
}

int SymMatMul() {
  // Matrix-Vector
  int n = 3;
  double Adata[9] = {9,-3,-6,-3,17,-10,-6,-10,38};
  double bdata[3] = {6,22,28};
  double xdata[3] = {3,3,2};
  Matrix A = {n, n, Adata};
  Matrix bans = {n, 1, bdata};
  Matrix x = {n, 1, xdata};
  Matrix b = NewMatrix(n, 1);
  MatrixSymmetricMultiply(&A, &x, &b, 1.0, 0.0);
  mu_assert(MatrixNormedDifference(&b, &bans) < 1e-6);

  // Matrix-matrix
  double Xdata[6] = {3,3,2,1,1,1};
  double Bdata[6] = {6,22,28, 0, 4, 22};
  Matrix X = {n, 2, Xdata};
  Matrix Bans = {n, 2, Bdata};
  Matrix B = NewMatrix(n, 2);
  MatrixSymmetricMultiply(&A, &X, &B, 1.0, 0.0);
  mu_assert(MatrixNormedDifference(&B, &Bans) < 1e-6);

  FreeMatrix(&b);
  FreeMatrix(&B);
  return 1;
}
 
int DiagonalCholesky() {
  int n = 5;
  Matrix A = NewMatrix(n, n);
  MatrixSetConst(&A, 0.0);
  for (int i = 0; i < n; ++i) {
    MatrixSetElement(&A, i, i, 9);
  }
  CholeskyInfo cholinfo = DefaultCholeskyInfo();
  MatrixCholeskyFactorizeWithInfo(&A, &cholinfo);
  for (int i = 0; i < n; ++i) {
    mu_assert(*MatrixGetElement(&A, i, i) == 3);
  }
  FreeMatrix(&A);
  return 1;
}

int DiagonalCholeskySolve() {
  int n = 5;
  Matrix A = NewMatrix(n, n);
  Matrix b = NewMatrix(n, 1);
  MatrixSetConst(&A, 0.0);
  for (int i = 0; i < n; ++i) {
    MatrixSetElement(&A, i, i, 9);
    MatrixSetElement(&b, i, 0, 9 * (i + 1));
  }
  CholeskyInfo cholinfo = DefaultCholeskyInfo();
  int res = MatrixCholeskyFactorizeWithInfo(&A, &cholinfo);
  mu_assert(res == 0);
  MatrixCholeskySolveWithInfo(&A, &b, &cholinfo);
  for (int i = 0; i < n; ++i) {
    mu_assert(*MatrixGetElement(&A, i, i) == 3);
    mu_assert(*MatrixGetElement(&b, i, 0) == (i + 1));
  }
  FreeMatrix(&A);
  FreeMatrix(&b);
  return 1;
}

int CholeskySolve3x3() {
  int n = 3;  
  double Adata[9] = {9,-3,-6,-3,17,-10,-6,-10,38};
  double bdata[3] = {6,22,28};
  double xdata[3] = {3,3,2};
  Matrix A = {n, n, Adata};
  Matrix b = {n, 1, bdata};
  Matrix x = {n, 1, xdata};
  CholeskyInfo cholinfo = DefaultCholeskyInfo();
  MatrixCholeskyFactorizeWithInfo(&A, &cholinfo);
  MatrixCholeskySolveWithInfo(&A, &b, &cholinfo);
  mu_assert(MatrixNormedDifference(&b, &x) < 1e-6);
  return 1;
}

void AllTests() {
  mu_run_test(DiagonalCholesky);
  mu_run_test(DiagonalCholeskySolve);
  mu_run_test(CholeskySolve3x3);
  mu_run_test(MatMul);
  mu_run_test(SymMatMul);
  MatrixPrintLinearAlgebraLibrary();
}

mu_test_main
