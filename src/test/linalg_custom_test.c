#include "matrix.h"
#include "test/minunit.h"
#include "linalg.h"
#include "linalg_custom.h"
#include "eigen_c/eigen_c.h"

int MatMul() {
  // Matrix-matrix
  Matrix A = NewMatrix(3,4);
  Matrix B = NewMatrix(4,5);
  Matrix C = NewMatrix(3,5);
  Matrix D = NewMatrix(4,3);
  Matrix Cans = NewMatrix(3,5);
  Matrix Bans = NewMatrix(4,5);
  Matrix Dans = NewMatrix(4,3);
  MatrixSetConst(&A, 4);
  MatrixSetConst(&B, 3);
  MatrixSetConst(&C, 2);
  MatrixSetConst(&Cans, 50);
  MatrixSetConst(&Bans, 3);
  MatrixSetConst(&Dans, 750);
  clap_MatrixMultiply(&A, &B, &C, 0, 0, 1.0, 1.0);
  mu_assert(MatrixNormedDifference(&C, &Cans) < 1e-6);

  clap_MatrixMultiply(&A, &C, &B, 1, 0, 1.0, -199);
  mu_assert(MatrixNormedDifference(&B, &Bans) < 1e-6);

  clap_MatrixMultiply(&B, &C, &D, 0, 1, 1.0, 0.0);
  mu_assert(MatrixNormedDifference(&D, &Dans) < 1e-6);

  // Matrix-vector
  double xdata[4] = {1,2,3,4};
  double bdata[3] = {40, 40, 40};
  Matrix x = {4, 1, xdata};
  Matrix bans = {3, 1, bdata};
  Matrix b = NewMatrix(3, 1);
  clap_MatrixMultiply(&A, &x, &b, 0, 0, 1.0, 0.0);
  mu_assert(MatrixNormedDifference(&b, &bans) < 1e-6);

  FreeMatrix(&A);
  FreeMatrix(&B);
  FreeMatrix(&C);
  FreeMatrix(&D);
  FreeMatrix(&Cans);
  FreeMatrix(&Bans);
  FreeMatrix(&Dans);
  FreeMatrix(&b);
  return 1;
}

int SymMatMulTest() {
  double Adata[9] = {1,2,3, 4,5,6, 7,8,9};
  double Bdata[6] = {2,4,6, 8,6,4};
  double Cdata[6] = {3,6,9, 12,11,10};
  double Ddata[6] = {34,72,102, 56,92,116};
  Matrix A = {3, 3, Adata};
  Matrix B = {3, 2, Bdata};
  Matrix C = {3, 2, Cdata};
  Matrix D = {3, 2, Ddata};
  clap_SymmetricMatrixMultiply(&A, &B, &C, 1.0, 2.0);
  mu_assert(MatrixNormedDifference(&C, &D) < 1e-6);
  return 1;
}

int MatAddTest() {
  double Adata[6] = {1,2,3, 4,5,6};
  double Bdata[6] = {2,4,6, 8,6,4};
  double Cdata[6] = {3,6,9, 12,11,10};
  double Ddata[6] = {1,2,3, 4,1,-2};
  Matrix A = {2, 3, Adata};
  Matrix B = {2, 3, Bdata};
  Matrix C = {2, 3, Cdata};
  Matrix D = {2, 3, Ddata};
  clap_MatrixAddition(&A, &B, 1.0);
  mu_assert(MatrixNormedDifference(&B, &C) < 1e-6);

  clap_MatrixAddition(&A, &C, -2);
  mu_assert(MatrixNormedDifference(&C, &D) < 1e-6);

  return 1;
}

int MatScale() {
  double Adata[6] = {1,2,3, 4,5,6};
  double Bdata[6] = {3,6,9, 12,15,18};
  Matrix A = {2, 3, Adata};
  Matrix B = {2, 3, Bdata};
  clap_MatrixScale(&A, 3);
  mu_assert(MatrixNormedDifference(&A, &B) < 1e-6);
  return 1;
}

int CholeskyFactorizeTest() {
  int n = 10;
  Matrix A1 = NewMatrix(n, n);
  Matrix A2 = NewMatrix(n, n);
  Matrix A = NewMatrix(n, n);
  Matrix Achol = NewMatrix(n, n);
  for (int i = 0; i < n*n; ++i) {
    A1.data[i] = (i-4)*(i+3)/6.0;
    A2.data[i] = A1.data[i];
  }
  clap_MatrixMultiply(&A1, &A2, &A, 1, 0, 1.0, 0.0);
  clap_AddDiagonal(&A, 1.0);
  MatrixCopy(&Achol, &A);
  int res = clap_CholeskyFactorize(&Achol);
  mu_assert(res == clap_kCholeskySuccess);
  
  // Check answer with Eigen
  void* fact;
  eigen_CholeskyFactorize(n, A.data, &fact);
  mu_assert(MatrixNormedDifference(&A, &Achol) < 1e-6);

  // Try to factorize an indefinite matrix
  clap_MatrixMultiply(&A1, &A2, &A, 1, 0, 1.0, 0.0);
  clap_AddDiagonal(&A, -1.0);
  MatrixCopy(&Achol, &A);
  res = clap_CholeskyFactorize(&Achol);
  mu_assert(res == clap_kCholeskyFail);

  FreeMatrix(&A1);
  FreeMatrix(&A2);
  FreeMatrix(&A);
  FreeMatrix(&Achol);
  return 1;
}

int TriBackSubTest() {
  int n = 3; 
  double Ldata[9] = {1, 2, 5, 0, 1, 6, 0, 0, 7};
  double bdata[3] = {-2, 3, 10};
  double ydata[3] = {-2.0, 7.0, -3.142857142857143};
  double xdata[3] = {-19.142857142857142, 9.693877551020408, -0.4489795918367347};

  Matrix L = {n, n, Ldata};
  Matrix b = {n, 1, bdata};
  Matrix y = {n, 1, ydata};
  Matrix x = {n, 1, xdata};
  clap_LowerTriBackSub(&L, &b, 0);
  mu_assert(MatrixNormedDifference(&b, &y) < 1e-6);

  clap_LowerTriBackSub(&L, &y, 1);
  mu_assert(MatrixNormedDifference(&x, &y) < 1e-6);

  return 1;
}

int CholeskySolveTest() {
  int n = 10;
  int m = 1;
  Matrix A1 = NewMatrix(n, n);
  Matrix A2 = NewMatrix(n, n);
  Matrix A = NewMatrix(n, n);
  Matrix Achol = NewMatrix(n, n);
  Matrix b = NewMatrix(n, m);
  Matrix x = NewMatrix(n, m);
  Matrix x_eigen = NewMatrix(n, m);
  for (int i = 0; i < n*n; ++i) {
    A1.data[i] = (i-4)*(i+3)/6.0;
    A2.data[i] = A1.data[i];
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      MatrixSetElement(&b, i, j, -i - (n+m) / 2 + j);
    }
  }

  clap_MatrixMultiply(&A1, &A2, &A, 1, 0, 1.0, 0.0);
  clap_AddDiagonal(&A, 1.0);
  MatrixCopy(&Achol, &A);
  MatrixCopy(&x, &b);
  MatrixCopy(&x_eigen, &b);

  clap_CholeskyFactorize(&Achol);
  clap_CholeskySolve(&Achol, &x);
  
  // Check answer with Eigen
  void* fact = NULL;
  eigen_CholeskyFactorize(n, A.data, &fact);
  eigen_CholeskySolve(n, m, fact, x_eigen.data);
  mu_assert(MatrixNormedDifference(&x, &x_eigen) < 1e-6);

  FreeMatrix(&A1);
  FreeMatrix(&A2);
  FreeMatrix(&A);
  FreeMatrix(&Achol);
  FreeMatrix(&b);
  FreeMatrix(&x);
  FreeMatrix(&x_eigen);
  eigen_FreeFactorization(fact);
  return 1;
}

void AllTests() {
  mu_run_test(MatMul);
  mu_run_test(MatAddTest);
  mu_run_test(MatScale);
  mu_run_test(CholeskyFactorizeTest);
  mu_run_test(TriBackSubTest);
  mu_run_test(CholeskySolveTest);
  mu_run_test(SymMatMulTest);
}

mu_test_main
