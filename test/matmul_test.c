#include <stdlib.h>
#include <stdio.h>

#include <omp.h>

#include "test/minunit.h"
#include "matrix.h"
#include "linalg.h"
#include "matmul.h"

mu_test_init

Matrix* A;
Matrix* B;
Matrix* C;

void InitMatrices(int N, int n) {
  A = (Matrix*) malloc(N * sizeof(Matrix));
  B = (Matrix*) malloc(N * sizeof(Matrix));
  C = (Matrix*) malloc(N * sizeof(Matrix));
  for (int i = 0; i < N; ++i) {
    A[i] = NewMatrix(n, n);
    B[i] = NewMatrix(n, n);
    C[i] = NewMatrix(n, n);
    for (int j = 0; j < n * n; ++j) {
      A[i].data[j] = cos(0.01 * j);
      B[i].data[j] = 2.1 * j  - 3.2 * j * j;
    }
    MatrixSetConst(&C[i], 0.0);
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

int Matrix4x4() {
  int N = 16 * 1000;
  int n = 4;
  InitMatrices(N, n);

  
  MatMul4x4_unrolled(A[0].data, B[0].data, C[0].data);
  MatMul4x4_SIMD(A[0].data, B[0].data, C[1].data);
  mu_assert(MatrixNormedDifference(C, C + 1) < 1e-10);

  double t_start = omp_get_wtime();
  for (int i = 0; i < N; ++i) {
    MatrixMultiply(A + i, B + i, C + i, 0, 0, 1.0, 0.0);
  }
  double t_blas = (omp_get_wtime() - t_start) * 1000.0;

  t_start = omp_get_wtime();
  for (int i = 0; i < N; ++i) {
    MatMulLoops(A + i, B + i, C + i);
  }
  double t_loops = (omp_get_wtime() - t_start) * 1000.0;

  t_start = omp_get_wtime();
  for (int i = 0; i < N; ++i) {
    MatMul4x4_SIMD(A[i].data, B[i].data, C[i].data);
  }
  double t_simd = (omp_get_wtime() - t_start) * 1000.0;

  t_start = omp_get_wtime();
  for (int i = 0; i < N; ++i) {
    MatMul4x4_unrolled(A[i].data, B[i].data, C[i].data);
  }
  double t_unroll = (omp_get_wtime() - t_start) * 1000.0;

  t_start = omp_get_wtime();
  for (int i = 0; i < N; ++i) {
    eigen_MatrixMultiply(n, n, n, A[i].data, B[i].data, C[i].data, 0, 0, 1.0, 0.0);
  }
  double t_eigen = (omp_get_wtime() - t_start) * 1000.0;
  printf("BLAS:   %f ms\n", t_blas);
  printf("Loops:  %f ms\n", t_loops);
  printf("SIMD:   %f ms\n", t_simd);
  printf("Unroll: %f ms\n", t_unroll);
  printf("Eigen:  %f ms\n", t_eigen);
  return 1;
}

int Matrix5x5() {
  double Adata[25] = {
    -11, -11, -10, -9, -8, 
     -7,  -5,  -5, -4, -3, 
     -2,  -1,   1,   1, 2, 
      3,   4,   5,   7, 7, 
      8,   9,  10,  11, 13
  };
  double Bdata[25];
  for (int i = 0; i < 25; ++i) {
    Bdata[i] = Adata[i] * 2;
  }

  Matrix A = {5, 5, Adata};
  Matrix B = {5, 5, Bdata};
  Matrix C = NewMatrix(5, 5);
  MatrixMultiply(&A, &B, &C, 0, 0, 1.0, 0.0);
  PrintMatrix(&C);
  printf("\nUnrolled\n");
  MatrixSetConst(&C, 0.0);
  MatMul5x5_unrolled(A.data, B.data, C.data);
  PrintMatrix(&C);

  FreeMatrix(&C);
  return 1;
}

int SIMD_Mult() {
  int n = 8;
  Matrix A = NewMatrix(n, n);
  Matrix B = NewMatrix(n, n);
  Matrix C1 = NewMatrix(n, n);
  Matrix C2 = NewMatrix(n, n);
  for (int i = 0; i < n * n; ++i) {
    A.data[i] = 2.1 * i + (i - n / 2) * (i - n / 2);
    B.data[i] = A.data[i] * 2;
  }
  MatrixSetConst(&C1, 0.0);
  MatrixSetConst(&C2, 0.0);
  MatrixMultiply(&A, &B, &C1, 0, 0, 1.0, 0.0);
  PrintMatrix(&C1);
  MatMulSIMD(n, A.data, B.data, C2.data);
  printf("\nSIMD:\n");
  PrintMatrix(&C2);
  printf("error = %e\n", MatrixNormedDifference(&C1, &C2));

  FreeMatrix(&A);
  FreeMatrix(&B);
  FreeMatrix(&C1);
  FreeMatrix(&C2);
  return 1;
}

int MatMul8x8() {
  int n = 8;
  Matrix A = NewMatrix(n, n);
  Matrix B = NewMatrix(n, n);
  Matrix C1 = NewMatrix(n, n);
  Matrix C2 = NewMatrix(n, n);
  for (int i = 0; i < n * n; ++i) {
    A.data[i] = 2.1 * i + (i - n / 2) * (i - n / 2);
    B.data[i] = A.data[i] * 2;
  }
  MatrixSetConst(&C1, 0.0);
  MatrixSetConst(&C2, 0.0);
  MatrixMultiply(&A, &B, &C1, 0, 0, 1.0, 0.0);
  PrintMatrix(&C1);

  printf("\nUnrolled:\n");
  MatMul8x8_unrolled(A.data, B.data, C2.data);
  PrintMatrix(&C2);
  printf("error = %e\n", MatrixNormedDifference(&C1, &C2));

  printf("\nEigen\n");
  eigen_MatrixMultiply(n, n, n, A.data, B.data, C2.data, 0, 0, 1.0, 0.0);
  printf("error = %e\n", MatrixNormedDifference(&C1, &C2));

  return 1;
}

void AllTests() {
  mu_run_test(Matrix4x4);
  // mu_run_test(Matrix5x5);
  // mu_run_test(SIMD_Mult);
  // mu_run_test(MatMul8x8);
}

mu_test_main
