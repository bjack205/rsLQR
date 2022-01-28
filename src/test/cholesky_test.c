#include <stdlib.h>
#include <stdio.h>

#include <omp.h>

#include "cholesky_factors.h"
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
      B[i].data[j] = cos(0.01 * j);
      A[i].data[j] = (2.1 * j  - 3.2 * j * j) / 100.0;
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
 
int Cholesky4x4() {
  int N = 16 * 1000;
  int n = 8;
  InitMatrices(N, n);
  MakePSD(N);

  MatrixCopy(A + 1, A + 0);
  MatrixCopy(B + 1, B + 0);
  Matrix* A0 = A + 2;
  Matrix* B0 = B + 2;
  MatrixCopy(A0, A + 0);
  MatrixCopy(B0, B + 0);
  MatrixCholeskyFactorize(A + 0);
  void* achol;
  eigen_CholeskyFactorize(n, A[1].data, &achol);
  mu_assert(MatrixNormedDifference(A + 0, A + 1) < 1e-10);

  MatrixCholeskySolve(A + 0, B + 0);
  eigen_CholeskySolve(n, n, achol, B[1].data);
  mu_assert(MatrixNormedDifference(B + 0, B + 1) < 1e-10);

  MatrixSymmetricMultiply(A0, B + 0, C + 0, 'L', 1.0, 0.0);
  mu_assert(MatrixNormedDifference(C + 0, B0) < 1e-10);
  mu_assert(MatrixNormedDifference(C + 0, B0) < 1e-10);
  eigen_SymmetricMatrixMultiply(n, n, A0->data, B[1].data, C[1].data);
  mu_assert(MatrixNormedDifference(C + 1, B0) < 1e-10);
  mu_assert(MatrixNormedDifference(C + 1, C + 0) < 1e-10);

  MatrixCopy(A + 0, A0);
  MatrixCopy(A + 1, A0);

  printf("Cholesky Factorization\n");
  double t_start = omp_get_wtime();
  for (int i = 0; i < N; ++i) {
    MatrixCholeskyFactorize(A + i);
  }
  double t_blas = (omp_get_wtime() - t_start) * 1000.0;
  printf("  BLAS:  %f ms\n", t_blas);

  void** achols = (void**) malloc(N * sizeof(void*));
  t_start = omp_get_wtime();
  for (int i = 0; i < N; ++i) {
    eigen_CholeskyFactorize(n, A[i].data, achols + i);
  }
  double t_eigen = (omp_get_wtime() - t_start) * 1000.0;
  printf("  Eigen: %f ms\n", t_eigen);

  // In parallel
  MakePSD(N);
  int nthreads = 4;
  int tasks_per_thread = N / nthreads;
  eigen_SetNumThreads(1);
  eigen_InitParallel();
  omp_set_num_threads(nthreads);
  #pragma omp parallel firstprivate(tasks_per_thread)
  {
    #pragma omp single
    {
      t_start = omp_get_wtime();
    }
    int id = omp_get_thread_num();
    int start = tasks_per_thread * id;
    int stop = tasks_per_thread * (id + 1);
    for (int i = start; i < stop; ++i) {
      MatrixCholeskyFactorize(A + i);
    }
    #pragma omp single
    {
      double t_blas_parallel = (omp_get_wtime() - t_start) * 1000.0;
      printf("  BLAS:  %f ms (speedup = %.2f)\n", t_blas_parallel, t_blas / t_blas_parallel);
    }

    #pragma omp single
    {
      t_start = omp_get_wtime();
    }
    for (int i = start; i < stop; ++i) {
      eigen_CholeskyFactorize(n, A[i].data, achols + i);
    }
    #pragma omp single
    {
      double t_eigen_parallel = (omp_get_wtime() - t_start) * 1000.0;
      printf("  Eigen: %f ms (speedup = %.2f)\n", t_eigen_parallel, t_eigen / t_eigen_parallel);
    }
  }

  printf("\nCholesky Solve\n");
  t_start = omp_get_wtime();
  for (int i = 0; i < N; ++i) {
    MatrixCholeskySolve(A + i, B + i);
  }
  t_blas = (omp_get_wtime() - t_start) * 1000.0;
  printf("  BLAS:  %f ms\n", t_blas);

  t_start = omp_get_wtime();
  for (int i = 0; i < N; ++i) {
    eigen_CholeskySolve(n, n, achols[i], B[i].data);
  }
  t_eigen = (omp_get_wtime() - t_start) * 1000.0;
  printf("  Eigen: %f ms\n", t_eigen);

  #pragma omp parallel firstprivate(tasks_per_thread)
  {
    #pragma omp single
    {
      t_start = omp_get_wtime();
    }
    int id = omp_get_thread_num();
    int start = tasks_per_thread * id;
    int stop = tasks_per_thread * (id + 1);
    for (int i = start; i < stop; ++i) {
      MatrixCholeskySolve(A + i, B + i);
    }
    #pragma omp single
    {
      double t_blas_parallel = (omp_get_wtime() - t_start) * 1000.0;
      printf("  BLAS:  %f ms (speedup = %.2f)\n", t_blas_parallel, t_blas / t_blas_parallel);
    }

    #pragma omp single
    {
      t_start = omp_get_wtime();
    }
    for (int i = start; i < stop; ++i) {
      eigen_CholeskySolve(n, n, achols[i], B[i].data);
    }
    #pragma omp single
    {
      double t_eigen_parallel = (omp_get_wtime() - t_start) * 1000.0;
      printf("  Eigen: %f ms (speedup = %.2f)\n", t_eigen_parallel, t_eigen / t_eigen_parallel);
    }
  }

  FreeMatrices(N);
  return 1;
}

int CholeskyFactors() {
  int n = 8;
  int N = 2;
  InitMatrices(N, n);
  MakePSD(N);

  NdLqrCholeskyFactors* cholfacts = ndlqr_NewCholeskyFactors(3, 8);
  CholeskyInfo* cholinfo = NULL;
  ndlqr_GetQFactorizon(cholfacts, 0, &cholinfo);
  mu_assert(cholinfo == cholfacts->cholinfo);
  MatrixCholeskyFactorizeWithInfo(A, cholinfo);
  MatrixCholeskySolveWithInfo(A, B, cholinfo);
  if (kUseEigen) {
    mu_assert(cholinfo->is_freed == false);
    mu_assert(cholinfo->lib == 'E');
  }

  ndlqr_GetRFactorizon(cholfacts, 0, &cholinfo);
  if (kUseEigen) {
    mu_assert(cholinfo->is_freed == true);
  }
  MatrixCholeskyFactorizeWithInfo(A + 1, cholinfo);
  if (kUseEigen) {
    mu_assert(cholinfo->is_freed == false);
  }

  int res = ndlqr_FreeCholeskyFactors(cholfacts);
  mu_assert(cholinfo->is_freed == true);
  FreeMatrices(N);
  mu_assert(res == 0);
  return 1;
}

void AllTests() {
  // mu_run_test(CholeskyFree);
  // mu_run_test(Cholesky4x4);
  mu_run_test(CholeskyFactors);
}

mu_test_main
