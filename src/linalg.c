#include "linalg.h"

#ifdef USE_MKL
#include "mkl_lapacke.h"
#include "mkl_cblas.h"
static const int kUseMKL = 1;
#else
static const int kUseMKL = 0;
#endif

#ifdef USE_BLAS
#include "lapacke.h"
#include "cblas.h"
#endif

#include <stdio.h>

#include "eigen_c/eigen_c.h"
#include "linalg_utils.h"
#include "linalg_custom.h"

CholeskyInfo DefaultCholeskyInfo() {
  CholeskyInfo cholinfo = {'\0', 0, '\0', NULL, 1};
  return cholinfo;
}

void FreeFactorization(CholeskyInfo* cholinfo) {
  if (cholinfo->lib == 'E' && !cholinfo->is_freed) {  // Eigen
    eigen_FreeFactorization(cholinfo->fact);
    cholinfo->is_freed = true;
  }
}

int MatrixAddition(Matrix* A, Matrix* B, double alpha) {
  if (!A || !B) return -1;
  int n = MatrixNumElements(A);
  switch (MatrixGetLinearAlgebraLibrary()) {
    case libEigen:
      eigen_MatrixAddition(n, A->data, B->data, alpha);
      break;

    default:
      clap_MatrixAddition(A, B, alpha);
      break;
  }
  return 0;
}

int MatrixCholeskyFactorizeWithInfo(Matrix* mat, CholeskyInfo* cholinfo) {
  MATRIX_LATIME_START;
  int out = 0;
  switch (MatrixGetLinearAlgebraLibrary()) {
    case libEigen:
      FreeFactorization(cholinfo);  // free the previous factorization since the code below allocates a new one
      out = eigen_CholeskyFactorize(mat->rows, mat->data, &cholinfo->fact);
      cholinfo->lib = 'E';
      cholinfo->is_freed = false;
      break;
    
#ifdef USE_BLAS
    case libBLAS:
      out = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', mat->rows, mat->data, mat->rows);
      cholinfo->lib = 'B';
      break;
#endif

    default: 
      cholinfo->lib = 'I';  // Internal
      cholinfo->is_freed = true;
      out = clap_CholeskyFactorize(mat);
      cholinfo->success = out == 0;
      break;
  }
  MATRIX_LATIME_STOP;
  cholinfo->success = out;
  cholinfo->uplo = 'L';
  return out;
}

int MatrixCholeskySolveWithInfo(Matrix* A, Matrix* b, CholeskyInfo* cholinfo) {
  MATRIX_LATIME_START;
  int out = 0;
  switch (MatrixGetLinearAlgebraLibrary()) {
    case libEigen:
      eigen_CholeskySolve(A->rows, b->cols, cholinfo->fact, b->data);
      break;

#ifdef USE_BLAS
    case libBLAS: 
      (void) cholinfo;
      out = LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', A->rows, b->cols, A->data, A->rows, b->data, b->rows);
      break;
#endif

    default:
      (void) cholinfo;
      clap_CholeskySolve(A, b);
      break;
  }
  MATRIX_LATIME_STOP;
  return out;
}

int MatrixCholeskyFactorize(Matrix* mat) {
  MATRIX_LATIME_START;
  int out = 0;
  switch (MatrixGetLinearAlgebraLibrary()) {

#ifdef USE_BLAS
    case libBLAS:
      out = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', mat->rows, mat->data, mat->rows);
      break;
#endif

    default:
      out = clap_CholeskyFactorize(mat);
      break;
  }
  MATRIX_LATIME_STOP;
  return out;
}

int MatrixCholeskySolve(Matrix* A, Matrix* b) {
  MATRIX_LATIME_START;
  int out = 0;
  switch (MatrixGetLinearAlgebraLibrary()) {

#ifdef USE_BLAS
    case libBLAS:
      out = LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', A->rows, b->cols, A->data, A->rows, b->data, b->rows);
      break;
#endif

    default:
      out = clap_CholeskySolve(A, b);
      break;

  }
  MATRIX_LATIME_STOP;
  return out;
}

void MatrixMultiply(Matrix* A, Matrix* B, Matrix* C, bool tA, bool tB, double alpha, double beta) {
  MATRIX_LATIME_START;
  // printf("Multiplying matrices of size (%d, %d), (%d, %d), (%d, %d)\n", A->rows, A->cols, B->rows, B->cols, C->rows, C->cols);
  switch (MatrixGetLinearAlgebraLibrary()) {
    case libEigen: {
      int m = tA ? A->cols : A->rows;
      int n = tA ? A->rows : A->cols;
      int k = tB ? B->rows : B->cols;
      eigen_MatrixMultiply(m, n, k, A->data, B->data, C->data, tA, tB, alpha, beta);
    } break;

#ifdef USE_BLAS
    case libBLAS: {
      CBLAS_TRANSPOSE transA = tA ? CblasTrans : CblasNoTrans; 
      CBLAS_TRANSPOSE transB = tB ? CblasTrans : CblasNoTrans; 
      if (B->cols == 1) {
        cblas_dgemv(CblasColMajor, transA, A->rows, A->cols, alpha, A->data, A->rows, B->data, 
          1, beta, C->data, 1);
      } else {
        int m = tA ? A->cols : A->rows;
        int n = tB ? B->rows : B->cols;
        int k = tA ? A->rows : A->cols;
        cblas_dgemm(CblasColMajor, transA, transB, m, n, k, alpha,
          A->data, A->rows, B->data, B->rows, beta, C->data, C->rows);
      }
    } break;
#endif

    default: {
      clap_MatrixMultiply(A, B, C, tA, tB, alpha, beta);
    } break;

  }
  MATRIX_LATIME_STOP;
}

void MatrixSymmetricMultiply(Matrix* Asym, Matrix* B, Matrix* C, double alpha, double beta) {
  MATRIX_LATIME_START;
  switch (MatrixGetLinearAlgebraLibrary()) {
    case libBLAS: {
#ifdef USE_BLAS
      CBLAS_SIDE blasside = CblasLeft;
      if (B->cols == 1) {
        cblas_dsymv(CblasColMajor, CblasLower, Asym->rows, alpha, Asym->data, Asym->rows, B->data, 1.0, beta, C->data, 1.0);
      } else {
        cblas_dsymm(CblasColMajor, blasside, CblasLower, Asym->rows, C->cols, alpha, Asym->data, 
          Asym->rows, B->data, B->rows, beta, C->data, C->rows);
      }
#endif
    } break;
    default: {
      clap_SymmetricMatrixMultiply(Asym, B, C, alpha, beta);
    } break;  
  }
  MATRIX_LATIME_STOP;
}

void MatrixCopyDiagonal(Matrix* dest, Matrix* src) {
  if (!dest || !src) return;
  MatrixSetConst(dest, 0.0);
  for (int i = 0; i < MatrixNumElements(src); ++i) {
    MatrixSetElement(dest, i, i, src->data[i]);
  }
}

enum MatrixLinearAlgebraLibrary MatrixGetLinearAlgebraLibrary() {
  if (kUseEigen) {
    return libEigen;
  } else if (kUseClap) {
    return libInternal;
  } else if (kUseMKL) {
    return libMKL;
  } else {
    return libBLAS;
  }
}

void MatrixPrintLinearAlgebraLibrary() {
  printf("Using ");
  switch (MatrixGetLinearAlgebraLibrary()) {
    case libEigen:
      printf("C++ Eigen Library.\n");
      break;
    case libInternal: 
      printf("Internal Linear Algebra Library.\n");
      break;
    case libMKL:
      printf("Intel MKL Library.\n");
      break;
    case libBLAS:
      printf("Open-source BLAS library.\n");
      break;
    default:
      printf("Not a recognized library.\n");
      break;
  }
}
