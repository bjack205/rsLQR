#pragma once

#include "matrix.h"

#ifdef USE_MKL
static const int kUseIntelMKL = 1;
#else
static const int kUseIntelMKL = 0;
#endif

#ifdef USE_EIGEN
static const int kUseEigen = 1;
#else
static const int kUseEigen = 0;
#endif

#ifdef USE_CLAP
static const int kUseClap = 1;
#else
static const int kUseClap = 0;
#endif

typedef struct {
  char uplo;    // 'L' or 'U' 
  int success;  // 0 if success, failure otherwise
  char lib;     // 'B' for BLAS, 'E' for eigen
  void* fact;   // pointer to Eigen data
  int is_freed;  // has the Eigen data been freed
} CholeskyInfo;

enum MatrixLinearAlgebraLibrary {
  libBLAS = 0,
  libMKL = 1,
  libEigen = 2,
  libInternal = 3,
};

void FreeFactorization(CholeskyInfo* cholinfo);

int MatrixAddition(Matrix* A, Matrix* B, double alpha);
int MatrixCholeskyFactorize(Matrix* mat);
int MatrixCholeskyFactorizeWithInfo(Matrix* mat, CholeskyInfo* cholinfo);
int MatrixCholeskySolve(Matrix* A, Matrix* b);
int MatrixCholeskySolveWithInfo(Matrix* A, Matrix* b, CholeskyInfo* cholinfo);
int MatrixCholeskyInverse(Matrix* A);
void MatrixMultiply(Matrix* A, Matrix* B, Matrix* C, bool tA, bool tB, double alpha, double beta);
void MatrixSymmetricMultiply(Matrix* Asym, Matrix* B, Matrix* C, char side, double alpha, double beta);

void MatrixCopyDiagonal(Matrix* dest, Matrix* src);

enum MatrixLinearAlgebraLibrary MatrixGetLinearAlgebraLibrary(); 

void MatrixPrintLinearAlgebraLibrary();
