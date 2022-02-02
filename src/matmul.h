#pragma once

#include <linalg.h>
#include <matrix.h>

#include "eigen_c/eigen_c.h"

int MatMulLoops(Matrix* A, Matrix* B, Matrix* C);
int MatMulBlas(int n, double* A, double* B, double* C);
int MatMulSIMD(int n, double* A, double* B, double* C);

int MatMul4x4_SIMD(double* A, double* B, double* C);
int MatMul4x4_unrolled(double* A, double* B, double* C);
int MatMul5x5_unrolled(double* A, double* B, double* C);
int MatMul8x8_unrolled(double* a, double* b, double* c);
