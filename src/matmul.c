#include "matmul.h"

#include <immintrin.h>
#include <math.h>
#include <omp.h>
#include <cblas.h>

int MatMulBlas(int n, double* A, double* B, double* C) {
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0,
    A, n, B, n, 0.0, C, n);
  return 0;
}

int MatMulLoops(Matrix* A, Matrix* B, Matrix* C) {
  int n = A->rows;
  int m = A->cols;
  int p = B->cols;
  for (int k = 0; k < m; ++k) {
    for (int i = 0; i < p; ++i) {
      double b = *MatrixGetElement(B, k, i);
      for (int j = 0; j < n; ++j) {
        double* c = MatrixGetElement(C, j, k);
        if (k == 0) *c = 0;
        *c += *MatrixGetElement(A, j, k) * b;
      }
    }
  }
  return 0;
}

int MatMulSIMD(int n, double* A, double* B, double* C) {
  static const int w = 4;  // SIMD width
  int nw = n / w;          // number of SIMD vectors per side

  __m256d a, b, c, d;
  // Loop over columns of A (each element of the sum)
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < nw; ++k) {
      a = _mm256_loadu_pd(A + k * w + i * n);

      // Loop over columns of the output
      for (int j = 0; j < n; ++j) {
        if (i == 0) {
          b = _mm256_set1_pd(B[i + j * n]);
          c = _mm256_mul_pd(a, b);
          _mm256_storeu_pd(C + k * w + j * n, c);
        } else {
          d = _mm256_loadu_pd(C + k * w + j * n);
          b = _mm256_set1_pd(B[i + j * n]);
          c = _mm256_fmadd_pd(a, b, d);
          _mm256_storeu_pd(C + k * w + j * n, c);
        }
      }
    }
  }
  return 0;
}

int MatMul4x4_SIMD(double* A, double* B, double* C) {
  static const int n = 4;

  __m256d a, b, c, d;
  for (int i = 0; i < 4; ++i) {
    if (i == 0) {
      a = _mm256_loadu_pd(A + i * n);
      for (int j = 0; j < 4; ++j) {
        b = _mm256_set1_pd(B[i + j * n]);
        c = _mm256_mul_pd(a, b);
        _mm256_storeu_pd(C + j * n, c);
      }
    } else {
      a = _mm256_loadu_pd(A + i * n);
      for (int j = 0; j < 4; ++j) {
        d = _mm256_loadu_pd(C + j * n);
        b = _mm256_set1_pd(B[i + j * n]);
        c = _mm256_fmadd_pd(a, b, d);
        _mm256_storeu_pd(C + j * n, c);
      }
    }
  }
  return 0;
}

int MatMul4x4_unrolled(double* A, double* B, double* C) {
  static const int n = 4;

  // First term in the sum (first col of A)
  __m256d a = _mm256_loadu_pd(A);
  __m256d b = _mm256_set1_pd(B[0]);
  __m256d c = _mm256_mul_pd(a, b);
  _mm256_storeu_pd(C, c);

  b = _mm256_set1_pd(B[n]);
  c = _mm256_mul_pd(a, b);
  _mm256_storeu_pd(C + n, c);

  b = _mm256_set1_pd(B[2 * n]);
  c = _mm256_mul_pd(a, b);
  _mm256_storeu_pd(C + 2 * n, c);

  b = _mm256_set1_pd(B[3 * n]);
  c = _mm256_mul_pd(a, b);
  _mm256_storeu_pd(C + 3 * n, c);

  // Second term in the sum
  a = _mm256_loadu_pd(A + n);      // 2nd column of A
  __m256d d = _mm256_loadu_pd(C);  // first column of C
  b = _mm256_set1_pd(B[1]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C, c);

  d = _mm256_loadu_pd(C + n);  // 2nd column of C
  b = _mm256_set1_pd(B[1 + n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + n, c);

  d = _mm256_loadu_pd(C + 2 * n);  // 3rd column of C
  b = _mm256_set1_pd(B[1 + 2 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 2 * n, c);

  d = _mm256_loadu_pd(C + 3 * n);  // 4th column of C
  b = _mm256_set1_pd(B[1 + 3 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 3 * n, c);

  // Third term in the sum (third column of A, third row of B)
  a = _mm256_loadu_pd(A + 2 * n);  // 3rd column of A
  d = _mm256_loadu_pd(C);          // first column of C
  b = _mm256_set1_pd(B[2]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C, c);

  d = _mm256_loadu_pd(C + n);
  b = _mm256_set1_pd(B[2 + n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + n, c);

  d = _mm256_loadu_pd(C + 2 * n);
  b = _mm256_set1_pd(B[2 + 2 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 2 * n, c);

  d = _mm256_loadu_pd(C + 3 * n);
  b = _mm256_set1_pd(B[2 + 3 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 3 * n, c);

  // Fourth term in the sum (fourth column of A, fourth row of B)
  a = _mm256_loadu_pd(A + 3 * n);
  d = _mm256_loadu_pd(C);  // first column of C
  b = _mm256_set1_pd(B[3]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C, c);

  d = _mm256_loadu_pd(C + n);
  b = _mm256_set1_pd(B[3 + n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + n, c);

  d = _mm256_loadu_pd(C + 2 * n);
  b = _mm256_set1_pd(B[3 + 2 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 2 * n, c);

  d = _mm256_loadu_pd(C + 3 * n);
  b = _mm256_set1_pd(B[3 + 3 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 3 * n, c);
  return 0;
}

int MatMul5x5_unrolled(double* A, double* B, double* C) {
  static const int n = 5;

  // Handle last element
  __m256d a =
      _mm256_setr_pd(A[4 + 0 * n], A[4 + 1 * n], A[4 + 2 * n], A[4 + 3 * n]);
  __m256d b = _mm256_load_pd(B + 4 * n);
  __m256d c = _mm256_mul_pd(a, b);
  __m128d c0 = _mm256_extractf128_pd(c, 0);
  __m128d c1 = _mm256_extractf128_pd(c, 1);
  __m128d sum1 = _mm_add_pd(c0, c1);
  __m128d sum2 = _mm_permute_pd(sum1, 1);
  __m128d sum = _mm_add_pd(sum1, sum2);
  _mm_storeu_pd(C, sum);
  C[4 + 4 * n] = C[0] + A[24] * B[24];

  // Handle last row
  a = _mm256_set1_pd(A[4 + 0 * n]);
  b = _mm256_setr_pd(B[0 + 0 * n], B[0 + 1 * n], B[0 + 2 * n], B[0 + 3 * n]);
  c = _mm256_mul_pd(a, b);

  a = _mm256_set1_pd(A[4 + 1 * n]);
  b = _mm256_setr_pd(B[1 + 0 * n], B[1 + 1 * n], B[1 + 2 * n], B[1 + 3 * n]);
  c = _mm256_fmadd_pd(a, b, c);

  a = _mm256_set1_pd(A[4 + 2 * n]);
  b = _mm256_setr_pd(B[2 + 0 * n], B[2 + 1 * n], B[2 + 2 * n], B[2 + 3 * n]);
  c = _mm256_fmadd_pd(a, b, c);

  a = _mm256_set1_pd(A[4 + 3 * n]);
  b = _mm256_setr_pd(B[3 + 0 * n], B[3 + 1 * n], B[3 + 2 * n], B[3 + 3 * n]);
  c = _mm256_fmadd_pd(a, b, c);

  a = _mm256_set1_pd(A[4 + 4 * n]);
  b = _mm256_setr_pd(B[4 + 0 * n], B[4 + 1 * n], B[4 + 2 * n], B[4 + 3 * n]);
  c = _mm256_fmadd_pd(a, b, c);
  _mm256_storeu_pd(C + 4, c);
  C[4 + 1 * n] = C[4 + 1];
  C[4 + 2 * n] = C[4 + 2];
  C[4 + 3 * n] = C[4 + 3];

  // First term in the sum (first col of A)
  a = _mm256_loadu_pd(A);
  b = _mm256_set1_pd(B[0]);
  c = _mm256_mul_pd(a, b);
  _mm256_storeu_pd(C, c);

  b = _mm256_set1_pd(B[n]);
  c = _mm256_mul_pd(a, b);
  _mm256_storeu_pd(C + n, c);

  b = _mm256_set1_pd(B[2 * n]);
  c = _mm256_mul_pd(a, b);
  _mm256_storeu_pd(C + 2 * n, c);

  b = _mm256_set1_pd(B[3 * n]);
  c = _mm256_mul_pd(a, b);
  _mm256_storeu_pd(C + 3 * n, c);

  b = _mm256_set1_pd(B[4 * n]);
  c = _mm256_mul_pd(a, b);
  _mm256_storeu_pd(C + 4 * n, c);

  // Second term in the sum
  a = _mm256_loadu_pd(A + n);      // 2nd column of A
  __m256d d = _mm256_loadu_pd(C);  // first column of C
  b = _mm256_set1_pd(B[1]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C, c);

  d = _mm256_loadu_pd(C + n);  // 2nd column of C
  b = _mm256_set1_pd(B[1 + n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + n, c);

  d = _mm256_loadu_pd(C + 2 * n);  // 3rd column of C
  b = _mm256_set1_pd(B[1 + 2 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 2 * n, c);

  d = _mm256_loadu_pd(C + 3 * n);  // 4th column of C
  b = _mm256_set1_pd(B[1 + 3 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 3 * n, c);

  d = _mm256_loadu_pd(C + 4 * n);  // 5th column of C
  b = _mm256_set1_pd(B[1 + 4 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 4 * n, c);

  // Third term in the sum (third column of A, third row of B)
  a = _mm256_loadu_pd(A + 2 * n);  // 3rd column of A
  d = _mm256_loadu_pd(C);          // first column of C
  b = _mm256_set1_pd(B[2]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C, c);

  d = _mm256_loadu_pd(C + n);
  b = _mm256_set1_pd(B[2 + n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + n, c);

  d = _mm256_loadu_pd(C + 2 * n);
  b = _mm256_set1_pd(B[2 + 2 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 2 * n, c);

  d = _mm256_loadu_pd(C + 3 * n);
  b = _mm256_set1_pd(B[2 + 3 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 3 * n, c);

  d = _mm256_loadu_pd(C + 4 * n);  // 5th column of C
  b = _mm256_set1_pd(B[2 + 4 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 4 * n, c);

  // Fourth term in the sum (fourth column of A, fourth row of B)
  a = _mm256_loadu_pd(A + 3 * n);
  d = _mm256_loadu_pd(C);  // first column of C
  b = _mm256_set1_pd(B[3]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C, c);

  d = _mm256_loadu_pd(C + n);
  b = _mm256_set1_pd(B[3 + n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + n, c);

  d = _mm256_loadu_pd(C + 2 * n);
  b = _mm256_set1_pd(B[3 + 2 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 2 * n, c);

  d = _mm256_loadu_pd(C + 3 * n);
  b = _mm256_set1_pd(B[3 + 3 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 3 * n, c);

  d = _mm256_loadu_pd(C + 4 * n);  // 5th column of C
  b = _mm256_set1_pd(B[3 + 4 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 4 * n, c);

  // Fifth term in the sum (fifth column of A, fifth row of B)
  a = _mm256_loadu_pd(A + 4 * n);
  d = _mm256_loadu_pd(C);  // first column of C
  b = _mm256_set1_pd(B[4]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C, c);

  d = _mm256_loadu_pd(C + n);
  b = _mm256_set1_pd(B[4 + n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + n, c);

  d = _mm256_loadu_pd(C + 2 * n);
  b = _mm256_set1_pd(B[4 + 2 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 2 * n, c);

  d = _mm256_loadu_pd(C + 3 * n);
  b = _mm256_set1_pd(B[4 + 3 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 3 * n, c);

  d = _mm256_loadu_pd(C + 4 * n);  // 5th column of C
  b = _mm256_set1_pd(B[4 + 4 * n]);
  c = _mm256_fmadd_pd(a, b, d);
  _mm256_storeu_pd(C + 4 * n, c);

  return 0;
}

int MatMul8x8_unrolled(double* a, double* b, double* c) {
  double a_0 = a[0];
  double a_1 = a[1];
  double a_2 = a[2];
  double a_3 = a[3];
  double a_4 = a[4];
  double a_5 = a[5];
  double a_6 = a[6];
  double a_7 = a[7];

  double tmp_0_0 = a_0 * b[0];
  double tmp_1_0 = a_1 * b[0];
  double tmp_2_0 = a_2 * b[0];
  double tmp_3_0 = a_3 * b[0];
  double tmp_4_0 = a_4 * b[0];
  double tmp_5_0 = a_5 * b[0];
  double tmp_6_0 = a_6 * b[0];
  double tmp_7_0 = a_7 * b[0];

  double tmp_0_1 = a_0 * b[8];
  double tmp_1_1 = a_1 * b[8];
  double tmp_2_1 = a_2 * b[8];
  double tmp_3_1 = a_3 * b[8];
  double tmp_4_1 = a_4 * b[8];
  double tmp_5_1 = a_5 * b[8];
  double tmp_6_1 = a_6 * b[8];
  double tmp_7_1 = a_7 * b[8];

  double tmp_0_2 = a_0 * b[16];
  double tmp_1_2 = a_1 * b[16];
  double tmp_2_2 = a_2 * b[16];
  double tmp_3_2 = a_3 * b[16];
  double tmp_4_2 = a_4 * b[16];
  double tmp_5_2 = a_5 * b[16];
  double tmp_6_2 = a_6 * b[16];
  double tmp_7_2 = a_7 * b[16];

  double tmp_0_3 = a_0 * b[24];
  double tmp_1_3 = a_1 * b[24];
  double tmp_2_3 = a_2 * b[24];
  double tmp_3_3 = a_3 * b[24];
  double tmp_4_3 = a_4 * b[24];
  double tmp_5_3 = a_5 * b[24];
  double tmp_6_3 = a_6 * b[24];
  double tmp_7_3 = a_7 * b[24];

  double tmp_0_4 = a_0 * b[32];
  double tmp_1_4 = a_1 * b[32];
  double tmp_2_4 = a_2 * b[32];
  double tmp_3_4 = a_3 * b[32];
  double tmp_4_4 = a_4 * b[32];
  double tmp_5_4 = a_5 * b[32];
  double tmp_6_4 = a_6 * b[32];
  double tmp_7_4 = a_7 * b[32];

  double tmp_0_5 = a_0 * b[40];
  double tmp_1_5 = a_1 * b[40];
  double tmp_2_5 = a_2 * b[40];
  double tmp_3_5 = a_3 * b[40];
  double tmp_4_5 = a_4 * b[40];
  double tmp_5_5 = a_5 * b[40];
  double tmp_6_5 = a_6 * b[40];
  double tmp_7_5 = a_7 * b[40];

  double tmp_0_6 = a_0 * b[48];
  double tmp_1_6 = a_1 * b[48];
  double tmp_2_6 = a_2 * b[48];
  double tmp_3_6 = a_3 * b[48];
  double tmp_4_6 = a_4 * b[48];
  double tmp_5_6 = a_5 * b[48];
  double tmp_6_6 = a_6 * b[48];
  double tmp_7_6 = a_7 * b[48];

  double tmp_0_7 = a_0 * b[56];
  double tmp_1_7 = a_1 * b[56];
  double tmp_2_7 = a_2 * b[56];
  double tmp_3_7 = a_3 * b[56];
  double tmp_4_7 = a_4 * b[56];
  double tmp_5_7 = a_5 * b[56];
  double tmp_6_7 = a_6 * b[56];
  double tmp_7_7 = a_7 * b[56];

  for (int j = 1; j < 8; ++j) {
    a_0 = a[0 + 8 * j];
    a_1 = a[1 + 8 * j];
    a_2 = a[2 + 8 * j];
    a_3 = a[3 + 8 * j];
    a_4 = a[4 + 8 * j];
    a_5 = a[5 + 8 * j];
    a_6 = a[6 + 8 * j];
    a_7 = a[7 + 8 * j];
    tmp_0_0 = fma(a_0, b[j + 0], tmp_0_0);
    tmp_1_0 = fma(a_1, b[j + 0], tmp_1_0);
    tmp_2_0 = fma(a_2, b[j + 0], tmp_2_0);
    tmp_3_0 = fma(a_3, b[j + 0], tmp_3_0);
    tmp_4_0 = fma(a_4, b[j + 0], tmp_4_0);
    tmp_5_0 = fma(a_5, b[j + 0], tmp_5_0);
    tmp_6_0 = fma(a_6, b[j + 0], tmp_6_0);
    tmp_7_0 = fma(a_7, b[j + 0], tmp_7_0);
    tmp_0_1 = fma(a_0, b[j + 8], tmp_0_1);
    tmp_1_1 = fma(a_1, b[j + 8], tmp_1_1);
    tmp_2_1 = fma(a_2, b[j + 8], tmp_2_1);
    tmp_3_1 = fma(a_3, b[j + 8], tmp_3_1);
    tmp_4_1 = fma(a_4, b[j + 8], tmp_4_1);
    tmp_5_1 = fma(a_5, b[j + 8], tmp_5_1);
    tmp_6_1 = fma(a_6, b[j + 8], tmp_6_1);
    tmp_7_1 = fma(a_7, b[j + 8], tmp_7_1);
    tmp_0_2 = fma(a_0, b[j + 16], tmp_0_2);
    tmp_1_2 = fma(a_1, b[j + 16], tmp_1_2);
    tmp_2_2 = fma(a_2, b[j + 16], tmp_2_2);
    tmp_3_2 = fma(a_3, b[j + 16], tmp_3_2);
    tmp_4_2 = fma(a_4, b[j + 16], tmp_4_2);
    tmp_5_2 = fma(a_5, b[j + 16], tmp_5_2);
    tmp_6_2 = fma(a_6, b[j + 16], tmp_6_2);
    tmp_7_2 = fma(a_7, b[j + 16], tmp_7_2);
    tmp_0_3 = fma(a_0, b[j + 24], tmp_0_3);
    tmp_1_3 = fma(a_1, b[j + 24], tmp_1_3);
    tmp_2_3 = fma(a_2, b[j + 24], tmp_2_3);
    tmp_3_3 = fma(a_3, b[j + 24], tmp_3_3);
    tmp_4_3 = fma(a_4, b[j + 24], tmp_4_3);
    tmp_5_3 = fma(a_5, b[j + 24], tmp_5_3);
    tmp_6_3 = fma(a_6, b[j + 24], tmp_6_3);
    tmp_7_3 = fma(a_7, b[j + 24], tmp_7_3);
    tmp_0_4 = fma(a_0, b[j + 32], tmp_0_4);
    tmp_1_4 = fma(a_1, b[j + 32], tmp_1_4);
    tmp_2_4 = fma(a_2, b[j + 32], tmp_2_4);
    tmp_3_4 = fma(a_3, b[j + 32], tmp_3_4);
    tmp_4_4 = fma(a_4, b[j + 32], tmp_4_4);
    tmp_5_4 = fma(a_5, b[j + 32], tmp_5_4);
    tmp_6_4 = fma(a_6, b[j + 32], tmp_6_4);
    tmp_7_4 = fma(a_7, b[j + 32], tmp_7_4);
    tmp_0_5 = fma(a_0, b[j + 40], tmp_0_5);
    tmp_1_5 = fma(a_1, b[j + 40], tmp_1_5);
    tmp_2_5 = fma(a_2, b[j + 40], tmp_2_5);
    tmp_3_5 = fma(a_3, b[j + 40], tmp_3_5);
    tmp_4_5 = fma(a_4, b[j + 40], tmp_4_5);
    tmp_5_5 = fma(a_5, b[j + 40], tmp_5_5);
    tmp_6_5 = fma(a_6, b[j + 40], tmp_6_5);
    tmp_7_5 = fma(a_7, b[j + 40], tmp_7_5);
    tmp_0_6 = fma(a_0, b[j + 48], tmp_0_6);
    tmp_1_6 = fma(a_1, b[j + 48], tmp_1_6);
    tmp_2_6 = fma(a_2, b[j + 48], tmp_2_6);
    tmp_3_6 = fma(a_3, b[j + 48], tmp_3_6);
    tmp_4_6 = fma(a_4, b[j + 48], tmp_4_6);
    tmp_5_6 = fma(a_5, b[j + 48], tmp_5_6);
    tmp_6_6 = fma(a_6, b[j + 48], tmp_6_6);
    tmp_7_6 = fma(a_7, b[j + 48], tmp_7_6);
    tmp_0_7 = fma(a_0, b[j + 56], tmp_0_7);
    tmp_1_7 = fma(a_1, b[j + 56], tmp_1_7);
    tmp_2_7 = fma(a_2, b[j + 56], tmp_2_7);
    tmp_3_7 = fma(a_3, b[j + 56], tmp_3_7);
    tmp_4_7 = fma(a_4, b[j + 56], tmp_4_7);
    tmp_5_7 = fma(a_5, b[j + 56], tmp_5_7);
    tmp_6_7 = fma(a_6, b[j + 56], tmp_6_7);
    tmp_7_7 = fma(a_7, b[j + 56], tmp_7_7);
  }
  c[0 + 0 * 8] = tmp_0_0;
  c[1 + 0 * 8] = tmp_1_0;
  c[2 + 0 * 8] = tmp_2_0;
  c[3 + 0 * 8] = tmp_3_0;
  c[4 + 0 * 8] = tmp_4_0;
  c[5 + 0 * 8] = tmp_5_0;
  c[6 + 0 * 8] = tmp_6_0;
  c[7 + 0 * 8] = tmp_7_0;

  c[0 + 1 * 8] = tmp_0_1;
  c[1 + 1 * 8] = tmp_1_1;
  c[2 + 1 * 8] = tmp_2_1;
  c[3 + 1 * 8] = tmp_3_1;
  c[4 + 1 * 8] = tmp_4_1;
  c[5 + 1 * 8] = tmp_5_1;
  c[6 + 1 * 8] = tmp_6_1;
  c[7 + 1 * 8] = tmp_7_1;

  c[0 + 2 * 8] = tmp_0_2;
  c[1 + 2 * 8] = tmp_1_2;
  c[2 + 2 * 8] = tmp_2_2;
  c[3 + 2 * 8] = tmp_3_2;
  c[4 + 2 * 8] = tmp_4_2;
  c[5 + 2 * 8] = tmp_5_2;
  c[6 + 2 * 8] = tmp_6_2;
  c[7 + 2 * 8] = tmp_7_2;

  c[0 + 3 * 8] = tmp_0_3;
  c[1 + 3 * 8] = tmp_1_3;
  c[2 + 3 * 8] = tmp_2_3;
  c[3 + 3 * 8] = tmp_3_3;
  c[4 + 3 * 8] = tmp_4_3;
  c[5 + 3 * 8] = tmp_5_3;
  c[6 + 3 * 8] = tmp_6_3;
  c[7 + 3 * 8] = tmp_7_3;

  c[0 + 4 * 8] = tmp_0_4;
  c[1 + 4 * 8] = tmp_1_4;
  c[2 + 4 * 8] = tmp_2_4;
  c[3 + 4 * 8] = tmp_3_4;
  c[4 + 4 * 8] = tmp_4_4;
  c[5 + 4 * 8] = tmp_5_4;
  c[6 + 4 * 8] = tmp_6_4;
  c[7 + 4 * 8] = tmp_7_4;

  c[0 + 5 * 8] = tmp_0_5;
  c[1 + 5 * 8] = tmp_1_5;
  c[2 + 5 * 8] = tmp_2_5;
  c[3 + 5 * 8] = tmp_3_5;
  c[4 + 5 * 8] = tmp_4_5;
  c[5 + 5 * 8] = tmp_5_5;
  c[6 + 5 * 8] = tmp_6_5;
  c[7 + 5 * 8] = tmp_7_5;

  c[0 + 6 * 8] = tmp_0_6;
  c[1 + 6 * 8] = tmp_1_6;
  c[2 + 6 * 8] = tmp_2_6;
  c[3 + 6 * 8] = tmp_3_6;
  c[4 + 6 * 8] = tmp_4_6;
  c[5 + 6 * 8] = tmp_5_6;
  c[6 + 6 * 8] = tmp_6_6;
  c[7 + 6 * 8] = tmp_7_6;

  c[0 + 7 * 8] = tmp_0_7;
  c[1 + 7 * 8] = tmp_1_7;
  c[2 + 7 * 8] = tmp_2_7;
  c[3 + 7 * 8] = tmp_3_7;
  c[4 + 7 * 8] = tmp_4_7;
  c[5 + 7 * 8] = tmp_5_7;
  c[6 + 7 * 8] = tmp_6_7;
  c[7 + 7 * 8] = tmp_7_7;
  return 0;
}
