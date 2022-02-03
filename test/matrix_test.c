#include "matrix.h"

#include "test/minunit.h"

mu_test_init

    int
    TestNewMatrix() {
  Matrix mat = NewMatrix(5, 4);
  mu_assert(mat.rows == 5);
  mu_assert(mat.cols == 4);
  mu_assert(MatrixNumElements(&mat) == 20);
  FreeMatrix(&mat);
  return 1;
}

int SetConst() {
  Matrix mat = NewMatrix(3, 4);
  MatrixSetConst(&mat, 5.0);
  for (int i = 0; i < 12; ++i) {
    mu_assert(mat.data[i] == 5);
  }
  MatrixSetConst(&mat, -4.2);
  for (int i = 0; i < 12; ++i) {
    mu_assert(mat.data[i] == -4.2);
  }
  FreeMatrix(&mat);
  return 1;
}

int GetIndex() {
  Matrix mat = NewMatrix(3, 4);
  for (int i = 0; i < 12; ++i) {
    mat.data[i] = i;
  }
  mu_assert(*MatrixGetElement(&mat, 0, 0) == 0);
  mu_assert(*MatrixGetElement(&mat, 1, 0) == 1);
  mu_assert(*MatrixGetElement(&mat, 0, 1) == 3);
  mu_assert(*MatrixGetElement(&mat, 2, 3) == 11);
  FreeMatrix(&mat);
  return 1;
}

int TestPrintMatrix() {
  Matrix mat = NewMatrix(3, 4);
  for (int i = 0; i < 12; ++i) {
    mat.data[i] = i;
  }
  PrintMatrix(&mat);
  FreeMatrix(&mat);
  return 1;
}

int CopyMatrix() {
  Matrix src = NewMatrix(10, 12);
  for (int i = 0; i < MatrixNumElements(&src); ++i) {
    src.data[i] = i * i - sqrt(i * 2.1);
  }
  Matrix dest = NewMatrix(10, 12);
  MatrixCopy(&dest, &src);
  for (int i = 0; i < MatrixNumElements(&src); ++i) {
    mu_assert(dest.data[i] == i * i - sqrt(i * 2.1));
  }
  FreeMatrix(&src);
  FreeMatrix(&dest);
  return 1;
}

int CopyTranspose() {
  Matrix src = NewMatrix(3, 4);
  for (int i = 0; i < MatrixNumElements(&src); ++i) {
    src.data[i] = i * i - sqrt(i * 2.1);
  }
  Matrix dest = NewMatrix(4, 3);
  MatrixCopyTranspose(&dest, &src);
  for (int i = 0; i < src.rows; ++i) {
    for (int j = 0; j < src.rows; ++j) {
      mu_assert(*MatrixGetElement(&src, i, j) == *MatrixGetElement(&dest, j, i));
    }
  }
  FreeMatrix(&src);
  FreeMatrix(&dest);
  return 1;
}

void AllTests() {
  mu_run_test(TestNewMatrix);
  mu_run_test(SetConst);
  mu_run_test(GetIndex);
  mu_run_test(TestPrintMatrix);
  mu_run_test(CopyMatrix);
  mu_run_test(CopyTranspose);
}

mu_test_main
