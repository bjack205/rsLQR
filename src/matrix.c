#include "matrix.h"

#ifndef PRECISION
#define PRECISION 5
#endif

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

Matrix NewMatrix(int rows, int cols) {
  double* data = (double*) malloc(rows * cols * sizeof(double));
  Matrix mat = {rows, cols, data};
  return mat; 
}

int MatrixSetConst(Matrix* mat, double val) {
  if (!mat) return -1;
  for (int i = 0; i < MatrixNumElements(mat); ++i) {
    mat->data[i] = val;
  }
  return 0;
}

int FreeMatrix(Matrix* mat) {
  if (mat) {
    if (mat->data) {
      free(mat->data);
      return 0;
    }
  }
  return -1;
}

int MatrixNumElements(const Matrix* mat) {
  if (!mat) return -1;
  return mat->rows * mat->cols;
}

int MatrixGetLinearIndex(const Matrix* mat, int row, int col) {
  if (!mat) return -1;
  return row + mat->rows * col;
}

double* MatrixGetElement(const Matrix* mat, int row, int col) {
  if (!mat) return NULL;
  return mat->data + MatrixGetLinearIndex(mat, row, col);
}

double* MatrixGetElementTranspose(const Matrix* mat, int row, int col, bool istranposed) {
  if (!istranposed) {
    return MatrixGetElement(mat, row, col);
  } else {
    return MatrixGetElement(mat, col, row);
  }
}

int MatrixSetElement(Matrix* mat, int row, int col, double val) {
  if (!mat) return -1;
  int linear_index = MatrixGetLinearIndex(mat, row, col);
  if (linear_index >= 0) {
    mat->data[linear_index] = val;
  } else {
    return -1;
  }
  return 1;
}

int MatrixCopy(Matrix* dest, Matrix* src) {
  if (!dest || !src) return -1;
  if ((dest->rows != src->rows) || (dest->cols != src->cols)) {
    fprintf(stderr, "Can't copy matrices of different sizes.\n");
    return -1;
  }
  memcpy(dest->data, src->data, MatrixNumElements(dest) * sizeof(double));
  return 0;
}

int MatrixCopyTranspose(Matrix* dest, Matrix* src) {
  if (!dest || !src) return -1;
  if ((dest->rows != src->cols) || (dest->cols != src->rows)) {
    fprintf(stderr, "Matrix sizes are not transposes of each other. Got (%d,%d) and (%d,%d).\n", 
        dest->rows, dest->cols, src->rows, src->cols);
    return -1;
  }
  for (int i = 0; i < dest->rows; ++i) {
    for (int j = 0; j < dest->cols; ++j) {
      int dest_index = MatrixGetLinearIndex(dest, i, j);
      int src_index = MatrixGetLinearIndex(src, j, i);
      dest->data[dest_index] = src->data[src_index];
    }
  }
  return 0;
}

int MatrixScaleByConst(Matrix* mat, double alpha) {
  if (!mat) return -1;
  for (int i = 0; i < MatrixNumElements(mat); ++i) {
    mat->data[i] *= alpha;
  }
  return 0;
}

double MatrixNormedDifference(Matrix* A, Matrix* B) {
  if (!A || !B) return INFINITY;
  if ((A->rows != B->rows) || (A->cols != B->cols)) {
    fprintf(stderr, "Can't compare matrices of different sizes. Got (%d,%d) and (%d,%d)\n",
        A->rows, A->cols, B->rows, B->cols);
    return INFINITY;
  }

  double diff = 0;
  for (int i = 0; i < MatrixNumElements(A); ++i) {
    double d = A->data[i] - B->data[i];
    diff += d * d;
  }
  return sqrt(diff);
}

int MatrixFlatten(Matrix* mat) {
  if (!mat) return -1;
  int size = MatrixNumElements(mat);
  mat->rows = size;
  mat->cols = 1;
  return 0;
}

int MatrixFlattenToRow(Matrix* mat) {
  if (!mat) return -1;
  int size = MatrixNumElements(mat);
  mat->rows = 1;
  mat->cols = size;
  return 0;
}

int PrintMatrix(const Matrix* mat) {
  if (!mat) return -1;
  for (int row = 0; row < mat->rows; ++row) {
    for (int col = 0; col < mat->cols; ++col) {
      printf("% 6.*g ", PRECISION, *MatrixGetElement(mat, row, col));
    }
    printf("\n");
  }
  return 0;
}

int PrintRowVector(const Matrix* mat) {
  if (!mat) return -1;
  printf("[ ");
  for (int i = 0; i < MatrixNumElements(mat); ++i) {
    printf("% 6.*g ", PRECISION, mat->data[i]);
  }
  printf("]\n");
  return 0;
}
