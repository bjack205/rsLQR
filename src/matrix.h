#pragma once

#include <stdbool.h>

typedef struct {
  int rows;
  int cols;
  double* data;
} Matrix;

Matrix NewMatrix(int rows, int cols);
int MatrixSetConst(Matrix* mat, double val);
int FreeMatrix(Matrix* mat);

int MatrixNumElements(const Matrix* mat);
int MatrixGetLinearIndex(const Matrix* mat, int row, int col);
double* MatrixGetElement(const Matrix* mat, int row, int col);
double* MatrixGetElementTranspose(const Matrix* mat, int row, int col, bool istranposed);
int MatrixSetElement(Matrix* mat, int row, int col, double val);
int MatrixCopy(Matrix* dest, Matrix* src);
int MatrixCopyTranspose(Matrix* dest, Matrix* src);

int MatrixScaleByConst(Matrix* mat, double alpha);
double MatrixNormedDifference(Matrix* A, Matrix* B);

int MatrixFlatten(Matrix* mat);
int MatrixFlattenToRow(Matrix* mat);

int PrintMatrix(const Matrix* mat);
int PrintRowVector(const Matrix* mat);

