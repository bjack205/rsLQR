#include "matrix/linalg_custom.h"

#include "math.h"
#include "stdio.h"

int clap_MatrixAddition(Matrix* A, Matrix* B, double alpha) {
  for (int i = 0; i < MatrixNumElements(A); ++i) {
    B->data[i] += alpha * A->data[i];
  }
  return 0;
}

int clap_MatrixScale(Matrix* A, double alpha) {
  for (int i = 0; i < MatrixNumElements(A); ++i) {
    A->data[i] *= alpha;
  }
  return 0;
}

int clap_MatrixMultiply(Matrix* A, Matrix* B, Matrix* C, bool tA, bool tB, double alpha,
                        double beta) {
  int n, m;
  if (tA) {
    n = A->cols;
    m = A->rows;
  } else {
    n = A->rows;
    m = A->cols;
  }
  int p = tB ? B->rows : B->cols;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      double* Cij = MatrixGetElement(C, i, j);
      *Cij *= beta;
      for (int k = 0; k < m; ++k) {
        double Aik = *MatrixGetElementTranspose(A, i, k, tA);
        double Bkj = *MatrixGetElementTranspose(B, k, j, tB);
        *Cij += alpha * Aik * Bkj;
      }
    }
  }
  return 0;
}

int clap_SymmetricMatrixMultiply(Matrix* Asym, Matrix* B, Matrix* C, double alpha,
                                 double beta) {
  int n, m;
  bool tA = false;
  bool tB = false;
  if (tA) {
    n = Asym->cols;
    m = Asym->rows;
  } else {
    n = Asym->rows;
    m = Asym->cols;
  }
  int p = tB ? B->rows : B->cols;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      double* Cij = MatrixGetElement(C, i, j);
      *Cij *= beta;
      for (int k = 0; k < m; ++k) {
        int row = i;
        int col = k;
        if (i < k) {
          row = k;
          col = i;
        }
        double Aik = *MatrixGetElement(Asym, row, col);
        double Bkj = *MatrixGetElement(B, k, j);
        *Cij += alpha * Aik * Bkj;
      }
    }
  }
  return 0;
  return 0;
}

int clap_AddDiagonal(Matrix* A, double alpha) {
  int n = A->rows;
  for (int i = 0; i < n; ++i) {
    double* Aii = MatrixGetElement(A, i, i);
    *Aii += alpha;
  }
  return 0;
}

int clap_CholeskyFactorize(Matrix* A) {
  int n = A->rows;
  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < j; ++k) {
      for (int i = j; i < n; ++i) {
        double* Aij = MatrixGetElement(A, i, j);
        double Aik = *MatrixGetElement(A, i, k);
        double Ajk = *MatrixGetElement(A, j, k);
        *Aij -= Aik * Ajk;
      }
    }
    double Ajj = *MatrixGetElement(A, j, j);
    if (Ajj <= 0) {
      return clap_kCholeskyFail;
    }
    double ajj = sqrt(Ajj);

    for (int i = j; i < n; ++i) {
      double* Aij = MatrixGetElement(A, i, j);
      *Aij /= ajj;
    }
  }
  return clap_kCholeskySuccess;
}

int clap_LowerTriBackSub(Matrix* L, Matrix* b, bool istransposed) {
  int n = b->rows;
  int m = b->cols;
  for (int j_ = 0; j_ < n; ++j_) {
    int j = istransposed ? n - j_ - 1 : j_;
    for (int k = 0; k < m; ++k) {
      double* xjk = MatrixGetElement(b, j, k);
      double Ljj = *MatrixGetElement(L, j, j);
      *xjk /= Ljj;

      for (int i_ = j_ + 1; i_ < n; ++i_) {
        int i = istransposed ? i_ - (j_ + 1) : i_;
        double* xik = MatrixGetElement(b, i, k);
        double Lij = *MatrixGetElementTranspose(L, i, j, istransposed);
        *xik -= Lij * (*xjk);
      }
    }
  }
  return 0;
}

int clap_CholeskySolve(Matrix* L, Matrix* b) {
  clap_LowerTriBackSub(L, b, 0);
  clap_LowerTriBackSub(L, b, 1);
  return 0;
}
