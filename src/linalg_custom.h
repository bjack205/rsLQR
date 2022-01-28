#include "matrix.h"

static const int clap_kCholeskySuccess = 0;
static const int clap_kCholeskyFail = -1;

int clap_MatrixAddition(Matrix* A, Matrix* B, double alpha);
int clap_MatrixScale(Matrix* A, double alpha);
int clap_MatrixMultiply(Matrix* A, Matrix* B, Matrix* C, bool tA, bool tB, double alpha, double beta);
int clap_MatrixTransposeMultiply(Matrix* A, Matrix* B, Matrix* C);
int clap_SymmetricMatrixMultiply(Matrix* Asym, Matrix* B, Matrix* C, double alpha, double beta);
int clap_AddDiagonal(Matrix* A, double alpha);

int clap_CholeskyFactorize(Matrix* A);
int clap_CholeskySolve(Matrix* A, Matrix* b);
int clap_LowerTriBackSub(Matrix* L, Matrix* b, bool istransposed);
