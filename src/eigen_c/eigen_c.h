#ifdef __cplusplus
extern "C" {
#endif

#define _eigen_MatrixMultiply(a, b, c)

void eigen_SetNumThreads(int n);
void eigen_InitParallel();

/**
 * @brief Add two vectors with scaling
 * 
 * `b = alpha * a + b`
 * 
 * Equivalent to the BLAS axpy routine.
 * 
 * @param n length of the inp
 * @param a Vector to add
 * @param b Destination vector
 * @param alpha Scaling on a
 */
void eigen_MatrixAddition(int n, double* a, double* b, double alpha);

void eigen_MatrixMultiply(int m, int n, int k, double* a, double* b, double* c,
                          bool tA, bool tB, double alpha, double beta);
void eigen_SymmetricMatrixMultiply(int n, int m, double* a, double* b,
                                   double* c);
void eigen_MatrixMultiply8x8(double* a, double* b, double* c);
void eigen_MatrixMultiply6x6(double* a, double* b, double* c);
void eigen_MatrixMultiply6x3(double* a, double* b, double* c);

int eigen_CholeskyFactorize(int n, double* a, void** fact);
void eigen_CholeskySolve(int n, int m, void* achol, double* b);
void eigen_FreeFactorization(void* achol);

#ifdef __cplusplus
}
#endif
