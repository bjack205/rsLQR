/**
 * @file linalg_custom.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief Simple linear algebra routines
 * @version 0.1
 * @date 2022-01-31
 *
 * @copyright Copyright (c) 2022
 *
 * @ingroup LinearAlgebra
 * @{
 *
 */
#include "matrix.h"

/**
 * @brief Flag for a successful Cholesky decomposition,
 *        i.e. the matrix is positive definite.
 *
 */
static const int clap_kCholeskySuccess = 0;

/**
 * @brief Flag if a Cholesky decomposition fails,
 *        i.e. the matrix is not positive definite.
 *
 */
static const int clap_kCholeskyFail = -1;

/**
 * @brief Add two matrices of the same size, storing the result in @p B
 *
 * Performs the following operation:
 *
 * \f[
 *  B = B + \alpha A
 * \f]
 *
 * @param[in]    A any matrix of size (m,n)
 * @param[inout] B any matrix of size (m,n)
 * @param[in]    alpha scalar factor on A
 * @return       0 if successful
 */
int clap_MatrixAddition(Matrix* A, Matrix* B, double alpha);

/**
 * @brief Scale a matrix by a constant
 *
 * @param A     any matrix of non-zero size
 * @param alpha scalar multiplier
 * @return      0 if successful
 */
int clap_MatrixScale(Matrix* A, double alpha);

/**
 * @brief Matrix multiplication
 *
 * \f[
 * C = \alpha A B + \beta C
 * \f]
 *
 * @param[in]    A     Matrix of size (m,n)
 * @param[in]    B     Matrix of size (n,p)
 * @param[inout] C     Output matrix of size (m,p)
 * @param[in]    tA    Should @p A be transposed
 * @param[in]    tB    Should @p B be transposed
 * @param[in]    alpha scalar on the \f$ A B \f$ term
 * @param[in]    beta  scalar on the \f$ C \f$ term. Set to zero for pure
 *                     matrix multiplication.
 */
int clap_MatrixMultiply(Matrix* A, Matrix* B, Matrix* C, bool tA, bool tB, double alpha,
                        double beta);

/**
 * @brief A shortcut to perform transposed matrix multiplication
 *
 * Calculates
 * \f[
 *  C = A^T B
 * \f]
 *
 * @param[in]  A any matrix of size (n,m)
 * @param[in]  B any matrix of size (n,p)
 * @param[out] C any matrix of size (m,p)
 * @return
 */
int clap_MatrixTransposeMultiply(Matrix* A, Matrix* B, Matrix* C);

/**
 * @brief Matrix multiplication with a symmetric matrix A
 *
 * Perform the following computation
 * \f[
 * C = \alpha A B + \beta C
 * \f]
 *
 * For a symmetric matrix \f$ A \f$.
 *
 *
 * @param[in]    Asym
 * @param[in]    B
 * @param[inout] C
 * @param[in]    alpha
 * @param[in]    beta
 */
int clap_SymmetricMatrixMultiply(Matrix* Asym, Matrix* B, Matrix* C, double alpha,
                                 double beta);

/**
 * @brief Add a constant value to the diagonal of a matrix
 *
 * @param[inout] A     a matrix of size (n,m) where n <= m
 * @param[in]    alpha scalar to add to the diagonal
 * @return 0 if successful
 */
int clap_AddDiagonal(Matrix* A, double alpha);

/**
 * @brief Perform a Cholesky decomposition
 *
 * Performs a Cholesky decomposition on the square matrix @p A, storing the result in the
 * lower triangular portion of @p A.
 *
 * @param  A a square symmetric matrix
 * @return clap_kCholeskySuccess if successful, and clap_kCholeskyFail if not.
 */
int clap_CholeskyFactorize(Matrix* A);

/**
 * @brief Solve a linear system of equation with a precomputed Cholesky decomposition.
 *
 * @param[in]    A A square matrix whose Cholesky decomposition is stored in the lower
 *               triangular portion of the matrix
 * @param[inout] b The right-hand-side vector. Stores the solution upon completion of the
 *               function.
 * @return 0 if successful
 */
int clap_CholeskySolve(Matrix* A, Matrix* b);

/**
 * @brief Solve a linear system of equation for a lower triangular matrix
 *
 * Uses backsubstitution to solve a system of equations of the following form:
 * \f[
 *  L x = b
 * \f]
 * for a lower-triangular matrix \f$ L \f$, or
 * \f[
 *  L^T x = b
 * \f]
 * if @p istransposed is true.
 *
 *
 * @param[in]          L A lower-triangular matrix
 * @param[inout]       b The right-hand-side vector. Stores the solution upon completion.
 * @param istransposed Should @p L be transposed when solving the system of equations.
 * @return 0 if successful
 */
int clap_LowerTriBackSub(Matrix* L, Matrix* b, bool istransposed);

/**@} */
