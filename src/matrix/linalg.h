/**
 * @file linalg.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief Defines core linear algebra routines needed by the solvers
 * @version 0.1
 * @date 2022-01-31
 *
 * @copyright Copyright (c) 2022
 *
 * @addtogroup LinearAlgebra Linear Algebra
 * @{
 */
#pragma once

#include "matrix/matrix.h"

#ifdef USE_MKL
static const int kUseMKL = 1;
#else
static const int kUseMKL = 0;
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

#ifdef USE_BLAS
static const int kUseBLAS = 1;
#else
static const int kUseBLAS = 0;
#endif

/**
 * @brief Stores info about a Cholesky decomposition
 *
 * Thie provides information about a Cholesky decomposition, including which
 * trianglular portion the data is stored in, if the decomposition was successful,
 * which library was used to compute the decomposition, as well as a pointer to any
 * storage that library needed in addition to the matrix itself.
 *
 * ## Methods
 * - DefaultCholeskyInfo()
 * - FreeFactorization()
 */
typedef struct {
  char uplo;     ///< 'L' or 'U'
  int success;   ///< 0 if success, failure otherwise
  char lib;      ///< 'B' for BLAS, 'E' for eigen, 'I' for internal
  void* fact;    ///< pointer to Eigen data
  int is_freed;  ///< has the Eigen data been freed
} CholeskyInfo;

/**
 * @brief Construct a default CholeskyInfo object
 *
 * @return A CholeskyInfo with default values
 */
CholeskyInfo DefaultCholeskyInfo();

/**
 * @brief Frees any data stored by the external library.
 *
 * This should be called prior to creating a new factorization, which will
 * allocate new memory for the factorization.
 *
 * @param cholinfo
 */
void FreeFactorization(CholeskyInfo* cholinfo);

/**
 * @brief List of supported linear algebra libraries
 *
 */
enum MatrixLinearAlgebraLibrary {
  libBLAS = 0,
  libMKL = 1,
  libEigen = 2,
  libInternal = 3,
};

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
int MatrixAddition(Matrix* A, Matrix* B, double alpha);

/**
 * @brief Compute the Cholesky decomposition on the matrix @p A
 *
 * Not supported by all libraries. Prefer to use MatrixCholeskyFactorizeWithInfo().
 *
 * @param  mat A square, positive-definite matrix
 * @return 0 if the call was successful
 */
int MatrixCholeskyFactorize(Matrix* mat);

/**
 * @brief Compute the Cholesky decomposition on the matrix @p A
 *
 * @param[inout] mat A square, positive-definite matrix
 * @param cholinfo   CholeskyInfo object for storing info about the factorization
 * @post             @p cholinfo.success can be checked to see if the factorization was
 *                   successful
 * @return           0 if successful
 */
int MatrixCholeskyFactorizeWithInfo(Matrix* mat, CholeskyInfo* cholinfo);

/**
 * @brief Solve a linear system using a precomputed Cholesky factorization
 *
 * Overwrite the input vector @p b.
 * Prefer to use the more robust MatrixCholeskySolveWithInfo().
 *
 * @param[in]    A A square matrix whose Cholesky decomposition has already been computed.
 * @param[inout] b The right-hand-side vector. Stores the solution vector.
 * @return       0 if successful
 */
int MatrixCholeskySolve(Matrix* A, Matrix* b);

/**
 * @brief Solve a linear system using a precomputed Cholesky factorization
 *
 * @param[in]      A A square matrix whose Cholesky decomposition has already been computed.
 * @param[inout]   b The right-hand-side vector. Stores the solution vector.
 * @param cholinfo Information about the precomputed Cholesky factorization in @p A.
 * @return
 */
int MatrixCholeskySolveWithInfo(Matrix* A, Matrix* b, CholeskyInfo* cholinfo);

/**
 * @brief Matrix multiplication with scaling
 *
 * Perform the computation
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
void MatrixMultiply(Matrix* A, Matrix* B, Matrix* C, bool tA, bool tB, double alpha,
                    double beta);

/**
 * @brief Matrix multiplication with a symmetric matrix A
 *
 * Perform the following computation
 * \f[
 *  C = \alpha A B + \beta C
 * \f]
 * For a symmetric matrix \f$ A \f$.
 *
 * @param[in]    Asym
 * @param[in]    B
 * @param[inout] C
 * @param[in]    alpha
 * @param[in]    beta
 */
void MatrixSymmetricMultiply(Matrix* Asym, Matrix* B, Matrix* C, double alpha, double beta);

/**
 * @brief Copy just the diagonal element of @p src to the diagonal of @p dest
 *
 * @param dest Destination matrix
 * @param src  Source matrix
 */
void MatrixCopyDiagonal(Matrix* dest, Matrix* src);

/**
 * @brief Get the linear algebra library currently being used.
 *
 * Its value is determined by the build system and cannot be changed at runtime.
 *
 * @return The linear algebra library being used by the system.
 */
enum MatrixLinearAlgebraLibrary MatrixGetLinearAlgebraLibrary();

/**
 * @brief Prints which linear algebra library is being used to stdout
 *
 */
void MatrixPrintLinearAlgebraLibrary();

/**@}*/
