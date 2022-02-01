/**
 * @file lqr_data.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief LQRData type
 * @version 0.1
 * @date 2022-01-31
 * 
 * @copyright Copyright (c) 2022
 * 
 * @addtogroup probdef 
 * @{ 
 */
#pragma once

#include "matrix.h"

/**
 * @brief Holds the data for a single time step of LQR
 * 
 * Stores the \f$ Q, R, q, r, c \f$ values for the cost function:
 * \f[
 * \frac{1}{2} x^T Q x + q^T x + \frac{1}{2} u^T R u + r^T r + c
 * \f]
 * 
 * and the \f$ A, B, d \f$ values for the dynamics:
 * \f[
 * x_{k+1} = A x_k + B u_k + d
 * \f]
 * 
 * ## Construction and destruction
 * A new LQRData object is constructed using ndlqr_NewLQRData(), which must be 
 * freed with a call to ndlqr_FreeLQRData().
 * 
 * ## Methods
 * - ndlqr_NewLQRData()
 * - ndlqr_FreeLQRData()
 * - ndlqr_InitializeLQRData()
 * - ndlqr_CopyLQRData()
 * - ndlqr_PrintLQRData()
 * 
 * ## Getters
 * The follow methods return a Matrix object wrapping the data from an LQRData object.
 * The user should NOT call FreeMatrix() on this data since it is owned by the LQRData
 * object.
 * - ndlqr_GetA()
 * - ndlqr_GetB()
 * - ndlqr_Getd()
 * - ndlqr_GetQ()
 * - ndlqr_GetR()
 * - ndlqr_Getq()
 * - ndlqr_Getr()
 * 
 */
typedef struct {
  int nstates;
  int ninputs;
  double* Q;
  double* R;
  double* q;
  double* r;
  double* c;
  double* A;
  double* B;
  double* d;
} LQRData;

/**
 * @brief Copy data into an initialized LQRData structure
 * 
 * Does not allocate any new memory.
 * 
 * @param lqrdata Initialized LQRData struct
 * @param Q       State cost Hessian
 * @param R       Control cost Hessian
 * @param q       State cost affine term
 * @param r       Control cost affine term
 * @param c       Constant cost term
 * @param A       Dynamics state matrix
 * @param B       Dynamics control matrix
 * @param d       Dynamics affine term
 * @return 0 if successful
 */
int ndlqr_InitializeLQRData(LQRData* lqrdata, double* Q, double* R, double* q, 
                            double* r, double c, double* A, double* B, double* d);

/**
 * @brief Allocate memory for a new LQRData structure
 * 
 * Must be paired with a single call to ndlqr_FreeLQRData().
 * 
 * @param nstates Length of the state vector
 * @param ninputs Number of control inputs
 * @return 0 if successful
 */
LQRData* ndlqr_NewLQRData(int nstates, int ninputs);

/**
 * @brief Free the memory for and LQRData object
 * 
 * @param lqrdata Initialized LQRData object
 * @post lqrdata = NULL
 * @return 0 if successful
 */
int ndlqr_FreeLQRData(LQRData* lqrdata);

/**
 * @brief Copies one LQRData object to another 
 * 
 * The two object must have equivalent dimensionality.
 * 
 * @param dest Copy destination
 * @param src  Source data
 * @return 0 if successful
 */
int ndlqr_CopyLQRData(LQRData* dest, LQRData* src);

Matrix ndlqr_GetA(LQRData* lqrdata); ///< @brief Get (n,n) state transition matrix
Matrix ndlqr_GetB(LQRData* lqrdata); ///< @brief Get (n,m) control input matrix
Matrix ndlqr_Getd(LQRData* lqrdata); ///< @brief Get (n,) affine dynamice term
Matrix ndlqr_GetQ(LQRData* lqrdata); ///< @brief Get state cost Hessian
Matrix ndlqr_GetR(LQRData* lqrdata); ///< @brief Get control cost Hessian
Matrix ndlqr_Getq(LQRData* lqrdata); ///< @brief Get affine state cost
Matrix ndlqr_Getr(LQRData* lqrdata); ///< @brief Get affine control cost

/**
 * @brief Prints the data contained in LQRData
 * 
 * Cost data is printed in rows and dynamics data is printed as normal matrices.
 * 
 * @param lqrdata 
 */
void ndlqr_PrintLQRData(LQRData* lqrdata);

/**@} */
