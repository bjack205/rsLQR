/**
 * @file nddata.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief Defines the core storage types used by the rsLQR solver
 * @version 0.1
 * @date 2022-01-30
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include "matrix.h"
#include "lqr_data.h"

/**
 * @brief A chunk of memory for a single time step
 *
 * Stores a matrix of size (2n+m,n), divided into blocks:
 *
 * \f[
 * \begin{bmatrix} \Lambda \\ X \\ U \end{bmatrix}
 * \f]
 * 
 * which correspond to the NdFactor.lambda, NdFactor.state, and NdFactor.input 
 * fields, which can also be extracted using the following methods:
 * - ndlqr_GetLambdaFactor()
 * - ndlqr_GetStateFactor()
 * - ndlqr_GetInputFactor()
 * 
 * Each of which return a Matrix of the corresponding size. Each of these blocks can have 
 * an arbitrary width, since these can represent either data from the KKT matrix or the 
 * right-hand-side vector(s).
 *
 * Internally, the solver stores arrays of these objects, which allow the solver
 * to extract chunks out of the original matrix data, by time step.
 */
typedef struct {
  Matrix lambda;  ///< (n,w) block for the dual variables
  Matrix state;   ///< (n,w) block for the state variables
  Matrix input;   ///< (m,w) block for the control input variables
} NdFactor;

Matrix ndlqr_GetLambdaFactor(NdFactor* factor);
Matrix ndlqr_GetStateFactor(NdFactor* factor);
Matrix ndlqr_GetInputFactor(NdFactor* factor);

/**
 * @brief Core storage container for the rsLQR solver
 *
 * Represents an array of memory blocks, arranged as follows:
 * 
 * \f[
 * \begin{bmatrix} 
 * F_1^{(1)} & X_1^{(2)} & \dots & X_1^{(K)} \\
 * F_2^{(1)} & X_2^{(2)} & \dots & X_2^{(K)} \\
 * \vdots    & \vdots    & \ddots & \vdots \\
 * F_{N-1}^{(1)} & X_{N-1}^{(2)} & \dots & X_{N-1}^{(K)} \\
 * \end{bmatrix}
 * \f]
 * 
 * Each \f$ F \f$ is a NdFactor further dividing this memory into chunks of size `(n,w)` or 
 * `(m,w)`, where the width `w` is equal to NdData.width. Each block is stored as an 
 * individual Matrix, which stores the data column-wise. This keeps the data for a single 
 * block together in one contiguous block of memory. The entire block of memory for all of 
 * the factors is allocated as one large block (with pointer NdData.data).
 *
 * In the solver, this is used to represent both the KKT matrix data and the right-hand-side
 * vector. When storing the matrix data, each column represents a level of the binary tree.
 * The current implementation only allows for a single right-hand-side vector, so that 
 * when `width` is passed to the initializer, it only creates a single column of factors.
 * Future modifications could alternatively make the right-hand-side vector the last column
 * in the matrix data, as suggested in the original paper.
 *
 * ## Methods
 * - ndlqr_NewNdData()
 * - ndlqr_FreeNdData()
 * - ndlqr_GetNdFactor()
 * - ndlqr_ResetNdFactor()
 */
typedef struct {
  int nstates;    ///< size of state vector
  int ninputs;    ///< number of control inputs
  int nsegments;  ///< number of segments, or one less than the length of the horizon
  int depth;      ///< number of columns of factors to store
  int width;      ///< width of each factor. Will be `n` for matrix data and typically 1 for the right-hand-side vector.
  double* data;       ///< pointer to entire chunk of allocated memory
  NdFactor* factors;  ///< (nsegments, depth) array of factors. Stored in column-order.
} NdData;

/**
 * @brief Initialize the NdData structure
 *
 * Note this allocates a large block of memory. This should be followed by a single 
 * call to ndlqr_FreeNdData().
 * 
 * @param nstates Number of variables in the state vector
 * @param ninputs Number of control inputs
 * @param nhorizon Length of the time horizon
 * @param width With of each factor. Should be `nstates` for KKT matrix data, 
                or 1 for the right-hand side vector.
 * @return The initialized NdData structure
 */
NdData* ndlqr_NewNdData(int nstates, int ninputs, int nhorizon, int width);

/**
 * @brief Frees the memory allocated in an NdData structure
 * 
 * @param nddata Initialized NdData structure
 * @return 0 if successful
 */
int ndlqr_FreeNdData(NdData* nddata);

/**
 * @brief Retrieve an individual NdFactor out of the NdData
 *
 * Retrieves a block of memory out of NdData, stored as an NdFactor. Typical usage
 * will look like:
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 * NdData* ndata = ndlqr_NewNdData(nstates, ninputs, nhorizon, nstates);
 * NdFactor* factor;
 * int index = 0;  // set to desired index
 * int level = 1;  // set to desired level
 * ndlqr_GetNdFactor(nddata, index, level, &factor);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 * 
 * @param nddata Storage location for the desired block of memory
 * @param index Time step of the factor to extract
 * @param level Level (or column in the NdData) of the desired factor 
 * @param factor Storage location for the factor. 
 * @return 0 if successful
 */
int ndlqr_GetNdFactor(NdData* nddata, int index, int level, NdFactor** factor);

/**
 * @brief Resets all of the memory for an NdData to zero.
 * 
 * @param nddata Initialized NdData structure
 */
void ndlqr_ResetNdData(NdData* nddata);
