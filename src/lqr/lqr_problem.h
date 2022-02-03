/**
 * @file lqr_problem.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief Defines the LQRProblem type
 * @version 0.1
 * @date 2022-01-30
 *
 * @copyright Copyright (c) 2022
 *
 * @addtogroup probdef Problem Definition
 * @{
 */
#pragma once

#include "lqr/lqr_data.h"

/**
 * @brief Describes an LQR problem with affine terms
 *
 * Internally, stores a vector of LQRData, one for each knot point, along with
 * the initial state and the horizon length.
 *
 * ## Construction and desctruction
 * The user can initialize and empty problem using ndlqr_NewLQRProblem(),
 * which again must be paired with a call to ndlqr_FreeLQRProblem().
 *
 * After the problem is initialized, it can be filled in from a vector LQRData using
 * ndlqr_InitializeLQRProblem().
 *
 */
typedef struct {
  int nhorizon;
  double* x0;
  LQRData** lqrdata;
} LQRProblem;

/**
 * @brief Initialize the problem with an initial state and the LQR data
 *
 * @param lqrproblem  An initialized LQRProblem
 * @param x0          Initial state vector. The data is copied into the problem.
 * @param lqrdata     A vector of LQR data. Each element is copied into the problem.
 * @return 0 if successful
 */
int ndlqr_InitializeLQRProblem(LQRProblem* lqrproblem, double* x0, LQRData** lqrdata);

/**
 * @brief Initialize a new LQRProblem data with unitialized data
 *
 * Must be paired with a call to ndlqr_FreeLQRProblem().
 *
 * @param nstates Length of the state vector
 * @param ninputs Number of control inputs
 * @param nhorizon Length of the horizon (i.e. number of knot points)
 * @return
 */
LQRProblem* ndlqr_NewLQRProblem(int nstates, int ninputs, int nhorizon);

/**
 * @brief Free the data stored by and LQRProblem
 *
 * Also frees the LQRProblem itself.
 *
 * @param lqrprob An initialized LQRProblem
 * @post lqrprob is NULL
 * @return 0 if successful
 */
int ndlqr_FreeLQRProblem(LQRProblem* lqrprob);

/**@} */
