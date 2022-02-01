/**
 * @file riccati_solve.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief Core methods for the Riccati solve
 * @version 0.1
 * @date 2022-01-30
 * 
 * @copyright Copyright (c) 2022
 * 
 * @addtogroup riccati Riccati Solver 
 * @{
 */
#pragma once

#include "riccati_solver.h"

#include "linalg.h"

/**
 * @brief Solve the LQR problem using Riccati recursion and a forward simulation of the 
 *        linear dynamics.
 * 
 * @param solver An initialized RiccatiSolver
 * @return 0 if successful
 */
int ndlqr_SolveRiccati(RiccatiSolver* solver);

/**
 * @brief Run the Riccati solver backward pass
 *
 * Uses backward Riccati recursion to solve for the feedback and feedforward LQR gains
 * along the trajectory. Also computes the quadratic cost-to-go and the expansions of the 
 * action-value function. All the data is stored in the solver.
 * 
 * @param solver An initialized RiccatiSolver
 * @return 0 if successful
 */
int ndlqr_BackwardPass(RiccatiSolver* solver);

/**
 * @brief Run the Riccati forward pass to solve for the solution vector 
 * 
 * Computes the solution vector by simulating the linear dynamics forward using the 
 * the feedback law \f$ u = -K x + d \f$.
 *
 * @pre The LQR gains must have been computed using ndlqr_BackwardPass()
 * @param solver An initialized RiccatiSolver
 * @return 0 if successful
 */
int ndlqr_ForwardPass(RiccatiSolver* solver);

/**@} */
