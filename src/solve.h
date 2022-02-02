/**
 * @file solve.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief Core methods for solving problems with nsLQR
 * @version 0.1
 * @date 2022-01-29
 * 
 * @copyright Copyright (c) 2022
 * 
 * @addtogroup rsLQR 
 * @{
 */
#include "solver.h"

/**
 * @brief Solve an LQR problem using rsLQR
 *
 * This is the main method of the package. Using the data already allocated in the 
 * solver, this first calculates the factorization of the matrix data and then 
 * solves for solution vector. Note that calling this method twice will not result
 * in the same data since it replaces the original right-hand-side data with the 
 * solution vector. 
 *
 * ## Extracting the solution
 * The solution can be retrieved after the solve by calling either ndlqr_GetSolution()
 * or ndlqr_CopySolution().
 *
 * ## Calling solve multiple times
 * If you want to call this method multiple times on the same data (e.g. when 
 * benchmarking solve times), use ndlqr_ResetNdData() between solves to reset the 
 * factorization data. Then re-initialize the solver with the problem data via
 * ndlqr_InitializeWithLQRProblem().
 * 
 * @param solver An nsLQR solver that has been initialized with the desired problem data.
 * @return 0 if successful.
 */
int ndlqr_Solve(NdLqrSolver* solver);

/**
 * @brief Return the solution vector
 * 
 * Returns the solution vector as a Matrix object, which is a simple wrapper around 
 * a raw pointer which points to the data actually stored by the solver. The user must not
 * free the data, as it is owned by the solver. To get a solution vector owned by the 
 * caller, use ndlqr_CopySolution() instead.
 *
 * ## Variable Ordering
 * The variables are ordered in the solution vector as follows:
 * 
 * \f[
 * \begin{bmatrix} 
 * \lambda_1^T & x_1^T & u_1^T & \lambda_2^T & \dots & x_{N-1}^T & u_{N-1}^T & \lambda_N^T & x_N^T 
 * \end{bmatrix}^T
 * \f]
 *
 * @pre ndlqr_Solve() has already been called
 * @param solver 
 * @return
 */
Matrix ndlqr_GetSolution(NdLqrSolver* solver);

/**
 * @brief Copies the solution vector to a user-supplied array
 * 
 * See ndlqr_GetSolution() for variable order.
 * 
 * @pre ndlqr_Solve() has already been called
 * @param solver The solver with the solution. 
 * @param soln Output location for the solution vector. Should have length at least equal to
 *             the output of ndlqr_GetNumVars().
 * @return
 */
int ndlqr_CopySolution(NdLqrSolver* solver, double* soln);

/**@} */
