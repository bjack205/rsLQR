/**
 * @file riccati_solver.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief Basic methods for creating and using the Riccati solver
 * @version 0.1
 * @date 2022-01-30
 *
 * @copyright Copyright (c) 2022
 *
 * @addtogroup riccati
 * @{
 */
#pragma once

#include "lqr/lqr_problem.h"
#include "matrix/matrix.h"

/**
 * @brief Solver that uses Riccati recursion to solve an LQR problem.
 *
 * Solves the generic LQR problem with affine terms using Riccati recursion and
 * a forward simulation of the linear dynamics. Assumes problems are of the following form:
 *
 * \f{align*}{
 * \underset{x_{1:N}, u_{1:N-1}}{\text{minimize}} &&& \frac{1}{2} x_N^T Q_N + x_N + q_N^T
 * x_N + \sum_{k-1}^{N-1} \frac{1}{2} x_k^T Q_k + x_k + q_k^T x_k + u_k^T R_k + u_k + r_k^T
 * u_k \\
 * \text{subject to} &&& x_{k+1} = A_k x_k + B_k u_k + f_k \\
 * &&& x_1 = x_\text{init}
 * \f}
 *
 * All the memory required by the solver
 * is initialized upon the creation of the solver to avoid any dynamic memory allocations
 * during the solve.
 *
 * ## Construction and destruction
 * Use ndlqr_NewRiccatiSolver() to initialize a new solver, which much be paired
 * with a single call to ndlqr_FreeRiccatiSolver() to free all of solver's memory.
 *
 * ## Typical Usage
 * Standard usage will typically look like the following:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * LQRProblem* lqrprob = ndlqr_ReadTestLQRProblem();  // your data here
 * RiccatiSolver* solver = ndlqr_NewRiccatiSolver(lqrprob);
 * ndlqr_SolveRiccati(solver);
 * ndlqr_PrintRiccatiSummary(solver);
 * double* soln = (double*) malloc(solver->nvars * sizeof(double));
 * ndlqr_CopyRiccatiSolution(solver, soln);
 * ndlqr_FreeRiccatiSolver();
 * ndlqr_FreeLQRProblem();
 * free(soln);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * ## Methods
 * - ndlqr_NewRiccatiSolver()
 * - ndlqr_FreeRiccatiSolver()
 * - ndlqr_PrintRiccatiSummary()
 * - ndlqr_GetRiccatiSolution()
 * - ndlqr_CopyRiccatiSolution()
 * - ndlqr_GetRiccatiSolveTimes()
 */
typedef struct {
  // clang-format off
  LQRProblem* prob;  ///< Problem data
  int nhorizon;  ///< length of the time horizon
  int nstates;   ///< size of state vector (n)
  int ninputs;   ///< number of control inputs (m)
  int nvars;     ///< total number of decision variables, including the dual variables
  double* data;  ///< pointer to the beginning of the single block of memory allocated by the solver
  Matrix* K;     ///< N-1 feedback gain matrices of size (m,n) 
  Matrix* d;     ///< N-1 feedforward gains of size (m,)
  Matrix* P;     ///< N cost-to-go Hessians of size (n,n)
  Matrix* p;     ///< N cost-to-go gradients of size (n,)
  Matrix* X;     ///< State trajectory. N vectors of size (n,)
  Matrix* U;     ///< Control trajectory. N-1 vectors of size (m,)
  Matrix* Y;     ///< Lagrange multipliers. N vectors of size (n,)
  Matrix* Qx;    ///< Gradient of the action-value function with respect to the state. N vectors of size (n,).
  Matrix* Qu;    ///< Gradient of the action-value function with respect to the control. N vectors of size (m,).
  Matrix* Qxx;   ///< Hessian of the action-value function with respect to the state. N vectors of size (n,n).
  Matrix* Qux;   ///< Cross-term Hessian of the action-value function. N vectors of size (m,n).
  Matrix* Quu;   ///< Hessian of the action-value function with respect to the control. N vectors of size (m,m).
  double t_solve_ms;          ///< Total solve time in milliseconds
  double t_backward_pass_ms;  ///< Time spent in the backward pass in milliseconds
  double t_forward_pass_ms;   ///< Time spent in the forward pass in milliseconds
  // clang-format on
} RiccatiSolver;

/**
 * @brief Initialize a new Riccati solver
 *
 * Create a new Riccati solver, provided the problem data given by lqrprob.
 *
 * @param lqrprob Contains all the data to describe the LQR problem to be solved.
 * @return An initialized Riccati solver.
 */
RiccatiSolver* ndlqr_NewRiccatiSolver(LQRProblem* lqrprob);

/**
 * @brief Free the memory for a Riccati solver
 *
 * @param solver Initialized Riccati solver.
 * @return 0 if successful
 * @post solver will be NULL
 */
int ndlqr_FreeRiccatiSolver(RiccatiSolver* solver);

/**
 * @brief Prints a summary of the solve
 *
 * ## Sample output
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.txt}
 * NDLQR Riccati Solve Summary
 * Solve time:    1.24 ms
 * Backward Pass: 1.13 ms (91.1 % of total)
 * Foward Pass:   0.11 ms (8.9 % of total)
 * Final error: 5.05696e-12
 * Final error after 2nd solve: 5.05696e-12
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * @pre ndlqr_SolveRiccati() has already been called
 * @param solver An initialized solver
 * @return 0 if successful
 */
int ndlqr_PrintRiccatiSummary(RiccatiSolver* solver);

/**
 * @brief Get the solution vector
 *
 * Returns the solution vector as a Matrix object, which is a simple wrapper around
 * a raw pointer which points to the data actually stored by the solver. The user must not
 * free the data, as it is owned by the solver. To get a solution vector owned by the
 * caller, use ndlqr_CopySolution() instead.
 *
 * ## Variable ordering
 * The variabled are ordered as follows:
 *
 * \f[
 * \begin{bmatrix}
 * \lambda_1^T & x_1^T & u_1^T & \lambda_2^T & \dots & x_{N-1}^T & u_{N-1}^T & \lambda_N^T &
 * x_N^T \end{bmatrix}^T \f]
 *
 * @pre ndlqr_SolveRiccati() has already been called
 * @param solver
 * @return
 */
Matrix ndlqr_GetRiccatiSolution(RiccatiSolver* solver);

int ndlqr_GetNumVarsRiccati(RiccatiSolver* solver);

/**
 * @brief Copies the solution to a user-supplied array
 *
 * See ndlqr_CopyRiccatiSolution() for variable ordering
 *
 * @pre ndlqr_SolveRiccati() has already been called
 * @param solver An initialized RiccatiSolver that has been solved but not freed
 * @param soln Destination for solution vector. Must have length at least equal to
 * [solver.nvars](@ref RiccatiSolver.nvars).
 * @return 0 if successful
 */
int ndlqr_CopyRiccatiSolution(RiccatiSolver* solver, double* soln);

/**
 * @brief Get the Riccati solve times
 *
 * Writes the solve times to the given pointers.
 *
 * @pre ndlqr_SolveRiccati() has already been called
 * @param solver An initialized RiccatiSolver
 * @param[out] t_solve Total solve time, in milliseconds
 * @param[out] t_bp    Backward pass time, in milliseconds
 * @param[out] t_fp    Forward pass time, in milliseconds
 * @return 0 if successful
 */
int ndlqr_GetRiccatiSolveTimes(RiccatiSolver* solver, double* t_solve, double* t_bp,
                               double* t_fp);

/**@} */
