/**
 * @file solver.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief High-level methods for the nsLQR solver.
 * @version 0.1
 * @date 2022-01-29
 *
 * @copyright Copyright (c) 2022
 *
 * @addtogroup rsLQR
//  * @{
 */
#pragma once

#include "rslqr/binary_tree.h"
#include "rslqr/cholesky_factors.h"
#include "linalg.h"
#include "lqr/lqr_problem.h"
#include "rslqr/nddata.h"

/**
 * @brief A struct describing how long each part of the solve took, in milliseconds.
 *
 * ## Methods
 * - ndlqr_NewNdLqrProfile()
 * - ndlqr_ResetProfile()
 * - ndlqr_CopyProfile()
 * - ndlqr_PrintProfile()
 * - ndlqr_CompareProfile()
 */
typedef struct {
  double t_total_ms;
  double t_leaves_ms;
  double t_products_ms;
  double t_cholesky_ms;
  double t_cholsolve_ms;
  double t_shur_ms;
  int num_threads;
} NdLqrProfile;

/**
 * @brief Create a profile initialized with zeros
 */
NdLqrProfile ndlqr_NewNdLqrProfile();

/**
 * @brief Reset the profile to its initialized state
 *
 * @param prof A profile
 */
void ndlqr_ResetProfile(NdLqrProfile* prof);

/**
 * @brief Copy the profile information to a new profile
 *
 * @param dest New location for data. Existing data will be overwritten.
 * @param src Data to be copied.
 */
void ndlqr_CopyProfile(NdLqrProfile* dest, NdLqrProfile* src);

/**
 * @brief Print a summary fo the profile
 *
 * @param profile
 */
void ndlqr_PrintProfile(NdLqrProfile* profile);

/**
 * @brief Compare two profiles, printing the comparison to stdout
 *
 * @param base The baseline profile
 * @param prof The "new" profile
 */
void ndlqr_CompareProfile(NdLqrProfile* base, NdLqrProfile* prof);

/**
 * @brief Main solver for rsLQR
 *
 * Core struct for solving problems with rsLQR. Allocates all the required memory
 * up front to avoid any dynamic memory allocations at runtime. Right now, the
 * horizon length is required to be a power of 2 (e.g. 32,64,128,256,etc.).
 *
 * ## Construction and destruction
 * Use ndlqr_NewNdLqrSolver() to initialize a new solver. This should always be
 * paired with a single call to ndlqr_FreeNdLqrSolver().
 *
 * ## Typical Usage
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * LQRProblem* lqrprob = ndlqr_ReadTestLQRProblem();  // your data here
 * int nstates = lqrprob->lqrdata[0]->nstates;
 * int ninputs = lqrprob->lqrdata[0]->ninputs;
 * int nhorizon = lqrprob->nhorizon;
 *
 * NdLqrSolver* solver = ndlqr_NewNdLqrSolver(nstates, ninputs, nhorizon);
 * ndlqr_InitializeWithLQRProblem(lqrprob, solver);
 * ndlqr_Solve(solver);
 * ndlqr_PrintSolveSummary();
 * ndlqr_FreeLQRProblem(lqrprob);
 * ndlqr_FreeNdLqrSolver(solver);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * ## Methods
 * - ndlqr_NewNdLqrSolver()
 * - ndlqr_FreeNdLqrSolver()
 * - ndlqr_InitializeWithLQRProblem()
 * - ndlqr_Solve()
 * - ndlqr_ResetSolver()
 * - ndlqr_GetNumVars()
 * - ndlqr_SetNumThreads()
 * - ndlqr_PrintSolveProfile()
 * - ndlqr_GetProfile()
 */
typedef struct {
  int nstates;   ///< size of state vector
  int ninputs;   ///< number of control inputs
  int nhorizon;  ///< length of the time horizon
  int depth;     ///< depth of the binary tree
  int nvars;     ///< number of decision variables (size of the linear system)
  OrderedBinaryTree tree;
  Matrix* diagonals;  ///< (nhorizon,2) array of diagonal blocks (Q,R)
  NdData* data;       ///< original matrix data
  NdData* fact;       ///< factorization
  NdData* soln;       ///< solution vector (also the initial RHS)
  NdLqrCholeskyFactors* cholfacts;
  double solve_time_ms;  ///< total solve time in milliseconds.
  double linalg_time_ms;
  NdLqrProfile profile;
  int num_threads;  ///< Number of threads used by the solver.
} NdLqrSolver;

/**
 * @brief Create a new solver, allocating all the required memory.
 *
 * Must be followed by a later call to ndlqr_FreeNdLqrSolver().
 *
 * @param nstates Number of elements in the state vector
 * @param ninputs Number of control inputs
 * @param nhorizon Length of the time horizon. Must be a power of 2.
 * @return A pointer to the new solver
 */
NdLqrSolver* ndlqr_NewNdLqrSolver(int nstates, int ninputs, int nhorizon);

/**
 * @brief Deallocates the memory for the solver.
 *
 * @param solver An initialized solver.
 * @return 0 if successful.
 * @post solver == NULL
 */
int ndlqr_FreeNdLqrSolver(NdLqrSolver* solver);

/**
 * @brief Initialize the solver with data from an LQR Problem.
 *
 * @pre Solver has already been initialized via ndlqr_NewNdLqrSolver()
 * @param lqrprob An initialized LQR problem with the data to be be solved.
 * @param solver An initialized solver.
 * @return 0 if successful
 */
int ndlqr_InitializeWithLQRProblem(const LQRProblem* lqrprob, NdLqrSolver* solver);

/**
 * @brief Resets the rsLQR solver
 *
 * Resets all of the data in the solver to how it was when it was first initialized.
 *
 * @param solver
 */
void ndlqr_ResetSolver(NdLqrSolver* solver);

/**
 * @brief Prints a summary of the solve
 *
 * Prints solve time, the residual norm, and the number of theads.
 *
 * @pre ndlqr_Solve() has already been called
 * @param solver
 */
void ndlqr_PrintSolveSummary(NdLqrSolver* solver);

/**
 * @brief Gets the total number of decision variables for the problem.
 *
 * @param solver
 */
int ndlqr_GetNumVars(NdLqrSolver* solver);

/**
 * @brief Set the number of threads to be used during the solve
 *
 * Does not guarantee that the specified number of threads will be used.
 * To query the actual number of threads used during the solve, use the
 * ndlqr_GetNumThreads() function after the solve.
 *
 * @param solver rsLQR solver
 * @param num_threads requested number of threads
 * @return 0 if successful
 */
int ndlqr_SetNumThreads(NdLqrSolver* solver, int num_threads);

/**
 * @brief Get the number of threads used during the rsLQR solve
 *
 * @param solver A solver which has already been initialized and solved
 * @return number of OpenMP threads used the by solver
 */
int ndlqr_GetNumThreads(NdLqrSolver* solver);

/**
 * @brief Prints a summary of how long individual components took
 *
 * @pre ndlqr_Solve() has already been called
 * @param solver A solver which has already been initialized and solved
 * @return 0 if successful
 */
int ndlqr_PrintSolveProfile(NdLqrSolver* solver);

/**
 * @brief Ge the internal profile data from a solve
 *
 * @param solver A solver which has already been initialized and solved
 * @return A profile object containing timing information about the solve
           See NdLqrProfile for more info.
 */
NdLqrProfile ndlqr_GetProfile(NdLqrSolver* solver);

/**@} */
