/**
 * @file nested_dissection.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief Core methods using during the rsLQR solve
 * @version 0.1
 * @date 2022-01-31
 * 
 * @copyright Copyright (c) 2022
 * 
 * @addtogroup rsLQR 
 * @{
 */
#pragma once

#include <stdlib.h>

#include "solver.h"

/**
 * @brief Solve all the equations for the lowest-level diagonal blocks, by timestep
 *
 * Called at the beginning of the solve, this method calculates the initial pieces 
 * of the factorization needed to start the recursion up through the rest of the levels.
 * For and LQR problem, it calculates
 *
 * \f[
 * Q_k^{-1} A_k^T
 * Q_k^{-1} q_k
 * R_k^{-1} B_k^T
 * R_k^{-1} r_k
 * \f]
 * 
 * for each index \f$ k \f$, with special-casing applied to the first and last indices.
 * 
 * @param solver An initialized rsLQR solver
 * @param index Knotpoint index
 * @return 0 if successful
 */
int ndlqr_SolveLeaf(NdLqrSolver* solver, int index);

int ndlqr_SolveLeaves(NdLqrSolver* solver);

/**
 * @brief Calculates one of the inner products needed at level @p data_level
 * 
 * Calculates the following:
 * 
 * \f[
 * \bar{\Lambda_{k+1}^{(p)}} = \bar{\Lambda_k^{(p)}}S + 
 *   \langle X_k^{(j)}, \bar{X}_k^{(p)} \rangle + 
 *   \langle U_k^{(j)}, \bar{U}_k^{(p)} \rangle +
 *   \langle X_{k+1}^{(j)}, \bar{X}_{k+1}^{(p)} \rangle +
 *   \langle U_{k+1}^{(j)}, \bar{U}_{k+1}^{(p)} \rangle
 * \f]
 * 
 * where \f$ j \f$ is @p data_level, \f$ p \f$ is @p fact_level, and \f$ k \f$ is
 * @p index.
 * 
 * @param data       The data for the original KKT matrix
 * @param fact       The current data for the factorization
 * @param index      Knot point index
 * @param data_level Level index for the original data, or equivalently the current level
 *                   being processed by the outer solve
 * @param fact_level Level index for the factorization data, or equivalently the parent or 
 *                   upper level. @p fact_level >= @p data_level.
 * @return 0 if successful
 */
int ndlqr_FactorInnerProduct(NdData* data, NdData* fact, int index, int data_level, int fact_level);

/**
 * @brief Use the precomputed Cholesky factorization to solve for y at each parent level
 * 
 * Solve the following linear system of equations, overwriting the right-hand-side.
 * 
 * \f[
 * \bar{\Lambda}_{k+1}^{(j)} y = \bar{\Lambda}_{k+1}^{(p)}
 * \f]
 * 
 * where \f$ j \f$ is @p level, \f$ p \f$ is @p upper_level, and \f$ k \f$ is
 * @p index.
 * 
 * This is the same as solving 
 * 
 * \f[
 * -\bar{B}_i^{(j)} y_i^{(j,p)} = \bar{b}_i^{(j,p)}
 * \f]
 * 
 * using the notation from the paper.
 * 
 * @param fact        Data for the factorization
 * @param cholinfo    Cached Cholesky factorization for \f$ \bar{\Lambda}_{k+1}^{(j)} \f$.
 * @param index       Knot point index. Should be calculated using ndlqr_GetIndexFromLeaf().
 * @param level       Level index for the level currently being processed by the 
 *                    upper-level solve.
 * @param upper_level Level index for the right-hand-side. @p upper_level > @p level.
 * @return 0 if successful
 */
int ndlqr_SolveCholeskyFactor(NdData* fact, CholeskyInfo* cholinfo, int index, int level, int upper_level);

/**
 * @brief Determines if the \f$ \Lambda \f$ should be updated during 
 *        ndlqr_UpdateSchurFactor()
 * 
 * Basically checks to see if the \f$ \Lambda \f$ for the given index is a $B_i^{(p)}$ for
 * any level greater than or equal to the current level.
 * 
 * @param tree  Binary tree with cached information about the structure of the problem
 * @param index Knot point index calculated using ndlqr_GetIndexAtLevel()
 * @param i     Knot point index being processed
 * @return 0 if successful
 */
bool ndlqr_ShouldCalcLambda(OrderedBinaryTree* tree, int index, int i);

/**
 * @brief Calculates \f$ x \f$ and \f$ z \f$ to complete the factorization at the current 
 *        level
 * 
 * Calculates the following
 * 
 * \f[
 * \bar{\Lambda}_k^{(p)} =  \bar{\Lambda}_k^{(p)} - \bar{\Lambda}_k^{(j)} \bar{\Lambda}_{k_\text{sep}}^{(p)}
 * \bar{X}_k^{(p)} =  \bar{X}_k^{(p)} - \bar{X}_k^{(j)} \bar{\Lambda}_{k_\text{sep}}^{(p)}
 * \bar{U}_k^{(p)} =  \bar{U}_k^{(p)} - \bar{U}_k^{(j)} \bar{\Lambda}_{k_\text{sep}}^{(p)}
 * \f]
 * 
 * where \f$ \bar{\Lambda}_{k_\text{sep}}^{(p)} \f$ is equivalent to $y_i^{(j,p)}$ from the 
 * paper and is the result of the ndlqr_SolveCholeskyFactor() for @p index, @p level, and 
 * @p upper_level.
 * 
 * @param fact        NdData for the factorization data
 * @param soln        NdData for the factorization data or solution vector
 * @param index       Knot point index of the "separator" at level @p level.
 *                    Should be calculated using ndlqr_GetIndexAtlevel().
 * @param i           Knot point index to be processed
 * @param level       Level index for the level currently being processed
 * @param upper_level Level index of the upper level. 
 *                    @p upper_level > @p level
 * @param calc_lambda Whether the \f$ \Lambda \f$ factor should be updated.
 *                    This should be calculated using ndlqr_ShouldCalcLambda().
 * @return
 */
int ndlqr_UpdateShurFactor(NdData* fact, NdData* soln, int index, int i,
                           int level, int upper_level, bool calc_lambda);

int ndlqr_ComputeShurCompliment(NdLqrSolver* solver, int index, int level, int upper_level);

/**@} */
