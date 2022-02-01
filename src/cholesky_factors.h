/**
 * @file cholesky_factors.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief A struct for storing all of the info for the Cholesky factorization of the 
 *        rsLQR solver. 
 * @version 0.1
 * @date 2022-01-31
 * 
 * @copyright Copyright (c) 2022
 * 
 * @addtogroup rsLQR 
 * @{
 */
#include "linalg.h"

/**
 * @brief Stores a list of CholeskyInfo structs for the rsLQR solver
 * 
 * ## Construction and deconstruction
 * Initialize a new struct using ndlqr_NewCholeskyFactors(), which must be paired 
 * with a call to ndlqr_FreeCholeskyFactors().
 * 
 * ## Methods
 * - ndlqr_NewCholeskyFactors()
 * - ndlqr_FreeCholeskyFactors()
 * - ndlqr_GetQFactorization()
 * - ndlqr_GetRFactorization()
 * - ndlqr_GetSFactorization()
 */
typedef struct {
  int depth;
  int nhorizon;
  CholeskyInfo* cholinfo;
  int numfacts;
} NdLqrCholeskyFactors;

/**
 * @brief Initialize a new NdLqrCholeskyFactors object
 * 
 * Must be paired with a call to ndlqr_FreeCholeskyFactors().
 * 
 * @param depth    Depth of the binary tree
 * @param nhorizon Length of the time horizon. @p log2(@depth) = @p nhorizon.
 * @return         An initialized NdLqrCholeskyFactors object.
 */
NdLqrCholeskyFactors* ndlqr_NewCholeskyFactors(int depth, int nhorizon);

/**
 * @brief Free the memory of a CholeskyFactors object
 * 
 * @param  cholfacts An initialized NdLqrCholeskyFactors object
 * @post   cholfacts = NULL
 * @return 0 if successful
 */
int ndlqr_FreeCholeskyFactors(NdLqrCholeskyFactors* cholfacts);

/**
 * @brief Get the CholeskyInfo for the matrix Q at index @p index
 * 
 * @param[in]  cholfacts All the stored info for the Cholesky solves
 * @param[in]  index     Knot point index 
 * @param[out] cholfact  Location for the retrieved CholeskyInfo
 * @return               0 if successful
 */
int ndlqr_GetQFactorizon(NdLqrCholeskyFactors* cholfacts, int index,
                         CholeskyInfo** cholfact);


/**
 * @brief Get the CholeskyInfo for the matrix R at index @p index
 * 
 * @param[in]  cholfacts All the stored info for the Cholesky solves
 * @param[in]  index     Knot point index 
 * @param[out] cholfact  Location for the retrieved CholeskyInfo
 * @return               0 if successful
 */
int ndlqr_GetRFactorizon(NdLqrCholeskyFactors* cholfacts, int index,
                         CholeskyInfo** cholfact);

/**
 * @brief Get the CholeskyInfo for the output of ndlqr_SolveCholeskyFactor(). 
 * 
 * Gets the CholeskyInfo for \f$ \bar{B}_i^{(j)} \f$.
 * 
 * @param[in]  cholfacts All the stored info for the Cholesky solves
 * @param[in]  leaf      The leaf index for the desired factor
 * @param[in]  level     The level of the binary tree for the factor
 * @param[out] cholfact  Location for the retrieved CholeskyInfo
 * @return               0 if successful
 */
int ndlqr_GetSFactorization(NdLqrCholeskyFactors* cholfacts, int leaf,
                            int level, CholeskyInfo** cholfact);

/**@} */