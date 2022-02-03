#include "rslqr/cholesky_factors.h"

#include "stdlib.h"
#include "rslqr/utils.h"

NdLqrCholeskyFactors* ndlqr_NewCholeskyFactors(int depth, int nhorizon) {
  if (depth <= 0) return NULL;
  if (nhorizon <= 0) return NULL;
  NdLqrCholeskyFactors* cholfacts =
      (NdLqrCholeskyFactors*)malloc(sizeof(NdLqrCholeskyFactors));
  if (!cholfacts) return NULL;
  int num_leaf_factors = 2 * nhorizon;
  int num_S_factors = 0;
  for (int level = 0; level < depth; ++level) {
    int numleaves = PowerOfTwo(depth - level - 1);
    num_S_factors += numleaves;
  }
  int numfacts = num_leaf_factors + num_S_factors;
  CholeskyInfo* cholinfo = (CholeskyInfo*)malloc(numfacts * sizeof(CholeskyInfo));
  if (!cholinfo) {
    free(cholfacts);
    return NULL;
  }
  for (int i = 0; i < numfacts; ++i) {
    cholinfo[i].fact = NULL;
    cholinfo[i].is_freed = true;
    cholinfo[i].lib = '\0';
    cholinfo[i].success = 1;
    cholinfo[i].uplo = '\0';
  }
  cholfacts->depth = depth;
  cholfacts->nhorizon = nhorizon;
  cholfacts->cholinfo = cholinfo;
  cholfacts->numfacts = numfacts;
  return cholfacts;
}

int ndlqr_FreeCholeskyFactors(NdLqrCholeskyFactors* cholfacts) {
  if (!cholfacts) return -1;
  for (int i = 0; i < cholfacts->numfacts; ++i) {
    FreeFactorization(cholfacts->cholinfo + i);
  }
  free(cholfacts->cholinfo);
  free(cholfacts);
  cholfacts = NULL;
  return 0;
}

int ndlqr_GetQFactorizon(NdLqrCholeskyFactors* cholfacts, int index,
                         CholeskyInfo** cholfact) {
  if (!cholfacts) return -1;
  if (index < 0 || index >= cholfacts->nhorizon) return -1;
  *cholfact = &cholfacts->cholinfo[2 * index];
  return 0;
}

int ndlqr_GetRFactorizon(NdLqrCholeskyFactors* cholfacts, int index,
                         CholeskyInfo** cholfact) {
  if (!cholfacts) return -1;
  if (index < 0 || index >= cholfacts->nhorizon - 1) return -1;
  *cholfact = &cholfacts->cholinfo[2 * index + 1];
  return 0;
}

int ndlqr_GetSFactorization(NdLqrCholeskyFactors* cholfacts, int leaf, int level,
                            CholeskyInfo** cholfact) {
  int numleaves = PowerOfTwo(cholfacts->depth - level - 1);
  int num_leaf_factors = 2 * cholfacts->nhorizon;

  if (!cholfacts) return -1;
  if (level < 0 || level >= cholfacts->depth) return -1;
  if (leaf < 0 || leaf >= numleaves) return -1;

  int leaf_index = 0;
  for (int lvl = 0; lvl < level; ++lvl) {
    int numleaves = PowerOfTwo(cholfacts->depth - lvl - 1);
    leaf_index += numleaves;
  }
  leaf_index += leaf;
  *cholfact = &cholfacts->cholinfo[num_leaf_factors + leaf_index];
  return 0;
}
