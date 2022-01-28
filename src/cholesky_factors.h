#include "linalg.h"

typedef struct {
  int depth;
  int nhorizon;
  CholeskyInfo* cholinfo;
  int numfacts;
} NdLqrCholeskyFactors;

NdLqrCholeskyFactors* ndlqr_NewCholeskyFactors(int depth, int nhorizon);
int ndlqr_FreeCholeskyFactors();
int ndlqr_GetQFactorizon(NdLqrCholeskyFactors* cholfacts, int index,
                         CholeskyInfo** cholfact);
int ndlqr_GetRFactorizon(NdLqrCholeskyFactors* cholfacts, int index,
                         CholeskyInfo** cholfact);
int ndlqr_GetSFactorization(NdLqrCholeskyFactors* cholfacts, int leaf,
                            int level, CholeskyInfo** cholfact);
