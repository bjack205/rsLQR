#pragma once

#include <stdlib.h>

#include "solver.h"

int ndlqr_SolveLeaf(NdLqrSolver* solver, int index);

int ndlqr_SolveLeaves(NdLqrSolver* solver);

int ndlqr_FactorInnerProduct(NdData* data, NdData* fact, int index, int data_level, int fact_level);

int ndlqr_SolveCholeskyFactor(NdData* fact, CholeskyInfo* cholinfo, int index, int level, int upper_level);

bool ndlqr_ShouldCalcLambda(OrderedBinaryTree* tree, int index, int i);

int ndlqr_UpdateShurFactor(NdData* fact, NdData* soln, int index, int i,
                           int level, int upper_level, bool calc_lambda);

int ndlqr_ComputeShurCompliment(NdLqrSolver* solver, int index, int level, int upper_level);
