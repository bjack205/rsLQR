#include "rslqr/nested_dissection.h"

#include <stdio.h>
#include <time.h>

#include "matrix/linalg.h"
#include "matrix/linalg_utils.h"
#include "utils.h"

int ndlqr_SolveLeaf(NdLqrSolver* solver, int index) {
  int nstates = solver->nstates;
  // int ninputs = solver->ninputs;
  int nhorizon = solver->nhorizon;

  NdFactor* C;
  NdFactor* F;
  NdFactor* z;
  Matrix* Q;
  Matrix* R;
  CholeskyInfo* Qchol = NULL;
  CholeskyInfo* Rchol = NULL;

  int k = index;
  if (index == 0) {
    ndlqr_GetNdFactor(solver->data, k, 0, &C);
    ndlqr_GetNdFactor(solver->fact, k, 0, &F);
    ndlqr_GetNdFactor(solver->soln, k, 0, &z);
    Q = &solver->diagonals[2 * k];
    R = &solver->diagonals[2 * k + 1];

    // Solve the block system of equations:
    // [   -I   ] [Fy]   [Cy]   [ 0 ]    [-A'    ]
    // [-I  Q   ] [Fx] = [Cx] = [ A'] => [ 0     ]
    // [      R ] [Fu]   [Cu]   [ B']    [ R \ B']
    MatrixCopy(&F->lambda, &C->state);
    MatrixScaleByConst(&F->lambda, -1.0);
    MatrixSetConst(&F->state, 0.0);
    MatrixCopy(&F->input, &C->input);
    ndlqr_GetRFactorizon(solver->cholfacts, 0, &Rchol);
    MatrixCholeskyFactorizeWithInfo(R, Rchol);
    MatrixCholeskySolveWithInfo(R, &F->input, Rchol);  // Fu = R \ Cu
    MatrixCholeskySolveWithInfo(R, &z->input, Rchol);  // zu = R \ zu

    // Solve the block system of equations (overwriting the rhs vector):
    // [   -I   ] [zy]   [zy]   [ -x0 ]    [ Qx0 + q ]   [-Q zy - zx ]
    // [-I  Q   ] [zx] = [zx] = [ -q  ] => [ x0      ] = [-zy        ]
    // [      R ] [zu]   [zu]   [ -r  ]    [-R \ r   ]   [ R \ zu    ]
    Matrix zy_temp = {nstates, 1,
                      C->lambda.data};  // grab an unused portion of the matrix data
    MatrixCopy(&zy_temp, &z->lambda);
    MatrixCopy(&z->lambda, &z->state);
    MatrixMultiply(Q, &zy_temp, &z->lambda, 0, 0, -1.0,
                   -1.0);  // zy = - Q * zy - zx

    MatrixCopy(&z->state, &zy_temp);
    MatrixScaleByConst(&z->state, -1.0);  // zx = -zy
    ndlqr_GetQFactorizon(solver->cholfacts, 0, &Qchol);
    MatrixCholeskyFactorizeWithInfo(Q, Qchol);

  } else {
    int level = 0;

    Q = &solver->diagonals[2 * k];
    ndlqr_GetQFactorizon(solver->cholfacts, k, &Qchol);
    MatrixCholeskyFactorizeWithInfo(Q, Qchol);

    ndlqr_GetNdFactor(solver->soln, k, 0, &z);

    // All the terms that don't apply at the last time step
    if (k < nhorizon - 1) {
      level = ndlqr_GetIndexLevel(&solver->tree, k);
      ndlqr_GetNdFactor(solver->data, k, level, &C);
      ndlqr_GetNdFactor(solver->fact, k, level, &F);

      R = &solver->diagonals[2 * k + 1];
      ndlqr_GetRFactorizon(solver->cholfacts, k, &Rchol);
      MatrixCholeskyFactorizeWithInfo(R, Rchol);

      MatrixCholeskySolveWithInfo(R, &z->input,
                                  Rchol);  // solve zu = R \ zu  (R \ -r)
      MatrixCopy(&F->state, &C->state);
      MatrixCholeskySolveWithInfo(Q, &F->state,
                                  Qchol);  // solve Fx = Q \ Cx  (Q \ A')
      MatrixCopy(&F->input, &C->input);
      MatrixCholeskySolveWithInfo(R, &F->input,
                                  Rchol);  // solve Fu = Q \ Cu  (R \ B')
    }
    // Only term at the last time step
    MatrixCholeskySolveWithInfo(Q, &z->state,
                                Qchol);  // solve zx = Q \ zx  (Q \ -q)

    // Solve for the terms from the dynamics of the previous time step
    // NOTE: This is -I on the state for explicit integration
    //       For implicit integrators we'd use the A2, B2 partials wrt the next
    //       state and control
    int prev_level = ndlqr_GetIndexLevel(&solver->tree, k - 1);
    ndlqr_GetNdFactor(solver->data, k, prev_level, &C);
    ndlqr_GetNdFactor(solver->fact, k, prev_level, &F);
    MatrixCopy(&F->state, &C->state);  // the -I matrix
    MatrixCholeskySolveWithInfo(Q, &F->state,
                                Qchol);  // solve Q \ -I from previous time step
    MatrixSetConst(&F->input, 0.0);      // Initialize the B2 matrix to zeros
  }
  return 0;
}

int ndlqr_SolveLeaves(NdLqrSolver* solver) {
  for (int k = 0; k < solver->nhorizon; ++k) {
    ndlqr_SolveLeaf(solver, k);
  }
  return 0;
}

int ndlqr_FactorInnerProduct(NdData* data, NdData* fact, int index, int data_level,
                             int fact_level) {
  NdFactor* C1;
  NdFactor* F1;
  NdFactor* C2;
  NdFactor* F2;
  ndlqr_GetNdFactor(data, index, data_level, &C1);
  ndlqr_GetNdFactor(fact, index, fact_level, &F1);
  ndlqr_GetNdFactor(data, index + 1, data_level, &C2);
  ndlqr_GetNdFactor(fact, index + 1, fact_level, &F2);
  Matrix S = ndlqr_GetLambdaFactor(F2);
  MatrixMultiply(&C1->state, &F1->state, &S, true, false, 1.0,
                 -1.0);  // S = C1x'F1x
  MatrixMultiply(&C1->input, &F1->input, &S, true, false, 1.0,
                 1.0);  // S = C1u'F1u + S
  MatrixMultiply(&C2->state, &F2->state, &S, true, false, 1.0,
                 1.0);  // S = C2x'F2x + S
  MatrixMultiply(&C2->input, &F2->input, &S, true, false, 1.0,
                 1.0);  // S = C2u'F2u + S
  return 0;
}

int ndlqr_SolveCholeskyFactor(NdData* fact, CholeskyInfo* cholinfo, int index, int level,
                              int upper_level) {
  if (!fact) return -1;
  if (upper_level <= level) {
    fprintf(stderr, "ERROR: `upper_level` must be greater than `level`.");
  }
  NdFactor* F;
  ndlqr_GetNdFactor(fact, index + 1, level, &F);
  Matrix Sbar = F->lambda;

  NdFactor* G;
  ndlqr_GetNdFactor(fact, index + 1, upper_level, &G);
  Matrix f = G->lambda;

  MatrixCholeskySolveWithInfo(&Sbar, &f, cholinfo);
  return 0;
}

int ndlqr_UpdateShurFactor(NdData* fact, NdData* soln, int index, int i, int level,
                           int upper_level, bool calc_lambda) {
  if (!fact || !soln) return -1;

  NdFactor* f_factor;
  NdFactor* g;
  NdFactor* F;
  ndlqr_GetNdFactor(soln, index + 1, upper_level, &f_factor);
  ndlqr_GetNdFactor(soln, i, upper_level, &g);
  ndlqr_GetNdFactor(fact, i, level, &F);
  Matrix* f = &f_factor->lambda;
  if (calc_lambda) {
    MatrixMultiply(&F->lambda, f, &g->lambda, 0, 0, -1.0, 1.0);
  }
  MatrixMultiply(&F->state, f, &g->state, 0, 0, -1.0, 1.0);
  MatrixMultiply(&F->input, f, &g->input, 0, 0, -1.0, 1.0);
  return 0;
}

bool ndlqr_ShouldCalcLambda(OrderedBinaryTree* tree, int index, int i) {
  BinaryNode* node = tree->node_list + index;
  bool is_start = i == node->left_inds.start || i == node->right_inds.start;
  return !is_start || i == 0;
}

int ndlqr_ComputeShurCompliment(NdLqrSolver* solver, int index, int level,
                                int upper_level) {
  BinaryNode* node = solver->tree.node_list + index;
  int left_start = node->left_inds.start;
  int right_stop = node->right_inds.stop;
  NdData* fact = solver->fact;
  NdData* soln = upper_level == 0 ? solver->soln : solver->fact;

  for (int i = left_start; i <= right_stop; ++i) {
    bool calc_lambda = ndlqr_ShouldCalcLambda(&solver->tree, index, i);
    ndlqr_UpdateShurFactor(fact, soln, index, i, level, upper_level, calc_lambda);
  }
  return 0;
}
