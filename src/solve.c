#include "solve.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "binary_tree.h"
#include "linalg.h"
#include "linalg_utils.h"
#include "nested_dissection.h"
#include "omp.h"
#include "utils.h"

#define ENABLE_PROFILER
#ifdef ENABLE_PROFILER
#define OMP_TICK \
  _Pragma("omp single") { t_start = omp_get_wtime(); }

#define OMP_TOC(t_elapsed) \
  _Pragma("omp single nowait") { t_elapsed += (omp_get_wtime() - t_start) * 1000.0; }
#else
#define OMP_TICK
#define OMP_TOC(t_elapsed)
#endif

UnitRange get_work(int total_work, int num_threads, int threadid) {
  int tasks_per_thread = total_work / num_threads;
  int start = tasks_per_thread * threadid;
  int stop = tasks_per_thread * (threadid + 1);
  if (threadid == num_threads - 1) {
    stop = total_work;
  }
  UnitRange rng = {start, stop};
  return rng;
}

int ndlqr_Solve(NdLqrSolver* solver) {
  // clock_t t_start_total = clock();
  double t_start_total = omp_get_wtime();
  double t_start = 0;
  (void)t_start;
  MatrixLinAlgTimeReset();

  int depth = solver->depth;
  int nhorizon = solver->nhorizon;

  omp_set_num_threads(solver->num_threads);

#pragma omp parallel
  {
#pragma omp single
    { solver->num_threads = omp_get_num_threads(); }
    // implicit barrier

    // Solve for independent diagonal blocks
    int num_threads = solver->num_threads;
    int threadid = omp_get_thread_num();
    UnitRange rng = get_work(solver->nhorizon, num_threads, threadid);
    OMP_TICK;
    for (int k = rng.start; k < rng.stop; ++k) {
      ndlqr_SolveLeaf(solver, k);
    }
    OMP_TOC(solver->profile.t_leaves_ms);
#pragma omp barrier

    // Solve factorization
    for (int level = 0; level < depth; ++level) {
      int numleaves = PowerOfTwo(depth - level - 1);

      // Calc Inner Products
      int cur_depth = depth - level;
      int num_products = numleaves * cur_depth;

      rng = get_work(num_products, num_threads, threadid);
      OMP_TICK;
      for (int i = rng.start; i < rng.stop; ++i) {
        int leaf = i / cur_depth;
        int upper_level = level + (i % cur_depth);
        int index = ndlqr_GetIndexFromLeaf(&solver->tree, leaf, level);
        ndlqr_FactorInnerProduct(solver->data, solver->fact, index, level, upper_level);
      }
      OMP_TOC(solver->profile.t_products_ms);
#pragma omp barrier

      // Cholesky factorization
      rng = get_work(numleaves, num_threads, threadid);
      OMP_TICK;
      for (int leaf = rng.start; leaf < rng.stop; ++leaf) {
        int index = ndlqr_GetIndexFromLeaf(&solver->tree, leaf, level);
        // Get the Sbar Matrix calculated above
        NdFactor* F;
        ndlqr_GetNdFactor(solver->fact, index + 1, level, &F);
        Matrix Sbar = F->lambda;
        CholeskyInfo* cholinfo;
        ndlqr_GetSFactorization(solver->cholfacts, leaf, level, &cholinfo);
        MatrixCholeskyFactorizeWithInfo(&Sbar, cholinfo);
      }
      OMP_TOC(solver->profile.t_cholesky_ms);
#pragma omp barrier

      // Solve with Cholesky factor for f
      int upper_levels = cur_depth - 1;
      int num_solves = numleaves * upper_levels;
      rng = get_work(num_solves, num_threads, threadid);
      OMP_TICK;
      for (int i = rng.start; i < rng.stop; ++i) {
        int leaf = i / upper_levels;
        int upper_level = level + 1 + (i % upper_levels);
        int index = ndlqr_GetIndexFromLeaf(&solver->tree, leaf, level);

        CholeskyInfo* cholinfo;
        ndlqr_GetSFactorization(solver->cholfacts, leaf, level, &cholinfo);
        ndlqr_SolveCholeskyFactor(solver->fact, cholinfo, index, level, upper_level);
      }
      OMP_TOC(solver->profile.t_cholsolve_ms);
#pragma omp barrier

      // Shur compliments
      int num_factors = nhorizon * upper_levels;
      rng = get_work(num_factors, num_threads, threadid);
      OMP_TICK;
      for (int i = rng.start; i < rng.stop; ++i) {
        int k = i / upper_levels;
        int upper_level = level + 1 + (i % upper_levels);

        int index = ndlqr_GetIndexAtLevel(&solver->tree, k, level);
        bool calc_lambda = ndlqr_ShouldCalcLambda(&solver->tree, index, k);
        ndlqr_UpdateShurFactor(solver->fact, solver->fact, index, k, level, upper_level,
                               calc_lambda);
      }
      OMP_TOC(solver->profile.t_shur_ms);
#pragma omp barrier
    }

    // Solve for solution vector using the cached factorization
    for (int level = 0; level < depth; ++level) {
      int numleaves = PowerOfTwo(depth - level - 1);

      // Calculate inner products with right-hand-side, with the factors
      // computed above
      rng = get_work(numleaves, num_threads, threadid);
      for (int leaf = rng.start; leaf < rng.stop; ++leaf) {
        int index = ndlqr_GetIndexFromLeaf(&solver->tree, leaf, level);

        // Calculate z = d - F'b1 - F2'b2
        ndlqr_FactorInnerProduct(solver->data, solver->soln, index, level, 0);
      }
#pragma omp barrier

      // Solve for separator variables with cached Cholesky decomposition
      rng = get_work(numleaves, num_threads, threadid);
      for (int leaf = rng.start; leaf < rng.stop; ++leaf) {
        int index = ndlqr_GetIndexFromLeaf(&solver->tree, leaf, level);

        // Get the Sbar Matrix calculated above
        NdFactor* F;
        NdFactor* z;
        ndlqr_GetNdFactor(solver->fact, index + 1, level, &F);
        ndlqr_GetNdFactor(solver->soln, index + 1, 0, &z);
        Matrix Sbar = F->lambda;
        Matrix zy = z->lambda;

        // Solve (S - C1'F1 - C2'F2)^{-1} (d - F1'b1 - F2'b2) -> Sbar \ z = zbar
        //                 |                       |
        //    reuse Cholesky factorization   Inner product calculated above
        CholeskyInfo* cholinfo;
        ndlqr_GetSFactorization(solver->cholfacts, leaf, level, &cholinfo);
        MatrixCholeskySolveWithInfo(&Sbar, &zy, cholinfo);
      }
#pragma omp barrier

      // Propagate information to solution vector
      //    y = y - F zbar
      rng = get_work(nhorizon, num_threads, threadid);
      for (int k = rng.start; k < rng.stop; ++k) {
        int index = ndlqr_GetIndexAtLevel(&solver->tree, k, level);
        bool calc_lambda = ndlqr_ShouldCalcLambda(&solver->tree, index, k);
        ndlqr_UpdateShurFactor(solver->fact, solver->soln, index, k, level, 0, calc_lambda);
      }
#pragma omp barrier
    }
  }
  double diff = omp_get_wtime() - t_start_total;
  solver->solve_time_ms = diff * 1000.0;
  solver->linalg_time_ms = MatrixGetLinAlgTimeMilliseconds();
  solver->profile.t_total_ms = solver->solve_time_ms;
  solver->profile.num_threads = solver->num_threads;
  return 0;
}

Matrix ndlqr_GetSolution(NdLqrSolver* solver) {
  Matrix soln = {solver->nvars, 1, solver->soln->data};
  return soln;
}

int ndlqr_CopySolution(NdLqrSolver* solver, double* soln) {
  if (!solver) return -1;
  memcpy(soln, solver->soln->data, solver->nvars * sizeof(double));
  return solver->nvars;
}
