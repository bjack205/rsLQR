#include <time.h>

#include "test/minunit.h"
#include "ndlqr.h"
#include "solve.h"
#include "nested_dissection.h"
#include "test/test_problem.h"
#include "linalg.h"
#include "utils.h"
#include "matmul.h"

#include "omp.h"

mu_test_init

#define kNumTests 5

void ParallelTiming(int (*func)(NdLqrSolver* solver, int i), NdLqrSolver* solver, int len) {
  printf("This problem has a horizon of %d\n", solver->nhorizon);

  for (int k = 0; k < len; ++k) {
    func(solver, k);
    // ndlqr_SolveLeaf(solver, k);
  }

  double t_start = omp_get_wtime();
  for (int k = 0; k < solver->nhorizon; ++k) {
    // ndlqr_SolveLeaf(solver, k);
    func(solver, k);
  }
  double t_elapsed = (omp_get_wtime() - t_start) * 1000.0; 
  double t_serial_ms = t_elapsed;
  printf("Took %f ms without threading\n", t_elapsed);

  int num_threads[kNumTests] = {1, 2, 4, 8, 16};
  double times_ms[kNumTests];
  for (int i = 0; i < kNumTests; ++i) {
    omp_set_num_threads(num_threads[i]);
    int tasks_per_thread = len / num_threads[i];
    double t_start, t_parallel;
    #pragma omp parallel shared(t_start, t_parallel)
    { 
      int id = omp_get_thread_num();
      int start = tasks_per_thread * id;
      int stop = tasks_per_thread * (id + 1);
      #pragma omp single
      {
        t_start = omp_get_wtime();
      }

      for (int k = start; k < stop; ++k) {
      // for (int k = 0; k < len; ++k) {
        func(solver, k);
        // ndlqr_SolveLeaf(solver, k);
      }

      #pragma omp single
      {
        t_parallel = (omp_get_wtime() - t_start) * 1000.0;
      }
    }
    times_ms[i] = t_parallel;
    double speedup = t_serial_ms / t_parallel;
    printf("Took %f ms (speedup of %.2fx) with %d threads\n", t_parallel, speedup, num_threads[i]);
  }
  (void) times_ms;
}

int SolveLeaves() {
  NdLqrSolver* solver = ndlqr_GenLongTestSolver();
  ParallelTiming(ndlqr_SolveLeaf, solver, solver->nhorizon);
  ndlqr_FreeNdLqrSolver(solver);
  return 1;
}

int InnerProductTask(NdLqrSolver* solver, int i) {
  int level = 0;
  int cur_depth = solver->depth - level;
  int leaf = i / cur_depth;
  int upper_level = level + (i % cur_depth);
  int index = ndlqr_GetIndexFromLeaf(&solver->tree, leaf, level);
  ndlqr_FactorInnerProduct(solver->data, solver->fact, index, level,
                            upper_level);
  return 0;
}

int InnerProducts() {
  NdLqrSolver* solver = ndlqr_GenLongTestSolver();
  int level = 0;
  int numleaves = PowerOfTwo(solver->depth - level - 1);
  int cur_depth = solver->depth - level;
  int num_products = numleaves * cur_depth;
  ParallelTiming(InnerProductTask, solver, num_products);
  ndlqr_FreeNdLqrSolver(solver);
  return 1;
}

int MatMul() {
  int N = 16 * 1000;
  int n = 8;
  eigen_SetNumThreads(1);
  eigen_InitParallel();
  Matrix* A = (Matrix*) malloc(N * sizeof(Matrix));
  Matrix* B = (Matrix*) malloc(N * sizeof(Matrix));
  Matrix* C = (Matrix*) malloc(N * sizeof(Matrix));
  for (int i = 0; i < N; ++i) {
    A[i] = NewMatrix(n, n);
    B[i] = NewMatrix(n, n);
    C[i] = NewMatrix(n, n);
    for (int j = 0; j < n * n; ++j) {
      A[i].data[j] = cos(2 * i);
      B[i].data[j] = 2.1 * i  - 3.2 * i * i;
    }
    MatrixSetConst(&C[i], 0.0);
  }

  double t_start = omp_get_wtime();
  for (int i = 0; i < N; ++i) {
    // MatrixMul(A + i, B + i, C + i);
    MatrixMultiply(A + i, B + i, C + i, 0, 0, 1.0, 0.0);
    // MatMul4x4_unrolled(A[i].data, B[i].data, C[i].data);
    // eigen_MatrixMultiply(n, n, n, A[i].data, B[i].data, C[i].data, 0, 0, 1.0, 0.0);
    // eigen_MatrixMultiply8x8(A[i].data, B[i].data, C[i].data);
    // MatMulSIMD(n, A[i].data, B[i].data, C[i].data);
  }
  double t_serial = (omp_get_wtime() - t_start) * 1000.0;
  printf("Serial time: %g ms\n", t_serial);

  // t_start = omp_get_wtime();
  // for (int i = 0; i < N; ++i) {
  //   MatMul8x8_unrolled(A[i].data, B[i].data, C[i].data);
  // }
  // double t_avx = (omp_get_wtime() - t_start) * 1000.0;
  // printf("Avx mul: %g ms\n", t_avx);


  int num_threads = 4; 
  omp_set_num_threads(num_threads);
  double t_parallel;
  (void) t_parallel;
  int tasks_per_thread = N / num_threads;
  int num_threads_actual = 0;
  printf("Total tasks: %d, tasks per thread = %d\n", N, tasks_per_thread);
  #pragma omp parallel default(none) shared(t_start, t_parallel, A, B, C, num_threads_actual) firstprivate(N, tasks_per_thread, n)
  {
    #pragma omp single
    {
      t_start = omp_get_wtime();
      num_threads_actual = omp_get_num_threads();
    }
    int id = omp_get_thread_num();
    int start = tasks_per_thread * id;
    int stop = tasks_per_thread * (id + 1);
    // printf("Id = %d, range = %d - %d\n", id, start, stop);
    for (int i = start; i < stop; ++i) {
    // #pragma omp for schedule(static, 1000)
    // for (int i = 0; i < N; ++i) {
      // MatrixMul(A + i, B + i, C + i);
      MatrixMultiply(A + i, B + i, C + i, 0, 0, 1.0, 0.0);

      // eigen_MatrixMultiply(n, n, n, A[i].data, B[i].data, C[i].data, 0, 0, 1.0, 0.0);
      // eigen_MatrixMultiply8x8(A[i].data, B[i].data, C[i].data);
      // MatMul4x4_unrolled(A[i].data, B[i].data, C[i].data);
      // MatMulSIMD(n, A[i].data, B[i].data, C[i].data);
    }
    #pragma omp single
    {
      t_parallel = (omp_get_wtime() - t_start) * 1000.0;
    }
  }
  printf("Actual number of threads = %d\n", num_threads_actual);
  printf("Parallel time = %g ms\n", t_parallel);
  printf("Expected / actual speedup: (%dx / %.2fx) = %.2f\n", num_threads, t_serial / t_parallel, t_serial / t_parallel / num_threads);

  for (int i = 0; i < N; ++i) {
    FreeMatrix(&A[i]);
    FreeMatrix(&B[i]);
    FreeMatrix(&C[i]);
  }

  free(A);
  free(B);
  free(C);
  return 1;
}

int SolveComp() {
  LQRProblem* lqrprob = ndlqr_ReadLongTestLQRProblem();
  NdLqrSolver* solver = ndlqr_GenLongTestSolver();
  solver->num_threads = 1;

  // Solve twice
  ndlqr_Solve(solver);
  ndlqr_ResetSolver(solver);
  ndlqr_InitializeWithLQRProblem(lqrprob, solver);
  ndlqr_Solve(solver);
  ndlqr_PrintSolveSummary(solver);
  NdLqrProfile serial;
  ndlqr_CopyProfile(&serial, &solver->profile);

  // Switch to Parallel
  solver->num_threads = 16;
  ndlqr_ResetSolver(solver);
  ndlqr_InitializeWithLQRProblem(lqrprob, solver);
  ndlqr_Solve(solver);
  ndlqr_ResetSolver(solver);
  ndlqr_InitializeWithLQRProblem(lqrprob, solver);
  ndlqr_Solve(solver);
  ndlqr_PrintSolveSummary(solver);
  NdLqrProfile parallel;
  ndlqr_CopyProfile(&parallel, &solver->profile);
  ndlqr_CompareProfile(&serial, &parallel);

  ndlqr_FreeNdLqrSolver(solver);
  ndlqr_FreeLQRProblem(lqrprob);
  return 1;
}

void AllTests() {
  // mu_run_test(SolveLeaves);
  // mu_run_test(InnerProducts);
  // mu_run_test(MatMul);
  mu_run_test(SolveComp);
}

mu_test_main
