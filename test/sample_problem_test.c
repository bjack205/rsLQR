#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ndlqr.h"
#include "riccati_solve.h"
#include "test/test_problem.h"

#ifdef FULLTEST
int kRunFullTest = FULLTEST;
#if FULLTEST
#define kNruns 100
#else
#define kNruns 10
#endif
#else
int kRunFullTest = 0;
#define kNruns 10
#endif

int compare_doubles(const void* a, const void* b) {
  int arg1 = *(const double*)a;
  int arg2 = *(const double*)b;
  if (arg1 < arg2) return -1;
  if (arg2 < arg1) return 1;
  return 0;
}

typedef struct {
  double mean;
  double std;
  double min;
  double median;
  double max;
} TimeStats;

void PrintStats(TimeStats stats) {
  printf("  Average = %.2g +- %.2g ms\n", stats.mean, stats.std);
  printf("  Mininum = %.2g\n", stats.min);
  printf("  Median  = %.2g\n", stats.median);
  printf("  Maximum = %.2g\n", stats.max);
}

TimeStats CalcStats(double* times, int len) {
  double time_total = 0.0;
  for (int i = 0; i < len; ++i) {
    time_total += times[i];
  }
  double time_mean = time_total / (double)len;
  double time_std = 0;
  double time_max = 0;
  double time_min = INFINITY;
  for (int i = 0; i < len; ++i) {
    double diff = times[i] - time_mean;
    time_std += diff * diff;
    time_max = fmax(time_max, times[i]);
    time_min = fmin(time_min, times[i]);
  }
  time_std = sqrt(time_std / (double)len);
  qsort(times, len, sizeof(double), compare_doubles);
  double time_median = times[len / 2];
  TimeStats stats = {time_mean, time_std, time_min, time_median, time_max};
  return stats;
}

int main(int argc, char* argv[]) {
  int len = 8;
  int valid_lengths[2] = {8, 256};
  int num_lengths = 2;
  if (argc == 2) {
    int tmp = len;
    if (sscanf(argv[1], "%d", &tmp) == 1) {
      bool is_input_valid = false;
      for (int i = 0; i < num_lengths; ++i) {
        if (tmp == valid_lengths[i]) {
          len = tmp;
          is_input_valid = true;
          break;
        }
      }
      if (!is_input_valid) {
        printf("Argument Error: %d is not a valid trajectory length. Please choose from (",
               tmp);
        int i;
        for (i = 0; i < num_lengths - 1; ++i) {
          printf("%d, ", valid_lengths[i]);
        }
        printf("%d)\n", valid_lengths[i]);
        return 0;
      }
    }
  }
  printf("Using a trajectory length of %d\n", len);

  LQRProblem* lqrprob = NULL;
  Matrix x_ans;
  switch (len) {
    case 8:
      lqrprob = ndlqr_ReadTestLQRProblem();
      x_ans = ReadMatrixJSONFile(LQRPROBFILE, "soln");
      break;

    case 256:
      lqrprob = ndlqr_ReadLongTestLQRProblem();
      x_ans = ReadMatrixJSONFile(LQRPROB256FILE, "soln");
      break;

    default:
      fprintf(stderr, "%d is not a known trajectory length.\n", len);
      break;
  }

  int nstates = lqrprob->lqrdata[0]->nstates;
  int ninputs = lqrprob->lqrdata[0]->ninputs;
  int nhorizon = lqrprob->nhorizon;

  NdLqrSolver* solver = ndlqr_NewNdLqrSolver(nstates, ninputs, nhorizon);
  ndlqr_InitializeWithLQRProblem(lqrprob, solver);
  solver->num_threads = 1;
  ndlqr_Solve(solver);

  RiccatiSolver* riccati = ndlqr_NewRiccatiSolver(lqrprob);

  double times[kNruns];
  double times_ric[kNruns];
  (void)times;

  Matrix x = ndlqr_GetSolution(solver);

  bool right_answer = true;
  bool same_answer = true;
  int nruns = kNruns;
  if (!kRunFullTest) {
    nruns = 1;
  }

  for (int i = 0; i < nruns; ++i) {
    ndlqr_InitializeWithLQRProblem(lqrprob, solver);

    clock_t start = clock();
    ndlqr_ResetNdData(solver->fact);
    ndlqr_Solve(solver);
    clock_t diff = clock() - start;

    double msec = diff * 1000.0 / (double)CLOCKS_PER_SEC;
    times[i] = msec;
    double err = MatrixNormedDifference(&x, &x_ans);
    right_answer &= err < 1e-6;

    ndlqr_SolveRiccati(riccati);
    times_ric[i] = riccati->t_solve_ms;
    Matrix x_ric = ndlqr_GetRiccatiSolution(riccati);
    err = MatrixNormedDifference(&x_ric, &x);
    same_answer &= err < 1e-6;
  }
  ndlqr_PrintSolveSummary(solver);
  double ric_time = 0.0;
  double ric_bp = 0.0;
  double ric_fp = 0.0;
  ndlqr_GetRiccatiSolveTimes(riccati, &ric_time, &ric_bp, &ric_fp);

  // Get time statistics
  TimeStats stats = CalcStats(times, kNruns);
  TimeStats stats_ric = CalcStats(times_ric, kNruns);

  // Print the summary
  PrintStats(stats);
  printf("Got the right answer? %d\n", right_answer);
  printf("Using Eigen? %d\n", kUseEigen);

  printf("\nRiccati Solver:\n");
  ndlqr_PrintRiccatiSummary(riccati);
  PrintStats(stats_ric);
  printf("Got the same answer? %d\n", same_answer);

  FreeMatrix(&x_ans);
  ndlqr_FreeLQRProblem(lqrprob);
  ndlqr_FreeNdLqrSolver(solver);
  ndlqr_FreeRiccatiSolver(riccati);
  return 0;
}
