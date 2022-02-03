#include "riccati/riccati_solver.h"

#include <stdlib.h>
#include <string.h>

#include "rslqr.h"
#include "riccati/riccati_solve.h"
#include "test/minunit.h"
#include "test/test_problem.h"

mu_test_init

#ifdef FULLTEST
    int kRunFullTest = FULLTEST;
#else
    int kRunFullTest = 0;
#endif

int RiccatiSolverTest() {
  LQRProblem* lqrprob = ndlqr_ReadTestLQRProblem();
  RiccatiSolver* solver = ndlqr_NewRiccatiSolver(lqrprob);
  int nhorizon = lqrprob->nhorizon;
  mu_assert(solver->Qx->rows == 6);
  mu_assert(solver->Qx->cols == 1);
  int dist = solver->Qx->data - solver->X[nhorizon - 1].data;
  mu_assert(dist == 6);
  ndlqr_FreeRiccatiSolver(solver);
  ndlqr_FreeLQRProblem(lqrprob);
  return 1;
}

int RiccatiStepTest() {
  LQRProblem* lqrprob = ndlqr_ReadTestLQRProblem();
  RiccatiSolver* solver = ndlqr_NewRiccatiSolver(lqrprob);

  int nhorizon = solver->prob->nhorizon;
  LQRData** lqrdata = solver->prob->lqrdata;

  int k = nhorizon - 1;
  Matrix Qd = ndlqr_GetQ(lqrdata[k]);
  Matrix q = ndlqr_Getq(lqrdata[k]);
  Matrix* Pn = solver->P + k;
  Matrix* pn = solver->p + k;
  MatrixCopyDiagonal(Pn, &Qd);
  MatrixCopy(pn, &q);

  for (--k; k >= 0; --k) {
    Pn = solver->P + k + 1;
    pn = solver->p + k + 1;

    Matrix A = ndlqr_GetA(lqrdata[k]);
    Matrix B = ndlqr_GetB(lqrdata[k]);
    Matrix f = ndlqr_Getd(lqrdata[k]);
    Matrix Qd = ndlqr_GetQ(lqrdata[k]);
    Matrix q = ndlqr_Getq(lqrdata[k]);
    Matrix Rd = ndlqr_GetR(lqrdata[k]);
    Matrix r = ndlqr_Getr(lqrdata[k]);

    Matrix* Qx = solver->Qx;
    Matrix* Qu = solver->Qu;
    Matrix* Qx_tmp = solver->Qx + 1;
    Matrix* Qu_tmp = solver->Qu + 1;
    MatrixCopy(Qx_tmp, pn);                          // Qx = p
    MatrixMultiply(Pn, &f, Qx_tmp, 0, 0, 1.0, 1.0);  // Qx = P * f + p

    MatrixMultiply(&B, Qx_tmp, Qu, 1, 0, 1.0, 0.0);  // Qu = B' * (P * f + p)
    MatrixMultiply(&A, Qx_tmp, Qx, 1, 0, 1.0, 0.0);  // Qx = A' * (P * f + p)
    MatrixAddition(&r, Qu, 1.0);                     // Qu = r + B' * (P * f + p)
    MatrixAddition(&q, Qx, 1.0);                     // Qx = q + A' * (P * f + p)

    Matrix* Qxx = solver->Qxx;
    Matrix* Qux = solver->Qux;
    Matrix* Quu = solver->Quu;
    Matrix* Qxx_tmp = solver->Qxx + 1;
    Matrix* Qux_tmp = solver->Qux + 1;
    Matrix* Quu_tmp = solver->Quu + 1;

    MatrixCopyDiagonal(Qxx, &Qd);
    MatrixCopyDiagonal(Quu, &Rd);

    MatrixMultiply(&A, Pn, Qxx_tmp, 1, 0, 1.0, 0.0);   // Qxx = A'P
    MatrixMultiply(&B, Pn, Qux_tmp, 1, 0, 1.0, 0.0);   // Qux = B'P
    MatrixMultiply(Qxx_tmp, &A, Qxx, 0, 0, 1.0, 1.0);  // Qxx = Q + A'P*A
    MatrixMultiply(Qux_tmp, &B, Quu, 0, 0, 1.0, 1.0);  // Quu = R + B'P*B
    MatrixMultiply(Qux_tmp, &A, Qux, 0, 0, 1.0, 0.0);  // Qux = B'P*A

    double Qx_data[6] = {-69.0, 0.5999999999999996, 70.2, 134.3, 210.3, 286.3};
    Matrix Qx_ans = {6, 1, Qx_data};
    double Qu_data[3] = {6.425000000000001, 20.145, 33.865};
    Matrix Qu_ans = {3, 1, Qu_data};
    double Qxx_data[36] = {11.0, 0.0,  0.0, 1.0,  0.0, 0.0,  0.0, 11.0, 0.0,
                           0.0,  1.0,  0.0, 0.0,  0.0, 11.0, 0.0, 0.0,  1.0,
                           1.0,  0.0,  0.0, 11.1, 0.0, 0.0,  0.0, 1.0,  0.0,
                           0.0,  11.1, 0.0, 0.0,  0.0, 1.0,  0.0, 0.0,  11.1};
    Matrix Qxx_ans = {6, 6, Qxx_data};
    double Quu_data[9] = {0.11025, 0.0, 0.0, 0.0, 0.11025, 0.0, 0.0, 0.0, 0.11025};
    Matrix Quu_ans = {3, 3, Quu_data};
    double Qux_data[18] = {0.05000000000000001,
                           0.0,
                           0.0,
                           0.0,
                           0.05000000000000001,
                           0.0,
                           0.0,
                           0.0,
                           0.05000000000000001,
                           1.005,
                           0.0,
                           0.0,
                           0.0,
                           1.005,
                           0.0,
                           0.0,
                           0.0,
                           1.005};
    Matrix Qux_ans = {3, 6, Qux_data};

    mu_assert(MatrixNormedDifference(&Qx_ans, Qx) < 1e-6);
    mu_assert(MatrixNormedDifference(&Qu_ans, Qu) < 1e-6);
    mu_assert(MatrixNormedDifference(&Qxx_ans, Qxx) < 1e-6);
    mu_assert(MatrixNormedDifference(&Quu_ans, Quu) < 1e-6);
    mu_assert(MatrixNormedDifference(&Qux_ans, Qux) < 1e-6);

    Matrix* K = solver->K + k;
    Matrix* d = solver->d + k;
    MatrixCopy(Quu_tmp, Quu);
    MatrixCopy(K, Qux);
    MatrixCopy(d, Qu);

    CholeskyInfo cholinfo = {'L', 0, 'E', NULL, 0};
    MatrixCholeskyFactorizeWithInfo(Quu_tmp, &cholinfo);
    MatrixCholeskySolveWithInfo(Quu_tmp, K, &cholinfo);
    MatrixCholeskySolveWithInfo(Quu_tmp, d, &cholinfo);
    MatrixScaleByConst(K, -1);
    MatrixScaleByConst(d, -1);
    FreeFactorization(&cholinfo);

    double K_data[18] = {-0.4535147392290251,
                         -0.0,
                         -0.0,
                         -0.0,
                         -0.4535147392290251,
                         -0.0,
                         -0.0,
                         -0.0,
                         -0.4535147392290251,
                         -9.1156462585034,
                         -0.0,
                         -0.0,
                         -0.0,
                         -9.1156462585034,
                         -0.0,
                         -0.0,
                         -0.0,
                         -9.1156462585034};
    Matrix K_ans = {3, 6, K_data};
    double d_data[3] = {-58.27664399092971, -182.72108843537413, -307.1655328798186};
    Matrix d_ans = {3, 1, d_data};

    mu_assert(MatrixNormedDifference(&K_ans, K) < 1e-6);
    mu_assert(MatrixNormedDifference(&d_ans, d) < 1e-6);

    Matrix* P = solver->P + k;
    Matrix* p = solver->p + k;

    MatrixCopy(P, Qxx);
    MatrixMultiply(Quu, K, Qux_tmp, 0, 0, 1.0, 0.0);  // Qux_tmp = Quu * K
    MatrixMultiply(K, Qux_tmp, P, 1, 0, 1.0, 1.0);    // P = Qxx + K'Quu*K
    MatrixMultiply(K, Qux, P, 1, 0, 1.0, 1.0);        // P = Quu + K'Quu*K + K'Qux
    MatrixMultiply(Qux, K, P, 1, 0, 1.0, 1.0);        // P = Quu + K'Quu*K + K'Qux + Qux'K

    MatrixCopy(p, Qx);
    MatrixMultiply(Quu, d, Qu_tmp, 0, 0, 1.0, 1.0);  // Qu_tmp = Quu * d
    MatrixMultiply(K, Qu_tmp, p, 1, 0, 1.0, 1.0);    // p = Qx + K'Quu*d
    MatrixMultiply(K, Qu, p, 1, 0, 1.0, 1.0);        // p = Qx + K'Quu*d + K'Qu
    MatrixMultiply(Qux, d, p, 1, 0, 1.0, 1.0);       // p = Qx + K'Quu*d + K'Qu + Qux'd

    double P_data[36] = {10.977324263038549,
                         0.0,
                         0.0,
                         0.5442176870748299,
                         0.0,
                         0.0,
                         0.0,
                         10.977324263038549,
                         0.0,
                         0.0,
                         0.5442176870748299,
                         0.0,
                         0.0,
                         0.0,
                         10.977324263038549,
                         0.0,
                         0.0,
                         0.5442176870748299,
                         0.5442176870748299,
                         0.0,
                         0.0,
                         1.9387755102040813,
                         0.0,
                         0.0,
                         0.0,
                         0.5442176870748299,
                         0.0,
                         0.0,
                         1.9387755102040813,
                         0.0,
                         0.0,
                         0.0,
                         0.5442176870748299,
                         0.0,
                         0.0,
                         1.9387755102040813};
    Matrix P_ans = {6, 6, P_data};
    double p_data[6] = {-71.91383219954649, -8.536054421768709, 54.84172335600907,
                        75.73197278911566,  26.66530612244904,  -22.401360544217596};
    Matrix p_ans = {6, 1, p_data};

    mu_assert(MatrixNormedDifference(&P_ans, P) < 1e-6);
    mu_assert(MatrixNormedDifference(&p_ans, p) < 1e-6);

    break;
  }

  ndlqr_FreeRiccatiSolver(solver);
  ndlqr_FreeLQRProblem(lqrprob);
  return 1;
}

int BackwardPassTest() {
  LQRProblem* lqrprob = ndlqr_ReadTestLQRProblem();
  RiccatiSolver* solver = ndlqr_NewRiccatiSolver(lqrprob);
  ndlqr_BackwardPass(solver);

  double P_data[36] = {11.774910288989293,
                       0.0,
                       0.0,
                       1.139328540729876,
                       0.0,
                       0.0,
                       0.0,
                       11.774910288989293,
                       0.0,
                       0.0,
                       1.139328540729876,
                       0.0,
                       0.0,
                       0.0,
                       11.774910288989293,
                       0.0,
                       0.0,
                       1.139328540729876,
                       1.139328540729876,
                       0.0,
                       0.0,
                       1.7402346445435521,
                       0.0,
                       0.0,
                       0.0,
                       1.139328540729876,
                       0.0,
                       0.0,
                       1.7402346445435521,
                       0.0,
                       0.0,
                       0.0,
                       1.139328540729876,
                       0.0,
                       0.0,
                       1.7402346445435521};
  Matrix P_ans = {6, 6, P_data};
  double p_data[6] = {109.00822409796677, 181.20262227329562, 253.3970204486244,
                      32.229649977292816, 26.00963298587046,  19.78961599444808};
  Matrix p_ans = {6, 1, p_data};
  double K_data[18] = {-6.005830262804116,
                       -0.0,
                       -0.0,
                       -0.0,
                       -6.005830262804116,
                       -0.0,
                       -0.0,
                       -0.0,
                       -6.005830262804116,
                       -6.832682175070581,
                       -0.0,
                       -0.0,
                       -0.0,
                       -6.832682175070581,
                       -0.0,
                       -0.0,
                       -0.0,
                       -6.832682175070581};
  Matrix K_ans = {3, 6, K_data};
  double d_data[3] = {-162.79238772394484, -156.8950187220568, -150.99764972016862};
  Matrix d_ans = {3, 1, d_data};

  mu_assert(MatrixNormedDifference(&P_ans, solver->P) < 1e-6);
  mu_assert(MatrixNormedDifference(&p_ans, solver->p) < 1e-6);
  mu_assert(MatrixNormedDifference(&K_ans, solver->K) < 1e-6);
  mu_assert(MatrixNormedDifference(&d_ans, solver->d) < 1e-6);

  ndlqr_FreeRiccatiSolver(solver);
  ndlqr_FreeLQRProblem(lqrprob);
  return 1;
}

int ForwardPassTest() {
  LQRProblem* lqrprob = ndlqr_ReadTestLQRProblem();
  RiccatiSolver* solver = ndlqr_NewRiccatiSolver(lqrprob);
  ndlqr_BackwardPass(solver);
  ndlqr_ForwardPass(solver);

  double Y_data[6] = {125.4177529215159,  166.85459270424298, 230.98171631063963,
                      56.151722488333256, -7.809331509236316, -74.33992381590357};
  double X_data[6] = {28.54177529215159,  26.285459270424298,  26.298171631063965,
                      2.4151722488333256, -10.380933150923632, -23.433992381590357};
  double U_data[3] = {75.7738986559091, -5.333981259758474, -72.09161999628498};

  Matrix Y_ans = {6, 1, Y_data};
  Matrix X_ans = {6, 1, X_data};
  Matrix U_ans = {3, 1, U_data};

  mu_assert(MatrixNormedDifference(&Y_ans, solver->Y + 7));
  mu_assert(MatrixNormedDifference(&X_ans, solver->X + 7));
  mu_assert(MatrixNormedDifference(&U_ans, solver->U + 6));

  ndlqr_FreeRiccatiSolver(solver);
  ndlqr_FreeLQRProblem(lqrprob);
  return 1;
}

int RiccatiSolveTest() {
  LQRProblem* lqrprob = ndlqr_ReadTestLQRProblem();
  RiccatiSolver* solver = ndlqr_NewRiccatiSolver(lqrprob);

  ndlqr_SolveRiccati(solver);
  ndlqr_PrintRiccatiSummary(solver);

  Matrix x_ans = ReadMatrixJSONFile(SAMPLEPROBFILE, "soln");
  Matrix x = ndlqr_GetRiccatiSolution(solver);
  double err = MatrixNormedDifference(&x, &x_ans);
  printf("  Final error: %g\n", err);
  mu_assert(err < 1e-10);
  FreeMatrix(&x_ans);

  ndlqr_FreeRiccatiSolver(solver);
  ndlqr_FreeLQRProblem(lqrprob);
  return 1;
}

int RiccatiSolveTwiceTest() {
  LQRProblem* lqrprob = ndlqr_ReadTestLQRProblem();
  RiccatiSolver* solver = ndlqr_NewRiccatiSolver(lqrprob);

  ndlqr_SolveRiccati(solver);
  ndlqr_SolveRiccati(solver);

  Matrix x_ans = ReadMatrixJSONFile(SAMPLEPROBFILE, "soln");
  double* soln = (double*)malloc(solver->nvars * sizeof(double));
  ndlqr_CopyRiccatiSolution(solver, soln);
  Matrix x = {solver->nvars, 1, soln};
  double err = MatrixNormedDifference(&x, &x_ans);
  printf("Final error after 2nd solve: %g\n", err);
  mu_assert(err < 1e-10);
  FreeMatrix(&x_ans);

  ndlqr_FreeRiccatiSolver(solver);
  ndlqr_FreeLQRProblem(lqrprob);
  return 1;
}

int SolveLongProblem() {
  LQRProblem* lqrprob = ndlqr_ReadLongTestLQRProblem();
  RiccatiSolver* riccati = ndlqr_NewRiccatiSolver(lqrprob);
  NdLqrSolver* solver = ndlqr_GenLongTestSolver();

  ndlqr_SolveRiccati(riccati);
  ndlqr_Solve(solver);

  Matrix x_soln = ReadMatrixJSONFile(LQRPROB256FILE, "soln");

  Matrix x_ndlqr = ndlqr_GetSolution(solver);
  Matrix x_ric = ndlqr_GetRiccatiSolution(riccati);
  double err_ndlqr = MatrixNormedDifference(&x_ndlqr, &x_soln);
  double err_ric = MatrixNormedDifference(&x_ric, &x_soln);
  double diff = MatrixNormedDifference(&x_ric, &x_ndlqr);
  printf("NDLQR Error:   %g\n", err_ndlqr);
  printf("Riccati Error: %g\n", err_ric);
  printf("Difference:    %g\n", diff);

  Matrix ndlqr_start = {10, 1, x_ndlqr.data};
  Matrix ric_start = {10, 1, x_ric.data};
  Matrix soln_start = {10, 1, x_soln.data};
  PrintRowVector(&ndlqr_start);
  PrintRowVector(&ric_start);
  PrintRowVector(&soln_start);

  FreeMatrix(&x_soln);
  ndlqr_FreeRiccatiSolver(riccati);
  ndlqr_FreeNdLqrSolver(solver);
  ndlqr_FreeLQRProblem(lqrprob);
  return 1;
}

void AllTests() {
  mu_run_test(RiccatiSolverTest);
  mu_run_test(RiccatiStepTest);
  mu_run_test(BackwardPassTest);
  mu_run_test(ForwardPassTest);
  mu_run_test(RiccatiSolveTest);
  mu_run_test(RiccatiSolveTwiceTest);
  if (kRunFullTest) {
    mu_run_test(SolveLongProblem);
  }
}

int main(int argc, char* argv[]) {
  if (argc > 1) {
    kRunFullTest = strcmp(argv[1], "full") == 0;
  }
  ResetTests();
  AllTests();
  PrintTestResult();
  return 0;
}
