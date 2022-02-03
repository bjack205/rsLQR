#include "riccati/riccati_solve.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int ndlqr_SolveRiccati(RiccatiSolver* solver) {
  if (!solver) return -1;
  clock_t t_start_total = clock();

  ndlqr_BackwardPass(solver);
  clock_t t_start_fp = clock();
  ndlqr_ForwardPass(solver);

  // Calculate timing
  clock_t t_stop = clock();
  clock_t diff_bp = t_start_fp - t_start_total;
  clock_t diff_fp = t_stop - t_start_fp;
  clock_t diff_total = t_stop - t_start_total;
  solver->t_solve_ms = diff_total * 1000.0 / (double)CLOCKS_PER_SEC;
  solver->t_backward_pass_ms = diff_bp * 1000.0 / (double)CLOCKS_PER_SEC;
  solver->t_forward_pass_ms = diff_fp * 1000.0 / (double)CLOCKS_PER_SEC;
  return 0;
}

int ndlqr_BackwardPass(RiccatiSolver* solver) {
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

    // Calculate gradient terms
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

    // Calculate Hessian terms
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

    // Calculate Gains
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

    // Calulate Cost-to-Go
    Matrix* P = solver->P + k;
    Matrix* p = solver->p + k;

    MatrixCopy(P, Qxx);
    MatrixMultiply(Quu, K, Qux_tmp, 0, 0, 1.0, 0.0);  // Qux_tmp = Quu * K
    MatrixMultiply(K, Qux_tmp, P, 1, 0, 1.0, 1.0);    // P = Qxx + K'Quu*K
    MatrixMultiply(K, Qux, P, 1, 0, 1.0, 1.0);        // P = Quu + K'Quu*K + K'Qux
    MatrixMultiply(Qux, K, P, 1, 0, 1.0, 1.0);        // P = Quu + K'Quu*K + K'Qux + Qux'K

    MatrixCopy(p, Qx);
    MatrixMultiply(Quu, d, Qu_tmp, 0, 0, 1.0, 0.0);  // Qu_tmp = Quu * d
    MatrixMultiply(K, Qu_tmp, p, 1, 0, 1.0, 1.0);    // p = Qx + K'Quu*d
    MatrixMultiply(K, Qu, p, 1, 0, 1.0, 1.0);        // p = Qx + K'Quu*d + K'Qu
    MatrixMultiply(Qux, d, p, 1, 0, 1.0, 1.0);       // p = Qx + K'Quu*d + K'Qu + Qux'd
  }
  return 0;
}

int ndlqr_ForwardPass(RiccatiSolver* solver) {
  if (!solver) return -1;
  int nhorizon = solver->prob->nhorizon;
  int nstates = solver->prob->lqrdata[0]->nstates;
  Matrix x0 = {nstates, 1, solver->prob->x0};

  MatrixCopy(solver->X + 0, &x0);
  int k;
  for (k = 0; k < nhorizon - 1; ++k) {
    Matrix A = ndlqr_GetA(solver->prob->lqrdata[k]);
    Matrix B = ndlqr_GetB(solver->prob->lqrdata[k]);
    Matrix f = ndlqr_Getd(solver->prob->lqrdata[k]);
    Matrix* Pk = solver->P + k;
    Matrix* pk = solver->p + k;
    Matrix* Kk = solver->K + k;
    Matrix* dk = solver->d + k;
    Matrix* xk = solver->X + k;
    Matrix* uk = solver->U + k;
    Matrix* yk = solver->Y + k;
    Matrix* xn = solver->X + k + 1;

    MatrixCopy(yk, pk);
    MatrixMultiply(Pk, xk, yk, 0, 0, 1.0, 1.0);  // y = P * x + p
    MatrixCopy(uk, dk);
    MatrixMultiply(Kk, xk, uk, 0, 0, 1.0, 1.0);  // un = K * x + d
    MatrixCopy(xn, &f);
    MatrixMultiply(&A, xk, xn, 0, 0, 1.0, 1.0);  // xn = A * x + f
    MatrixMultiply(&B, uk, xn, 0, 0, 1.0, 1.0);  // xn = A * x + B * u + f
  }
  Matrix* Pk = solver->P + k;
  Matrix* pk = solver->p + k;
  Matrix* xk = solver->X + k;
  Matrix* yk = solver->Y + k;
  MatrixCopy(yk, pk);
  MatrixMultiply(Pk, xk, yk, 0, 0, 1.0, 1.0);  // y = P * x + p
  return 0;
}
