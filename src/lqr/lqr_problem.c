#include "lqr/lqr_problem.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int ndlqr_InitializeLQRProblem(LQRProblem* lqrproblem, double* x0, LQRData** lqrdata) {
  if (!lqrproblem) return -1;
  for (int k = 0; k < lqrproblem->nhorizon; ++k) {
    ndlqr_CopyLQRData(lqrproblem->lqrdata[k], lqrdata[k]);
  }
  memcpy(lqrproblem->x0, x0, lqrproblem->lqrdata[0]->nstates * sizeof(double));
  return 0;
}

LQRProblem* ndlqr_NewLQRProblem(int nstates, int ninputs, int nhorizon) {
  if (nhorizon <= 0) {
    fprintf(stderr, "ERROR: Horizon must be positive.\n");
    return NULL;
  }

  LQRData** lqrdata = (LQRData**)malloc(nhorizon * sizeof(LQRData*));
  if (lqrdata == NULL) {
    fprintf(stderr, "ERROR: Couldn't allocate memory for LQRProblem.\n");
    return NULL;
  }
  for (int k = 0; k < nhorizon; ++k) {
    lqrdata[k] = ndlqr_NewLQRData(nstates, ninputs);
  }
  double* x0 = (double*)malloc(nstates * sizeof(double));
  LQRProblem* lqrproblem = (LQRProblem*)malloc(sizeof(LQRProblem));
  lqrproblem->nhorizon = nhorizon;
  lqrproblem->x0 = x0;
  lqrproblem->lqrdata = lqrdata;
  return lqrproblem;
}

int ndlqr_FreeLQRProblem(LQRProblem* lqrprob) {
  if (!lqrprob) return -1;
  if (lqrprob->lqrdata) {
    for (int k = 0; k < lqrprob->nhorizon; ++k) {
      ndlqr_FreeLQRData(lqrprob->lqrdata[k]);
    }
  }
  if (lqrprob->x0) {
    free(lqrprob->x0);
  }
  free(lqrprob->lqrdata);
  free(lqrprob);
  lqrprob = NULL;
  return 0;
}
