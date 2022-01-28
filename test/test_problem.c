#include "test/test_problem.h"

#include <stdlib.h>

LQRData* ndlqr_ReadTestLQRData() {
  const char* filename = LQRDATAFILE;
  LQRData* lqrdata = ndlqr_ReadLQRDataJSONFile(filename);
  return lqrdata;
}

LQRProblem* ndlqr_ReadTestLQRProblem() {
  const char* filename = LQRPROBFILE;
  LQRProblem* lqrprob = ndlqr_ReadLQRProblemJSONFile(filename);
  return lqrprob;
}

NdLqrSolver* ndlqr_GenTestSolver() {
  LQRProblem* lqrprob = ndlqr_ReadTestLQRProblem();
  int nstates = lqrprob->lqrdata[0]->nstates;
  int ninputs = lqrprob->lqrdata[0]->ninputs;
  int nhorizon = lqrprob->nhorizon;
  NdLqrSolver* solver = ndlqr_NewNdLqrSolver(nstates, ninputs, nhorizon);
  ndlqr_InitializeWithLQRProblem(lqrprob, solver);
  ndlqr_FreeLQRProblem(lqrprob);
  return solver;
}

LQRProblem* ndlqr_ReadLongTestLQRProblem() {
  const char* filename = LQRPROB256FILE;
  LQRProblem* lqrprob = ndlqr_ReadLQRProblemJSONFile(filename);
  return lqrprob;
}

NdLqrSolver* ndlqr_GenLongTestSolver() {
  LQRProblem* lqrprob = ndlqr_ReadLongTestLQRProblem();
  int nstates = lqrprob->lqrdata[0]->nstates;
  int ninputs = lqrprob->lqrdata[0]->ninputs;
  int nhorizon = lqrprob->nhorizon;
  NdLqrSolver* solver = ndlqr_NewNdLqrSolver(nstates, ninputs, nhorizon);
  ndlqr_InitializeWithLQRProblem(lqrprob, solver);
  ndlqr_FreeLQRProblem(lqrprob);
  return solver;
}
