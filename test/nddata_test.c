#include "nddata.h"

#include <stdlib.h>
#include <string.h>

#include "test/minunit.h"

mu_test_init

    int
    NewNdDataTest() {
  int nstates = 6;
  int ninputs = 3;
  int nsegments = 7;
  NdData* nddata_bad = ndlqr_NewNdData(nstates, ninputs, nsegments, nstates);
  mu_assert(nddata_bad == NULL);
  nddata_bad = ndlqr_NewNdData(nstates * 0, ninputs, nsegments + 1, nstates);
  mu_assert(nddata_bad == NULL);
  nddata_bad = ndlqr_NewNdData(nstates, ninputs * 0, nsegments + 1, nstates);
  mu_assert(nddata_bad == NULL);
  nddata_bad = ndlqr_NewNdData(nstates, ninputs, nsegments * 1, nstates);
  mu_assert(nddata_bad == NULL);

  NdData* nddata = ndlqr_NewNdData(nstates, ninputs, nsegments + 1, nstates);
  mu_assert(nddata->nstates == nstates);
  mu_assert(nddata->ninputs == ninputs);
  mu_assert(nddata->nsegments == nsegments);
  mu_assert(nddata->depth == 3);
  int res = ndlqr_FreeNdData(nddata);
  mu_assert(res == 0);
  return 1;
}

void SetNdDataBlock(NdData* nddata, int off, double valy, double valx, double valu) {
  int nstates = nddata->nstates;
  int ninputs = nddata->ninputs;
  int width = nddata->width;
  for (int i = 0; i < nstates * width; ++i) {
    nddata->data[i + off] = valy;
    nddata->data[i + off + nstates * width] = valx;
  }
  for (int i = 0; i < width * ninputs; ++i) {
    nddata->data[i + off + 2 * (nstates * width)] = valu;
  }
}

int CheckFactors(NdFactor* factor, double valy, double valx, double valu) {
  Matrix Cy = ndlqr_GetLambdaFactor(factor);
  Matrix Cx = ndlqr_GetStateFactor(factor);
  Matrix Cu = ndlqr_GetInputFactor(factor);
  for (int i = 0; i < MatrixNumElements(&Cy); ++i) {
    mu_assert(Cy.data[i] == valy);
    mu_assert(Cx.data[i] == valx);
  }
  for (int i = 0; i < MatrixNumElements(&Cu); ++i) {
    mu_assert(Cu.data[i] == valu);
  }
  return 1;
}

int SetFactors() {
  int nstates = 6;
  int ninputs = 3;
  int nsegments = 7;
  int nhorizon = nsegments + 1;
  NdData* nddata = ndlqr_NewNdData(nstates, ninputs, nsegments + 1, nstates);
  int depth = nddata->depth;
  int factorsize = (2 * nstates + ninputs) * nstates;
  int levelsize = factorsize * nhorizon;  // number of doubles in a level
  int totalsize = levelsize * depth;      // number of doubles in NdData
  NdFactor* factor;

  // Check First Factor
  SetNdDataBlock(nddata, 0, 1.0, 2.1, 3.3);
  ndlqr_GetNdFactor(nddata, 0, 0, &factor);
  mu_assert(CheckFactors(factor, 1.0, 2.1, 3.3) == 1);

  // Check Second Factor
  SetNdDataBlock(nddata, factorsize, 1.1, 2.2, 3.4);
  ndlqr_GetNdFactor(nddata, 1, 0, &factor);
  mu_assert(nddata->data + factorsize == factor->lambda.data);
  mu_assert(CheckFactors(factor, 1.1, 2.2, 3.4) == 1);

  // Check Second Factor on next level
  SetNdDataBlock(nddata, factorsize + levelsize, 1.2, 2.3, 3.5);
  ndlqr_GetNdFactor(nddata, 1, 1, &factor);
  mu_assert(CheckFactors(factor, 1.2, 2.3, 3.5) == 1);

  // Write to the last element in the memory block
  nddata->data[totalsize - 1] = 101.23;
  ndlqr_GetNdFactor(nddata, nhorizon - 1, depth - 1, &factor);
  Matrix Cu = ndlqr_GetInputFactor(factor);
  double lastelement = *MatrixGetElement(&Cu, ninputs - 1, nstates - 1);
  mu_assert(lastelement == 101.23);

  ndlqr_FreeNdData(nddata);
  return 1;
}

int SetSolutionFactors() {
  int nstates = 6;
  int ninputs = 3;
  int nsegments = 7;
  int width = 1;  // the solution vector
  NdData* nddata = ndlqr_NewNdData(nstates, ninputs, nsegments + 1, width);
  int factorsize = (2 * nstates + ninputs) * width;
  NdFactor* factor;

  SetNdDataBlock(nddata, 0, 1.1, 2.1, 3.1);
  ndlqr_GetNdFactor(nddata, 0, 0, &factor);
  Matrix Cy = ndlqr_GetLambdaFactor(factor);
  Matrix Cx = ndlqr_GetStateFactor(factor);
  mu_assert(Cy.cols == 1);
  mu_assert(Cy.rows == nstates);
  mu_assert(Cx.cols == 1);
  mu_assert(Cx.rows == nstates);
  mu_assert(CheckFactors(factor, 1.1, 2.1, 3.1) == 1);

  // Check Second Factor
  SetNdDataBlock(nddata, factorsize, 1.3, 2.3, 3.3);
  ndlqr_GetNdFactor(nddata, 1, 0, &factor);
  mu_assert(CheckFactors(factor, 1.3, 2.3, 3.3) == 1);

  ndlqr_FreeNdData(nddata);
  return 1;
}

void AllTests() {
  mu_run_test(NewNdDataTest);
  mu_run_test(SetFactors);
  mu_run_test(SetSolutionFactors);
}

mu_test_main
