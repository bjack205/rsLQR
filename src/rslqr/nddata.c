#include "rslqr/nddata.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

Matrix ndlqr_GetLambdaFactor(NdFactor* factor) { return factor->lambda; }

Matrix ndlqr_GetStateFactor(NdFactor* factor) { return factor->state; }
Matrix ndlqr_GetInputFactor(NdFactor* factor) { return factor->input; }

NdData* ndlqr_NewNdData(int nstates, int ninputs, int nhorizon, int width) {
  int nsegments = nhorizon - 1;
  if (nstates <= 0 || ninputs <= 0 || nsegments <= 0) return NULL;
  if (!IsPowerOfTwo(nhorizon)) {
    fprintf(stderr, "ERROR: Number of segments must be one less than a power of 2.\n");
    return NULL;
  }

  // A little hacky, but set depth to 1 for the rhs vector
  int depth;
  if (width == 1) {
    depth = 1;
  } else {
    depth = LogOfTwo(nhorizon);
  }

  // Allocate one large block of memory for the data
  int numfactors = nhorizon * depth;
  int factorsize = (2 * nstates + ninputs) * width;
  double* data = (double*)calloc(numfactors * factorsize, sizeof(double));
  if (data == NULL) {
    fprintf(stderr, "ERROR: Failed to allocate memory for NdData.\n");
    return NULL;
  }

  // Create the factors using the allocated memory
  NdFactor* factors = (NdFactor*)malloc(numfactors * sizeof(NdFactor));
  for (int i = 0; i < numfactors; ++i) {
    double* factordata = data + i * factorsize;
    factors[i].lambda.rows = nstates;
    factors[i].lambda.cols = width;
    factors[i].lambda.data = factordata;
    factors[i].state.rows = nstates;
    factors[i].state.cols = width;
    factors[i].state.data = factordata + nstates * width;
    factors[i].input.rows = ninputs;
    factors[i].input.cols = width;
    factors[i].input.data = factordata + 2 * (nstates * width);
  }

  // Create the NdData struct
  NdData* nddata = (NdData*)malloc(sizeof(NdData));
  nddata->nstates = nstates;
  nddata->ninputs = ninputs;
  nddata->nsegments = nsegments;
  nddata->depth = depth;
  nddata->width = width;
  nddata->data = data;
  nddata->factors = factors;
  return nddata;
}

void ndlqr_ResetNdData(NdData* nddata) {
  int nhorizon = nddata->nsegments + 1;
  int numfactors = nhorizon * nddata->depth;
  int factorsize = (2 * nddata->nstates + nddata->ninputs) * nddata->width;
  memset(nddata->data, 0, numfactors * factorsize * sizeof(double));
}

int ndlqr_FreeNdData(NdData* nddata) {
  if (!nddata) return -1;
  free(nddata->factors);
  free(nddata->data);
  free(nddata);
  return 0;
}

int ndlqr_GetNdFactor(NdData* nddata, int index, int level, NdFactor** factor) {
  if (index < 0 || index > nddata->nsegments) {
    fprintf(stderr, "Invalid index. Must be between %d and %d, got %d.\n", 0,
            nddata->nsegments, index);
    return -1;
  }
  if (level < 0 || level >= nddata->depth) {
    fprintf(stderr, "Invalid level. Must be between %d and %d, got %d.\n", 0,
            nddata->depth - 1, level);
    return -1;
  }
  int linear_index = index + (nddata->nsegments + 1) * level;
  *factor = nddata->factors + linear_index;
  return 0;
}
