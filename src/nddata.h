#pragma once

#include "matrix.h"
#include "lqr_data.h"

typedef struct {
  Matrix lambda;
  Matrix state;
  Matrix input;
} NdFactor;

Matrix ndlqr_GetLambdaFactor(NdFactor* factor);
Matrix ndlqr_GetStateFactor(NdFactor* factor);
Matrix ndlqr_GetInputFactor(NdFactor* factor);

typedef struct {
  int nstates;
  int ninputs;
  int nsegments;
  int depth;
  int width;
  double* data;       // pointer to entire chunk of allocated memory
  NdFactor* factors;  // (nsegments, depth) array of factors. Stored in column-order.
} NdData;

NdData* ndlqr_NewNdData(int nstates, int ninputs, int nsegments, int width);
int ndlqr_FreeNdData(NdData* nddata);
int ndlqr_GetNdFactor(NdData* nddata, int index, int level, NdFactor** factor);
void ndlqr_ResetNdData(NdData* nddata);
