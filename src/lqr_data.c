#include "lqr_data.h"

#include <cjson/cJSON.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

int ndlqr_InitializeLQRData(LQRData* lqrdata, double* Q, double* R, double* q,
                            double* r, double c, double* A, double* B,
                            double* d) {
  if (!lqrdata) return -1;
  memcpy(lqrdata->Q, Q, lqrdata->nstates * sizeof(double));
  memcpy(lqrdata->R, R, lqrdata->ninputs * sizeof(double));
  memcpy(lqrdata->q, q, lqrdata->nstates * sizeof(double));
  memcpy(lqrdata->r, r, lqrdata->ninputs * sizeof(double));
  *lqrdata->c = c;
  memcpy(lqrdata->A, A, lqrdata->nstates * lqrdata->nstates * sizeof(double));
  memcpy(lqrdata->B, B, lqrdata->nstates * lqrdata->ninputs * sizeof(double));
  memcpy(lqrdata->d, d, lqrdata->nstates * sizeof(double));
  return 0;
}

LQRData* ndlqr_NewLQRData(int nstates, int ninputs) {
  int cost_size = 2 * nstates + 2 * ninputs + 1;  // Q,R,q,r,c
  int dynamics_size = nstates * nstates + nstates * ninputs + nstates;  // A,B,d
  int total_size = cost_size + dynamics_size;
  double* data = (double*)malloc(total_size * sizeof(double));
  double* Q = data;
  double* R = data + nstates;
  double* q = data + nstates + ninputs;
  double* r = data + 2 * nstates + ninputs;
  double* c = data + 2 * nstates + 2 * ninputs;
  double* A = data + cost_size;
  double* B = data + cost_size + nstates * nstates;
  double* d = data + cost_size + nstates * nstates + nstates * ninputs;
  LQRData* lqrdata = (LQRData*)malloc(sizeof(LQRData));
  lqrdata->nstates = nstates;
  lqrdata->ninputs = ninputs;
  lqrdata->Q = Q;
  lqrdata->R = R;
  lqrdata->q = q;
  lqrdata->r = r;
  lqrdata->c = c;
  lqrdata->A = A;
  lqrdata->B = B;
  lqrdata->d = d;
  return lqrdata;
}

int ndlqr_FreeLQRData(LQRData* lqrdata) {
  if (!lqrdata) return -1;
  if (lqrdata->Q) {
    free(
        lqrdata
            ->Q);  // This points to the beginning of the allocated memory block
  }
  free(lqrdata);
  return 0;
}

int ndlqr_CopyLQRData(LQRData* dest, LQRData* src) {
  if (dest->nstates != src->nstates || dest->ninputs != src->ninputs) {
    fprintf(stderr,
            "Can't copy LQRData of different sizes: (%d,%d) and (%d,%d).\n",
            dest->nstates, dest->ninputs, src->nstates, src->ninputs);
    return -1;
  }
  int nstates = dest->nstates;
  int ninputs = dest->ninputs;
  int cost_size = 2 * nstates + 2 * ninputs + 1;  // Q,R,q,r,c
  int dynamics_size = nstates * nstates + nstates * ninputs + nstates;  // A,B,d
  int total_size = cost_size + dynamics_size;
  memcpy(dest->Q, src->Q, total_size * sizeof(double));
  return 0;
}

Matrix ndlqr_GetA(LQRData* lqrdata) {
  Matrix mat = {lqrdata->nstates, lqrdata->nstates, lqrdata->A};
  return mat;
}

Matrix ndlqr_GetB(LQRData* lqrdata) {
  Matrix mat = {lqrdata->nstates, lqrdata->ninputs, lqrdata->B};
  return mat;
}

Matrix ndlqr_Getd(LQRData* lqrdata) {
  Matrix mat = {lqrdata->nstates, 1, lqrdata->d};
  return mat;
}

Matrix ndlqr_GetQ(LQRData* lqrdata) {
  Matrix mat = {lqrdata->nstates, 1, lqrdata->Q};
  return mat;
}

Matrix ndlqr_Getq(LQRData* lqrdata) {
  Matrix mat = {lqrdata->nstates, 1, lqrdata->q};
  return mat;
}

Matrix ndlqr_GetR(LQRData* lqrdata) {
  Matrix mat = {lqrdata->ninputs, 1, lqrdata->R};
  return mat;
}

Matrix ndlqr_Getr(LQRData* lqrdata) {
  Matrix mat = {lqrdata->ninputs, 1, lqrdata->r};
  return mat;
}

void PrintAsRow(Matrix* mat) {
  printf("[");
  for (int i = 0; i < MatrixNumElements(mat); ++i) {
    printf("%6.2f ", mat->data[i]);
  }
  printf("]\n");
}

void ndlqr_PrintLQRData(LQRData* lqrdata) {
  // clang-format off
  printf("LQR Data with n=%d, m=%d:\n", lqrdata->nstates, lqrdata->ninputs);
  Matrix mat = ndlqr_GetQ(lqrdata);
  printf("Q = "); PrintAsRow(&mat);
  mat = ndlqr_GetR(lqrdata);
  printf("R = "); PrintAsRow(&mat);
  mat = ndlqr_Getq(lqrdata);
  printf("q = "); PrintAsRow(&mat);
  mat = ndlqr_Getr(lqrdata);
  printf("r = "); PrintAsRow(&mat);
  printf("c = %f\n", *lqrdata->c);
  printf("A:\n");
  mat = ndlqr_GetA(lqrdata);
  PrintMatrix(&mat);
  printf("B:\n");
  mat = ndlqr_GetB(lqrdata);
  PrintMatrix(&mat);
  mat = ndlqr_Getd(lqrdata);
  printf("d = "); PrintAsRow(&mat);
  // clang-format on
}
