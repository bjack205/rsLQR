#pragma once

#include "matrix.h"

typedef struct {
  int nstates;
  int ninputs;
  double* Q;
  double* R;
  double* q;
  double* r;
  double* c;
  double* A;
  double* B;
  double* d;
} LQRData;

int ndlqr_InitializeLQRData(LQRData* lqrdata, double* Q, double* R, double* q, 
                                 double* r, double c, double* A, double* B, double* d);
LQRData* ndlqr_NewLQRData(int nstates, int ninputs);
int ndlqr_FreeLQRData(LQRData* lqrdata);
int ndlqr_CopyLQRData(LQRData* dest, LQRData* src);

Matrix ndlqr_GetA(LQRData* lqrdata);
Matrix ndlqr_GetB(LQRData* lqrdata);
Matrix ndlqr_Getd(LQRData* lqrdata);
Matrix ndlqr_GetQ(LQRData* lqrdata);
Matrix ndlqr_GetR(LQRData* lqrdata);
Matrix ndlqr_Getq(LQRData* lqrdata);
Matrix ndlqr_Getr(LQRData* lqrdata);

void ndlqr_PrintLQRData(LQRData* lqrdata);
