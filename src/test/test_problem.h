#pragma once

#include "lqr_problem.h"
#include "json_utils.h"
#include "solver.h"

#ifndef LQRDATAFILE
#define LQRDATAFILE "../lqr_data.json"
#endif

#ifndef LQRPROBFILE
#define LQRPROBFILE "../lqr_data.json"
#endif

#ifndef LQRPROB256FILE
#define LQRPROBFILE "../lqr_data_256.json"
#endif

LQRData* ndlqr_ReadTestLQRData();

LQRProblem* ndlqr_ReadTestLQRProblem();

NdLqrSolver* ndlqr_GenTestSolver();

LQRProblem* ndlqr_ReadLongTestLQRProblem();

NdLqrSolver* ndlqr_GenLongTestSolver();
