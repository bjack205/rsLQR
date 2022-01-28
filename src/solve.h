#include "solver.h"

int ndlqr_Solve(NdLqrSolver* solver);

Matrix ndlqr_GetSolution(NdLqrSolver* solver);

int ndlqr_CopySolution(NdLqrSolver* solver, double* soln);
