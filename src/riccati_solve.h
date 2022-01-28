#pragma once

#include "riccati_solver.h"

#include "linalg.h"

int ndlqr_SolveRiccati(RiccatiSolver* solver);
int ndlqr_BackwardPass(RiccatiSolver* solver);
int ndlqr_ForwardPass(RiccatiSolver* solver);
