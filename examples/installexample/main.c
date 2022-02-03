#include <stdio.h>
#include "rslqr.h"

int main() {
    printf("This is the example using an installed rsLQR package!.\n\n");
    // Read an LQRProblem from a json file
    const char* filename = LQRPROBFILE;
    LQRProblem* lqrprob = ndlqr_ReadLQRProblemJSONFile(filename);
    int nstates = lqrprob->lqrdata[0]->nstates;
    int ninputs = lqrprob->lqrdata[0]->ninputs;
    int nhorizon = lqrprob->nhorizon;

    // Initialize the solver and solve 
    NdLqrSolver* solver = ndlqr_NewNdLqrSolver(nstates, ninputs, nhorizon);
    if (ndlqr_InitializeWithLQRProblem(lqrprob, solver) == 0) {
        ndlqr_Solve(solver);
    }

    // Print the solve summary
    ndlqr_PrintSolveSummary(solver);

    // Free the problem data and the solver
    ndlqr_FreeLQRProblem(lqrprob);
    ndlqr_FreeNdLqrSolver(solver);
    return 0;
}
