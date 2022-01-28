#pragma once

#include "lqr_problem.h"
#include "matrix.h"

LQRData* ndlqr_ReadLQRDataJSONFile(const char* filename);

LQRProblem* ndlqr_ReadLQRProblemJSONFile(const char* filename);

/**
 * @brief Read a Matrix from a JSON file.
 * 
 * @param filename  Name of the JSON file
 * @param name      Name of the JSON field where the data is stored 
 * @return          Transfers ownership. User must call `FreeMatrix`.
 */
Matrix ReadMatrixJSONFile(const char* filename, const char* name);
