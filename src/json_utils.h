/**
 * @file json_utils.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief Some basic utility function for working with json files. 
 *        Uses the cJSON library.
 * @version 0.1
 * @date 2022-01-31
 * 
 * @copyright Copyright (c) 2022
 * 
 * @addtogroup probdef 
 * @{
 */
#pragma once

#include "lqr_problem.h"
#include "matrix.h"

/**
 * @brief Read an LQRData structure from JSON data
 * 
 * The user is responsible for calling ndlqr_FreeLQRData() after calling this function.
 * 
 * ## Expected format
 * The json data is expected to be structured as follows:
 * 
 * ~~~~~{.json}
 * {
 *  "index": <integer>,
 *  "nstates": <integer>,
 *  "ninputs": <integer>,
 *  "Q": <array>,
 *  "R": <array>,
 *  "q": <array>,
 *  "r": <array>,
 *  "c": <double>,
 *  "A": <2D array, stored columnwise>,
 *  "B": <2D array, stored columnwise>,
 *  "d": <double>,
 * }
 * ~~~~~
 * 
 * @param filename path to the json file
 * @return An initialized LQRData structure. NULL if unsuccessful.
 */
LQRData* ndlqr_ReadLQRDataJSONFile(const char* filename);

/**
 * @brief Read and LQRProblem structure from JSON data
 * 
 * The user is expected to call ndlqr_FreeLQRProblem() after calling this function.
 * 
 * ## Expected Format
 * The json data is expected to be structured as follows:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.json}
 * {
 *   "nhorizon": <integer>,
 *   "x0": <array>,
 *   "lqrdata": <array of LQRData>,
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * where each element of `lqrdata` is of the form specified in ndlqr_ReadLQRDataJSONFile().
 * 
 * @param filename Path to the json file
 * @return         An initialized LQRProblem struct. NULL if unsuccessful.
 */
LQRProblem* ndlqr_ReadLQRProblemJSONFile(const char* filename);

/**
 * @brief Read a Matrix from a JSON file.
 * 
 * @param filename  Name of the JSON file
 * @param name      Name of the JSON field where the data is stored 
 * @return          Transfers ownership. User must call `FreeMatrix`.
 */
Matrix ReadMatrixJSONFile(const char* filename, const char* name);

/**@} */
