#include "json_utils.h"

#include <cjson/cJSON.h>
#include <stdio.h>
#include <stdlib.h>

#include "utils.h"

/**
 * @brief Read a JSON array into a vector of doubles.
 *
 * @param json JSON object containing the array
 * @param name Name of the JSON array
 * @param buf  Storage location for the data
 * @param len  Expected length of the array
 * @return int 0 if successful, -1 otherwise.
 */
int ReadJSONArray(cJSON* json, const char* name, double* buf, int len) {
  cJSON* array = cJSON_GetObjectItemCaseSensitive(json, name);
  cJSON* number;
  int i = 0;
  if (cJSON_IsArray(array)) {
    cJSON_ArrayForEach(number, array) {
      if (cJSON_IsNumber(number) && i < len) {
        buf[i] = number->valuedouble;
      }
      ++i;
    }
    if (i == len) {
      return 0;
    } else if (i > len) {
      fprintf(stderr, "JSON array was longer than expected (%d instead of %d).\n", i, len);
    } else {
      fprintf(stderr, "JSON array was shorter than expected (%d instead of %d).\n", i, len);
    }
  } else {
    fprintf(stderr, "Couldn't find an array of name %s.\n", name);
  }
  return -1;
}

int GetJSONMatrixSize(cJSON* json, const char* name, int* rows, int* cols) {
  cJSON* array = cJSON_GetObjectItemCaseSensitive(json, name);
  cJSON* column;
  cJSON* number;
  if (cJSON_IsArray(array)) {
    // Loop over columns
    int j = 0;
    cJSON_ArrayForEach(column, array) {
      if (cJSON_IsArray(column)) {
        // Loop over rows
        int i = 0;
        cJSON_ArrayForEach(number, column) { ++i; }
        if (j == 0) {
          *rows = i;
        } else if (*rows != i) {
          fprintf(stderr,
                  "ERROR: The number of rows changed during parsing. Failed to read as a "
                  "2D array.\n");
          rows = NULL;
          cols = NULL;
          return -1;
        }
        ++j;
      }
    }
    *cols = j;
  } else {
    rows = NULL;
    cols = NULL;
    fprintf(stderr, "ERROR: Unable to parse the field %s as a JSON array.\n", name);
    return -1;
  }
  return 0;
}

/**
 * @brief Read a column-major JSON array into a vector of doubles.
 *
 * @param json JSON object containing the array
 * @param name Name of the JSON array
 * @param buf  Storage location for the data
 * @param rows Expected number of rows in the JSON array
 * @param cols Expected number of columns in the JSON array
 * @return int 0 if successful, -1 otherwise.
 */
int ReadJSONMatrix(cJSON* json, const char* name, double* buf, int rows, int cols) {
  cJSON* array = cJSON_GetObjectItemCaseSensitive(json, name);
  cJSON* column;
  cJSON* number;
  int j = 0;
  int status = 0;
  if (cJSON_IsArray(array)) {
    // Loop over columns
    cJSON_ArrayForEach(column, array) {
      if (cJSON_IsArray(column) && j < cols) {
        // Loop over rows
        int i = 0;
        cJSON_ArrayForEach(number, column) {
          if (cJSON_IsNumber(number) && i < rows) {
            buf[i + j * rows] = number->valuedouble;
          }
          ++i;
        }
        // Check that we got the expected number of rows
        if (i != rows) {
          status = -1;
          fprintf(stderr,
                  "Got unexpected length of JSON column number %d (%d instead of %d).\n", j,
                  i, rows);
        }
        ++j;
      }
    }
    // Check that we got the expected number of columns
    if (j != cols) {
      status = -1;
      fprintf(stderr, "Got an unexpected number of JSON columns (%d instead of %d).\n", j,
              cols);
    }
  } else {
    status = -1;
    fprintf(stderr, "Couldn't find an array of name %s.\n", name);
  }
  return status;
}

/**
 * @brief Parse a JSON object containing LQRData into an LQRData type
 *
 * @param json        JSON data containing the LQR data to be parsed.
 * @param lqrdata_out Address of LQRData pointer. Needs to be freed using
 * `ndlqr_FreeLQRData`.
 * @return int        0 if successful, -1 otherwise.
 */
int ReadLQRDataJSON(cJSON* json, LQRData** lqrdata_out) {
  LQRData* lqrdata = *lqrdata_out;
  cJSON* item;
  int nstates = 0;
  item = cJSON_GetObjectItemCaseSensitive(json, "nstates");
  if (cJSON_IsNumber(item)) {
    nstates = item->valueint;
  }

  int ninputs = 0;
  item = cJSON_GetObjectItemCaseSensitive(json, "ninputs");
  if (cJSON_IsNumber(item)) {
    ninputs = item->valueint;
  }
  if (nstates != lqrdata->nstates || ninputs != lqrdata->ninputs) {
    fprintf(stderr,
            "ERROR: The state and input dimensions in the JSON file didn't match the "
            "expected dimensions\n");
    return -1;
  }

  // Read the cost constant
  item = cJSON_GetObjectItemCaseSensitive(json, "c");
  if (cJSON_IsNumber(item)) {
    *(lqrdata->c) = item->valuedouble;
  } else {
    return -1;
  }

  // Read all the array fields
  int status = 0;
  status += ReadJSONArray(json, "Q", lqrdata->Q, nstates);
  status += ReadJSONArray(json, "R", lqrdata->R, ninputs);
  status += ReadJSONArray(json, "q", lqrdata->q, nstates);
  status += ReadJSONArray(json, "r", lqrdata->r, ninputs);
  status += ReadJSONArray(json, "d", lqrdata->d, nstates);
  status += ReadJSONMatrix(json, "A", lqrdata->A, nstates, nstates);
  status += ReadJSONMatrix(json, "B", lqrdata->B, nstates, ninputs);

  if (status == 0) {
    *lqrdata_out = lqrdata;
    return 0;
  } else {
    fprintf(stderr,
            "ERROR: The LQR data file wasn't successfully parsed. Returning an empty "
            "struct.\n");
    return -1;
  }
}

LQRProblem* ndlqr_ReadLQRProblemJSONFile(const char* filename) {
  // Read the JSON data from the file into a long string
  char* jsondata = NULL;
  int len;
  if (0 != ReadFile(filename, &jsondata, &len)) {
    fprintf(stderr, "ERROR: Reading LQR Problem file failed.\n");
    return NULL;
  }

  // Parse the JSON string
  cJSON* json = cJSON_Parse(jsondata);
  if (json == NULL) {
    const char* error_ptr = cJSON_GetErrorPtr();
    if (error_ptr != NULL) {
      fprintf(stderr, "ERROR: Error parsing JSON file: %s\n", error_ptr);
    }
    cJSON_Delete(json);
    free(jsondata);
    return NULL;
  }

  cJSON* item;
  int nhorizon = 0;
  item = cJSON_GetObjectItemCaseSensitive(json, "nhorizon");
  if (cJSON_IsNumber(item)) {
    nhorizon = item->valueint;
  }

  // Get nstates and ninputs by reading the first knotpoint
  int nstates = 0;
  int ninputs = 0;
  cJSON* data_array = cJSON_GetObjectItemCaseSensitive(json, "lqrdata");
  cJSON* json_lqrdata = cJSON_GetArrayItem(data_array, 0);
  item = cJSON_GetObjectItemCaseSensitive(json_lqrdata, "nstates");
  if (cJSON_IsNumber(item)) {
    nstates = item->valueint;
  }
  item = cJSON_GetObjectItemCaseSensitive(json_lqrdata, "ninputs");
  if (cJSON_IsNumber(item)) {
    ninputs = item->valueint;
  }

  // Create an LQR Problem. Must be free later with `ndlqr_FreeLQRProblem`
  LQRProblem* lqrprob = ndlqr_NewLQRProblem(nstates, ninputs, nhorizon);

  // Loop over knot points and extract out the JSON data into `LQRData`
  if (cJSON_IsArray(data_array)) {
    cJSON_ArrayForEach(json_lqrdata, data_array) {
      cJSON* json_index = cJSON_GetObjectItemCaseSensitive(json_lqrdata, "index");
      int index = -1;
      if (cJSON_IsNumber(json_index)) {
        index = json_index->valueint - 1;  // Since Julia is 1-based indexed
      }

      LQRData* lqrdata = lqrprob->lqrdata[index];
      if (0 != ReadLQRDataJSON(json_lqrdata, &lqrdata)) {
        fprintf(stderr, "WARNING: Failed to parse the LQR JSON data at index %d\n", index);
      }
    }
  }

  // Get initial state
  int out = ReadJSONArray(json, "x0", lqrprob->x0, nstates);

  cJSON_Delete(json);
  free(jsondata);

  if (out != 0) {
    ndlqr_FreeLQRProblem(lqrprob);
    return NULL;
  }

  return lqrprob;
}

LQRData* ndlqr_ReadLQRDataJSONFile(const char* filename) {
  char* jsondata = NULL;
  int len;
  if (0 != ReadFile(filename, &jsondata, &len)) {
    fprintf(stderr, "ERROR: Reading LQR file failed.\n");
    return NULL;
  }

  cJSON* json = cJSON_Parse(jsondata);
  free(jsondata);
  if (json == NULL) {
    const char* error_ptr = cJSON_GetErrorPtr();
    if (error_ptr != NULL) {
      fprintf(stderr, "ERROR: Error parsing JSON file: %s\n", error_ptr);
    }
    return NULL;
  }

  // Read state and control dimension
  cJSON* item;
  int nstates = 0;
  item = cJSON_GetObjectItemCaseSensitive(json, "nstates");
  if (cJSON_IsNumber(item)) {
    nstates = item->valueint;
  }

  int ninputs = 0;
  item = cJSON_GetObjectItemCaseSensitive(json, "ninputs");
  if (cJSON_IsNumber(item)) {
    ninputs = item->valueint;
  }
  if (ninputs <= 0 || nstates <= 0) {
    fprintf(stderr,
            "ERROR: Couldn't get a valid state and control dimension from the LQR Data "
            "file: %s\n",
            filename);
    return NULL;
  }

  LQRData* lqrdata = ndlqr_NewLQRData(nstates, ninputs);
  int out = ReadLQRDataJSON(json, &lqrdata);
  cJSON_Delete(json);
  if (out == 0) {
    return lqrdata;
  } else {
    fprintf(stderr, "ERROR: Error parsing LQR JSON data.\n");
    return NULL;
  }
}

Matrix ReadMatrixJSONFile(const char* filename, const char* name) {
  // TODO: allow this to read 1D arrays as well
  Matrix nullmat = {0, 0, NULL};

  char* jsondata = NULL;
  int len;
  if (0 != ReadFile(filename, &jsondata, &len)) {
    fprintf(stderr, "ERROR: Reading LQR file failed.\n");
    return nullmat;
  }

  cJSON* json = cJSON_Parse(jsondata);
  free(jsondata);
  if (json == NULL) {
    const char* error_ptr = cJSON_GetErrorPtr();
    if (error_ptr != NULL) {
      fprintf(stderr, "ERROR: Error parsing JSON file: %s\n", error_ptr);
    }
    return nullmat;
  }

  int rows, cols;
  if (0 != GetJSONMatrixSize(json, name, &rows, &cols)) {
    cJSON_Delete(json);
    fprintf(stderr, "ERROR: Couldn't get the size of the JSON Matrix.\n");
    return nullmat;
  }
  // printf("Rows = %d, cols = %d\n", rows, cols);
  Matrix mat = NewMatrix(rows, cols);
  if (0 != ReadJSONMatrix(json, name, mat.data, rows, cols)) {
    cJSON_Delete(json);
    FreeMatrix(&mat);
    fprintf(stderr, "ERROR: Wasn't able to parse the JSON Matrix.\n");
    return nullmat;
  }
  cJSON_Delete(json);
  return mat;
}
