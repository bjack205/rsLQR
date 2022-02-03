#include "utils.h"

#include <stdlib.h>
#include <string.h>

#include "lqr/json_utils.h"
#include "matrix/matrix.h"
#include "test/minunit.h"

#ifndef LQRDATAFILE
#define LQRDATAFILE "../lqr_data.json"
#endif

mu_test_init

    // clang-format off
int ReadFileTest() {
  // clang-format on
  const char* filename = LQRDATAFILE;
  char* data = NULL;
  int len;
  int out = ReadFile(filename, &data, &len);
  mu_assert((int)strlen(data) == len);
  mu_assert(data[len] == '\0');
  mu_assert(out == 0);
  mu_assert(len == 473);
  mu_assert(data[0] == '{');
  mu_assert(data[1] == '"');
  mu_assert(data[2] == 'i');
  free(data);

  out = ReadFile(LQRPROBFILE, &data, &len);
  mu_assert(out == 0);
  // mu_assert(len == 12974);
  mu_assert(data[0] == '{');
  mu_assert(data[1] == '\n');
  mu_assert(data[5] == 'n');
  free(data);
  return 1;
}

int ReadMatrixFromJSON() {
  Matrix mat = ReadMatrixJSONFile(SAMPLEPROBFILE, "test");
  for (int i = 0; i < 12; ++i) {
    mu_assert(mat.data[i] == i + 1);
  }
  FreeMatrix(&mat);
  return 1;
}

void AllTests() {
  mu_run_test(ReadFileTest);
  mu_run_test(ReadMatrixFromJSON);
}

mu_test_main
