#pragma once

#include <stdbool.h>

bool IsPowerOfTwo(int x);

static inline int PowerOfTwo(int x) {
  return 1 << x;
}

void ndlqr_LinAlgTimeStart();
void ndlqr_LinAlgTimeStop();

int LogOfTwo(int x);

/**
 * @brief Read the contents of a file into a heap-allocated `char` array.
 * 
 * It is the user's responsibility to call `free` on `out` after the string data is 
 * allocated.
 * 
 * @param filename Name of the file to read.
 * @param out      Pointer to the array (pointer) to the heap-allocated string data.
 * @param len      Length of the string data.
 * @return int     0 if successful, -1 otherwise.
 */
int ReadFile(const char* filename, char** out, int* len);
