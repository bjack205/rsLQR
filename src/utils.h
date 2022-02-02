/**
 * @file utils.h
 * @author Brian Jackson (bjack205@gmail.com)
 * @brief Provides basic functions like powers of 2 and reading files to a string.
 * @version 0.1
 * @date 2022-01-30
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <stdbool.h>

/**
 * @brief Determines if the input integer is a power of 2
 */
bool IsPowerOfTwo(int x);

/**
 * @brief Calculate 2^x efficiently for integers
 */
static inline int PowerOfTwo(int x) { return 1 << x; }

void ndlqr_LinAlgTimeStart();
void ndlqr_LinAlgTimeStop();

/**
 * @brief Efficient computation of log2(x) for integers
 */
int LogOfTwo(int x);

/**
 * @brief Read the contents of a file into a heap-allocated `char` array.
 *
 * It is the user's responsibility to call `free` on `out` after the string data is
 * allocated.
 *
 * # Example
 * The following example reads
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * const char* filename = "mydata.txt"
 * char* data = NULL;
 * int len;
 * int out = ReadFile(filename, &data, &len);
 * free(data);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * @param filename Name of the file to read.
 * @param out      Pointer to the array (pointer) to the heap-allocated string data.
 * @param len      Length of the string data.
 * @return int     0 if successful, -1 otherwise.
 */
int ReadFile(const char* filename, char** out, int* len);
