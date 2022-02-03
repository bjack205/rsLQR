#include "utils.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

bool IsPowerOfTwo(int x) { return (x != 0) && ((x & (x - 1)) == 0); }

int LogOfTwo(int x) {
  int shift = 0;
  while (((x >> shift) & 1) != 1) {
    shift++;
  }
  return shift;
}

int ReadFile(const char* filename, char** out, int* len) {
  FILE* fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Couldn't open file\n");
    return -1;
  }

  fseek(fp, 0L, SEEK_END);
  int size = ftell(fp);
  rewind(fp);

  char* buf = (char*)malloc(size + 1);
  if (!buf) {
    fclose(fp);
    fprintf(stderr, "Couldn't allocate memory for the file contents.");
    return -1;
  }

  if (1 != fread(buf, size, 1, fp)) {
    fprintf(stderr, "Failed to read the entire file.");
    free(buf);
    fclose(fp);
    return -1;
  }

  // Add null-termination
  buf[size] = '\0';

  *out = buf;
  *len = size;
  fclose(fp);
  return 0;
}
