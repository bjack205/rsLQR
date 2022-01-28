#include "test/minunit.h"

extern int tests_run;
extern int tests_passed;

bool NearlyEqual(double a, double b) {
  if (a == b) {
    return true;
  }
  if (a > b) {
    int c = a;
    a = b;
    b = c;
  }
  if (b <= nextafter(a, b)) {
    return true;
  } else {
    return false;
  }
}

void ResetTests() {
  tests_run = 0;
  tests_passed = 0;
}

void PrintTestResult() {
  printf("Passed %d / %d tests\n", tests_passed, tests_run);
  if (tests_run == tests_passed) {
    printf("ALL TESTS PASSED!\n");
  }
}

int TestResult() {
  return !(tests_run == tests_passed);
}
