#include "test/minunit.h"

// extern int tests_run;
// extern int tests_passed;

int tests_run = 0;
int tests_passed = 0;

void TestFail() { ++tests_passed; }
void TestPass() {
  ++tests_passed;
  ++tests_run;
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

int TestResult() { return !(tests_run == tests_passed); }
