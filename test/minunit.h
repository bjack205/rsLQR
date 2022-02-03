#pragma once

#include <math.h>
#include <stdbool.h>
#include <stdio.h>

void TestFail();
void TestPass();

// #define mu_assert_eq(a, b) do { if (a != b) return "Expected equality of "; } while (0)
// #define mu_assert_lt(a, b) do { if (a < b) return "Expected " + ; } while (0)
#define mu_assert(test)                                                                 \
  do {                                                                                  \
    if (!(test)) {                                                                      \
      printf("TEST FAILED! %s\n%s%s:%d\n", #test, "             ", __FILE__, __LINE__); \
      TestFail();                                                                       \
      return 0;                                                                         \
    } else {                                                                            \
      TestPass();                                                                       \
    }                                                                                   \
  } while (0)
#define mu_run_test(test) \
  do {                    \
    test();               \
  } while (0)
#define mu_test_main     \
  int main() {           \
    ResetTests();        \
    AllTests();          \
    PrintTestResult();   \
    return TestResult(); \
  }
#define mu_test_init

void ResetTests();

void PrintTestResult();

int TestResult();
