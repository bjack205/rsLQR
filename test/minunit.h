#pragma once

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

// #define mu_assert_eq(a, b) do { if (a != b) return "Expected equality of "; } while (0)
// #define mu_assert_lt(a, b) do { if (a < b) return "Expected " + ; } while (0)
#define mu_assert(test) do { if (!(test)) { printf("TEST FAILED! %s\n%s%s:%d\n", #test, \
                                            "             ", __FILE__, __LINE__); \
                                      return 0; }} while (0)
#define mu_run_test(test) do { tests_passed += test(); tests_run++; } while (0)
#define mu_test_main int main() {ResetTests(); AllTests(); PrintTestResult(); return TestResult(); }
#define mu_test_init extern int tests_run; extern int tests_passed;
static int tests_run = 0;
static int tests_passed = 0;

bool NearlyEqual(double a, double b);

void ResetTests();

void PrintTestResult();

int TestResult();
