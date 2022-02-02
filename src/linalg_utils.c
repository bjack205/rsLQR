#include <time.h>

typedef struct {
  clock_t t_start;
  double ms_total;
} LinAlgTiming;
static LinAlgTiming la_timing = {0, 0.0};

void MatrixLinAlgTimeStart() { la_timing.t_start = clock(); }
void MatrixLinAlgTimeStop() {
  clock_t diff = clock() - la_timing.t_start;
  double millisec = diff * 1000.0 / (double)CLOCKS_PER_SEC;
  la_timing.ms_total += millisec;
}
void MatrixLinAlgTimeReset() { la_timing.ms_total = 0; }
double MatrixGetLinAlgTimeMilliseconds() { return la_timing.ms_total; }
