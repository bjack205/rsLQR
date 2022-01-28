
#ifdef MATRIX_LATIME_ENABLE
#define MATRIX_LATIME_START MatrixLinAlgTimeStart()
#define MATRIX_LATIME_STOP MatrixLinAlgTimeStop()
static const int kMatrixLinearAlgebraTimingEnabled = 1;
#else
#define MATRIX_LATIME_START
#define MATRIX_LATIME_STOP
static const int kMatrixLinearAlgebraTimingEnabled = 0;
#endif

void MatrixLinAlgTimeStart();
void MatrixLinAlgTimeStop();
void MatrixLinAlgTimeReset();
double MatrixGetLinAlgTimeMilliseconds();
