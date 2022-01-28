#include <Eigen/Dense>

using MapMatrixXd = Eigen::Map<Eigen::MatrixXd>;

extern "C" {

void eigen_SetNumThreads(int n) { Eigen::setNbThreads(n); }

void eigen_InitParallel() { Eigen::initParallel(); }

void eigen_MatrixAddition(int n, double* a, double* b, double alpha) {
  MapMatrixXd A(a, n, 1);
  MapMatrixXd B(b, n, 1);
  B = alpha * A + B;
}

void eigen_MatrixMultiply(int m, int n, int k, double* a, double* b, double* c,
                          bool tA, bool tB, double alpha, double beta) {
  int rowA = tA ? n : m;
  int colA = tA ? m : n;
  int rowB = tB ? k : n;
  int colB = tB ? n : k;
  Eigen::Map<Eigen::MatrixXd> A(a, rowA, colA);
  Eigen::Map<Eigen::MatrixXd> B(b, rowB, colB);
  Eigen::Map<Eigen::MatrixXd> C(c, m, k);
  if (!tA && !tB) {
    C = (A * B) * alpha + beta * C;
  } else if (tA && !tB) {
    C = (A.transpose() * B) * alpha + beta * C;
  } else if (tA && tB) {
    C = (A.transpose() * B.transpose()) * alpha + beta * C;
  } else if (!tA && tB) {
    C = (A * B.transpose()) * alpha + beta * C;
  }
}

int eigen_CholeskyFactorize(int n, double* a, void** fact) {
  using LLT = Eigen::LLT<Eigen::Ref<MapMatrixXd>>;
  MapMatrixXd A(a, n, n);
  LLT* llt = new LLT(A);
  *fact = static_cast<void*>(llt);
  int info = static_cast<int>(llt->info());
  return info;  // 0 is success
}

void eigen_FreeFactorization(void* achol) {
  using LLT = Eigen::LLT<Eigen::Ref<MapMatrixXd>>;
  LLT* llt = static_cast<LLT*>(achol);
  delete llt;
}

void eigen_CholeskySolve(int n, int m, void* achol, double* b) {
  using LLT = Eigen::LLT<Eigen::Ref<MapMatrixXd>>;
  MapMatrixXd B(b, n, m);

  LLT* llt = static_cast<LLT*>(achol);
  llt->solveInPlace(B);
}

void eigen_SymmetricMatrixMultiply(int n, int m, double* a, double* b,
                                   double* c) {
  Eigen::Map<Eigen::MatrixXd> A(a, n, m);
  Eigen::Map<Eigen::MatrixXd> B(b, n, n);
  Eigen::Map<Eigen::MatrixXd> C(c, n, m);
  C = A.selfadjointView<Eigen::Lower>() * B;
}

void eigen_MatrixMultiply8x8(double* a, double* b, double* c) {
  Eigen::Map<Eigen::Matrix<double, 8, 8>> A(a);
  Eigen::Map<Eigen::Matrix<double, 8, 8>> B(b);
  Eigen::Map<Eigen::Matrix<double, 8, 8>> C(c);
  C = A * B;
}

void eigen_MatrixMultiply6x6(double* a, double* b, double* c) {
  Eigen::Map<Eigen::Matrix<double, 6, 6>> A(a);
  Eigen::Map<Eigen::Matrix<double, 6, 6>> B(b);
  Eigen::Map<Eigen::Matrix<double, 6, 6>> C(c);
  C = A * B;
}

void eigen_MatrixMultiply6x3(double* a, double* b, double* c) {
  Eigen::Map<Eigen::Matrix<double, 6, 6>> A(a);
  Eigen::Map<Eigen::Matrix<double, 6, 6>> B(b);
  Eigen::Map<Eigen::Matrix<double, 6, 6>> C(c);
  C = A * B;
}

}  // extern "C"