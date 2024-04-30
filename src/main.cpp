#include "Matrix.hpp"

template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  for (const auto &elem : v)
    os << " " << elem << " ";
  os << std::endl;
  return os;
}

int main(int argc, char **argv) {
  algebra::Matrix<double, RowWise> M;

  // M(1, 1) = 2.;
  // M(2, 2) = 3.;

  // M(1, 1) = 4.;
  // M(1, 2) = 5;

  // std::cout << M << std::endl;

  // // std::vector<double> v1{1, 1};
  // // // std::cout << "1st vector:\n" << v1 << std::endl;

  // // // std::cout << "Result of the mul:\n" << M * v1 << std::endl;

  // // std::vector<double> v2{1, 1, 1};
  // // // std::cout << "2nd vector:\n" << v2 << std::endl;

  // // // std::cout << "Result of the mul:\n" << M * v2 << std::endl;

  // M.compress();
  // std::cout << M << std::endl;

  // // std::cout << "1st vector:\n" << v1 << std::endl;

  // // std::cout << "Result of the mul:\n" << M * v1 << std::endl;

  // // std::cout << "2nd vector:\n" << v2 << std::endl;

  // // std::cout << "Result of the mul:\n" << M * v2 << std::endl;

  // // M(1, 1) = 2.;

  // // std::cout << M << std::endl;

  // M.uncompress();
  // std::cout << M << std::endl;

  std::string filename{"lnsp_131.mtx"};

  M.read_MatrixMarket(filename);

  std::cout << M << std::endl;

  std::vector<double> v_file(M.cols(), 1);

  std::cout << "Result of the mul:\n" << M * v_file << std::endl;

  return 0;
}