#include "Matrix.hpp"

template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  for (const auto &elem : v)
    os << " " << elem << " ";
  os << std::endl;
  return os;
}

void func1(auto &M) {
  // Create the  following matrix by dynamically assigning new values
  //[4, 0, 0,
  // 0, 7, 5,
  // 0, 6, 0]
  // Then prints it, compress it, print it again, uncompress and print one last
  // time to check the conversions between the two format are ok
  std::cout << "\nTest 1" << std::endl;

  M(0, 0) = 4.;
  M(1, 1) = 7.;
  M(1, 2) = 5.;
  M(2, 1) = 6.;

  std::cout << M << std::endl;

  M.compress();
  std::cout << M << std::endl;

  M.uncompress();
  std::cout << M << std::endl;
}

void func2(const auto &M) {
  // Tries to access a 0 and a non-0 element in a dynamic format matrix passed
  // through a const reference
  std::cout << "\nTest 2" << std::endl;

  std::cout << "(0, 0): " << M(0, 0) << std::endl;
  std::cout << "(1, 1): " << M(1, 1) << std::endl;
}

void func3(const auto &M) {
  // Tries to access a 0 and a non-0 element in a compressed format matrix
  // passed through a const reference It is in fact the same test as before but
  // since the matrix is const the compression cannot be done here
  std::cout << "\nTest 3" << std::endl;

  std::cout << "(0, 0): " << M(0, 0) << std::endl;
  std::cout << "(1, 1): " << M(1, 1) << std::endl;
}

void func4(auto M) {
  // Tries to change an element and dynamically add a new one to a dynamic
  // format matrix passed through a copy
  std::cout << "\nTest 4" << std::endl;

  M(0, 0) = 2.;
  std::cout << "Element (0, 0) changed" << std::endl;
  M(1, 1) = 1.;
  std::cout << "Element (1, 1) added" << std::endl;

  std::cout << M << std::endl;
}

void func5(auto &M) {
  // Tries to change an element and dynamically add a new one to a compressed
  // format matrix passed through a copy It is in fact the same test as before,
  // but this time adding a new one should fail
  std::cout << "\nTest 5" << std::endl;

  M(0, 0) = 2.;
  std::cout << "Element (0, 0) changed" << std::endl;
  M(1, 1) = 1.;
  std::cout << "Element (1, 1) added" << std::endl;

  std::cout << M << std::endl;
}

int main(int argc, char **argv) {
  algebra::Matrix<double, RowWise> M1;

  func1(M1);

  algebra::Matrix<double, RowWise> M2(2, 2);
  M2(0, 0) = 1.;

  func2(M2);

  M2.compress();
  func3(M2);

  algebra::Matrix<double, RowWise> M3(1, 1);
  M3(0, 0) = 100.;

  func4(M3);

  //! ATTENTION: This test is set up to throw an exception
  //! Comment it away to proceed with the other tests
  M3.compress();
  func5(M3);

  // // std::vector<double> v1{1, 1};
  // // // std::cout << "1st vector:\n" << v1 << std::endl;

  // // // std::cout << "Result of the mul:\n" << M * v1 << std::endl;

  // // std::vector<double> v2{1, 1, 1};
  // // // std::cout << "2nd vector:\n" << v2 << std::endl;

  // // // std::cout << "Result of the mul:\n" << M * v2 << std::endl;

  // M.compress();
  // M(1, 1) = 8.;
  // std::cout << M << std::endl;

  // M(2, 1) = 8.;
  // std::cout << M << std::endl;

  // // std::cout << "1st vector:\n" << v1 << std::endl;

  // // std::cout << "Result of the mul:\n" << M * v1 << std::endl;

  // // std::cout << "2nd vector:\n" << v2 << std::endl;

  // // std::cout << "Result of the mul:\n" << M * v2 << std::endl;

  // M(1, 1) = 2.;

  // std::cout << M << std::endl;

  // M.uncompress();
  // std::cout << M << std::endl;

  // std::string filename{"lnsp_131.mtx"};

  // algebra::Matrix<double, RowWise> N;

  // N.read_MatrixMarket(filename);

  // std::cout << N << std::endl;

  // std::vector<double> v_file(M.cols(), 1);

  // std::cout << "Result of the mul:\n" << M * v_file << std::endl;

  return 0;
}