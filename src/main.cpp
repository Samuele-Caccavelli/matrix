#include "Matrix.hpp"
#include "chrono.hpp"
#include <chrono>
#include <cxxabi.h>
#include <typeinfo>

template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  for (const auto &elem : v)
    os << " " << elem << " ";
  os << std::endl;
  return os;
}

// Create the following matrix by dynamically assigning new values
// |4, 0, 0|
// |0, 7, 5|
// |0, 6, 0|
// Then prints it, compress it, print it again, uncompress and print one last
// time to check the conversions between the two format are ok
void test1(auto &M);

// Tries to access a 0 and a non-0 element in a matrix passed through a const
// reference
void test2(const auto &M);

// Tries to access a 0 and a non-0 element in a matrix passed through a
// non-const reference
void test3(auto &M);

// Tries to change an element and dynamically add a new one to a matrix passed
// through a copy
void test4(auto M);

// Perform a matrix x vector operation
void test5(const auto &M);

// Perform a matrix x vector operation without cluttering the terminal with the
// matrix
void test5_noprint(const auto &M);

// Perform a matrix x complex vector operation
//! The code to print the type is taken online
//! It should work for both GCC and Clang compilers
void test5_complex(const auto &M);

// Computes the 3 different norms of a matrix
void test6(const auto &M);

int main(int argc, char **argv) {
  //---------------------------------------------------------------------------
  // TEST FOR ROW-WISE MATRICES
  //---------------------------------------------------------------------------
  std::cout << "-----TESTS FOR ROW-WISE MATRICES-----" << std::endl;

  std::cout << "\nFORMAT CHANGE TEST" << std::endl;
  algebra::Matrix<double, algebra::RowWise> M1(3, 3);
  M1(0, 0) = 4.;
  M1(1, 1) = 7.;
  M1(1, 2) = 5.;
  M1(2, 1) = 6.;

  test1(M1);

  std::cout << "\nCONST ACCESS TEST" << std::endl;
  algebra::Matrix<double, algebra::RowWise> M2(2, 2);
  M2(0, 0) = 1.;
  std::cout << "Dynamic format" << std::endl;
  test2(M2);
  std::cout << "Compressed format" << std::endl;
  M2.compress();
  test2(M2);
  M2.uncompress();

  std::cout << "\nNON-CONST ACCESS TEST" << std::endl;
  std::cout << "Dynamic format" << std::endl;
  test3(M2);
  std::cout << "Compressed format" << std::endl;
  M2.compress();
  test3(M2);
  M2.uncompress();

  std::cout << "\nADDING AND CHANGING AN ELEMENT TEST" << std::endl;
  algebra::Matrix<double, algebra::RowWise> M3(1, 1);
  M3(0, 0) = 100.;
  std::cout << "Dynamic format" << std::endl;
  test4(M3);
  // //! ATTENTION: This test is set up to throw an exception
  // //! Comment it away to proceed with the other tests
  // std::cout << "Compressed format" << std::endl;
  // M3.compress();
  // test4(M3);

  std::cout << "\nMULTIPLICATION TEST" << std::endl;
  algebra::Matrix<double, algebra::RowWise> M4;
  M4(0, 0) = 4.;
  M4(1, 1) = 7.;
  M4(1, 2) = 5.;
  M4(2, 1) = 6.;
  // The matrix is
  // |4, 0, 0|
  // |0, 7, 5|
  // |0, 6, 0|
  // While the vector is
  // |1, 2, 3|'
  // So the expected value is
  // |4, 29, 12|'
  std::cout << "Dynamic format" << std::endl;
  test5(M4);
  std::cout << "The expected value is: [4, 29, 12]\n" << std::endl;
  std::cout << "Compressed format" << std::endl;
  M4.compress();
  test5(M4);
  std::cout << "The expected value is: [4, 29, 12]" << std::endl;

  std::cout << "\nBIG MULTIPLICATION TEST" << std::endl;
  std::string filename{"lnsp_131.mtx"};
  algebra::Matrix<double, algebra::RowWise> M5;
  M5.read_MatrixMarket(filename);
  std::cout << "Dynamic format" << std::endl;
  test5_noprint(M5);
  std::cout << "Compressed format" << std::endl;
  M5.compress();
  test5_noprint(M5);

  std::cout << "\nNORM TEST" << std::endl;
  // The matrix is
  // |4, 0, 0|
  // |0, 7, 5|
  // |0, 6, 0|
  // So its norms are respectively
  // Infinity   = 12
  // One        = 13
  // Frobenius ~= 11.225
  std::cout << "Dynamic format" << std::endl;
  M4.uncompress();
  test6(M4);
  std::cout << "The expected value is: [12, 13, 11.225]\n" << std::endl;
  std::cout << "Compressed format" << std::endl;
  M4.compress();
  test6(M4);
  std::cout << "The expected value is: [12, 13, 11.225]" << std::endl;

  std::cout << "\nBIG NORM TEST" << std::endl;
  std::cout << "Dynamic format" << std::endl;
  M5.uncompress();
  test6(M5);
  std::cout << "Compressed format" << std::endl;
  M5.compress();
  test6(M5);

  //---------------------------------------------------------------------------
  // TEST FOR ROW-WISE COMPLEX MATRICES
  //---------------------------------------------------------------------------
  std::cout << "\n\n-----TESTS FOR ROW-WISE COMPLEX MATRICES-----" << std::endl;
  using namespace std::complex_literals;

  std::cout << "\nMULTIPLICATION TEST" << std::endl;
  algebra::Matrix<std::complex<double>, algebra::RowWise> Mc1;
  Mc1(0, 0) = 1.0 + 1i;
  Mc1(1, 1) = 2. + 2i;
  // The matrix is
  // |1+i,  0  |
  // | 0 , 2+2i|
  // While the vector is
  // |1, i|'
  // So the expected value is
  // |1+i, -2+2i|'
  std::cout << "Dynamic format" << std::endl;
  test5_complex(Mc1);
  std::cout << "The expected value is: [1+i, -2+2i]\n" << std::endl;
  std::cout << "Compressed format" << std::endl;
  Mc1.compress();
  test5_complex(Mc1);
  std::cout << "The expected value is: [1+i, -2+2i]" << std::endl;

  std::cout << "\nNORM TEST" << std::endl;
  Mc1.uncompress();
  Mc1(0, 1) = 1i;
  // The matrix is
  // |1+i,    i|
  // | 0 , 2+2i|
  // So its norms are respectively
  // Infinity  ~= 2.8284
  // One       ~= 3.8284
  // Frobenius ~= 3.3166
  std::cout << "Dynamic format" << std::endl;
  test6(Mc1);
  std::cout << "The expected value is: [2.82843, 3.82843, 3.31662]\n"
            << std::endl;
  std::cout << "Compressed format" << std::endl;
  Mc1.compress();
  test6(Mc1);
  std::cout << "The expected value is: [2.82843, 3.82843, 3.31662]"
            << std::endl;

  return 0;
}

// int main(int argc, char **argv) {
//   //---------------------------------------------------------------------------
//   // TEST FOR COLUMN-WISE MATRICES
//   //---------------------------------------------------------------------------
//   std::cout << "-----TESTS FOR COLUMN-WISE MATRICES-----" << std::endl;

//   std::cout << "\nFORMAT CHANGE TEST" << std::endl;
//   algebra::Matrix<double, algebra::ColumnWise> M1(3, 3);
//   M1(0, 0) = 4.;
//   M1(1, 1) = 7.;
//   M1(1, 2) = 5.;
//   M1(2, 1) = 6.;

//   test1(M1);

//   std::cout << "\nCONST ACCESS TEST" << std::endl;
//   algebra::Matrix<double, algebra::ColumnWise> M2(2, 2);
//   M2(0, 0) = 1.;
//   std::cout << "Dynamic format" << std::endl;
//   test2(M2);
//   std::cout << "Compressed format" << std::endl;
//   M2.compress();
//   test2(M2);
//   M2.uncompress();

//   std::cout << "\nNON-CONST ACCESS TEST" << std::endl;
//   std::cout << "Dynamic format" << std::endl;
//   test3(M2);
//   std::cout << "Compressed format" << std::endl;
//   M2.compress();
//   test3(M2);
//   M2.uncompress();

//   std::cout << "\nADDING AND CHANGING AN ELEMENT TEST" << std::endl;
//   algebra::Matrix<double, algebra::ColumnWise> M3(1, 1);
//   M3(0, 0) = 100.;
//   std::cout << "Dynamic format" << std::endl;
//   test4(M3);
//   // //! ATTENTION: This test is set up to throw an exception
//   // //! Comment it away to proceed with the other tests
//   // std::cout << "Compressed format" << std::endl;
//   // M3.compress();
//   // test4(M3);

//   std::cout << "\nMULTIPLICATION TEST" << std::endl;
//   algebra::Matrix<double, algebra::ColumnWise> M4;
//   M4(0, 0) = 4.;
//   M4(1, 1) = 7.;
//   M4(1, 2) = 5.;
//   M4(2, 1) = 6.;
//   // The matrix is
//   // |4, 0, 0|
//   // |0, 7, 5|
//   // |0, 6, 0|
//   // While the vector is
//   // |1, 2, 3|'
//   // So the expected value is
//   // |4, 29, 12|'
//   std::cout << "Dynamic format" << std::endl;
//   test5(M4);
//   std::cout << "The expected value is: [4, 29, 12]\n" << std::endl;
//   std::cout << "Compressed format" << std::endl;
//   M4.compress();
//   test5(M4);
//   std::cout << "The expected value is: [4, 29, 12]" << std::endl;

//   std::cout << "\nBIG MULTIPLICATION TEST" << std::endl;
//   std::string filename{"lnsp_131.mtx"};
//   algebra::Matrix<double, algebra::ColumnWise> M5;
//   M5.read_MatrixMarket(filename);
//   std::cout << "Dynamic format" << std::endl;
//   test5_noprint(M5);
//   std::cout << "Compressed format" << std::endl;
//   M5.compress();
//   test5_noprint(M5);

//   std::cout << "\nNORM TEST" << std::endl;
//   // The matrix is
//   // |4, 0, 0|
//   // |0, 7, 5|
//   // |0, 6, 0|
//   // So its norms are respectively
//   // Infinity   = 12
//   // One        = 13
//   // Frobenius ~= 11.225
//   std::cout << "Dynamic format" << std::endl;
//   M4.uncompress();
//   test6(M4);
//   std::cout << "The expected value is: [12, 13, 11.225]\n" << std::endl;
//   std::cout << "Compressed format" << std::endl;
//   M4.compress();
//   test6(M4);
//   std::cout << "The expected value is: [12, 13, 11.225]" << std::endl;

//   std::cout << "\nBIG NORM TEST" << std::endl;
//   std::cout << "Dynamic format" << std::endl;
//   M5.uncompress();
//   test6(M5);
//   std::cout << "Compressed format" << std::endl;
//   M5.compress();
//   test6(M5);

//   //---------------------------------------------------------------------------
//   // TEST FOR COLUMN-WISE COMPLEX MATRICES
//   //---------------------------------------------------------------------------
//   std::cout << "\n\n-----TESTS FOR COLUMN-WISE COMPLEX MATRICES-----"
//             << std::endl;
//   using namespace std::complex_literals;

//   std::cout << "\nMULTIPLICATION TEST" << std::endl;
//   algebra::Matrix<std::complex<double>, algebra::ColumnWise> Mc1;
//   Mc1(0, 0) = 1.0 + 1i;
//   Mc1(1, 1) = 2. + 2i;
//   // The matrix is
//   // |1+i,  0  |
//   // | 0 , 2+2i|
//   // While the vector is
//   // |1, i|'
//   // So the expected value is
//   // |1+i, -2+2i|'
//   std::cout << "Dynamic format" << std::endl;
//   test5_complex(Mc1);
//   std::cout << "The expected value is: [1+i, -2+2i]\n" << std::endl;
//   std::cout << "Compressed format" << std::endl;
//   Mc1.compress();
//   test5_complex(Mc1);
//   std::cout << "The expected value is: [1+i, -2+2i]" << std::endl;

//   std::cout << "\nNORM TEST" << std::endl;
//   Mc1.uncompress();
//   Mc1(0, 1) = 1i;
//   // The matrix is
//   // |1+i,    i|
//   // | 0 , 2+2i|
//   // So its norms are respectively
//   // Infinity  ~= 2.8284
//   // One       ~= 3.8284
//   // Frobenius ~= 3.3166
//   std::cout << "Dynamic format" << std::endl;
//   test6(Mc1);
//   std::cout << "The expected value is: [2.82843, 3.82843, 3.31662]\n"
//             << std::endl;
//   std::cout << "Compressed format" << std::endl;
//   Mc1.compress();
//   test6(Mc1);
//   std::cout << "The expected value is: [2.82843, 3.82843, 3.31662]"
//             << std::endl;

//   return 0;
// }

void test1(auto &M) {
  std::cout << M << std::endl;

  M.compress();
  std::cout << M << std::endl;

  M.uncompress();
  std::cout << M << std::endl;
}

void test2(const auto &M) {
  std::cout << "(0, 0): " << M(0, 0) << std::endl;
  std::cout << "(1, 1): " << M(1, 1) << std::endl;
}

void test3(auto &M) {
  std::cout << "(0, 0): " << M(0, 0) << std::endl;
  std::cout << "(1, 1): " << M(1, 1) << std::endl;
}

void test4(auto M) {
  M(0, 0) = 2.;
  std::cout << "Element (0, 0) changed" << std::endl;
  M(1, 1) = 1.;
  std::cout << "Element (1, 1) added" << std::endl;

  std::cout << M << std::endl;
}

void test5(const auto &M) {
  using namespace Timings;

  std::vector<double> v(M.cols(), 1);
  for (double k = 1; k < v.size(); ++k) {
    v[k] += k;
  }

  std::cout << "The vector is:" << v << std::endl;

  Chrono time;

  time.start();

  auto res = M * v;

  time.stop();

  std::cout << "Their product is:" << res
            << "It has been evaluated in: " << time.wallTime()
            << " microseconds" << std::endl;
}

void test5_noprint(const auto &M) {
  using namespace Timings;

  std::vector<double> v(M.cols(), 1);
  for (double k = 1; k < v.size(); ++k) {
    v[k] += k;
  }

  Chrono time;

  time.start();

  auto res = M * v;

  time.stop();

  std::cout << "Their product has been evaluated in: " << time.wallTime()
            << " microseconds" << std::endl;
}

void test5_complex(const auto &M) {
  using namespace std::complex_literals;
  using namespace Timings;

  std::vector<std::complex<double>> v(M.cols(), 1);
  for (double k = 1; k < v.size(); ++k) {
    v[k] *= 1i;
  }

  std::cout << "The matrix type is: "
            << abi::__cxa_demangle(typeid(M).name(), 0, 0, 0) << std::endl;

  std::cout << "The vector type is: "
            << abi::__cxa_demangle(typeid(v).name(), 0, 0, 0) << std::endl;
  std::cout << "The vector is:" << v << std::endl;

  Chrono time;

  time.start();

  auto res = M * v;

  time.stop();

  std::cout << "The result type is: "
            << abi::__cxa_demangle(typeid(res).name(), 0, 0, 0) << std::endl;
  std::cout << "Their product is:" << res
            << "It has been evaluated in: " << time.wallTime()
            << " microseconds" << std::endl;
}

void test6(const auto &M) {
  using namespace Timings;
  Chrono time;

  time.start();
  auto inf_norm = M.template norm<algebra::Infinity>();
  time.stop();
  std::cout << "Infinity norm: " << inf_norm
            << " has been evaluated in: " << time.wallTime() << " microseconds"
            << std::endl;

  time.start();
  auto one_norm = M.template norm<algebra::One>();
  time.stop();
  std::cout << "One norm: " << one_norm
            << " has been evaluated in: " << time.wallTime() << " microseconds"
            << std::endl;

  time.start();
  auto frob_norm = M.template norm<algebra::Frobenius>();
  time.stop();
  std::cout << "Frobenius norm: " << frob_norm
            << " has been evaluated in: " << time.wallTime() << " microseconds"
            << std::endl;
}