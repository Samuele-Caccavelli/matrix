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
  // format matrix passed through a copy
  // It is in fact the same test as before,
  // but this time adding a new one should fail
  std::cout << "\nTest 5" << std::endl;

  M(0, 0) = 2.;
  std::cout << "Element (0, 0) changed" << std::endl;
  M(1, 1) = 1.;
  std::cout << "Element (1, 1) added" << std::endl;

  std::cout << M << std::endl;
}

void func6(const auto &M) {
  // Perform a dynamic format matrix x vector operation
  using namespace Timings;

  std::cout << "\nTest 6" << std::endl;

  std::vector<double> v(M.cols(), 1);
  for (double k = 0; k < v.size(); ++k) {
    v[k] += k;
  }

  std::cout << "The matrix is:" << M << std::endl;
  std::cout << "The vector is:" << v << std::endl;

  Chrono time;

  time.start();

  auto res = M * v;

  time.stop();

  std::cout << "Their product is:" << res
            << "It has been evaluated in: " << time.wallTime()
            << " microseconds" << std::endl;
}

// TODO delete this since we first compress the matrix in the main, then pass it
void func7(const auto &M) {
  // Perform a compressed format matrix x vector operation
  std::cout << "\nTest 7" << std::endl;

  M.compress();

  std::vector<double> v(M.cols(), 1);
  for (double k = 0; k < v.size(); ++k) {
    v[k] += k;
  }

  std::cout << "The matrix is:" << M << std::endl;
  std::cout << "The vector is:" << v << std::endl;

  auto start_time = std::chrono::high_resolution_clock::now();

  auto res = M * v;

  auto end_time = std::chrono::high_resolution_clock::now();
  auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(
      end_time - start_time);

  std::cout << "Their product is:" << res
            << "It has been evaluated in: " << elapsed_time.count()
            << " microseconds" << std::endl;
}

// Func 8-9 are very similar to the previous function, just without all the
// prints on the terminal to avoid cluttering it too much since they are
// intended to be used a big matrix
void func8(const auto &M) {
  // Perform a dynamic format matrix x vector operation
  std::cout << "\nTest 8" << std::endl;

  std::vector<double> v(M.cols(), 1);
  for (double k = 0; k < v.size(); ++k) {
    v[k] += k;
  }

  auto start_time = std::chrono::high_resolution_clock::now();

  auto res = M * v;

  auto end_time = std::chrono::high_resolution_clock::now();
  auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(
      end_time - start_time);

  std::cout << "Their product has been evaluated in: " << elapsed_time.count()
            << " microseconds" << std::endl;
}

void func9(auto &M) {
  // Perform a compressed format matrix x vector operation
  std::cout << "\nTest 9" << std::endl;

  M.compress();

  std::vector<double> v(M.cols(), 1);
  for (double k = 0; k < v.size(); ++k) {
    v[k] += k;
  }

  auto start_time = std::chrono::high_resolution_clock::now();

  auto res = M * v;

  auto end_time = std::chrono::high_resolution_clock::now();
  auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(
      end_time - start_time);

  std::cout << "Their product has been evaluated in: " << elapsed_time.count()
            << " microseconds" << std::endl;
}

void func10() {
  // Perform a RowWise complex matrix times a complex vector operation
  //|1+i,  0  | * |1| = | 1+ i|
  //| 0 , 2+2i|   |i|   |-2+2i|
  // The code to print the type is taken online
  // It should work for both GCC and Clang compilers
  std::cout << "\nTest 10" << std::endl;

  using namespace std::complex_literals;
  algebra::Matrix<std::complex<double>, algebra::RowWise> Mc;
  Mc(0, 0) = 1.0 + 1i;
  Mc(1, 1) = 2. + 2i;
  std::cout << "The matrix type is: "
            << abi::__cxa_demangle(typeid(Mc).name(), 0, 0, 0) << std::endl;
  std::cout << Mc << std::endl;

  std::vector<std::complex<double>> vec{1, 1i};
  std::cout << "The vector type is: "
            << abi::__cxa_demangle(typeid(vec).name(), 0, 0, 0) << std::endl;
  std::cout << vec << std::endl;

  auto res = Mc * vec;
  std::cout << "The result type is: "
            << abi::__cxa_demangle(typeid(res).name(), 0, 0, 0) << std::endl;
  std::cout << Mc * vec << std::endl;
}

void func11() {
  // Perform a ColumnWise complex matrix times a complex vector operation
  //|1+i,  0  | * |1| = | 1+ i|
  //| 0 , 2+2i|   |i|   |-2+2i|
  // The code to print the type is taken online
  // It should work for both GCC and Clang compilers
  std::cout << "\nTest 11" << std::endl;

  using namespace std::complex_literals;
  algebra::Matrix<std::complex<double>, algebra::ColumnWise> Mc;
  Mc(0, 0) = 1.0 + 1i;
  Mc(1, 1) = 2. + 2i;
  std::cout << "The matrix type is: "
            << abi::__cxa_demangle(typeid(Mc).name(), 0, 0, 0) << std::endl;
  std::cout << Mc << std::endl;

  std::vector<std::complex<double>> vec{1, 1i};
  std::cout << "The vector type is: "
            << abi::__cxa_demangle(typeid(vec).name(), 0, 0, 0) << std::endl;
  std::cout << vec << std::endl;

  auto res = Mc * vec;
  std::cout << "The result type is: "
            << abi::__cxa_demangle(typeid(res).name(), 0, 0, 0) << std::endl;
  std::cout << Mc * vec << std::endl;
}

void func12(const auto &M) {
  // Computes the 3 different norms of a matrix
  std::cout << "\nTest 12" << std::endl;

  auto inf_norm = M.template norm<algebra::Infinity>();
  std::cout << "Infinity norm: " << inf_norm << std::endl;

  auto one_norm = M.template norm<algebra::One>();
  std::cout << "One norm: " << one_norm << std::endl;

  auto frob_norm = M.template norm<algebra::Frobenius>();
  std::cout << "Frobenius norm: " << frob_norm << std::endl;
}

int main(int argc, char **argv) {
  /*
  //---------------------------------------------------------------------------
  //TEST FOR ROW-WISE MATRICES
  //---------------------------------------------------------------------------
  algebra::Matrix<double, algebra::RowWise> M1;

  func1(M1);

  algebra::Matrix<double, algebra::RowWise> M2(2, 2);
  M2(0, 0) = 1.;

  func2(M2);

  M2.compress();
  func3(M2);

  algebra::Matrix<double, algebra::RowWise> M3(1, 1);
  M3(0, 0) = 100.;

  func4(M3);

  // //! ATTENTION: This test is set up to throw an exception
  // //! Comment it away to proceed with the other tests
  // M3.compress();
  // func5(M3);


  algebra::Matrix<double, algebra::RowWise> M;
  M(0, 0) = 4.;
  M(1, 1) = 7.;
  M(1, 2) = 5.;
  M(2, 1) = 6.;

  // The matrix is
  //[4, 0, 0,
  // 0, 7, 5,
  // 0, 6, 0]
  // While the vector is
  //[1, 2, 3]'
  // So the expected value is
  //[4, 29, 12]'
  func6(M);
  M.compress();
  func6(M);
  M.uncompress();

  // The matrix is
  //[4, 0, 0,
  // 0, 7, 5,
  // 0, 6, 0]
  // So its norms should be respectively
  // Infinity   = 12
  // One        = 13
  // Frobenius ~= 11.22
  func12(M);
  M.compress();
  func12(M);
  */

  using namespace std::complex_literals;
  algebra::Matrix<std::complex<double>, algebra::RowWise> Mc;
  Mc(0, 0) = 1.0 + 1i;
  Mc(1, 1) = 2. + 2i;
  Mc(0, 1) = 1i;

  func12(Mc);
  Mc.compress();
  func12(Mc);

  /*
  //---------------------------------------------------------------------------
  //TEST FOR COLUMN-WISE MATRICES
  //---------------------------------------------------------------------------
  algebra::Matrix<double, algebra::ColumnWise> M4;

  func1(M4);

  algebra::Matrix<double, algebra::ColumnWise> M5(2, 2);
  M5(0, 0) = 1.;

  func2(M5);

  M5.compress();
  func3(M5);

  algebra::Matrix<double, algebra::ColumnWise> M6(1, 1);
  M6(0, 0) = 100.;

  func4(M6);

  // //! ATTENTION: This test is set up to throw an exception
  // //! Comment it away to proceed with the other tests
  // M6.compress();
  // func5(M6);

  algebra::Matrix<double, algebra::ColumnWise> M;
  M(0, 0) = 4.;
  M(1, 1) = 7.;
  M(1, 2) = 5.;
  M(2, 1) = 6.;

  // The matrix is
  //[4, 0, 0,
  // 0, 7, 5,
  // 0, 6, 0]
  // While the vector is
  //[1, 2, 3]'
  // So the expected value is
  //[4, 29, 12]'
  func6(M);
  func7(M);
  */

  std::string filename{"lnsp_131.mtx"};
  /*
  //---------------------------------------------------------------------------
  //TEST FOR EFFICIENCY OF * OPERATOR USING A ROW-WISE BIG MATRIX
  //(DYNAMIC VS COMPRESS FORMAT)
  //---------------------------------------------------------------------------
  algebra::Matrix<double, algebra::RowWise> N1;
  N1.read_MatrixMarket(filename);
  func8(N1);
  func9(N1);
  */

  /*
  //---------------------------------------------------------------------------
  //TEST FOR EFFICIENCY OF * OPERATOR USING A COLUMN-WISE BIG MATRIX
  //(DYNAMIC VS COMPRESS FORMAT)
  //---------------------------------------------------------------------------
  algebra::Matrix<double, algebra::ColumnWise> N2;
  N2.read_MatrixMarket(filename);
  func8(N2);
  func9(N2);
  */

  /*
  //---------------------------------------------------------------------------
  // TEST FOR * OPERATOR USING COMPLEX NUMBERS
  //---------------------------------------------------------------------------
  func10();
  func11();
  */

  return 0;
}