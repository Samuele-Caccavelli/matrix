#include "Matrix.hpp"

int main(int argc, char **argv) {
  algebra::Matrix<double, RowWise> M(2, 2);

  std::cout << "M(1, 1): " << M(1, 1) << "\n"
            << "M(1, 2): " << M(1, 2) << "\n"
            << "M(2, 1): " << M(2, 1) << "\n"
            << "M(2, 2): " << M(2, 2) << std::endl;

  std::cout << "Now let's change some value" << std::endl;

  M(1, 1) = 2.;
  M(2, 2) = 3.;

  std::cout << "M(1, 1): " << M(1, 1) << "\n"
            << "M(1, 2): " << M(1, 2) << "\n"
            << "M(2, 1): " << M(2, 1) << "\n"
            << "M(2, 2): " << M(2, 2) << std::endl;

  std::cout << "Let's try to access some value out of bound" << std::endl;

  std::cout << "M(1, 3): " << M(3, 1) << std::endl;

  return 0;
}