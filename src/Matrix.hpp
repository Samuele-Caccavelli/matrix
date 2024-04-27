#ifndef MATRIX_H
#define MATRIX_H

//! just to avoid squiggles for size_t
#include <cstddef>

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
// #include <functional>

enum StorageOrder {RowWise = 0, ColumnWise = 1};

namespace algebra {
using KeyType = std::array<std::size_t, 2>;

// return lhs < rhs using RowWise ordering
struct less_col {
  bool operator()(const KeyType lhs, const KeyType rhs) const {
    if (lhs[1] < rhs[1])
      return true;
    else if (lhs[1] > rhs[1])
      return false;
    return lhs[0] < rhs[0];
  }
};

template <class T, StorageOrder ORDERING> class Matrix {
private:
  std::map<KeyType, T> data;
  bool compressed = false;
  std::size_t n_rows;
  std::size_t n_cols;

  // getters
  std::size_t rows() const;
  std::size_t cols() const;

  // setters
  void rows(size_t r);
  void cols(size_t c);

public:
  // true if the matrix is in the compressed format
  bool is_compressed();

  // call operator
  auto operator()(std::size_t r, std::size_t c) const;
  auto &operator()(std::size_t r, std::size_t c);

  // default constructor
  Matrix(std::size_t r = 0, std::size_t c = 0) : n_rows{r}, n_cols{c} {
    std::cout << "I'm a " << n_rows << " x " << n_cols << " matrix"
              << std::endl;
  }

  void upper_bound() const;
  void lower_bound() const;
};

template <class T, StorageOrder ORDERING>
bool Matrix<T, ORDERING>::is_compressed() {
  return compressed;
}

template <class T, StorageOrder ORDERING>
std::size_t Matrix<T, ORDERING>::rows() const {
  return n_rows;
}

template <class T, StorageOrder ORDERING>
std::size_t Matrix<T, ORDERING>::cols() const {
  return n_cols;
}

template <class T, StorageOrder ORDERING>
void Matrix<T, ORDERING>::rows(std::size_t r) {
  n_rows = r;
}

template <class T, StorageOrder ORDERING>
void Matrix<T, ORDERING>::cols(std::size_t c) {
  n_cols = c;
}

template <class T, StorageOrder ORDERING>
auto Matrix<T, ORDERING>::operator()(std::size_t r, std::size_t c) const {

  // std::cout << "Call of the const version" << std::endl;

  if (r > n_rows || c > n_cols) {
    std::cerr << "You are trying to access an element out of bound"
              << std::endl;
    return std::numeric_limits<T>::quiet_NaN();
  }

  KeyType key{r, c};
  // find returns the iterator to the found object if it exists
  auto elem = data.find(key);

  if (elem != data.end())
    return elem->second;
  return 0.;
}

// here to allow assigning new elements even if they are "out of bound", there
// is no check on the dimensions of the matrix, and if not respected, the
// dimensions are redefined
template <class T, StorageOrder ORDERING>
auto &Matrix<T, ORDERING>::operator()(std::size_t r, std::size_t c) {

  // std::cout << "Call of the non-const version" << std::endl;

  if (r > n_rows)
    n_rows = r;

  if (c > n_cols)
    n_cols = c;

  KeyType key{r, c};
  return data[key];
}

template <class T, StorageOrder ORDERING>
void Matrix<T, ORDERING>::upper_bound() const {
  KeyType second_row{2, 1};

  std::cout << "The first row of the matrix is: " << std::endl;

  auto it = data.begin();

  while (it->first < second_row) {
    std::cout << "M(" << it->first[0] << ", " << it->first[1]
              << "): " << it->second << std::endl;
    ++it;
  }
}

template <class T, StorageOrder ORDERING>
void Matrix<T, ORDERING>::lower_bound() const {
  return;
}

// // class specialization for ColumnWise ordering
// template <class T> class Matrix<T, ColumnWise> {
//   private:
//   std::map<KeyType, T, less_col> data;
//   bool compressed = false;
//   std::size_t n_rows;
//   std::size_t n_cols;

//   // getters
//   std::size_t rows() const;
//   std::size_t cols() const;

//   // setters
//   void rows(size_t r);
//   void cols(size_t c);

// public:
//   // true if the matrix is in the compressed format
//   bool is_compressed();

//   // call operator
//   auto operator()(std::size_t r, std::size_t c) const;
//   auto &operator()(std::size_t r, std::size_t c);

//   // default constructor
//   Matrix(std::size_t r = 0, std::size_t c = 0) : n_rows{r}, n_cols{c} {
//     std::cout << "I'm a " << n_rows << " x " << n_cols << " matrix"
//               << std::endl;
//   }

//   void upper_bound() const;
//   void lower_bound() const;
// };

// template <class T>
// void Matrix<T, ColumnWise>::upper_bound() const {
//   KeyType second_col{1, 2};

//   std::cout << "The first column of the matrix is: " << std::endl;

//   auto it = data.begin();

//   while (it->first < second_col)
//   {
//     std::cout << "M(" << it->first[0] << ", " << it->first[1] << "): " <<
//     it->second << std::endl;
//     ++it;
//   }
// }

} // namespace algebra

#endif /* MATRIX_H */