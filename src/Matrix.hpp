#ifndef MATRIX_H
#define MATRIX_H

//! just to avoid squiggles for size_t
#include <cstddef>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

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
  // dynamic storage technique
  std::map<KeyType, T> data;

  // compressed storage technique
  std::vector<std::size_t> row_idx{};
  std::vector<std::size_t> col_idx{};
  std::vector<T> data_compressed{};

  bool compressed = false;
  std::size_t n_rows;
  std::size_t n_cols;

  // setters
  void rows(size_t r);
  void cols(size_t c);

public:
  // true if the matrix is in the compressed format
  bool is_compressed();

  // getters
  std::size_t rows() const;
  std::size_t cols() const;

  // call operator
  auto operator()(std::size_t r, std::size_t c) const;
  auto &operator()(std::size_t r, std::size_t c);

  // default constructor
  Matrix(std::size_t r = 0, std::size_t c = 0) : n_rows{r}, n_cols{c} {
    std::cout << "I'm a " << n_rows << " x " << n_cols << " matrix"
              << std::endl;
  }

  // compress
  void compress();

  // uncompress
  void uncompress();

  // print
  template <class U, StorageOrder ORDER>
  friend std::ostream &operator<<(std::ostream &os, const Matrix<U, ORDER> &M);

  // multiplication
  template <class U, StorageOrder ORDER>
  friend std::vector<U> operator*(const Matrix<U, ORDER> &M,
                                  const std::vector<U> &v);

  // resize
  void resize(std::size_t r = 0, std::size_t c = 0);

  // read from a stream
  void read_MatrixMarket(std::string &filename);
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

  // here we add an element to the map, that is empty since the matrix is in its
  // compressed state
  // TODO it will be removed when we try to go back to a dynamic state
  //! Non ho provato a usare una exception perchè non volevo un try ... catch
  //! nel main
  if (compressed) {
    std::cerr << "!!! TRYING TO ADD AN ELEMENT IN A COMPRESSED STATE !!!"
              << std::endl;
    std::cerr << "The operation will be ignored\n" << std::endl;
    KeyType key{0, 0};
    return data[key];
  }

  if (r > n_rows)
    n_rows = r;

  if (c > n_cols)
    n_cols = c;

  KeyType key{r, c};
  return data[key];
}

// template <class T, StorageOrder ORDERING>
// void Matrix<T, ORDERING>::upper_bound() const {
//   KeyType second_row{2, 1};

//   std::cout << "The first row of the matrix is: " << std::endl;

//   auto it = data.begin();

//   while (it->first < second_row) {
//     std::cout << "M(" << it->first[0] << ", " << it->first[1]
//               << "): " << it->second << std::endl;
//     ++it;
//   }
// }

template <class T, StorageOrder ORDERING> void Matrix<T, ORDERING>::compress() {
  if (is_compressed())
    return;

  std::size_t non_zero = data.size();

  row_idx.resize(n_rows + 1, 0);
  col_idx.resize(non_zero, 0);
  data_compressed.resize(non_zero, 0);

  std::size_t idx = 0;

  for (auto iter = data.cbegin(); iter != data.cend(); ++iter) {
    row_idx[iter->first[0]] += 1;
    col_idx[idx] = iter->first[1];
    data_compressed[idx] = iter->second;
    ++idx;
  }

  for (std::size_t idx = 1; idx < n_rows + 1; ++idx) {
    row_idx[idx] += row_idx[idx - 1];
  }

  data.clear();

  compressed = true;
}

template <class T, StorageOrder ORDERING>
void Matrix<T, ORDERING>::uncompress() {
  if (!is_compressed())
    return;

  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = row_idx[i]; j < row_idx[i + 1]; ++j) {
      KeyType key{i + 1, col_idx[j]};
      data.emplace(std::make_pair(key, data_compressed[j]));
    }
  }

  data_compressed.clear();
  row_idx.clear();
  col_idx.clear();

  compressed = false;
}

template <class T, StorageOrder ORDERING>
std::ostream &operator<<(std::ostream &os, const Matrix<T, ORDERING> &M) {
  if (M.compressed) {
    os << "Compressed format" << std::endl;

    os << "Values: ";
    for (const auto &elem : M.data_compressed)
      os << " " << elem << " ";
    os << std::endl;

    os << "Cols: ";
    for (const auto &elem : M.col_idx)
      os << " " << elem << " ";
    os << std::endl;

    os << "Rows: ";
    for (const auto &elem : M.row_idx)
      os << " " << elem << " ";
    os << std::endl;

    return os;
  }

  os << "Dynamic format" << std::endl;
  for (const auto &elem : M.data)
    os << "(" << elem.first[0] << ", " << elem.first[1] << "): " << elem.second
       << std::endl;
  return os;
}

template <class T, StorageOrder ORDERING>
void Matrix<T, ORDERING>::resize(std::size_t r, std::size_t c) {
  if (compressed) {
    std::cerr
        << "The matrix is in a compressed state, it's not possible to resize it"
        << std::endl;
    return;
  }

  if (r > n_rows)
    n_rows = r;

  if (c > n_cols)
    n_cols = c;
}

template <class T, StorageOrder ORDERING>
void Matrix<T, ORDERING>::read_MatrixMarket(std::string &filename) {
  // Try to open the file
  std::ifstream input(filename);
  if (!input.is_open())
    throw std::runtime_error("ERROR: problems while opening the file\n");

  // Read header line
  std::string line;
  std::getline(input, line);
  if (line.substr(0, 21) != "%%MatrixMarket matrix")
    throw std::runtime_error("ERROR: invalid Matrix Market file format\n");

  // Read matrix dimensions
  std::size_t non_zero{0};
  std::getline(input, line);
  std::stringstream ss(line);
  ss >> n_rows >> n_cols >> non_zero;

  // Insert non-zero elemets
  size_t i{0};
  size_t j{0};
  T value{};
  for (size_t k = 0; k < non_zero; ++k) {
    std::getline(input, line);
    std::stringstream ss(line);
    ss >> i >> j >> value;
    KeyType key{i, j};
    data.emplace(std::make_pair(key, value));
  }

  std::cout << "File read correctly" << std::endl;
}

template <class T, StorageOrder ORDERING>
std::vector<T> operator*(const Matrix<T, ORDERING> &M,
                         const std::vector<T> &v) {

  if (M.n_cols != v.size()) {
    std::cerr << "Dimensions for the multiplication are not compatible"
              << std::endl;
    std::vector<T> res{};
    return res;
  }

  std::vector<T> res(M.n_rows, 0);

  if (!M.compressed) {
    for (size_t i = 0; i < M.n_rows; ++i) {
      for (size_t j = 0; j < M.n_cols; ++j) {
        res[i] += M(i + 1, j + 1) * v[j];
      }
    }
    return res;
  }

  for (size_t i = 0; i < M.n_rows; ++i) {
    for (size_t j = M.row_idx[i]; j < M.row_idx[i + 1]; ++j) {
      res[i] += M.data_compressed[j] * v[M.col_idx[j] - 1];
    }
  }
  return res;
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