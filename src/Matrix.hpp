#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace algebra {
using KeyType = std::array<std::size_t, 2>;

enum StorageOrder { RowWise = 0, ColumnWise = 1 };

enum NormType { Infinity = 0, One = 1, Frobenius = 2 };

// Return lhs < rhs using RowWise ordering
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
  // Dynamic storage technique
  std::map<KeyType, T> data;

  // Compressed storage technique
  std::vector<std::size_t> row_idx{};
  std::vector<std::size_t> col_idx{};
  std::vector<T> data_compressed{};

  // Boolean that records the storage technique utilized
  bool compressed = false;

  // Dimensions of the matrix
  std::size_t n_rows;
  std::size_t n_cols;

public:
  // Default constructor
  Matrix(std::size_t r = 0, std::size_t c = 0) : n_rows{r}, n_cols{c} {}

  // True if the matrix is in the compressed format
  bool is_compressed() const;

  // Getters
  std::size_t rows() const;
  std::size_t cols() const;

  // Convert the matrix to CSR (or CSC)
  void compress();

  // Convert the matrix to COOmap
  void uncompress();

  // Resize
  void resize(std::size_t r = 0, std::size_t c = 0);

  // Read from a file in MatrixMarket format
  void read_MatrixMarket(const std::string &filename);

  // Calculate the norm of the matrix
  template <NormType TYPE> double norm() const;

  // Call operators
  auto operator()(std::size_t r, std::size_t c) const;
  auto &operator()(std::size_t r, std::size_t c);

  // Friend operator to print the matrix on screen
  // It handles both storage solutions
  template <class U, StorageOrder ORDER>
  friend std::ostream &operator<<(std::ostream &os, const Matrix<U, ORDER> &M);

  // Multiplication friend operator that handles Matrix * vector
  template <class U, StorageOrder ORDER>
  friend std::vector<U> operator*(const Matrix<U, ORDER> &M,
                                  const std::vector<U> &v);
};

template <class T, StorageOrder ORDERING>
bool Matrix<T, ORDERING>::is_compressed() const {
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

template <class T, StorageOrder ORDERING> void Matrix<T, ORDERING>::compress() {
  // If already compressed, do nothing
  if (compressed)
    return;

  std::size_t non_zero = data.size();

  row_idx.resize(n_rows + 1, 0);
  col_idx.resize(non_zero, 0);
  data_compressed.resize(non_zero, 0);

  std::size_t idx = 0;

  // Add element in the new storage
  for (auto iter = data.cbegin(); iter != data.cend(); ++iter) {
    row_idx[iter->first[0] + 1] += 1;
    col_idx[idx] = iter->first[1];
    data_compressed[idx] = iter->second;
    ++idx;
  }

  // As it was before, row_idx contained the number of non-zero elements on that
  // row
  for (std::size_t idx = 1; idx < n_rows + 1; ++idx) {
    row_idx[idx] += row_idx[idx - 1];
  }

  // Clear the old storage
  data.clear();

  compressed = true;
}

template <class T, StorageOrder ORDERING>
void Matrix<T, ORDERING>::uncompress() {
  // If already uncompressed, do nothing
  if (!compressed)
    return;

  // Add element in the new storage
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = row_idx[i]; j < row_idx[i + 1]; ++j) {
      KeyType key{i, col_idx[j]};
      data.emplace(std::make_pair(key, data_compressed[j]));
    }
  }

  // Clear the old storage
  data_compressed.clear();
  row_idx.clear();
  col_idx.clear();

  compressed = false;
}

template <class T, StorageOrder ORDERING>
void Matrix<T, ORDERING>::resize(std::size_t r, std::size_t c) {
  if (compressed) {
    throw std::runtime_error("ERROR: the matrix is in a compressed state, it's "
                             "not possible to resize it\n");
  }

  if (r > n_rows)
    n_rows = r;

  if (c > n_cols)
    n_cols = c;
}

template <class T, StorageOrder ORDERING>
void Matrix<T, ORDERING>::read_MatrixMarket(const std::string &filename) {
  if (!data.empty() || !data_compressed.empty())
    throw std::runtime_error(
        "ERROR: you can only read from a file in a matrix without any data, "
        "otherwise old data could be overwritten\n");
  if (compressed)
    throw std::runtime_error("ERROR: you can only read from a file in a "
                             "uncompressed format matrix\n");

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

  std::cout << "I'm reading a " << n_rows << " x " << n_cols << " matrix"
            << std::endl;

  // Insert non-zero elemets
  size_t i{0};
  size_t j{0};
  T value{};
  for (size_t k = 0; k < non_zero; ++k) {
    std::getline(input, line);
    std::stringstream ss(line);
    ss >> i >> j >> value;
    KeyType key{i - 1, j - 1};
    data.emplace(std::make_pair(key, value));
  }

  std::cout << "File read correctly" << std::endl;
}

template <class T, algebra::StorageOrder ORDERING>
template <algebra::NormType TYPE>
double algebra::Matrix<T, ORDERING>::norm() const {
  // The use of std::max_element instead of std::max is preferred for
  // consistency since std::max cannot be used for the One case

  // Infinity norm case
  if constexpr (TYPE == Infinity) {
    // In both cases each element of par_sum will have the absolute sum of the
    // elements on the corresponding row

    // -----------------
    // COMPRESSED FORMAT
    // -----------------
    if (compressed) {
      std::vector<double> par_sum(n_rows, 0);
      for (size_t i = 0; i < n_rows; ++i) {
        for (size_t j = row_idx[i]; j < row_idx[i + 1]; ++j) {
          par_sum[i] += std::abs(data_compressed[j]);
        }
      }
      auto max_iter = std::max_element(par_sum.begin(), par_sum.end());
      return *max_iter;
    }
    // -----------------
    // DYNAMIC FORMAT
    // -----------------

    // For dynamic format we have to access each element
    std::vector<double> par_sum(n_rows, 0);
    for (size_t i = 0; i < n_rows; ++i) {
      for (size_t j = 0; j < n_cols; ++j) {
        par_sum[i] += std::abs(this->operator()(i, j));
      }
    }
    auto max_iter = std::max_element(par_sum.begin(), par_sum.end());
    return *max_iter;
  }

  // One norm case
  if constexpr (TYPE == One) {
    // In both cases each element of par_sum will have the absolute sum of the
    // elements on the corresponding column

    // -----------------
    // COMPRESSED FORMAT
    // -----------------
    if (compressed) {
      std::vector<double> par_sum(n_cols, 0);
      for (size_t j = 0; j < col_idx.size(); ++j) {
        par_sum[col_idx[j]] += std::abs(data_compressed[j]);
      }
      auto max_iter = std::max_element(par_sum.begin(), par_sum.end());
      return *max_iter;
    }
    // -----------------
    // DYNAMIC FORMAT
    // -----------------

    // For dynamic format we have to access each element
    std::vector<double> par_sum(n_cols, 0);
    for (size_t i = 0; i < n_rows; ++i) {
      for (size_t j = 0; j < n_cols; ++j) {
        par_sum[j] += std::abs(this->operator()(i, j));
      }
    }
    auto max_iter = std::max_element(par_sum.begin(), par_sum.end());
    return *max_iter;
  }

  // Frobenius norm case
  if constexpr (TYPE == Frobenius) {
    // -----------------
    // COMPRESSED FORMAT
    // -----------------
    if (compressed) {
      double par_sum{0};
      for (const auto &elem : data_compressed) {
        par_sum += std::norm(elem);
      }
      return std::sqrt(par_sum);
    }
    // -----------------
    // DYNAMIC FORMAT
    // -----------------

    // For dynamic format we still can access only non-zero elements
    double par_sum{0};
    for (const auto &elem : data) {
      par_sum += std::norm(elem.second);
    }
    return std::sqrt(par_sum);
  }
}

template <class T, StorageOrder ORDERING>
auto Matrix<T, ORDERING>::operator()(std::size_t r, std::size_t c) const {
  if (r >= n_rows || c >= n_cols) {
    throw std::runtime_error("ERROR: you are trying to access an element out "
                             "of bound in a const function\n");
  }

  // -----------------
  // COMPRESSED FORMAT
  // -----------------
  if (compressed) {
    for (size_t i = row_idx[r]; i < row_idx[r + 1]; ++i) {
      if (col_idx[i] == c)
        return data_compressed[i];
    }
    // If the element is not found
    T ret{};
    return ret;
  }
  // -----------------
  // DYNAMIC FORMAT
  // -----------------
  KeyType key{r, c};
  auto elem = data.find(key);

  if (elem != data.end())
    return elem->second;

  // If the element is not found
  T ret{};
  return ret;
}

template <class T, StorageOrder ORDERING>
auto &Matrix<T, ORDERING>::operator()(std::size_t r, std::size_t c) {
  // To allow dynamic construction, check on matrix dimensions in only done for
  // compressed format

  // -----------------
  // COMPRESSED FORMAT
  // -----------------
  if (compressed) {
    // Check on the dimensions
    if (r >= n_cols || c >= n_cols)
      throw std::runtime_error(
          "ERROR: trying to add an element with the matrix "
          "in compressed format\n");

    // Check if the dimensions are correct but the element is not already
    // present
    for (size_t i = row_idx[r]; i < row_idx[r + 1]; ++i) {
      if (col_idx[i] == c)
        return data_compressed[i];
    }
    throw std::runtime_error("ERROR: trying to add an element with the matrix "
                             "in compressed format\n");
  }
  // -----------------
  // DYNAMIC FORMAT
  // -----------------
  if (r >= n_rows)
    n_rows = r + 1;

  if (c >= n_cols)
    n_cols = c + 1;

  KeyType key{r, c};
  return data[key];
}

template <class T, StorageOrder ORDERING>
std::ostream &operator<<(std::ostream &os, const Matrix<T, ORDERING> &M) {
  os << "This is a " << M.n_rows << " x " << M.n_cols << " matrix" << std::endl;

  // -----------------
  // COMPRESSED FORMAT
  // -----------------
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
  // -----------------
  // DYNAMIC FORMAT
  // -----------------
  os << "Dynamic format" << std::endl;
  for (const auto &elem : M.data)
    os << "(" << elem.first[0] << ", " << elem.first[1] << "): " << elem.second
       << std::endl;

  return os;
}

template <class T, StorageOrder ORDERING>
std::vector<T> operator*(const Matrix<T, ORDERING> &M,
                         const std::vector<T> &v) {
  // Check for the dimensions
  if (M.n_cols != v.size()) {
    throw std::runtime_error(
        "ERROR: dimensions for the multiplication are not compatible\n");
  }

  std::vector<T> res(M.n_rows, 0);

  // -----------------
  // COMPRESSED FORMAT
  // -----------------
  if (M.compressed) {
    for (size_t i = 0; i < M.n_rows; ++i) {
      for (size_t j = M.row_idx[i]; j < M.row_idx[i + 1]; ++j) {
        res[i] += M.data_compressed[j] * v[M.col_idx[j]];
      }
    }
    return res;
  }
  // -----------------
  // DYNAMIC FORMAT
  // -----------------
  for (size_t i = 0; i < M.n_rows; ++i) {
    for (size_t j = 0; j < M.n_cols; ++j) {
      res[i] += M(i, j) * v[j];
    }
  }
  return res;
}

//---------------------------------------------------------------------------
// COLUMNWISE ORDERING SPECIALIZATION
//---------------------------------------------------------------------------
template <class T> class Matrix<T, ColumnWise> {
private:
  // Dynamic storage technique
  std::map<KeyType, T, less_col> data;

  // Compressed storage technique
  std::vector<std::size_t> row_idx{};
  std::vector<std::size_t> col_idx{};
  std::vector<T> data_compressed{};

  // Boolean that records the storage technique utilized
  bool compressed = false;

  // Dimensions of the matrix
  std::size_t n_rows;
  std::size_t n_cols;

public:
  // Default constructor
  Matrix(std::size_t r = 0, std::size_t c = 0) : n_rows{r}, n_cols{c} {}

  // True if the matrix is in the compressed format
  bool is_compressed() const;

  // Getters
  std::size_t rows() const;
  std::size_t cols() const;

  // Convert the matrix to CSR (or CSC)
  void compress();

  // Convert the matrix to COOmap
  void uncompress();

  // Resize
  void resize(std::size_t r = 0, std::size_t c = 0);

  // Read from a file in MatrixMarket format
  void read_MatrixMarket(const std::string &filename);

  // Calculate the norm of the matrix
  template <NormType TYPE> double norm() const;

  // Call operators
  auto operator()(std::size_t r, std::size_t c) const;
  auto &operator()(std::size_t r, std::size_t c);

  // Friend operator to print the matrix on screen
  // It handles both storage solutions
  template <class U>
  friend std::ostream &operator<<(std::ostream &os,
                                  const Matrix<U, ColumnWise> &M);

  // Multiplication friend operator that handles Matrix * vector
  template <class U>
  friend std::vector<U> operator*(const Matrix<U, ColumnWise> &M,
                                  const std::vector<U> &v);
};

template <class T> bool Matrix<T, ColumnWise>::is_compressed() const {
  return compressed;
}

template <class T> std::size_t Matrix<T, ColumnWise>::rows() const {
  return n_rows;
}

template <class T> std::size_t Matrix<T, ColumnWise>::cols() const {
  return n_cols;
}

template <class T> void Matrix<T, ColumnWise>::compress() {
  // If already compressed, do nothing
  if (compressed)
    return;

  std::size_t non_zero = data.size();

  col_idx.resize(n_cols + 1, 0);
  row_idx.resize(non_zero, 0);
  data_compressed.resize(non_zero, 0);

  std::size_t idx = 0;

  // Add element in the new storage
  for (auto iter = data.cbegin(); iter != data.cend(); ++iter) {
    col_idx[iter->first[1] + 1] += 1;
    row_idx[idx] = iter->first[0];
    data_compressed[idx] = iter->second;
    ++idx;
  }

  // As it was before, col_idx contained the number of non-zero elements on that
  // column
  for (std::size_t idx = 1; idx < n_rows + 1; ++idx) {
    col_idx[idx] += col_idx[idx - 1];
  }

  // Clear the old storage
  data.clear();

  compressed = true;
}

template <class T> void Matrix<T, ColumnWise>::uncompress() {
  // If already uncompressed, do nothing
  if (!compressed)
    return;

  // Add element in the new storage
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = col_idx[i]; j < col_idx[i + 1]; ++j) {
      KeyType key{row_idx[j], i};
      data.emplace(std::make_pair(key, data_compressed[j]));
    }
  }

  // Clear the old storage
  data_compressed.clear();
  row_idx.clear();
  col_idx.clear();

  compressed = false;
}

template <class T>
void Matrix<T, ColumnWise>::resize(std::size_t r, std::size_t c) {
  if (compressed) {
    throw std::runtime_error("ERROR: the matrix is in a compressed state, it's "
                             "not possible to resize it\n");
  }

  if (r > n_rows)
    n_rows = r;

  if (c > n_cols)
    n_cols = c;
}

template <class T>
void Matrix<T, ColumnWise>::read_MatrixMarket(const std::string &filename) {
  if (!data.empty() || !data_compressed.empty())
    throw std::runtime_error(
        "ERROR: you can only read from a file in a matrix without any data, "
        "otherwise old data could be overwritten\n");
  if (compressed)
    throw std::runtime_error("ERROR: you can only read from a file in a "
                             "uncompressed format matrix\n");

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

  std::cout << "I'm reading a " << n_rows << " x " << n_cols << " matrix"
            << std::endl;

  // Insert non-zero elemets
  size_t i{0};
  size_t j{0};
  T value{};
  for (size_t k = 0; k < non_zero; ++k) {
    std::getline(input, line);
    std::stringstream ss(line);
    ss >> i >> j >> value;
    KeyType key{i - 1, j - 1};
    data.emplace(std::make_pair(key, value));
  }

  std::cout << "File read correctly" << std::endl;
}

template <class T>
template <algebra::NormType TYPE>
double algebra::Matrix<T, ColumnWise>::norm() const {
  // The use of std::max_element instead of std::max is preferred for
  // consistency since std::max cannot be used for the Infinity case

  // Infinity norm case
  if constexpr (TYPE == Infinity) {
    // In both cases each element of par_sum will have the absolute sum of the
    // elements on the corresponding row

    // -----------------
    // COMPRESSED FORMAT
    // -----------------
    if (compressed) {
      std::vector<double> par_sum(n_rows, 0);
      for (size_t i = 0; i < row_idx.size(); ++i) {
        par_sum[row_idx[i]] += std::abs(data_compressed[i]);
      }
      auto max_iter = std::max_element(par_sum.begin(), par_sum.end());
      return *max_iter;
    }
    // -----------------
    // DYNAMIC FORMAT
    // -----------------

    // For dynamic format we have to access each element
    std::vector<double> par_sum(n_rows, 0);
    for (size_t i = 0; i < n_rows; ++i) {
      for (size_t j = 0; j < n_cols; ++j) {
        par_sum[i] += std::abs(this->operator()(i, j));
      }
    }
    auto max_iter = std::max_element(par_sum.begin(), par_sum.end());
    return *max_iter;
  }

  // One norm case
  if constexpr (TYPE == One) {
    // In both cases each element of par_sum will have the absolute sum of the
    // elements on the corresponding column

    // -----------------
    // COMPRESSED FORMAT
    // -----------------
    if (compressed) {
      std::vector<double> par_sum(n_cols, 0);
      for (size_t j = 0; j < n_cols; ++j) {
        for (size_t i = col_idx[j]; i < col_idx[j + 1]; ++i) {
          par_sum[j] += std::abs(data_compressed[i]);
        }
      }
      auto max_iter = std::max_element(par_sum.begin(), par_sum.end());
      return *max_iter;
    }
    // -----------------
    // DYNAMIC FORMAT
    // -----------------

    // For dynamic format we have to access each element
    std::vector<double> par_sum(n_cols, 0);
    for (size_t i = 0; i < n_rows; ++i) {
      for (size_t j = 0; j < n_cols; ++j) {
        par_sum[j] += std::abs(this->operator()(i, j));
      }
    }
    auto max_iter = std::max_element(par_sum.begin(), par_sum.end());
    return *max_iter;
  }

  // Frobenius norm case
  if constexpr (TYPE == Frobenius) {
    // -----------------
    // COMPRESSED FORMAT
    // -----------------
    if (compressed) {
      double par_sum{0};
      for (const auto &elem : data_compressed) {
        par_sum += std::norm(elem);
      }
      return std::sqrt(par_sum);
    }
    // -----------------
    // DYNAMIC FORMAT
    // -----------------

    // For dynamic format we still can access only non-zero elements
    double par_sum{0};
    for (const auto &elem : data) {
      par_sum += std::norm(elem.second);
    }
    return std::sqrt(par_sum);
  }
}

template <class T>
auto Matrix<T, ColumnWise>::operator()(std::size_t r, std::size_t c) const {
  if (r >= n_rows || c >= n_cols) {
    throw std::runtime_error("ERROR: you are trying to access an element out "
                             "of bound in a const function\n");
  }

  // -----------------
  // COMPRESSED FORMAT
  // -----------------
  if (compressed) {
    for (size_t j = col_idx[c]; j < col_idx[c + 1]; ++j) {
      if (row_idx[j] == r)
        return data_compressed[j];
    }
    // If the element is not found
    T ret{};
    return ret;
  }
  // -----------------
  // DYNAMIC FORMAT
  // -----------------
  KeyType key{r, c};
  auto elem = data.find(key);

  if (elem != data.end())
    return elem->second;

  // If the element is not found
  T ret{};
  return ret;
}

template <class T>
auto &Matrix<T, ColumnWise>::operator()(std::size_t r, std::size_t c) {
  // To allow dynamic construction, check on matrix dimensions in only done for
  // compressed format

  // -----------------
  // COMPRESSED FORMAT
  // -----------------
  if (compressed) {
    // Check on the dimensions
    if (r >= n_cols || c >= n_cols)
      throw std::runtime_error(
          "ERROR: trying to add an element with the matrix "
          "in compressed format\n");

    // Check if the dimensions are correct but the element is not already
    // present
    for (size_t j = col_idx[c]; j < col_idx[c + 1]; ++j) {
      if (row_idx[j] == r)
        return data_compressed[j];
    }
    throw std::runtime_error("ERROR: trying to add an element with the matrix "
                             "in compressed format\n");
  }
  // -----------------
  // DYNAMIC FORMAT
  // -----------------
  if (r >= n_rows)
    n_rows = r + 1;

  if (c >= n_cols)
    n_cols = c + 1;

  KeyType key{r, c};
  return data[key];
}

template <class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T, ColumnWise> &M) {
  os << "This is a " << M.n_rows << " x " << M.n_cols << " matrix" << std::endl;

  // -----------------
  // COMPRESSED FORMAT
  // -----------------
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
  // -----------------
  // DYNAMIC FORMAT
  // -----------------
  os << "Dynamic format" << std::endl;
  for (const auto &elem : M.data)
    os << "(" << elem.first[0] << ", " << elem.first[1] << "): " << elem.second
       << std::endl;

  return os;
}

template <class T>
std::vector<T> operator*(const Matrix<T, ColumnWise> &M,
                         const std::vector<T> &v) {
  // Check for the dimensions
  if (M.n_cols != v.size()) {
    throw std::runtime_error(
        "ERROR: dimensions for the multiplication are not compatible\n");
  }

  std::vector<T> res(M.n_rows, 0);

  // -----------------
  // COMPRESSED FORMAT
  // -----------------
  if (M.compressed) {
    for (size_t j = 0; j < M.n_cols; ++j) {
      for (size_t i = M.col_idx[j]; i < M.col_idx[j + 1]; ++i) {
        res[M.row_idx[i]] += M.data_compressed[i] * v[j];
      }
    }
    return res;
  }
  // -----------------
  // DYNAMIC FORMAT
  // -----------------
  for (size_t i = 0; i < M.n_rows; ++i) {
    for (size_t j = 0; j < M.n_cols; ++j) {
      res[i] += M(i, j) * v[j];
    }
  }
  return res;
}

} // namespace algebra

#endif /* MATRIX_H */