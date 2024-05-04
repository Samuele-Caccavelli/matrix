#ifndef MATRIX_H
#define MATRIX_H

#include <array>
#include <cmath>
#include <cstddef>
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

  // Call operators
  auto operator()(std::size_t r, std::size_t c) const;
  auto &operator()(std::size_t r, std::size_t c);

  // Insertion operator that handles both storage solutions
  template <class U, StorageOrder ORDER>
  friend std::ostream &operator<<(std::ostream &os, const Matrix<U, ORDER> &M);

  // Multiplication operator that handles Matrix * vector
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
  if (compressed)
    return;

  std::size_t non_zero = data.size();

  row_idx.resize(n_rows + 1, 0);
  col_idx.resize(non_zero, 0);
  data_compressed.resize(non_zero, 0);

  std::size_t idx = 0;

  for (auto iter = data.cbegin(); iter != data.cend(); ++iter) {
    row_idx[iter->first[0] + 1] += 1;
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
  if (!compressed)
    return;

  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = row_idx[i]; j < row_idx[i + 1]; ++j) {
      KeyType key{i, col_idx[j]};
      data.emplace(std::make_pair(key, data_compressed[j]));
    }
  }

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
    // std::cerr
    //     << "The matrix is in a compressed state, it's not possible to resize
    //     it"
    //     << std::endl;
    // return;
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

template <class T, StorageOrder ORDERING>
auto Matrix<T, ORDERING>::operator()(std::size_t r, std::size_t c) const {
  if (r >= n_rows || c >= n_cols) {
    throw std::runtime_error("ERROR: you are trying to access an element out "
                             "of bound in a const function\n");
    // std::cerr << "You are trying to access an element out of bound"
    //           << std::endl;
    // return std::numeric_limits<T>::quiet_NaN();
  }

  if (compressed) {
    for (size_t j = row_idx[r]; j < row_idx[r + 1]; ++j) {
      if (col_idx[j] == c)
        return data_compressed[j];
    }
    // needed for avoiding errors with complex type
    T ret;
    return ret;
  }

  // find returns the iterator to the found object if it exists
  KeyType key{r, c};
  auto elem = data.find(key);

  if (elem != data.end())
    return elem->second;
  // needed for avoiding errors with complex type
  T ret;
  return ret;
}

// here to allow assigning new elements even if they are "out of bound", there
// is no check on the dimensions of the matrix, and if not respected, the
// dimensions are redefined
template <class T, StorageOrder ORDERING>
auto &Matrix<T, ORDERING>::operator()(std::size_t r, std::size_t c) {
  if (compressed) {
    for (size_t j = row_idx[r]; j < row_idx[r + 1]; ++j) {
      if (col_idx[j] == c)
        return data_compressed[j];
    }
    throw std::runtime_error("ERROR: trying to add an element with the matrix "
                             "in compressed format\n");
    // std::cerr << "!!! TRYING TO ADD AN ELEMENT IN A COMPRESSED STATE !!!"
    //           << std::endl;
    // std::cerr << "The operation will be ignored\n" << std::endl;
    // KeyType key{0, 0};
    // return data[key];
  }

  if ((r + 1) > n_rows)
    n_rows = r + 1;

  if ((c + 1) > n_cols)
    n_cols = c + 1;

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

template <class T, StorageOrder ORDERING>
std::ostream &operator<<(std::ostream &os, const Matrix<T, ORDERING> &M) {
  os << "This is a " << M.n_rows << " x " << M.n_cols << " matrix" << std::endl;
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
std::vector<T> operator*(const Matrix<T, ORDERING> &M,
                         const std::vector<T> &v) {

  if (M.n_cols != v.size()) {
    throw std::runtime_error(
        "ERROR: dimensions for the multiplication are not compatible\n");
    // std::cerr << "Dimensions for the multiplication are not compatible"
    //           << std::endl;
    // std::vector<T> res{};
    // return res;
  }

  std::vector<T> res(M.n_rows, 0);

  if (!M.compressed) {
    for (size_t i = 0; i < M.n_rows; ++i) {
      for (size_t j = 0; j < M.n_cols; ++j) {
        res[i] += M(i, j) * v[j];
      }
    }
    return res;
  }

  for (size_t i = 0; i < M.n_rows; ++i) {
    for (size_t j = M.row_idx[i]; j < M.row_idx[i + 1]; ++j) {
      res[i] += M.data_compressed[j] * v[M.col_idx[j]];
    }
  }
  return res;
}

// class specialization for ColumnWise ordering
template <class T> class Matrix<T, ColumnWise> {
private:
  // Dynamic storage technique
  std::map<KeyType, T, less_col> data;

  // Compressed storage technique
  std::vector<std::size_t> row_idx{};
  std::vector<std::size_t> col_idx{};
  std::vector<T> data_compressed{};

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

  // Call operators
  auto operator()(std::size_t r, std::size_t c) const;
  auto &operator()(std::size_t r, std::size_t c);

  // Insertion operator that handles both storage solutions
  template <class U>
  friend std::ostream &operator<<(std::ostream &os,
                                  const Matrix<U, ColumnWise> &M);

  // Multiplication operator that handles Matrix * vector
  template <class U>
  friend std::vector<U> operator*(const Matrix<U, ColumnWise> &M,
                                  const std::vector<U> &v);
};

// template <class T>
// void Matrix<T, ColumnWise>::upper_bound()

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
  if (compressed)
    return;

  std::size_t non_zero = data.size();

  col_idx.resize(n_cols + 1, 0);
  row_idx.resize(non_zero, 0);
  data_compressed.resize(non_zero, 0);

  std::size_t idx = 0;

  for (auto iter = data.cbegin(); iter != data.cend(); ++iter) {
    col_idx[iter->first[1] + 1] += 1;
    row_idx[idx] = iter->first[0];
    data_compressed[idx] = iter->second;
    ++idx;
  }

  for (std::size_t idx = 1; idx < n_rows + 1; ++idx) {
    col_idx[idx] += col_idx[idx - 1];
  }

  data.clear();

  compressed = true;
}

template <class T> void Matrix<T, ColumnWise>::uncompress() {
  if (!compressed)
    return;

  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = col_idx[i]; j < col_idx[i + 1]; ++j) {
      KeyType key{row_idx[j], i};
      data.emplace(std::make_pair(key, data_compressed[j]));
    }
  }

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
    // std::cerr
    //     << "The matrix is in a compressed state, it's not possible to resize
    //     it"
    //     << std::endl;
    // return;
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
auto Matrix<T, ColumnWise>::operator()(std::size_t r, std::size_t c) const {
  if (r >= n_rows || c >= n_cols) {
    throw std::runtime_error("ERROR: you are trying to access an element out "
                             "of bound in a const function\n");
    // std::cerr << "You are trying to access an element out of bound"
    //           << std::endl;
    // return std::numeric_limits<T>::quiet_NaN();
  }

  if (compressed) {
    for (size_t j = col_idx[c]; j < col_idx[c + 1]; ++j) {
      if (row_idx[j] == r)
        return data_compressed[j];
    }
    // needed for avoiding errors with complex type
    T ret;
    return ret;
  }

  // find returns the iterator to the found object if it exists
  KeyType key{r, c};
  auto elem = data.find(key);

  if (elem != data.end())
    return elem->second;
  // needed for avoiding errors with complex type
  T ret;
  return ret;
}

// here to allow assigning new elements even if they are "out of bound", there
// is no check on the dimensions of the matrix, and if not respected, the
// dimensions are redefined
template <class T>
auto &Matrix<T, ColumnWise>::operator()(std::size_t r, std::size_t c) {
  if (compressed) {
    for (size_t j = col_idx[c]; j < col_idx[c + 1]; ++j) {
      if (row_idx[j] == r)
        return data_compressed[j];
    }
    throw std::runtime_error("ERROR: trying to add an element with the matrix "
                             "in compressed format\n");
    // std::cerr << "!!! TRYING TO ADD AN ELEMENT IN A COMPRESSED STATE !!!"
    //           << std::endl;
    // std::cerr << "The operation will be ignored\n" << std::endl;
    // KeyType key{0, 0};
    // return data[key];
  }

  if ((r + 1) > n_rows)
    n_rows = r + 1;

  if ((c + 1) > n_cols)
    n_cols = c + 1;

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

template <class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T, ColumnWise> &M) {
  os << "This is a " << M.n_rows << " x " << M.n_cols << " matrix" << std::endl;
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

template <class T>
std::vector<T> operator*(const Matrix<T, ColumnWise> &M,
                         const std::vector<T> &v) {

  if (M.n_cols != v.size()) {
    throw std::runtime_error(
        "ERROR: dimensions for the multiplication are not compatible\n");
    // std::cerr << "Dimensions for the multiplication are not compatible"
    //           << std::endl;
    // std::vector<T> res{};
    // return res;
  }

  std::vector<T> res(M.n_rows, 0);

  if (!M.compressed) {
    for (size_t i = 0; i < M.n_rows; ++i) {
      for (size_t j = 0; j < M.n_cols; ++j) {
        res[i] += M(i, j) * v[j];
      }
    }
    return res;
  }

  for (size_t i = 0; i < M.n_cols; ++i) {
    for (size_t j = M.col_idx[i]; j < M.col_idx[i + 1]; ++j) {
      res[M.row_idx[j]] += M.data_compressed[j] * v[i];
    }
  }
  return res;
}

} // namespace algebra

#endif /* MATRIX_H */