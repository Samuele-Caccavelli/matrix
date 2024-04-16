#ifndef MATRIX_H
#define MATRIX_H

enum StorageOrder {RowWise = 0, ColumnWise = 1};

namespace algebra {
template <class T, StorageOrder ORDERING> class Matrix {
public:
  using KeyType = std::array<std::size_t, 2>;
  using DataType = T;

private:
  std::map<KeyType, DataType> data;
  bool compressed;

public:
  bool is_compressed();
};

    template<class T, StorageOrder ORDERING>
    bool
    Matrix<T, ORDERING>::is_compressed() {
        return compressed;
    }

} // namespace algebra

#endif /* MATRIX_H */