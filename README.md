# MATRIX

The code implement a header only library for sparse matrices that support both dynamic and compressed storage and both row-wise and column-wise storage ordering.

## MAKEFILE

The executable can be made by typing `make`.

In the Makefile there also are three utilities functions:
1. `clean`: remove all the object file.
2. `distclean`: remove all object file and executable file.
3. `clean_matrix`: remove the matrix utilized for the tests.

The executable only needs the lnsp_131.mtx file that is automatically downloaded by the Makefile, so it can simply be called by  `./main`.

## MAIN

The main present lots of test to be run to check the implementation of the class.

As it is, only the tests for row-wise matrices are run to avoid clutter the terminal too much. Another main is already given to run the tests for column-wise matrices. To do that just comment lines 50 to 198 and uncomment lines 200 to 349.

**ATTENTION:** At line 92 and 242 the test that can be run are set up to fail and throw an exception that will stop the program. Make sure they are commented to continue with the tests.

All the tests are run by some helpers functions that whose description is reported here (other then having it near their declaration):
- `void test1(auto &M)`:
```
Create the following matrix by dynamically assigning new values
 |4, 0, 0|
 |0, 7, 5|
 |0, 6, 0|
Then prints it, compress it, print it again, uncompress and print one last time to check the conversions between the two format are ok
```
- `void test2(const auto &M)`:
```
Tries to access a 0 and a non-0 element in a matrix passed through a const reference
```
- `void test3(auto &M)`:
```
Tries to access a 0 and a non-0 element in a matrix passed through a non-const reference
```
- `void test4(auto M)`:
```
Tries to change an element and dynamically add a new one to a matrix passed through a copy
```
- `void test5(const auto &M)`:
```
Perform a matrix x vector operation
```
- `void test5_noprint(const auto &M)`:
```
Perform a matrix x vector operation without cluttering the terminal with the matrix
```
- `void test5_complex(const auto &M)`:
```
Perform a matrix x complex vector operation
```
- `void test6(const auto &M)`:
```
Computes the 3 different norms of a matrix
```

## MATRIX CLASS
Here are the functions that form the public interface of the class:
- `Matrix(std::size_t r = 0, std::size_t c = 0)`
- `bool is_compressed() const`
- 
    ```
    std::size_t rows() const
    std::size_t cols() const
    ```
- `void compress()`
- `void resize(std::size_t r = 0, std::size_t c = 0)`
- `void read_MatrixMarket(const std::string &filename)`
- `template <NormType TYPE> double norm() const`
- 
    ```
    auto operator()(std::size_t r, std::size_t c) const
    auto &operator()(std::size_t r, std::size_t c)
    ```
- 
    ```
    template <class U, StorageOrder ORDER>
    friend std::ostream
    &operator<<(std::ostream &os, const Matrix<U, ORDER> &M)
    ```
- 
    ```
    template <class U, StorageOrder ORDER>
    friend std::vector<U> 
    operator*(const Matrix<U, ORDER> &M, const std::vector<U> &v)
    ```

