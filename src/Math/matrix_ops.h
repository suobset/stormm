// -*-c++-*-
#ifndef STORMM_MATRIX_OPS_H
#define STORMM_MATRIX_OPS_H

#include <algorithm>
#include <cmath>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "Accelerator/hybrid.h"
#include "FileManagement/file_enumerators.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "UnitTesting/file_snapshot.h"
#include "matrix.h"

namespace stormm {
namespace math {

using constants::ExceptionResponse;
using card::Hybrid;
using card::HybridTargetLevel;
using diskutil::PrintSituation;
using parse::NumberFormat;
using parse::PolyNumeric;
using parse::realToString;
using testing::writeSnapshot;

/// \brief A twin pointer, useful for performing partial-pivoting Gaussian elimination on
///        row-major matrices.
template <typename T> struct TwinPointer {
  T* ptr_a;  ///< Pointer to the first stretch of data
  T* ptr_b;  ///< Pointer to the second stretch of data
};

/// \brief Enumerate the transpose states of a matrix, for some basic matrix operations
enum class TransposeState {
  AS_IS,     ///< The matrix shall be taken as it is found
  TRANSPOSE  ///< The matrix shall be handled by first taking its transpose (the original matrix
             ///<   will not be disturbed)
};

/// \brief Multiply two matrices, with or without transposition of either, after multiplying each
///        by scalar prefactors and any pre-existing result matrix by its own prefactor.  This
///        routine is designed to emulate BLAS functionality, with a slightly different API.
///
/// Overloaded:
///   - Operate on raw pointers to C-style arrays, with trusted row and column bounds
///   - Take HpcMatrices as inputs, with desired transposes
///
/// \param a        Hereafter "matrix A", the first of the operand matrices
/// \param row_a    Number of rows in matrix A
/// \param col_a    Number of columns in matrix A
/// \param b        Hereafter "matrix B", the second of the operand matrices
/// \param row_b    Number of rows in matrix B
/// \param col_b    Number of columns in matrix B
/// \param c        Hereafter "matrix C", the third of the operand matrices
/// \param scale_a  Scalar prefactor for matrix A
/// \param scale_b  Scalar prefactor for matrix B
/// \param scale_c  Scalar prefactor for matrix C
/// \param x_a      Transpose instruction for matrix A
/// \param x_b      Transpose instruction for matrix B
/// \{
template <typename T>
void matrixMultiply(const T* a, size_t row_a, size_t col_a, const T* b, size_t row_b, size_t col_b,
                    T* c, T scale_a = 1.0, T scale_b = 1.0, T scale_c = 1.0,
                    TransposeState x_a = TransposeState::AS_IS,
                    TransposeState x_b = TransposeState::AS_IS);

template <typename T>
void matrixMultiply(const HpcMatrix<T> &a, const HpcMatrix<T> &b, HpcMatrix<T> *c, T scale_a = 1.0,
                    T scale_b = 1.0, T scale_c = 1.0, TransposeState x_a = TransposeState::AS_IS,
                    TransposeState x_b = TransposeState::AS_IS);
/// \}

/// \brief Multiply a vector by a matrix, to arrive at another vector.  The constants for scaling
///        the matrix, input vector, or output vector (prior to adding them matrix-vector product)
///        are carried over from the more general matrix-matrix multiply routine.  Transposing the
///        input matrix prior to the matrix-vector multiplication is also feasible.
///
/// Overloaded:
///   - Operate on raw pointers to C-style arrays, with trusted row and column bounds
///   - Take an HpcMatrix, with either desired transpose, plus Hybrid object "vectors" as inputs 
///
/// \param a        Hereafter "matrix A", the matrix in b = scale_b * b + scale_a * scale_x * (A x)
/// \param x        The left-hand side vector to use in multiplication, hereafter "vector X"
/// \param b        The right-hand side vector, hereafter, "vector B"
/// \param row_a    The number of rows in matrix A (taken to be the number of elements in vector B,
///                 unless A is transposed in the calculation, in which case it is taken as the
///                 number of elements in vector X))
/// \param col_a    The number of columns in matrix A (taken to be the number of elements in
///                 vector X unless A is transposed, when it would be the number of elements in
///                 vector B)
/// \param scale_a  Scaling factor for matrix A
/// \param scale_x  Scaling factor for vector X
/// \param scale_b  Scaling factor for vector B
/// \param x_a      Indicator of whether to do a virtual transpose of A prior to the computation
/// \{
template <typename T>
void matrixVectorMultiply(const T* a, const T* x, T* b, size_t row_a, size_t col_a,
                          T scale_a = 1.0, T scale_x = 1.0, T scale_b = 1.0,
                          TransposeState x_a = TransposeState::AS_IS);

template <typename T>
void matrixVectorMultiply(const HpcMatrix<T> &a, const Hybrid<T> &x, Hybrid<T> *b, T scale_a = 1.0,
                          T scale_x = 1.0, T scale_b = 1.0,
                          TransposeState x_a = TransposeState::AS_IS);
/// \}

/// \brief Invert a square matrix.  Simple function to avoid the heavy lift and compilation of
///        BLAS calls.  Templated for single- and double-precision functionality.
///
/// Overloaded:
///   - Preserve the original matrix
///   - Destroy the original matrix 
///
/// \param matrix   The matrix to invert, in column-mjaor (contiguous columns) format
/// \param inverse  The inverse matrix result, pre-allocated and returned in column-major format
/// \{
template <typename T> void invertSquareMatrix(const T* matrix, T* inverse, size_t rank);
template <typename T> void invertSquareMatrix(T* matrix, T* inverse, size_t rank);
/// \}

/// \brief Check the dimensions of matrices against the determined rank.  This function serves the
///        jacobiEigensolver routine.
///
/// \param actual_rank  The presumed rank of the matrix
/// \param s_a          Size of the array allocated for matrix A
/// \param s_v          Size of the array allocated for eigenvector matrix V (or some other
///                     matrix)
/// \param policy       Procedure in the event that the matrices are not of the right size
void checkArraySizeToMatrixRank(size_t actual_rank, size_t rank, size_t s_a, size_t s_v,
                                ExceptionResponse policy);
  
/// \brief Compute the eigenvalues and eigenvectors of a real, symmetric matrix.
///
///  Overloaded:
///    - Take templated type pointers for the matrices and vectors, with an explicitly stated rank
///    - Take std::vectors for the matrices and vector, with the rank to be inferred if it is not
///      given explicitly
///    - Take Hybrid objects for the matrices and vector, with the rank to be inferred if it is not
///      explicitly stated
///
/// \param a       The dense matrix of interest
/// \param v       Pre-allocated matrix which will eventually hold the eiegnvectors
/// \param d       Pre-allocated vector which will eventually hold the eigenvalues
/// \param rank    Rank of the matrix
/// \param policy  Procedure in the event that the iterative solver does not coverge
/// \{
template <typename T>
void jacobiEigensolver(T* a, T* v, T* d, size_t rank,
                       ExceptionResponse policy = ExceptionResponse::WARN);

template <typename T>
void jacobiEigensolver(std::vector<T> *a, std::vector<T> *v, std::vector<T> *d, size_t rank = 0,
                       ExceptionResponse policy = ExceptionResponse::WARN);

template <typename T>
void jacobiEigensolver(Hybrid<T> *a, Hybrid<T> *v, Hybrid<T> *d, size_t rank = 0,
                       ExceptionResponse policy = ExceptionResponse::WARN);

template <typename T>
void jacobiEigensolver(HpcMatrix<T> *a, HpcMatrix<T> *v, Hybrid<T> *d,
                       ExceptionResponse policy = ExceptionResponse::WARN);
/// \}

/// \brief Compute the upper-triangular form of a system of equations using the QR decomposition.
///        The matrix shall be presented in column major format

  
/// \brief Compute the box space transformation matrix given a sequence of six real numbers.
///        Templated for single- and double-precision real forms.
///
/// Overloaded:
///   - Take each dimension as a separate value
///   - Take all dimensions as a const std::vector reference
///   - Take all dimensions as a const C-style array pointer
///
/// \param lx      Length along the first principal axis
/// \param ly      Length along the second principal axis
/// \param lz      Length along the third principal axis
/// \param alpha   First box angle, in radians
/// \param beta    Second box angle, in radians
/// \param gamma   Third box angle, in radians
/// \param matrix  Pre-allocated space for the matrix 
/// \param dims    Vector of all dimensions in the order {lx, ly, lz, alpha, beta, gamma}
/// \{
template <typename T>
void computeBoxTransform(T lx, T ly, T lz, T alpha, T beta, T gamma, T* umat, T* invu);

template <typename T> void computeBoxTransform(const std::vector<T> &dims, T* umat, T* invu);

template <typename T> void computeBoxTransform(const T* dims, T* umat, T* invu);
/// \}

/// \brief Pull the box dimensions back out of the inverse transformation matrix.
///
/// \param lx     The first box dimension (returned)
/// \param ly     The second box dimension (returned)
/// \param lz     The third box dimension (returned)
/// \param alpha  The first box angle (returned)
/// \param beta   The second box angle (returned)
/// \param gamma  The third box angle (returned)
/// \param invu   Inverse transformation matrix (takes fractional coordinates into real space)
template <typename T>
void extractBoxDimensions(T *lx, T *ly, T *lz, T *alpha, T *beta, T *gamma, const T* invu);

/// \brief Print a matrix, either to the screen (if no file name is supplied) or to a file.  This
///        makes use of the writeSnapshot() function in 
///
/// Overloaded:
///   - Print a matrix based on a raw pointer to scalar data
///   - Print a matrix based on a std::vector with scalar data
///   - Print a matrix based on a Hybrid object with scalar data
///
/// \param matrix       The matrix to print
/// \param rows         Number of rows in the matrix
/// \param cols         Number of columns in the matrix
/// \param varname      Variable name to give the matrix when printed
/// \param file_name    File name to print the matrix into (stdout if this argument is empty)
/// \param expectation  Indicator of wheter to expect a file of the given name, and what to do
/// \{
template <typename T>
void printMatrix(const T* matrix, const int rows, const int cols, const std::string &varname,
                 const std::string &file_name = std::string(""),
                 const PrintSituation expectation = PrintSituation::APPEND);

template <typename T>
void printMatrix(const std::vector<T> &matrix, const int rows, const int cols,
                 const std::string &varname, const std::string &file_name = std::string(""),
                 const PrintSituation expectation = PrintSituation::APPEND);

template <typename T>
void printMatrix(const Hybrid<T> &matrix, const int rows, const int cols,
                 const std::string &varname, const std::string &file_name = std::string(""),
                 const PrintSituation expectation = PrintSituation::APPEND);

template <typename T>
void printMatrix(const HpcMatrix<T> &matrix, const std::string &varname,
                 const std::string &file_name = std::string(""),
                 const PrintSituation expectation = PrintSituation::APPEND);
/// \}

} // namespace math
} // namespace stormm

#include "matrix_ops.tpp"

#endif
