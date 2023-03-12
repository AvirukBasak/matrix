#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef enum {
    // possible exceptions
    ERR_0ROWS,       // matrix can't have 0 rows
    ERR_0COLS,       // matrix can't have 0 columns
    ERR_ROUTB,       // row index out of bounds
    ERR_COUTB,       // column index out of bounds
    ERR_INCMP,       // incompatible matrices
    ERR_NOSQR,       // not a square matrix
    ERR_DETR0,       // during inversion, determinant is 0
    // more general exceptions
    ERR_NULLPTR,     // null pointer exception
} MatrixError;

/**
 * Throw an exception with value = exception code.
 * @error MatrixError
 */
void matrix_errndie(MatrixError ex, const char *msg);

typedef struct {
    /**
     * Rows of this matrix.
     */
    int rows;
    /**
     * Columns of this matrix.
     */
    int cols;
    /**
     * Double pointer to matrix data location.
     */
    double **m;
} Matrix;

/**
 * Create a new matrix object.
 * @param rows Row size of DDA
 * @param cols Column size of DDA
 * @param val Default initial value
 * @return Matrix
 * @error MatrixError Matrix can't have 0 rows - ERR_0ROWS
 * @error MatrixError Matrix can't have 0 columns - ERR_0COLS
 */
Matrix matrix_new(int rows, int cols, double val);

/**
 * Clears the array of the matrix instance on scope exit.
 * @param this A reference to the new matrix
 */
void matrix_free(Matrix *this);

/**
 * Get an element of the matrix from an index.
 * @param m1 The current matrix
 * @param i Row wise position of element
 * @param j Column wise position of element
 * @return double The value at index i, j
 * @error MatrixError Row index out of bounds - ERR_ROUTB
 * @error MatrixError Column index out of bounds - ERR_COUTB
 */
double matrix_get(Matrix m1, int i, int j);

/**
 * Set an element of the matrix to an index.
 * @param m1 The current matrix
 * @param i Row wise position of element
 * @param j Column wise position of element
 * @param val Value to be set
 * @error MatrixError Row index out of bounds - ERR_ROUTB
 * @error MatrixError Column index out of bounds - ERR_COUTB
 */
void matrix_set(Matrix m1, int i, int j, double val);

/*
 * Compares two matrices for equality.
 * @param m1 The current matrix
 * @param m2 The matrix to compare to
 * @return boolean True if equal
 */
bool matrix_equals(Matrix m1, Matrix m2);

/**
 * Add or subtract two compatible matrices.
 * @param m1 The current matrix
 * @param m2 The matrix to add
 * @return Matrix The matrix of sums/diff
 * @error MatrixError If matrices aren't compatible - ERR_INCMP
 */
Matrix matrix_addorsub(Matrix m1, Matrix m2, bool sub);

/**
 * Add two compatible matrices.
 * @param m1 The current matrix
 * @param m2 The matrix to add
 * @return Matrix The matrix of sums
 * @error MatrixError If matrices aren't compatible - ERR_INCMP
 */
Matrix matrix_add(Matrix m1, Matrix m2, bool sub);

/**
 * Subtract 2nd from 1st matrix.
 * @param m1 The current matrix
 * @param m2 The matrix to subtract
 * @return Matrix The matrix of differences
 * @error MatrixError If matrices aren't compatible - ERR_INCMP
 */
Matrix matrix_subtract(Matrix m1, Matrix m2);

/**
 * Multiply a matrix by a scalar.
 * @param m1 The current matrix
 * @param scalar Scalar to multiply by
 * @return matrix The matrix of products
 */
Matrix matrix_scale(Matrix m1, double scalar);

/**
 * Multiply two compatible matrices.
 * @param m1 The current matrix
 * @param m2 The matrix to multiply by
 * @return Matrix The matrix after multiplication
 * @error MatrixError If matrices aren't compatible - ERR_INCMP
 */
Matrix matrix_multiply(Matrix m1, Matrix m2);

/**
 * Calculate matrix to the power of +ve integer.
 * @param m1 The current matrix
 * @param index Power of matrix
 * @return Matrix The resulting matrix
 * @error MatrixError Same as MatrixErrors of matrix::multiply method
 */
Matrix matrix_power(Matrix m1, int index);

/**
 * Excludes a row and a column and generates a sub matrix.
 * Useful for calculating determinants and cofactor matrices.
 * @param m1 The current matrix
 * @param row The row to exclude
 * @param col The column to exclude
 * @return Matrix The sub matrix
 * @error MatrixError Row index out of bounds - ERR_ROUTB
 * @error MatrixError Column index out of bounds - ERR_COUTB
 */
Matrix matrix_mksubmatrix(Matrix m1, int row, int col);

/**
 * Calculate determinant of this matrix.
 * @param m1 The current matrix
 * @return double The determinant
 * @throw MatrixError If matrix isn't a square matrix - ERR_NOSQR
 */
double matrix_determinant(Matrix m1);

/**
 * Calculate the transpose of this matrix.
 * @param m1 The current matrix
 * @return Matrix The transpose
 */
Matrix matrix_transpose(Matrix m1);

/**
 * Calculate the cofactor matrix of this matrix.
 * @param m1 The current matrix
 * @return Matrix The cofactor matrix
 * @error MatrixError If matrix isn't a square matrix - ERR_NOSQR
 */
Matrix matrix_cofactor(Matrix m1);

/**
 * Calculate the adjoint of this matrix.
 * @param m1 The current matrix
 * @return Matrix The adjoint
 * @error MatrixError If matrix isn't a square matrix - ERR_NOSQR
 */
Matrix matrix_adjoint(Matrix m1);

/**
 * Calculate the inverse of this matrix.
 * @param m1 The current matrix
 * @return Matrix The inverse
 * @error MatrixError If matrix isn't a square matrix - ERR_NOSQR
 * @error MatrixError If determinant is 0 - ERR_DETR0
 */
Matrix matrix_inverse(Matrix m1);

/**
 * Display the matrix.
 * @param m1 The current matrix
 */
void matrix_print(Matrix m1);

/**
 * Create a null matrix of given size.
 * @param rows Rows of matrix
 * @param cols Cols of matrix
 * @return Matrix A null matrix
 * @error MatrixError Row index out of bounds - ERR_ROUTB
 * @error MatrixError Column index out of bounds - ERR_COUTB
 */
Matrix matrix_mknull(int rows, int cols);

/**
 * Create a unit matrix of given size.
 * @param n Size of matrix
 * @return Matrix A unit matrix
 * @error MatrixError Row index out of bounds - ERR_ROUTB
 * @error MatrixError Column index out of bounds - ERR_COUTB
 */
Matrix matrix_mkunit(int n);

#endif
