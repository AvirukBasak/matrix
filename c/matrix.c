#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "matrix.h"

/**
 * Throw an exception with value = exception code.
 * @error MatrixError
 */
void matrix_errndie(MatrixError ex, const char *msg)
{
#ifdef DEBUG
    printf("MatrixError: %s\n", msg);
#endif
    exit(100 + ex);
}

/**
 * Create a new matrix object.
 * @param rows Row size of DDA
 * @param cols Column size of DDA
 * @param val Default initial value
 * @return Matrix
 * @error MatrixError Matrix can't have 0 rows - ERR_0ROWS
 * @error MatrixError Matrix can't have 0 columns - ERR_0COLS
 */
Matrix matrix_new(int rows, int cols, double val)
{
    if (rows < 1)
        matrix_errndie(ERR_0ROWS, "matrix can't have 0 rows");
    if (cols < 1)
        matrix_errndie(ERR_0COLS, "matrix can't have 0 columns");
    Matrix m1 = {
        .rows = rows,
        .cols = cols,
        .m = malloc(rows * sizeof(double*)),
    };
    for (int i = 0; i < m1.rows; i++) {
        m1.m[i] = malloc(cols * sizeof(double));
        for (int j = 0; j < m1.cols; j++) 
            m1.m[i][j] = val;
    }
    return m1;
}

/**
 * Clears the array of the matrix instance on scope exit.
 * @param this A reference to the new matrix
 */
void matrix_free(Matrix *this)
{
    for (int i = 0; i < this->rows; i++)
        free(this->m[i]);
    free(this->m);
    this->rows = this->cols = 0;
    this->m = NULL;
}

/**
 * Get an element of the matrix from an index.
 * @param m1 The current matrix
 * @param i Row wise position of element
 * @param j Column wise position of element
 * @return double The value at index i, j
 * @error MatrixError Row index out of bounds - ERR_ROUTB
 * @error MatrixError Column index out of bounds - ERR_COUTB
 */
double matrix_get(Matrix m1, int i, int j)
{
    if (i < 0 || i > m1.rows)
        matrix_errndie(ERR_ROUTB, "row index out of bounds");
    if (j < 0 || j > m1.rows)
        matrix_errndie(ERR_COUTB, "column index out of bounds");
    return m1.m[i][j];
}

/**
 * Set an element of the matrix to an index.
 * @param m1 The current matrix
 * @param i Row wise position of element
 * @param j Column wise position of element
 * @param val Value to be set
 * @error MatrixError Row index out of bounds - ERR_ROUTB
 * @error MatrixError Column index out of bounds - ERR_COUTB
 */
void matrix_set(Matrix m1, int i, int j, double val)
{
    if (i < 0 || i > m1.rows)
        matrix_errndie(ERR_ROUTB, "row index out of bounds");
    if (j < 0 || j > m1.cols)
        matrix_errndie(ERR_COUTB, "column index out of bounds");
    m1.m[i][j] = val;
}

/*
 * Compares two matrices for equality.
 * @param m1 The current matrix
 * @param m2 The matrix to compare to
 * @return boolean True if equal
 */
bool matrix_equals(Matrix m1, Matrix m2)
{
    if (m1.rows != m2.rows || m1.cols != m2.cols)
        return false;
    for (int i = 0; i < m1.rows; i++)
        for (int j = 0; j < m1.cols; j++)
            if (m1.m[i][j] != m2.m[i][j])
                return false;
    return true;
}

/**
 * Add or subtract two compatible matrices.
 * @param m1 The current matrix
 * @param m2 The matrix to add
 * @return Matrix The matrix of sums/diff
 * @error MatrixError If matrices aren't compatible - ERR_INCMP
 */
Matrix matrix_addorsub(Matrix m1, Matrix m2, bool sub)
{
    if (m2.rows != m1.rows || m2.cols != m1.cols)
        matrix_errndie(ERR_INCMP, "incompatible matrices");
    Matrix nm = matrix_new(m1.rows, m1.cols, 0);
    for (int i = 0; i < m1.rows; i++)
        for (int j = 0; j < m1.cols; j++) {
            double rslt = m1.m[i][j] + (sub ? (-1) * m2.m[i][j] : m2.m[i][j]);
            nm.m[i][j] = rslt;
        }
    return nm;
}

/**
 * Add two compatible matrices.
 * @param m1 The current matrix
 * @param m2 The matrix to add
 * @return Matrix The matrix of sums
 * @error MatrixError If matrices aren't compatible - ERR_INCMP
 */
Matrix matrix_add(Matrix m1, Matrix m2, bool sub)
{
    if (m2.rows != m1.rows || m2.cols != m1.cols)
        matrix_errndie(ERR_INCMP, "incompatible matrices for addition");
    return matrix_addorsub(m1, m2, false);
}

/**
 * Subtract 2nd from 1st matrix.
 * @param m1 The current matrix
 * @param m2 The matrix to subtract
 * @return Matrix The matrix of differences
 * @error MatrixError If matrices aren't compatible - ERR_INCMP
 */
Matrix matrix_subtract(Matrix m1, Matrix m2)
{
    if (m2.rows != m1.rows || m2.cols != m1.cols)
        matrix_errndie(ERR_INCMP, "incompatible matrices for subtraction");
    return matrix_addorsub(m1, m2, true);
}

/**
 * Multiply a matrix by a scalar.
 * @param m1 The current matrix
 * @param scalar Scalar to multiply by
 * @return matrix The matrix of products
 */
Matrix matrix_scale(Matrix m1, double scalar)
{
    Matrix nm = matrix_new(m1.rows, m1.cols, 0);
    for (int i = 0; i < m1.rows; i++)
        for (int j = 0; j < m1.cols; j++) {
            double rslt = scalar * m1.m[i][j];
            nm.m[i][j] = rslt;
        }
    return nm;
}

/**
 * Multiply two compatible matrices.
 * @param m1 The current matrix
 * @param m2 The matrix to multiply by
 * @return Matrix The matrix after multiplication
 * @error MatrixError If matrices aren't compatible - ERR_INCMP
 */
Matrix matrix_multiply(Matrix m1, Matrix m2)
{
    if (m1.cols != m2.rows)
        matrix_errndie(ERR_INCMP, "incompatible matrices for multiplication");
    int m = m1.rows;
    int n = m1.cols;     // same
    n = m2.rows;         // same
    int o = m2.cols;
    Matrix nm = matrix_new(m, o, 0);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < o; j++) {
            double sum = 0;
            for (int k = 0; k < n; k++)
                sum += m1.m[i][k] * m2.m[k][j];
            nm.m[i][j] = sum;
        }
    return nm;
}

/**
 * Calculate matrix to the power of +ve integer.
 * @param m1 The current matrix
 * @param index Power of matrix
 * @return Matrix The resulting matrix
 * @error MatrixError Same as MatrixErrors of matrix::multiply method
 */
Matrix matrix_power(Matrix m1, int index)
{
    Matrix nm = m1;
    for (int i = 0; i < index - 1; i++) {
        Matrix tmp = matrix_multiply(nm, m1);
        // free all temp matrices, but not this
        if (nm.m != m1.m) matrix_free(&nm);
        nm = tmp;
    }
    return nm;
}

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
Matrix matrix_mksubmatrix(Matrix m1, int row, int col)
{
    if (row < 0 || row > m1.rows)
        matrix_errndie(ERR_ROUTB, "row index out of bounds");
    if (col < 0 || col > m1.rows)
        matrix_errndie(ERR_COUTB, "column index out of bounds");
    Matrix subm = matrix_new(m1.rows - 1, m1.cols - 1, 0);
    bool skip_row = false;
    for (int j = 0; j < m1.rows - 1; j++) {
        bool skip_col = false;
        int j_self = j;
        if (j == row)
            skip_row = true;
        if (skip_row)
            j_self++;
        for (int k = 0; k < m1.cols - 1; k++) {
            int k_self = k;
            if (k == col)
                skip_col = true;
            if (skip_col)
                k_self++;
            subm.m[j][k] = m1.m[j_self][k_self];
        }
    }
    return subm;
}

/**
 * Calculate determinant of this matrix.
 * @param m1 The current matrix
 * @return double The determinant
 * @throw MatrixError If matrix isn't a square matrix - ERR_NOSQR
 */
double matrix_determinant(Matrix m1)
{
    if (m1.rows != m1.cols)
        matrix_errndie(ERR_NOSQR, "not a square matrix");
    int n = m1.rows;
    double det = 0;
    if (n == 1)
        return m1.m[0][0];
    else if (n == 2)
        return m1.m[0][0] * m1.m[1][1] - m1.m[0][1] * m1.m[1][0];
    else
        for (int i = 0; i < n; i++) {
            double coeff = pow(-1, i) * m1.m[0][i];
            Matrix subm = matrix_mksubmatrix(m1, 0, i);
            double subdet = matrix_determinant(subm);
            matrix_free(&subm);
            double term = coeff * subdet;
            det += term;
        }
    return det;
}

/**
 * Calculate the transpose of this matrix.
 * @param m1 The current matrix
 * @return Matrix The transpose
 */
Matrix matrix_transpose(Matrix m1)
{
    Matrix nm = matrix_new(m1.cols, m1.rows, 0);
    for (int i = 0; i < m1.rows; i++)
        for (int j = 0; j < m1.cols; j++)
            nm.m[j][i] = m1.m[i][j];
    return nm;
}

/**
 * Calculate the cofactor matrix of this matrix.
 * @param m1 The current matrix
 * @return Matrix The cofactor matrix
 * @error MatrixError If matrix isn't a square matrix - ERR_NOSQR
 */
Matrix matrix_cofactor(Matrix m1)
{
    if (m1.rows != m1.cols)
        matrix_errndie(ERR_NOSQR, "not a square matrix");
    Matrix nm = matrix_new(m1.cols, m1.rows, 0);
    for (int i = 0; i < m1.rows; i++)
        for (int j = 0; j < m1.cols; j++) {
            double coeff = pow(-1, i + 1 + j + 1);
            Matrix subm = matrix_mksubmatrix(m1, i, j);
            nm.m[i][j] = coeff * matrix_determinant(subm);
            matrix_free(&subm);
        }
    return nm;
}

/**
 * Calculate the adjoint of this matrix.
 * @param m1 The current matrix
 * @return Matrix The adjoint
 * @error MatrixError If matrix isn't a square matrix - ERR_NOSQR
 */
Matrix matrix_adjoint(Matrix m1)
{
    if (m1.rows != m1.cols)
        matrix_errndie(ERR_NOSQR, "not a square matrix");
    Matrix cof = matrix_cofactor(m1);
    Matrix adj = matrix_transpose(cof);
    matrix_free(&cof);
    return adj;
}

/**
 * Calculate the inverse of this matrix.
 * @param m1 The current matrix
 * @return Matrix The inverse
 * @error MatrixError If matrix isn't a square matrix - ERR_NOSQR
 * @error MatrixError If determinant is 0 - ERR_DETR0
 */
Matrix matrix_inverse(Matrix m1)
{
    if (m1.rows != m1.cols)
        matrix_errndie(ERR_NOSQR, "not a square matrix");
    double determinant = matrix_determinant(m1);
    if (determinant == 0)
        matrix_errndie(ERR_DETR0, "determinant is 0");
    Matrix adj = matrix_adjoint(m1);
    Matrix inv = matrix_scale(adj, 1/determinant);
    matrix_free(&adj);
    return inv;
}

/**
 * Display the matrix.
 * @param m1 The current matrix
 */
void matrix_print(Matrix m1)
{
    for (int i = 0; i < m1.rows; i++) {
        for (int j = 0; j < m1.cols; j++)
            printf("%.3lf%s", m1.m[i][j], ((j < m1.cols -1)? ", " : ""));
        printf("\n");
    }
}

/**
 * Create a null matrix of given size.
 * @param rows Rows of matrix
 * @param cols Cols of matrix
 * @return Matrix A null matrix
 * @error MatrixError Row index out of bounds - ERR_ROUTB
 * @error MatrixError Column index out of bounds - ERR_COUTB
 */
Matrix matrix_mknull(int rows, int cols)
{
    if (cols == 0) cols = rows;
    Matrix nm = matrix_new(rows, cols, 0);
    return nm;
}

/**
 * Create a unit matrix of given size.
 * @param n Size of matrix
 * @return Matrix A unit matrix
 * @error MatrixError Row index out of bounds - ERR_ROUTB
 * @error MatrixError Column index out of bounds - ERR_COUTB
 */
Matrix matrix_mkunit(int n)
{
    Matrix nm = matrix_new(n, n, 0);
    for (int i = 0; i < nm.rows; i++)
        nm.m[i][i] = 1;
    return nm;
}
