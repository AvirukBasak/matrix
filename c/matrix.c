/**
 * @author Aviruk Basak
 * @date March 14th, 2023
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "libmatrix/matrix.h"

void matrix_errndie(MatrixError ex, const char *msg)
{
#ifdef DEBUG
    printf("MatrixError: %s\n", msg);
    ex;
#endif
    exit(100 + ex);
    msg;
}

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

Matrix matrix_from(int rows, int cols, double dda[rows][cols])
{
    Matrix nm = matrix_new(rows, cols, 0);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            nm.m[i][j] = dda[i][j];
    return nm;
}

Matrix matrix_clone(const Matrix m1)
{
    Matrix nm = matrix_new(m1.rows, m1.cols, 0);
    for (int i = 0; i < m1.rows; i++)
        for (int j = 0; j < m1.cols; j++)
            nm.m[i][j] = m1.m[i][j];
    return nm;
}

void matrix_free(Matrix *this)
{
    for (int i = 0; i < this->rows; i++)
        free(this->m[i]);
    free(this->m);
    this->rows = this->cols = 0;
    this->m = NULL;
}

double matrix_get(const Matrix m1, int i, int j)
{
    if (i < 0 || i > m1.rows)
        matrix_errndie(ERR_ROUTB, "row index out of bounds");
    if (j < 0 || j > m1.rows)
        matrix_errndie(ERR_COUTB, "column index out of bounds");
    return m1.m[i][j];
}

void matrix_set(const Matrix m1, int i, int j, double val)
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
bool matrix_equals(const Matrix m1, const Matrix m2)
{
    if (m1.rows != m2.rows || m1.cols != m2.cols)
        return false;
    for (int i = 0; i < m1.rows; i++)
        for (int j = 0; j < m1.cols; j++)
            if (m1.m[i][j] != m2.m[i][j])
                return false;
    return true;
}

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
    matrix_free(&m1);
    matrix_free(&m2);
    return nm;
}

Matrix matrix_add(Matrix m1, Matrix m2)
{
    if (m2.rows != m1.rows || m2.cols != m1.cols)
        matrix_errndie(ERR_INCMP, "incompatible matrices for addition");
    return matrix_addorsub(m1, m2, false);
}

Matrix matrix_subtract(Matrix m1, Matrix m2)
{
    if (m2.rows != m1.rows || m2.cols != m1.cols)
        matrix_errndie(ERR_INCMP, "incompatible matrices for subtraction");
    return matrix_addorsub(m1, m2, true);
}

Matrix matrix_scale(Matrix m1, double scalar)
{
    Matrix nm = matrix_new(m1.rows, m1.cols, 0);
    for (int i = 0; i < m1.rows; i++)
        for (int j = 0; j < m1.cols; j++) {
            double rslt = scalar * m1.m[i][j];
            nm.m[i][j] = rslt;
        }
    matrix_free(&m1);
    return nm;
}

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
    matrix_free(&m1);
    matrix_free(&m2);
    return nm;
}

Matrix matrix_power(Matrix m1, int index)
{
    Matrix nm = matrix_clone(m1);
    for (int i = 0; i < index - 1; i++)
        nm = matrix_multiply(nm, matrix_clone(m1));
    matrix_free(&m1);
    return nm;
}

Matrix matrix_mksubmatrix(const Matrix m1, int row, int col)
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

double matrix_determinant(const Matrix m1)
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

Matrix matrix_transpose(Matrix m1)
{
    Matrix nm = matrix_new(m1.cols, m1.rows, 0);
    for (int i = 0; i < m1.rows; i++)
        for (int j = 0; j < m1.cols; j++)
            nm.m[j][i] = m1.m[i][j];
    matrix_free(&m1);
    return nm;
}

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
    matrix_free(&m1);
    return nm;
}

Matrix matrix_adjoint(Matrix m1)
{
    if (m1.rows != m1.cols)
        matrix_errndie(ERR_NOSQR, "not a square matrix");
    return matrix_transpose(matrix_cofactor(m1));
}

Matrix matrix_inverse(Matrix m1)
{
    if (m1.rows != m1.cols)
        matrix_errndie(ERR_NOSQR, "not a square matrix");
    double determinant = matrix_determinant(m1);
    if (determinant == 0)
        matrix_errndie(ERR_DETR0, "determinant is 0");
    return matrix_scale(matrix_adjoint(m1), 1/determinant);
}

void matrix_print(const Matrix m1)
{
    for (int i = 0; i < m1.rows; i++) {
        for (int j = 0; j < m1.cols; j++)
            printf("%.3lf%s", m1.m[i][j], ((j < m1.cols -1)? ", " : ""));
        printf("\n");
    }
}

Matrix matrix_mknull(int rows, int cols)
{
    if (cols == 0) cols = rows;
    Matrix nm = matrix_new(rows, cols, 0);
    return nm;
}

Matrix matrix_mkunit(int n)
{
    Matrix nm = matrix_new(n, n, 0);
    for (int i = 0; i < nm.rows; i++)
        nm.m[i][i] = 1;
    return nm;
}
