public class Matrix
{
    final String FLOAT_PRECISION = "%.2f";

    final String E_0ROWS = "E_0ROWS";    // matrix can't have 0 rows
    final String E_0COLS = "E_0COLS";    // matrix can't have 0 columns
    final String E_ROUTB = "E_ROUTB";    // row index out of bounds
    final String E_COUTB = "E_COUTB";    // column index out of bounds
    final String E_INCMP = "E_INCMP";    // incompatible matrices
    final String E_NOSQR = "E_NOSQR";    // not a square matrix
    final String E_DETR0 = "E_DETR0";    // during inversion, determinant is 0

    public boolean coderr;

    public int rows;
    public int cols;
    public double[][] matrix;

    /**
     * Default constructor not allowed
     */
    private Matrix() { }

    /**
     * Create a new Matrix object
     *
     * @param matrix If DDA is known, pass the DDA
     *
     * @throws Exception matrix can't have 0 rows - E_0ROWS
     * @throws Exception matrix can't have 0 columns - E_0COLS
     */
    public Matrix (double[][] matrix) throws Exception
    {
        this.rows = matrix.length;
        this.cols = matrix[0].length;
        if (this.rows < 1) {
            this.throwException (this.E_0ROWS, "matrix can't have 0 rows");
        }
        if (this.cols < 1) {
            this.throwException (this.E_0COLS, "matrix can't have 0 rows");
        }
        this.matrix = matrix;
        this.coderr = true;
    }

    /**
     * Create a new Matrix object
     *
     * @param rows If DDA is unknown, pass no. of rows
     * @param cols If DDA is unknown, pass no. of cols
     *
     * @throws Exception matrix can't have 0 rows - E_0ROWS
     * @throws Exception matrix can't have 0 columns - E_0COLS
     */
    public Matrix (int rows, int cols) throws Exception
    {
        if (rows < 1) {
            this.throwException (this.E_0ROWS, "matrix can't have 0 rows");
        }
        if (cols < 1) {
            this.throwException (this.E_0COLS, "matrix can't have 0 rows");
        }
        this.rows = rows;
        this.cols = cols;
        this.coderr = true;
        this.matrix = new double [rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                this.matrix [i][j] = 0;
            }
        }
    }

    private void throwException (String errCode, String errMsg) throws Exception
    {
        if (this.coderr) {
            throw new Exception (errCode);
        } else {
            throw new Exception ("Matrix: " + errMsg);
        }
    }

    /**
     * Create a null matrix of given size
     *
     * @param n size of matrix
     *
     * @return Matrix a null matrix
     *
     * @throws Exception row index out of bounds - E_ROUTB
     * @throws Exception column index out of bounds - E_COUTB
     */
    public static Matrix O (int n) throws Exception
    {
        Matrix om = new Matrix (n, n);
        return om;
    }

    /**
     * Create a null matrix of given size
     *
     * @param rows rows of matrix
     * @param cols cols of matrix
     *
     * @return Matrix a null matrix
     *
     * @throws Exception row index out of bounds - E_ROUTB
     * @throws Exception column index out of bounds - E_COUTB
     */
    public static Matrix O (int rows, int cols) throws Exception
    {
        Matrix om = new Matrix (rows, cols);
        return om;
    }

    /**
     * Create a unit matrix of given size
     *
     * @param n size of matrix
     *
     * @return Matrix a unit matrix
     *
     * @throws Exception row index out of bounds - E_ROUTB
     * @throws Exception column index out of bounds - E_COUTB
     */
    public static Matrix I (int n) throws Exception
    {
        Matrix im = new Matrix (n, n);
        for (int i = 0; i < n; i++) {
            im.set (i + 1, i + 1, 1);
        }
        return im;
    }

    /**
     * Get an element of the matrix from an index
     *
     * @param i row wise position of element (starts from 1)
     * @param j column wise position of element (starts from 1)
     *
     * @return double The value at index i, j
     *
     * @throws Exception row index out of bounds - E_ROUTB
     * @throws Exception column index out of bounds - E_COUTB
     */
    public double get (int i, int j) throws Exception
    {
        i--;
        j--;
        if (i < 0 || i > this.rows) {
            this.throwException (this.E_ROUTB, "row index out of bounds");
        }
        if (j < 0 || j > this.cols) {
            this.throwException (this.E_COUTB, "column index out of bounds");
        }
        return this.matrix [i][j];
    }

    /**
     * Set an element of the matrix to an index
     *
     * @param i row wise position of element (starts from 1)
     * @param j column wise position of element (starts from 1)
     *
     * @throws Exception row index out of bounds - E_ROUTB
     * @throws Exception column index out of bounds - E_COUTB
     */
    public void set (int i, int j,  double value) throws Exception
    {
        i--;
        j--;
        if (i < 0 || i > this.rows) {
            this.throwException (this.E_ROUTB, "row index out of bounds");
        }
        if (j < 0 || j > this.cols) {
            this.throwException (this.E_COUTB, "column index out of bounds");
        }
        this.matrix [i][j] = value;
    }

    /*
     * Compares two matrices for equality
     *
     * @param mtb The matrix to compare to
     *
     * @return boolean true if equal
     */
    public boolean equals (Matrix mtb)
    {
        if (this.rows != mtb.rows || this.cols != mtb.cols) {
             return false;
        }
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.cols; j++) {
                if (this.matrix [i][j] != mtb.matrix [i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Add two compatible matrices
     *
     * @param mtb The matrix to add
     *
     * @return Matrix The matrix of sums
     *
     * @throws Exception If matrices aren't compatible - E_INCMP
     */
    public Matrix add (Matrix mtb) throws Exception
    {
        if (this.rows != mtb.rows || this.cols != mtb.cols) {
            this.throwException (this.E_INCMP, "incompatible matrices for addition");
        }
        Matrix rslt = new Matrix (this.rows, this.cols);
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.cols; j++) {
                rslt.matrix [i][j] = this.matrix [i][j] + mtb.matrix [i][j];
            }
        }
        return rslt;
    }

    /**
     * Subtract two compatible matrices
     *
     * @param mtb The matrix to subtract
     *
     * @return Matrix The matrix of differences
     *
     * @throws Exception If matrices aren't compatible - E_INCMP
     */
    public Matrix subtract (Matrix mtb) throws Exception
    {
        if (this.rows != mtb.rows || this.cols != mtb.cols) {
            this.throwException (this.E_INCMP, "incompatible matrices for subtraction");
        }
        Matrix rslt = new Matrix (this.rows, this.cols);
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.cols; j++) {
                rslt.matrix [i][j] = this.matrix [i][j] - mtb.matrix [i][j];
            }
        }
        return rslt;
    }

    /**
     * Multiplies a matrix by a scalar
     *
     * @param scalar Scalar to multiply by
     *
     * @return Matrix The matrix of products
     */
    public Matrix scale (double scalar) throws Exception
    {
        Matrix rslt = new Matrix (this.rows, this.cols);
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.cols; j++) {
                rslt.matrix [i][j] = scalar * this.matrix [i][j];
            }
        }
        return rslt;
    }

    /**
     * Multiplies two compatible matrices
     *
     * @param mtb The matrix to multiply by
     *
     * @return Matrix The matrix after multiplication
     *
     * @throws Exception If matrices aren't compatible - E_INCMP
     */
    public Matrix multiply (Matrix mtb) throws Exception
    {
        if (this.cols != mtb.rows) {
            this.throwException (this.E_INCMP, "incompatible matrices for multiplication");
        }
        int m = this.rows;
        int n = this.cols; // same
        n = mtb.rows;  // same
        int o = mtb.cols;
        Matrix rslt = new Matrix (m, o);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < o; j++) {
                // zeroing the sum without changing datatypes
                double sum = this.matrix [0][0] - this.matrix [0][0];
                for (int k = 0; k < n; k++) {
                    sum += this.matrix [i][k] * mtb.matrix [k][j];
                }
                rslt.matrix [i][j] = sum;
            }
        }
        return rslt;
    }

    /**
     * Calculate matrix to the power of index
     *
     * @param index Power of matrix
     *
     * @return Matrix The resulting matrix
     *
     * @throws Exception same as errors of Matrix.multiply method
     */
    public Matrix power (int index) throws Exception
    {
        Matrix result = this;
        for (int i = 0; i < index - 1; i++) {
            result = result.multiply (this);
        }
        return result;
    }

    /**
     * Excludes a row and s column and generates a sub matrix. Useful for calculating determinants and cofactor matrices.
     *
     * @param row The row to exclude (starts from 1)
     * @param col The column to exclude (starts from 1)
     *
     * @return Matrix The sub matrix
     *
     * @throws Exception row index out of bounds - E_ROUTB
     * @throws Exception column index out of bounds - E_COUTB
     */
    private Matrix excludeRowCol (int row, int col) throws Exception
    {
        row--;
        col--;
        if (row < 0 || row > this.rows) {
            this.throwException (this.E_ROUTB, "row index out of bounds");
        }
        if (col < 0 || col > this.cols) {
            this.throwException (this.E_COUTB, "column index out of bounds");
        }
        Matrix subMatrix = new Matrix (this.rows - 1, this.cols - 1);
        boolean skipRow = false;
        for (int j = 0; j < this.rows - 1; j++) {
            boolean skipCol = false;
            int jSelf = j;
            if (j == row) {
                skipRow = true;
            }
            if (skipRow) {
                jSelf++;
            }
            for (int k = 0; k < this.cols - 1; k++) {
                int kSelf = k;
                if (k == col) {
                    skipCol = true;
                }
                if (skipCol) {
                    kSelf++;
                }
                subMatrix.matrix [j][k] = this.matrix [jSelf][kSelf];
            }
        }
        return subMatrix;
    }

    /**
     * Calculate determinant of matrix
     *
     * @return double The determinant
     *
     * @throw Exception If matrix isn't a square matrix - E_NOSQR
     */
    public double determinant() throws Exception
    {
        if (this.rows != this.cols) {
            this.throwException (this.E_NOSQR, "not a square matrix");
        }
        int n = this.rows;
        double determinant = 0;
        if (n == 1) {
            return this.matrix [0][0];
        } else if (n == 2) {
            return this.matrix [0][0] * this.matrix [1][1] - this.matrix [0][1] * this.matrix [1][0];
        } else {
            for (int i = 0; i < n; i++) {
                double coefficient = Math.pow (-1, i) * this.matrix [0][i];
                Matrix subMatrix = this.excludeRowCol (1, i + 1);
                double subDeterminant = subMatrix.determinant();
                double term = coefficient * subDeterminant;
                determinant += term;
            }
        }
        return determinant;
    }

    /**
     * Calculate the transpose of this matrix
     *
     * @return Matrix The transpose
     */
    public Matrix transpose() throws Exception
    {
        Matrix rslt = new Matrix (this.cols, this.rows);
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.cols; j++) {
                rslt.matrix [j][i] = this.matrix [i][j];
            }
        }
        return rslt;
    }

    /**
     * Calculate the cofactor matrix of this matrix
     *
     * @return Matrix The cofactor matrix
     *
     * @throws Exception If matrix isn't a square matrix - E_NOSQR
     */
    public Matrix cofactor() throws Exception
    {
        if (this.rows != this.cols) {
            this.throwException (this.E_NOSQR, "not a square matrix");
        }
        Matrix rslt = new Matrix (this.rows, this.cols);
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.cols; j++) {
                double coefficient = Math.pow (-1, i + 1 + j + 1);
                Matrix subMatrix = this.excludeRowCol (i + 1, j + 1);
                rslt.matrix [i][j] = coefficient * subMatrix.determinant();
            }
        }
        return rslt;
    }

    /**
     * Calculate the adjoint of this matrix
     *
     * @return Matrix The adjoint
     *
     * @throws Exception If matrix isn't a square matrix - E_NOSQR
     */
    public Matrix adjoint() throws Exception
    {
        if (this.rows != this.cols) {
            this.throwException (this.E_NOSQR, "not a square matrix");
        }
        return this.cofactor().transpose();
    }

    /**
     * Calculate the inverse of this matrix
     *
     * @return Matrix The inverse
     *
     * @throws Exception If matrix isn't a square matrix - E_NOSQR
     * @throws Exception If determinant is 0 - E_DETR0
     */
    public Matrix inverse() throws Exception
    {
        if (this.rows != this.cols) {
            this.throwException (this.E_NOSQR, "not a square matrix");
        }
        double det = this.determinant();
        if (det == 0) {
            this.throwException (this.E_DETR0, "determinant is 0");
        }
        Matrix rslt = this.adjoint().scale (1 / this.determinant());
        return rslt;
    }

    /**
     * Display the matrix
     */
    public void print() throws Exception
    {
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.cols; j++) {
                System.out.printf (this.FLOAT_PRECISION + " \t", this.matrix [i][j]);
            }
            System.out.println();
        }
    }

    /**
     * Display the matrix
     *
     * @param message The message to print
     */
    public void print (String message) throws Exception
    {
        System.out.println (message);
        this.print();
    }

    /**
     * Demo function, tests the matrix class
     * Expected output should be:
     * --------------------------
     * Test 1:
     * true
     *
     * Test 2:
     *  0.50   1.50  -2.00 	
     * -1.00   1.00  -1.00 	
     *  0.50  -1.50   2.00
     */
    public static void test() throws Exception
    {
        // Test 1: Given a matrix A
        Matrix A = new Matrix (new double [][] {
            {1, 2, 2},
            {2, 1, 2},
            {2, 2, 1}
        });

        // Prove: A^2 - 4A - 5I(3) = O
        System.out.println ("Test 1:\n" + A.power(2).subtract(A.scale(4)).subtract(Matrix.I(3).scale(5)).equals(Matrix.O(3)));

        // Test 2: Given a matrix B
        Matrix B = new Matrix (new double [][] {
            {1, 0, 1},
            {3, 4, 5},
            {2, 3, 4}
        });

        // find inverse of B
        B.inverse().print("\nTest 2:");
    }
}
