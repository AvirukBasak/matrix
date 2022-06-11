class Matrix:

    FLOAT_PRECISION = "%.2f"

    E_0ROWS = "E_0ROWS"    # matrix can't have 0 rows
    E_0COLS = "E_0COLS"    # matrix can't have 0 columns
    E_ROUTB = "E_ROUTB"    # row index out of bounds
    E_COUTB = "E_COUTB"    # column index out of bounds
    E_INCMP = "E_INCMP"    # incompatible matrices
    E_NOSQR = "E_NOSQR"    # not a square matrix
    E_DETR0 = "E_DETR0"    # during inversion, determinant is 0

    def __init__ (self, matrix, rows = 0, cols = 0):
        '''
        Create a new Matrix object
        
        Parameters:
            matrix (int [][]): A DDA to initialize, pass None of unknown
            matrix (float [][]): A DDA to initialize, pass None of unknown
            rows (int): If DDA is unknown, pass no. of rows
            cols (int): If DDA is unknown, pass no. of cols

        Raises:
            Error: matrix can't have 0 rows - E_0ROWS
            Error: matrix can't have 0 columns - E_0COLS
        '''
        self.coderr = True
        if rows != 0 and cols != 0:
            self.rows = rows
            self.cols = cols
            self.matrix = [[0 for x in range (cols)] for y in range (rows)]
            return
        self.rows = len (matrix)
        if self.rows == 1:
            matrix = [ matrix ]
        self.cols = len (matrix [0])
        if self.rows < 1:
            self.raise_exception (self.E_0ROWS, "matrix can't have 0 rows")
        if self.cols < 1:
            self.raise_exception (self.E_0COLS, "matrix can't have 0 rows")
        self.matrix = matrix

    def raise_exception (self, err_code, err_msg):
        if self.coderr:
            raise Exception (err_code)
        else:
            raise Exception ("Matrix: " + err_msg)

    def O (n, cols = 0):
        '''
        Create a null matrix of given size

        Parameters:
            n (int): size of matrix
            cols (int): cols of null matrix (if rows != cols)

        Returns:
            Matrix: a null matrix

        Raises:
            Error: row index out of bounds - E_ROUTB
            Error: column index out of bounds - E_COUTB
        '''
        if cols == 0:
            cols = n
        om = Matrix (None, n, cols)
        return om

    def I (n):
        '''
        Create a unit matrix of given size

        Parameters:
            n (int): size of matrix

        Returns:
            Matrix: a unit matrix

        Raises:
            Error: row index out of bounds - E_ROUTB
            Error: column index out of bounds - E_COUTB
        '''
        im = Matrix (None, n, n)
        for i in range (n):
            im.set (i + 1, i + 1, 1)
        return im

    def get (self, i, j):
        '''
        Get an element of the matrix from an index
        
        Parameters:
            i (int): row wise position of element (starts from 1)
            j (int): column wise position of element (starts from 1)
        
        Returns:
            any: The element

        Raises:
            Error: row index out of bounds - E_ROUTB
            Error: column index out of bounds - E_COUTB
        '''
        i -= 1
        j -= 1
        if i < 0 or i > self.rows:
            self.raise_exception (self.E_ROUTB, "row index out of bounds")
        if j < 0 or j > self.cols:
            self.raise_exception (self.E_COUTB, "column index out of bounds")
        return self.matrix [i][j]

    def set (self, i, j, value):
        '''
        Set an element of the matrix to an index
        
        Parameters:
            i (int): row wise position of element (starts from 1)
            j (int): column wise position of element (starts from 1)

        Raises:
            Error: row index out of bounds - E_ROUTB
            Error: column index out of bounds - E_COUTB
        '''
        i -= 1
        j -= 1
        if i < 0 or i > self.rows:
            self.raise_exception (self.E_ROUTB, "row index out of bounds")
        if j < 0 or j > self.cols:
            self.raise_exception (self.E_COUTB, "column index out of bounds")
        self.matrix [i][j] = value

    def equals (self, mtb):
        '''
        Compares two matrices for equality

        Parameters:
            mtb (Matrix): The matrix to compare to

        Returns:
            boolean: True if equal
        '''
        if self.rows != mtb.rows or self.cols != mtb.cols:
             return False
        for i in range (self.rows):
            for j in range (self.cols):
                if self.matrix [i][j] != mtb.matrix [i][j]:
                    return False
        return True

    def add (self, mtb):
        '''
        Adds two compatible matrices

        Parameters:
            mtb (Matrix): The matrix to add

        Returns:
            Matrix: The matrix of sums

        Raises:
            Error: If matrices aren't compatible - E_INCMP
        '''
        if self.rows != mtb.rows or self.cols != mtb.cols:
            self.raise_exception (self.E_INCMP, "incompatible matrices for addition")
        rslt = Matrix (None, self.rows, self.cols)
        for i in range (self.rows):
            for j in range (self.cols):
                rslt.matrix [i][j] = self.matrix [i][j] + mtb.matrix [i][j]
        return rslt

    def subtract (self, mtb):
        '''
        Subtracts two compatible matrices

        Parameters:
            mtb (Matrix): The matrix to subtract

        Returns:
            Matrix: The matrix of differences

        Raises:
            Error: If matrices aren't compatible - E_INCMP
        '''
        if self.rows != mtb.rows or self.cols != mtb.cols:
            self.raise_exception (self.E_INCMP, "incompatible matrices for subtraction")
        rslt = Matrix (None, self.rows, self.cols)
        for i in range (self.rows):
            for j in range (self.cols):
                rslt.matrix [i][j] = self.matrix [i][j] - mtb.matrix [i][j]
        return rslt

    def scale (self, scalar):
        '''
        Multiplies a matrix by a scalar

        Parameters:
            scalar (int): Scalar to multiply by
            scalar (float): Scalar to multiply by

        Returns:
            Matrix: The matrix of products
        '''
        rslt = Matrix (None, self.rows, self.cols)
        for i in range (self.rows):
            for j in range (self.cols):
                rslt.matrix [i][j] = scalar * self.matrix [i][j]
        return rslt

    def multiply (self, mtb):
        '''
        Multiplies two compatible matrices

        Parameters:
            mtb (Matrix): The matrix to multiply by

        Returns:
            Matrix: The matrix after multiplication

        Raises:
            Error: If matrices aren't compatible - E_INCMP
        '''
        if self.cols != mtb.rows:
            self.raise_exception (self.E_INCMP, "incompatible matrices for multiplication")
        m = self.rows
        n = self.cols # same
        n = mtb.rows  # same
        o = mtb.cols
        rslt = Matrix (None, m, o)
        for i in range (m):
            for j in range (o):
                # zeroing the sum without changing datatypes
                sum = self.matrix [0][0] - self.matrix [0][0]
                for k in range (n):
                    sum += self.matrix [i][k] * mtb.matrix [k][j]
                rslt.matrix [i][j] = sum
        return rslt

    def power (self, index):
        '''
        Calculate matrix to the power of index

        Parameters:
            index (int): Power of matrix

        Returns:
            Matrix: The resulting matrix

        Raises:
            Error: same as errors of Matrix.multiply method
        '''
        result = self
        for i in range (index - 1):
            result = result.multiply (self)
        return result

    def exclude_row_col (self, row, col):
        '''
        Excludes a row and s column and generates a sub matrix. Useful for calculating determinants and cofactor matrices.

        Parameters:
            row (int): The row to exclude (starts from 1)
            col (int): The column to exclude (starts from 1)

        Returns:
            Matrix: The sub matrix

        Raises:
            Error: row index out of bounds - E_ROUTB
            Error: column index out of bounds - E_COUTB
        '''
        row -= 1
        col -= 1
        if row < 0 or row > self.rows:
            self.raise_exception (self.E_ROUTB, "row index out of bounds")
        if col < 0 or col > self.cols:
            self.raise_exception (self.E_COUTB, "column index out of bounds")
        sub_matrix = Matrix (None, self.rows - 1, self.cols - 1)
        skip_row = False
        for j in range (self.rows - 1):
            skip_col = False
            j_self = j
            if j == row:
                skip_row = True
            if skip_row:
                j_self += 1
            for k in range (self.cols - 1):
                k_self = k
                if k == col:
                    skip_col = True
                if skip_col:
                    k_self += 1
                sub_matrix.matrix [j][k] = self.matrix [j_self][k_self]
        return sub_matrix

    def determinant (self):
        '''
        Calculate determinant of matrix

        Returns:
            int: The determinant
            float: The determinant

        Raises:
            Error: If matrix isn't a square matrix - E_NOSQR
        '''
        if self.rows != self.cols:
            self.raise_exception (self.E_NOSQR, "not a square matrix")
        n = self.rows
        determinant = self.matrix [0][0] - self.matrix [0][0]
        if n == 1:
            return self.matrix [0][0]
        elif n == 2:
            return self.matrix [0][0] * self.matrix [1][1] - self.matrix [0][1] * self.matrix [1][0]
        else:
            for i in range (n):
                coefficient = (-1) ** i * self.matrix [0][i]
                sub_matrix = self.exclude_row_col (1, i + 1)
                sub_determinant = sub_matrix.determinant()
                term = coefficient * sub_determinant
                determinant += term
        return determinant

    def transpose (self):
        '''
        Calculate the transpose of this matrix

        Returns:
            Matrix: The transpose
        '''
        rslt = Matrix (None, self.cols, self.rows)
        for i in range (self.rows):
            for j in range (self.cols):
                rslt.matrix [j][i] = self.matrix [i][j]
        return rslt

    def cofactor (self):
        '''
        Calculate the cofactor matrix of this matrix

        Returns:
            Matrix: The cofactor matrix

        Raises:
            Error: If matrix isn't a square matrix - E_NOSQR
        '''
        if self.rows != self.cols:
            self.raise_exception (self.E_NOSQR, "not a square matrix")
        rslt = Matrix (None, rows = self.rows, cols = self.cols)
        for i in range (self.rows):
            for j in range (self.cols):
                coefficient = (-1) ** (i + 1 + j + 1)
                sub_matrix = self.exclude_row_col (i + 1, j + 1)
                rslt.matrix [i][j] = coefficient * sub_matrix.determinant()
        return rslt

    def adjoint (self):
        '''
        Calculate the adjoint of this matrix

        Returns:
            Matrix: The adjoint

        Raises:
            Error: If matrix isn't a square matrix - E_NOSQR
        '''
        if self.rows != self.cols:
            self.raise_exception (self.E_NOSQR, "not a square matrix")
        return self.cofactor().transpose()

    def inverse (self):
        '''
        Calculate the inverse of this matrix

        Returns:
            Matrix: The inverse

        Raises:
            Error: If matrix isn't a square matrix - E_NOSQR
            Error: If determinant is 0 - E_DETR0
        '''
        if self.rows != self.cols:
            self.raise_exception (self.E_NOSQR, "not a square matrix")
        det = self.determinant()
        if det == 0:
            self.raise_exception (self.E_DETR0, "determinant is 0")
        rslt = self.adjoint().scale (1 / self.determinant())
        return rslt

    def print (self, message = None):
        '''
        Display the matrix
        '''
        if message != None:
            print (message)
        for i in range (self.rows):
            for j in range (self.cols):
                print (self.FLOAT_PRECISION % self.matrix [i][j], end = " \t")
            print()

    def test():
        '''
        Demo function, tests the matrix class
        Expected output should be:
        --------------------------
        Test 1:
        True

        Test 2:
         0.50   1.50  -2.00 	
        -1.00   1.00  -1.00 	
         0.50  -1.50   2.00
        '''
        # Test 1: Given a matrix A
        A = Matrix ([
            [1, 2, 2],
            [2, 1, 2],
            [2, 2, 1]
        ])

        # Prove: A^2 - 4A - 5I(3) = O
        print ("Test 1:\n%s" % A.power(2).subtract(A.scale(4)).subtract(Matrix.I(3).scale(5)).equals(Matrix.O(3)))

        # Test 2: Given a matrix B
        B = Matrix ([
            [1, 0, 1],
            [3, 4, 5],
            [2, 3, 4]
        ])

        # find inverse of B
        B.inverse().print("\nTest 2:")
