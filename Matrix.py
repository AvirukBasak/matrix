class Matrix:

    FLOAT_PRECISION = "%.2f"

    E_0ROWS = "E_0ROWS"    # matrix can't have 0 rows
    E_0COLS = "E_0COLS"    # matrix can't have 0 columns
    E_ROUTB = "E_ROUTB"    # row index out of bounds
    E_COUTB = "E_COUTB"    # column index out of bounds
    E_INCMP = "E_INCMP"    # incompatible matrices
    E_NOSQR = "E_NOSQR"    # not a square matrix
    E_DETR0 = "E_DETR0"    # during inversion, determinant is 0

    def __init__(self, matrix, rows = 0, cols = 0):
        '''
        Create a new Matrix object
        
        Parameters:
            matrix(int[][]): A DDA to initialize, pass None of unknown
            matrix(float[][]): A DDA to initialize, pass None of unknown
            rows(int): If DDA is unknown, pass no. of rows
            cols(int): If DDA is unknown, pass no. of cols

        Raises:
            Error: matrix can't have 0 rows - E_0ROWS
            Error: matrix can't have 0 columns - E_0COLS
        '''
        self.coderr = True
        if rows != 0 and cols != 0:
            self.rows = rows
            self.cols = cols
            self.matrix = [[0 for x in range(cols)] for y in range(rows)]
            return
        self.rows = len(matrix)
        if self.rows == 1:
            matrix = [ matrix ]
        self.cols = len(matrix[0])
        if self.rows < 1:
            self.raise_exception(self.E_0ROWS, "matrix can't have 0 rows")
        if self.cols < 1:
            self.raise_exception(self.E_0COLS, "matrix can't have 0 rows")
        self.matrix = matrix

    def __len__(self):
        return (self.rows, self.cols)

    def __getitem__(self, i):
        return self.matrix[i]

    def __setitem__(self, i, value):
        self.matrix[i] = value

    def raise_exception(self, err_code, err_msg):
        if self.coderr:
            raise Exception(err_code)
        else:
            raise Exception("Matrix: " + err_msg)

    def O(n, cols = 0):
        '''
        Create a null matrix of given size

        Parameters:
            n(int): size of matrix
            cols(int): cols of null matrix(if rows != cols)

        Returns:
            Matrix: a null matrix

        Raises:
            Error: row index out of bounds - E_ROUTB
            Error: column index out of bounds - E_COUTB
        '''
        if cols == 0:
            cols = n
        om = Matrix(None, n, cols)
        return om

    def I(n):
        '''
        Create a unit matrix of given size

        Parameters:
            n(int): size of matrix

        Returns:
            Matrix: a unit matrix

        Raises:
            Error: row index out of bounds - E_ROUTB
            Error: column index out of bounds - E_COUTB
        '''
        im = Matrix(None, n, n)
        for i in range(n):
            im.set(i + 1, i + 1, 1)
        return im

    def get(self, i, j):
        '''
        Get an element of the matrix from an index
        
        Parameters:
            i(int): row wise position of element(starts from 1)
            j(int): column wise position of element(starts from 1)
        
        Returns:
            any: The element

        Raises:
            Error: row index out of bounds - E_ROUTB
            Error: column index out of bounds - E_COUTB
        '''
        i -= 1
        j -= 1
        if i < 0 or i > self.rows:
            self.raise_exception(self.E_ROUTB, "row index out of bounds")
        if j < 0 or j > self.cols:
            self.raise_exception(self.E_COUTB, "column index out of bounds")
        return self[i][j]

    def set(self, i, j, value):
        '''
        Set an element of the matrix to an index
        
        Parameters:
            i(int): row wise position of element(starts from 1)
            j(int): column wise position of element(starts from 1)

        Raises:
            Error: row index out of bounds - E_ROUTB
            Error: column index out of bounds - E_COUTB
        '''
        i -= 1
        j -= 1
        if i < 0 or i > self.rows:
            self.raise_exception(self.E_ROUTB, "row index out of bounds")
        if j < 0 or j > self.cols:
            self.raise_exception(self.E_COUTB, "column index out of bounds")
        self[i][j] = value

    def equals(self, mtb):
        '''
        Compares two matrices for equality

        Parameters:
            mtb(Matrix): The matrix to compare to

        Returns:
            boolean: True if equal
        '''
        if self.rows != mtb.rows or self.cols != mtb.cols:
             return False
        for i in range(self.rows):
            for j in range(self.cols):
                if self[i][j] != mtb[i][j]:
                    return False
        return True

    def __eq__(self, mtb):
        return self.equals(mtb)

    def __contains__(self, el):
        for row in self.m:
            return el in row
        return False

    def add(self, mtb):
        '''
        Adds two compatible matrices

        Parameters:
            mtb(Matrix): The matrix to add

        Returns:
            Matrix: The matrix of sums

        Raises:
            Error: If matrices aren't compatible - E_INCMP
        '''
        if self.rows != mtb.rows or self.cols != mtb.cols:
            self.raise_exception(self.E_INCMP, "incompatible matrices for addition")
        rslt = Matrix(None, self.rows, self.cols)
        for i in range(self.rows):
            for j in range(self.cols):
                rslt[i][j] = self[i][j] + mtb[i][j]
        return rslt

    def __add__(self, mtb):
        return self.add(mtb)

    def subtract(self, mtb):
        '''
        Subtracts two compatible matrices

        Parameters:
            mtb(Matrix): The matrix to subtract

        Returns:
            Matrix: The matrix of differences

        Raises:
            Error: If matrices aren't compatible - E_INCMP
        '''
        if self.rows != mtb.rows or self.cols != mtb.cols:
            self.raise_exception(self.E_INCMP, "incompatible matrices for subtraction")
        rslt = Matrix(None, self.rows, self.cols)
        for i in range(self.rows):
            for j in range(self.cols):
                rslt[i][j] = self[i][j] - mtb[i][j]
        return rslt

    def __sub__(self, mtb):
        return self.subtract(mtb)

    def scale(self, scalar):
        '''
        Multiplies a matrix by a scalar

        Parameters:
            scalar(int): Scalar to multiply by
            scalar(float): Scalar to multiply by

        Returns:
            Matrix: The matrix of products
        '''
        rslt = Matrix(None, self.rows, self.cols)
        for i in range(self.rows):
            for j in range(self.cols):
                rslt[i][j] = scalar * self[i][j]
        return rslt

    def __pos__(self):
        return self

    def __neg__(self):
        return self.scale(-1)

    def multiply(self, mtb):
        '''
        Multiplies two compatible matrices

        Parameters:
            mtb(Matrix): The matrix to multiply by

        Returns:
            Matrix: The matrix after multiplication

        Raises:
            Error: If matrices aren't compatible - E_INCMP
        '''
        if self.cols != mtb.rows:
            self.raise_exception(self.E_INCMP, "incompatible matrices for multiplication")
        m = self.rows
        n = self.cols # same
        n = mtb.rows  # same
        o = mtb.cols
        rslt = Matrix(None, m, o)
        for i in range(m):
            for j in range(o):
                # zeroing the sum without changing datatypes
                sum = self[0][0] - self[0][0]
                for k in range(n):
                    sum += self[i][k] * mtb[k][j]
                rslt[i][j] = sum
        return rslt

    def __mul__(self, nm):
        if isinstance(nm, Matrix):
            return self.multiply(nm)
        elif isinstance(nm, int) or isinstance(nm, float):
            return self.scale(nm)
        else:
            self.raise_exception(self.E_INCMP, "incompatible matrices for multiplication")

    def __rmul__(self, nm):
        return self * nm

    def power(self, index):
        '''
        Calculate matrix to the power of index

        Parameters:
            index(int): Power of matrix

        Returns:
            Matrix: The resulting matrix

        Raises:
            Error: same as errors of Matrix.multiply method
        '''
        if index > 0: tmp = self
        elif index < 0: tmp = self.inverse()
        else: tmp = Matrix.I( min(self.rows, self.cols) )
        result, index = tmp, abs(index)
        for i in range(index - 1):
            result = result.multiply(tmp)
        return result

    def __pow__(self, num):
        return self.power(num)

    def exclude_row_col(self, row, col):
        '''
        Excludes a row and s column and generates a sub matrix. Useful for calculating determinants and cofactor matrices.

        Parameters:
            row(int): The row to exclude(starts from 1)
            col(int): The column to exclude(starts from 1)

        Returns:
            Matrix: The sub matrix

        Raises:
            Error: row index out of bounds - E_ROUTB
            Error: column index out of bounds - E_COUTB
        '''
        row -= 1
        col -= 1
        if row < 0 or row > self.rows:
            self.raise_exception(self.E_ROUTB, "row index out of bounds")
        if col < 0 or col > self.cols:
            self.raise_exception(self.E_COUTB, "column index out of bounds")
        sub_matrix = Matrix(None, self.rows - 1, self.cols - 1)
        skip_row = False
        for j in range(self.rows - 1):
            skip_col = False
            j_self = j
            if j == row:
                skip_row = True
            if skip_row:
                j_self += 1
            for k in range(self.cols - 1):
                k_self = k
                if k == col:
                    skip_col = True
                if skip_col:
                    k_self += 1
                sub_matrix[j][k] = self[j_self][k_self]
        return sub_matrix

    def determinant(self):
        '''
        Calculate determinant of matrix

        Returns:
            int: The determinant
            float: The determinant

        Raises:
            Error: If matrix isn't a square matrix - E_NOSQR
        '''
        if self.rows != self.cols:
            self.raise_exception(self.E_NOSQR, "not a square matrix")
        n = self.rows
        determinant = self[0][0] - self[0][0]
        if n == 1:
            return self[0][0]
        elif n == 2:
            return self[0][0] * self[1][1] - self[0][1] * self[1][0]
        else:
            for i in range(n):
                coefficient = (-1) ** i * self[0][i]
                sub_matrix = self.exclude_row_col(1, i + 1)
                sub_determinant = abs(sub_matrix)
                term = coefficient * sub_determinant
                determinant += term
        return determinant

    def __abs__(self):
        return self.determinant()

    def transpose(self):
        '''
        Calculate the transpose of this matrix

        Returns:
            Matrix: The transpose
        '''
        rslt = Matrix(None, self.cols, self.rows)
        for i in range(self.rows):
            for j in range(self.cols):
                rslt[j][i] = self[i][j]
        return rslt

    def cofactor(self):
        '''
        Calculate the cofactor matrix of this matrix

        Returns:
            Matrix: The cofactor matrix

        Raises:
            Error: If matrix isn't a square matrix - E_NOSQR
        '''
        if self.rows != self.cols:
            self.raise_exception(self.E_NOSQR, "not a square matrix")
        rslt = Matrix(None, rows = self.rows, cols = self.cols)
        for i in range(self.rows):
            for j in range(self.cols):
                coefficient = (-1) ** (i + 1 + j + 1)
                sub_matrix = self.exclude_row_col(i + 1, j + 1)
                rslt[i][j] = coefficient * abs(sub_matrix)
        return rslt

    def adjoint(self):
        '''
        Calculate the adjoint of this matrix

        Returns:
            Matrix: The adjoint

        Raises:
            Error: If matrix isn't a square matrix - E_NOSQR
        '''
        if self.rows != self.cols:
            self.raise_exception(self.E_NOSQR, "not a square matrix")
        return self.cofactor().transpose()

    def inverse(self):
        '''
        Calculate the inverse of this matrix

        Returns:
            Matrix: The inverse

        Raises:
            Error: If matrix isn't a square matrix - E_NOSQR
            Error: If determinant is 0 - E_DETR0
        '''
        if self.rows != self.cols:
            self.raise_exception(self.E_NOSQR, "not a square matrix")
        det = abs(self)
        if det == 0:
            self.raise_exception(self.E_DETR0, "determinant is 0")
        rslt =  1 / det * self.adjoint()
        return rslt

    def __invert__(self):
        return self ** -1

    def __str__(self):
        str = ''
        for i in range(self.rows):
            for j in range(self.cols):
                str += ( self.FLOAT_PRECISION % self[i][j] ) + " \t"
            str += "\n"
        return str

    def print(self, message=None):
        '''
        Display the matrix
        '''
        if message != None:
            print(message)
        print(String(self))

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
        A = Matrix([
            [1, 2, 2],
            [2, 1, 2],
            [2, 2, 1]
        ])

        # Prove: A^2 - 4A - 5I(3) = O
        print("Test 1:\n%s" % ( A**2 - 4*A - 5*Matrix.I(3) == Matrix.O(3) ))

        # Test 2: Given a matrix B
        B = Matrix([
            [1, 0, 1],
            [3, 4, 5],
            [2, 3, 4]
        ])

        # find inverse of B
        print("\nTest 2:\n", ~B)

if __name__ == '__main__':
    Matrix.test()
