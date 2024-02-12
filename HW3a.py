def is_symmetric(matrix):
    """
    Checks if a matrix is symmetric.

    Parameters:
        matrix (list of lists): The matrix to be checked.

    Returns:
        bool: True if the matrix is symmetric, False otherwise.
    """
    n = len(matrix)
    for i in range(n):
        for j in range(i+1, n):
            if matrix[i][j] != matrix[j][i]:
                return False
    return True

def cholesky_decomposition(matrix):
    """
    Computes the Cholesky decomposition of a symmetric positive definite matrix.

    Parameters:
        matrix (list of lists): The matrix to be decomposed.

    Returns:
        list of lists: The lower triangular Cholesky factor L.
    """
    n = len(matrix)
    L = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i+1):
            temp_sum = sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                L[i][j] = (matrix[i][i] - temp_sum)**0.5
            else:
                L[i][j] = (1.0 / L[j][j] * (matrix[i][j] - temp_sum))

    return L

def forward_substitution(L, b):
    """
    Solves a lower triangular linear system Lx = b using forward substitution.

    Parameters:
        L (list of lists): The lower triangular matrix.
        b (list): The vector of constants.

    Returns:
        list: The solution vector.
    """
    n = len(L)
    x = [0] * n

    for i in range(n):
        x[i] = b[i]
        for j in range(i):
            x[i] -= L[i][j] * x[j]
        x[i] /= L[i][i]

    return x

def backward_substitution(U, b):
    """
    Solves an upper triangular linear system Ux = b using backward substitution.

    Parameters:
        U (list of lists): The upper triangular matrix.
        b (list): The vector of constants.

    Returns:
        list: The solution vector.
    """
    n = len(U)
    x = [0] * n

    for i in range(n - 1, -1, -1):
        x[i] = b[i]
        for j in range(i + 1, n):
            x[i] -= U[i][j] * x[j]
        x[i] /= U[i][i]

    return x

def doolittle_decomposition(matrix):
    """
    Computes the Doolittle LU decomposition of a square matrix.

    Parameters:
        matrix (list of lists): The matrix to be decomposed.

    Returns:
        tuple: A tuple containing the lower triangular matrix L and the upper triangular matrix U.
    """
    n = len(matrix)
    L = [[0] * n for _ in range(n)]
    U = [[0] * n for _ in range(n)]

    for i in range(n):
        L[i][i] = 1

    for i in range(n):
        for j in range(i, n):
            U[i][j] = matrix[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
        for j in range(i + 1, n):
            L[j][i] = (matrix[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

    return L, U

def solve_matrix_equation(matrix, b):
    """
    Solves a matrix equation Ax = b using Cholesky decomposition if positive definite,
    and Doolittle decomposition otherwise.

    Parameters:
        matrix (list of lists): The coefficient matrix A.
        b (list): The vector of constants.

    Returns:
        list: The solution vector x.
    """
    if is_symmetric(matrix):
        try:
            L = cholesky_decomposition(matrix)
            y = forward_substitution(L, b)
            x = backward_substitution(L, y)
            return x
        except ValueError:
            pass
    # If Cholesky decomposition fails or matrix is not symmetric, use Doolittle method
    L, U = doolittle_decomposition(matrix)
    y = forward_substitution(L, b)
    x = backward_substitution(U, y)
    return x

# First problem
A1 = [[1, -1, 3, 2], [-1, 5, -5, -2], [3, -5, 19, 3], [2, -2, 3, 21]]
b1 = [15, -35, 94, 1]

if is_symmetric(A1):
    print("The matrix is symmetric.")
else:
    print("The matrix is not symmetric.")

solution1 = solve_matrix_equation(A1, b1)
print("The solution vector for the first problem is:", solution1)

# Second problem
A2 = [[4, 2, 4, 0], [2, 2, 3, 2], [4, 3, 6, 3], [2, 3, 9, 1]]
b2 = [20, 36, 60, 122]

if is_symmetric(A2):
    print("The matrix is symmetric.")
else:
    print("The matrix is not symmetric.")

solution2 = solve_matrix_equation(A2, b2)
print("The solution vector for the second problem is:", solution2)




