# choleskyDecompose
Cholesky Decomposition and Matrix Inversion

This code decomposes a symmetric positive definate square matrix into a lower triangle matrix such that
A = LxL'
where A is the original symmetric matrix, L is the output lower triangular matrix and L' is the transpose of the lower triangular matrix

The code has two functions: 

cholesky_decompose() which decomposes an input square matrix A of size n and outputs matrix L

cholesky_inverse() which inverts an square matrix A (doesn't have to be symmetric) and outputs A_inv
Time complexity is O(N^3) which beats conventional matrix inversion approach that uses recursive determinants and has O(N!) complexity.

The repository also contains a separate implementation file with some matrix operations used by the inversion function. 
