#pragma once
#include <cmath>
#include <iostream>

void gjelim(double** lhs, double** rhs, long nrows, long ncolsrhs);
void tred2(double** a, size_t n, double* d, double* e);
void tqli(double* d, double* e, int n, double** z);



/**
	Solves A·X = B for a vector X, where A is specified by the arrays u[1..m][1..n] , w[1..n],
	v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for
	square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
	No input quantities are destroyed, so the routine may be called sequentially with different b’s.
*/
void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);

/**
	Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
	U ·W ·V T . The matrix U replaces a on output. The diagonal matrix of singular values W is out-
	put as a vector w[1..n]. The matrix V (not the transpose V T ) is output as v[1..n][1..n] .
*/
bool svdcmp(double **a, int m, int n, double w[], double **v);


