#include "tumult/math/DoubleMatrix.h"
#include "tumult/math/matrix_algorithms.h"

Matrix<double> DoubleMatrix::identity(size_t size){
	Matrix<double> r(size,size);
	for(size_t i = 0; i < size; ++i){
		r.at(i,i)=1.0;
	};
	return r;
};
void DoubleMatrix::invert(Matrix<double>& m){
	assert(m.rows()==m.columns());
	Matrix<double> i=identity(m.rows());
	//gjelim(m.raw(),i.raw(),m.rows(),i.columns());
	gaussj(m,i);
};
//m wird ersetzt!
void DoubleMatrix::eigenVectors(Matrix<double>& m){
	assert(m.rows()==m.columns());
	double*	diagonal=new double[m.rows()];
	double*	offDiagonal=new double[m.rows()];
	tred2(m.raw(), m.rows(), diagonal, offDiagonal);
	tqli(diagonal, offDiagonal, m.rows(), m.raw());
	delete[]	diagonal;
	delete[]	offDiagonal;
};
void DoubleMatrix::eigenSystem(Matrix<double>& eigenValues, Matrix<double>& _eigenVectors,const Matrix<double>& m){
	assert(m.rows()==m.columns());
	_eigenVectors = m;
	eigenVectors(_eigenVectors);
	eigenValues=_eigenVectors.transpose()*m*_eigenVectors;
};

bool DoubleMatrix::singular_value_decomposition(Matrix<double>& U, Matrix<double>& W, Matrix<double>& V,const Matrix<double>& A){
	U=A;
	V=Matrix<double>(U.columns(),U.columns());
	W=Matrix<double>(U.columns(),U.columns());
	assert(U.columns());
	double*	w=new double[U.columns()];
	/**
		Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
		U ·W ·V T . The matrix U replaces a on output. The diagonal matrix of singular values W is out-
		put as a vector w[1..n]. The matrix V (not the transpose V T ) is output as v[1..n][1..n] .
	*/
	bool ok = svdcmp(U.raw(), U.rows(), U.columns(), w, V.raw());
	for(size_t i=0; i < U.columns(); ++i){
		for(size_t j=0; j < U.columns(); ++j){
			if(i==j){	W(i,j)=w[i];	}
			else{		W(i,j)=0.0;		};
		};
	};
	//std::cerr<<ok<<"blib\n";
	if(w){
		delete[] w;
	};
	//std::cerr<<ok<<"blub\n";
	return ok;
};

void DoubleMatrix::solve(Matrix<double>& x, const Matrix<double>& A, const Matrix<double>& b){
	/**
		Solves A·X = B for a vector X, where A is specified by the arrays u[1..m][1..n] , w[1..n],
		v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for
		square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
		No input quantities are destroyed, so the routine may be called sequentially with different b’s.
	*/
	
	x=Matrix<double>(A.columns(),1);
	Matrix<double>	U(A);
	Matrix<double>	V=Matrix<double>(U.columns(),U.columns());
	double*	w=new double[U.columns()];
	svdcmp(U.raw(), U.rows(), U.columns(), w, V.raw());
	
	double wmax=0.0;
	for(size_t i=0;i<U.columns();++i){
		if(w[i]>wmax){wmax=w[i];};
	};
	double wmin=1e-6*wmax;
	for(size_t i=0;i<U.columns();++i){
		if(w[i]<wmin){w[i]=0.0;};
	};
	
	double*	x_=new double[x.rows()];
	double*	b_=new double[b.rows()];
	
	for(size_t i = 0; i < x.rows(); ++i){ x_[i]=x(i,0);	};
	for(size_t i = 0; i < b.rows(); ++i){ b_[i]=b(i,0);	};
	
	svbksb(U.raw(),w, V.raw(), U.rows(), U.columns(), b_, x_);
	
	for(size_t i = 0; i < x.rows(); ++i){ x(i,0)=x_[i];	};
	
	if(w){	delete[] w;};
	if(x_){ delete[] x_;};
	if(b_){ delete[] b_;};
};
