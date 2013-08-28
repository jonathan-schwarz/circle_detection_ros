#pragma once
#include <list>
#include "circle_detection/tumult/math/Matrix.h"
#include "circle_detection/tumult/tools/Vertex.h"
#include <algorithm>

class DoubleMatrix
{
public:
	static Matrix<double> identity(size_t size);
	static void invert(Matrix<double>& m);
	static void eigenVectors(Matrix<double>& m);
	static void eigenSystem(Matrix<double>& eigenValues, Matrix<double>& _eigenVectors,const Matrix<double>& m);
	
	template<unsigned int SIZE>
	static Matrix<double>	covariance(const std::list<Vertex<double,SIZE> >& centralizedSamples);
	
	static void pseudo_inverse_left(Matrix<double>& pil,const Matrix<double>& A){
		Matrix<double> At=A.transpose();
		Matrix<double> B=At*A;
		invert(B);
		pil=B*At;
	};
	
	static void pseudo_inverse_right(Matrix<double>& pir,const Matrix<double>& A){
		Matrix<double> At=A.transpose();
		Matrix<double> B=A*At;
		DoubleMatrix::invert(B);
		pir=At*B;
	};
	
	static bool singular_value_decomposition(Matrix<double>& U, Matrix<double>& W, Matrix<double>& V,const Matrix<double>& A);
	static void solve(Matrix<double>& x, const Matrix<double>& A, const Matrix<double>& b);
	
	static Matrix<double> inverse(const Matrix<double>& M){
		Matrix<double> r(M);
		invert(r);
		return M;
	};
	
	
	static bool pseudo_inverse(Matrix<double>& pi,const Matrix<double>& A, double null_threshold=1e-12){
		Matrix<double> U;
		Matrix<double> W;
		Matrix<double> V;
		bool ok = singular_value_decomposition(U,W,V,A);
		for(size_t i = 0; i<W.rows(); ++i){
			for(size_t j = 0; j<W.rows(); ++j){
				if(fabs(W(i,j))<null_threshold){
					W(i,j)=0.0;
				}
				else{
					W(i,j)=1.0/W(i,j);
				};
			};
		};
		pi=V*W*U.transpose();
		return ok;
	};
	
	/**
		numerical recepies, gaussj
	*/
	static void gaussj(Matrix<double>& A, Matrix<double>& B){
		size_t	n=A.rows();
		size_t	m=B.columns();
		int*	indxc=new int[n];
		int*	indxr=new int[n];
		int*	ipiv=new int[n];
		
		for(size_t j=0; j<n; ++j){	ipiv[j]=0;	};
		for(size_t i=0; i<n; ++i){
			size_t	icol=0;
			size_t	irow=0;
			double big=0.0;
			for(size_t j=0; j<n; ++j){
				if(ipiv[j]!=1){
					for(size_t k=0; k<n; ++k){
						if(ipiv[k]==0){
							if(fabs(A(j,k))>=big){
								irow=j;
								icol=k;
							}
						}
					};
				}
			};
			++(ipiv[icol]);
			if(irow!=icol){
				for(size_t l=0;l<n;++l){std::swap(A(irow,l),A(icol,l));}
				for(size_t l=0;l<m;++l){std::swap(B(irow,l),B(icol,l));}
			};
			indxr[i]=irow;
			indxc[i]=icol;
			if(A(icol,icol)==0.0){throw("gaussj: Singular Matrix");};
			double pivinv=1.0/A(icol,icol);
			A(icol,icol)=1.0;
			for(size_t l=0;l<n;++l){A(icol,l)*=pivinv;};
			for(size_t l=0;l<m;++l){B(icol,l)*=pivinv;};
			for(size_t ll=0;ll<n;++ll){
				if(ll!=icol){
					double dum=A(ll,icol);
					A(ll,icol)=0.0;
					for(size_t l=0;l<n;++l){A(ll,l)-=(A(icol,l)*dum);};
					for(size_t l=0;l<m;++l){B(ll,l)-=(B(icol,l)*dum);};
				};
			};
		};
		
		//for(int l=n-1;l>=0;--l){
		for(size_t l=n-1;l<n;--l){ //check auf l<n wegen ueberlauf bei 0!
			if(indxr[l]!=indxc[l]){
				for(size_t k=0;k<n;++k){
					std::swap(A(k,indxr[l]),A(k,indxc[l]));
				};
			};
		};
		
		delete[]	indxc;
		delete[]	indxr;
		delete[]	ipiv;
	};
	
};

template<unsigned int SIZE>
inline Matrix<double>	DoubleMatrix::covariance(const std::list<Vertex<double,SIZE> >& centralizedSamples){
	Matrix<double> result(SIZE,SIZE);
	size_t	cnt	=	0;
	for(typename std::list<Vertex<double,SIZE> >::const_iterator i = centralizedSamples.begin(); i!=centralizedSamples.end(); ++i){
		for(unsigned int j = 0; j < SIZE; ++j){
			double	j_val=i->get(j);
			for(unsigned int k = 0; k < SIZE; ++k){
				result.at(j,k)	+=	(j_val*(i->get(k)));
			};
		};
		++cnt;
	};
	return	result/((double)(cnt-1));
};












