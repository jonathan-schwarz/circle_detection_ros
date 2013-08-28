#include "tumult/math/matrix_algorithms.h"
#include <assert.h>

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
////	inversion-stuff
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

#include <cmath>
#include <cassert>
#include <iostream>

template<typename T>
struct Array{
	T*		data;
	size_t	size;

	Array(size_t _size)
		:	size(_size)
	{
		data = new T[size];
	};
	~Array(){
		delete[] data;
	};
	T&	operator[](size_t i){
		assert((i-1)<size);
		return data[i-1];
	};
	T&	operator[](int i){
		assert(i>0);
		assert((i-1)<size);
		return data[i-1];
	};
};
template<typename T>
struct EArray{
	T*		data;
	size_t	size;

	EArray(T* _data=0, size_t _size=0)
		:	data(_data),
			size(_size)
	{};
	T&	operator[](size_t i){
		assert((i-1)<size);
		return data[i-1];
	};
	T&	operator[](int i){
		assert(i>0);
		assert((i-1)<size);
		return data[i-1];
	};
};

template<typename T>
struct EMatrix{
	EArray<T>*	data;
	size_t		size;

	EMatrix(T** _data, size_t rows, size_t cols)
		:	size(rows)
	{
		data=new EArray<T>[size];
		for(size_t i = 0; i < size; ++i){
			data[i].data=_data[i];
			data[i].size=cols;
		};
	};
	EArray<T>&	operator[](size_t i){
		assert((i-1)<size);
		return data[i-1];
	};
	EArray<T>&	operator[](int i){
		assert(i>0);
		assert((i-1)<size);
		return data[i-1];
	};
};

inline int		IMIN(int a, int b){return (a < b ? a : b);}
inline double	FMAX(double a, double b){	return	(a > b ?a : b);}
inline double	SIGN(double a, double b){ return ((b >= 0.0) ? fabs(a) : -fabs(a));}
inline double	SQR(double x){return x*x;}
double pythag(double a, double b)/*Computes (a2 + b2 )1/2 without destructive underflow or overflow.*/
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}


/**Solves AX = B for a vector X, where A is specified by the arrays u[1..m][1..n] , w[1..n],
v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for
square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
No input quantities are destroyed, so the routine may be called sequentially with different b?s.*/
void svbksb(double **_u, double _w[], double **_v, int m, int n, double _b[], double _x[]){
	int jj,j,i;
	double s;//,*tmp;
	EMatrix<double>	u(_u,m,n);
	EArray<double>	w(_w,n);
	EMatrix<double>	v(_v,n,n);
	EArray<double>	b(_b,m);
	EArray<double>	x(_x,n);

	Array<double> tmp(n);
	//tmp=vector(1,n);
	for (j=1;j<=n;j++) {/*Calculate U T B.*/
		s=0.0;
		if (w[j]) {/*Nonzero result only if wj is nonzero.*/
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];/*This is the divide by wj .*/
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {/*Matrix multiply by V to get answer.*/
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	//free_vector(tmp,1,n);
}


///beispiel
#if 0
double wmax,wmin,**a,**u,*w,**v,*b,*x;
int i,j;
...
for(i=1;i<=N;i++)/*Copy a into u if you don?t want it to be destroyed.*/
	for j=1;j<=N;j++)
		u[i][j]=a[i][j];
svdcmp(u,N,N,w,v);/*SVD the square matrix a.*/
wmax=0.0;/*Will be the maximum singular value obtained.*/
for(j=1;j<=N;j++) if (w[j] > wmax) wmax=w[j];
/*This is where we set the threshold for singular values allowed to be nonzero. The constant
is typical, but not universal. You have to experiment with your own application.*/
wmin=wmax*1.0e-6;
for(j=1;j<=N;j++) if (w[j] < wmin) w[j]=0.0;
svbksb(u,w,v,N,N,b,x);/*Now we can backsubstitute.*/
#endif


/**Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
U W V^T. The matrix U replaces a on output. The diagonal matrix of singular values W is out-
put as a vector w[1..n]. The matrix V (not the transpose V T ) is output as v[1..n][1..n] .*/
bool svdcmp(double **_a, int m, int n, double _w[], double **_v)
{
	EMatrix<double>	a(_a,m,n);
	EArray<double>	w(_w,n);
	EMatrix<double>	v(_v,n,n);
	
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z;//*rv1;
	//rv1=vector(1,n);
	Array<double> rv1(n);
	g=scale=anorm=0.0;/*Householder reduction to bidiagonal form.*/
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {/*Accumulation of right-hand transformations.*/
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)/*Double division to avoid possible underflow.*/
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {/*Accumulation of left-hand transformations.*/
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {/*Diagonalization of the bidiagonal form: Loop over singular values, and over allowed iterations.*/
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {/*Test for splitting.*/
				nm=l-1;/*Note that rv1[1] is always zero.*/
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if(nm==0){std::cerr<<"hi there\n"; return false;};
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;/*Cancellation of rv1[l], if l > 1.*/
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {/*Convergence.*/
				if (z < 0.0) {/*Singular value is made nonnegative.*/
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30){
				//nrerror("no convergence in 30 svdcmp iterations");
				return false;
			}
			x=w[l];/*Shift from bottom 2-by-2 minor.*/
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;/*Next QR transformation:*/
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;/*Rotation can be arbitrary if z = 0.*/
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	//free_vector(rv1,1,n);
	return true;
}



//        gjelim
void swaprows(double** arr, long row0, long row1){
	double* temp;
	temp=arr[row0];
	arr[row0]=arr[row1];
	arr[row1]=temp;
}

void gjelim(double** lhs, double** rhs, long nrows, long ncolsrhs) {
	//        augment lhs array with rhs array and store in arr2
	assert(nrows);
	assert(nrows+ncolsrhs);
	double** arr2=new double*[nrows];
	for(long row=0; row<nrows; ++row)
		arr2[row]=new double[nrows+ncolsrhs];
	
	for(long row=0; row<nrows; ++row){
		for(long col=0; col<nrows; ++col){
			arr2[row][col]=lhs[row][col];
		}
		for(long col=nrows; col<nrows+ncolsrhs; ++col){
			arr2[row][col]=rhs[row][col-nrows];
		}
	}
	//        perform forward elimination to get arr2 in row-echelon form
	for(long dindex=0; dindex<nrows; ++dindex){
		//        run along diagonal, swapping rows to move zeros in working position
		//        (along the diagonal) downwards
		if( (dindex==(nrows-1)) && (arr2[dindex][dindex]==0)){
			return; //  no solution
		} else if(arr2[dindex][dindex]==0){
			swaprows(arr2, dindex, dindex+1);
		}
		//        divide working row by value of working position to get a 1 on the
		//        diagonal
		if(arr2[dindex][dindex] == 0.0){
			return;
		}else{
			double tempval=arr2[dindex][dindex];
			for(long col=0; col<nrows+ncolsrhs; ++col){
				arr2[dindex][col]/=tempval;
			}
		}
		
		//        eliminate value below working position by subtracting a multiple of
		//        the current row
		for(long row=dindex+1; row<nrows; ++row){
			double wval=arr2[row][dindex];
			for(long col=0; col<nrows+ncolsrhs; ++col){
				arr2[row][col]-=wval*arr2[dindex][col];
			}
		}
	}
	
	//        backward substitution steps
	for(long dindex=nrows-1; dindex>=0; --dindex){
		//        eliminate value above working position by subtracting a multiple of
		//        the current row
		for(long row=dindex-1; row>=0; --row){
			double wval=arr2[row][dindex];
			for(long col=0; col<nrows+ncolsrhs; ++col){
				arr2[row][col]-=wval*arr2[dindex][col];
			}
		}
	}
	
	//        assign result to replace rhs
	for(long row=0; row<nrows; ++row){
		for(long col=0; col<ncolsrhs; ++col){
			rhs[row][col]=arr2[row][col+nrows];
		}
	}
	
	for(long row=0; row<nrows; ++row)
		delete[] arr2[row];
	delete[] arr2;
}

/**  Reduce a real, symmetric matrix to a symmetric, tridiag. matrix. */
/* Householder reduction of matrix a to tridiagonal form.
   Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
   Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
        Springer-Verlag, 1976, pp. 489-494.
        W H Press et al., Numerical Recipes in C, Cambridge U P,
        1988, pp. 373-374.  */
void tred2(double** a, size_t n, double* d, double* e){
	int		l, k, j, i;
	double	scale, hh, h, g, f;
	for (i = n; i >= 2; i--){
		l = i - 1;
		h = scale = 0.0;
		if (l > 1){
			for (k = 1; k <= l; k++)
				scale += fabs(a[i-1][k-1]);
			if (scale == 0.0)
				e[i-1] = a[i-1][l-1];
			else{
				for (k = 1; k <= l; k++){
					a[i-1][k-1] /= scale;
					h += a[i-1][k-1] * a[i-1][k-1];
				}
				f = a[i-1][l-1];
				g = f>0 ? -sqrt(h) : sqrt(h);
				e[i-1] = scale * g;
				h -= f * g;
				a[i-1][l-1] = f - g;
				f = 0.0;
				for (j = 1; j <= l; j++){
					a[j-1][i-1] = a[i-1][j-1]/h;
					g = 0.0;
					for (k = 1; k <= j; k++)
						g += a[j-1][k-1] * a[i-1][k-1];
					for (k = j+1; k <= l; k++)
						g += a[k-1][j-1] * a[i-1][k-1];
					e[j-1] = g / h;
					f += e[j-1] * a[i-1][j-1];
				}
				hh = f / (h + h);
				for (j = 1; j <= l; j++){
					f = a[i-1][j-1];
					e[j-1] = g = e[j-1] - hh * f;
					for (k = 1; k <= j; k++)
						a[j-1][k-1] -= (f * e[k-1] + g * a[i-1][k-1]);
				}
			}
		}
		else
			e[i-1] = a[i-1][l-1];
		d[i-1] = h;
	}
	d[1-1] = 0.0;
	e[1-1] = 0.0;
	for (i = 1; i <= n; i++){
		l = i - 1;
		if (d[i-1]){
			for (j = 1; j <= l; j++){
				g = 0.0;
				for (k = 1; k <= l; k++)
					g += a[i-1][k-1] * a[k-1][j-1];
				for (k = 1; k <= l; k++)
					a[k-1][j-1] -= g * a[k-1][i-1];
			}
		}
		d[i-1] = a[i-1][i-1];
		a[i-1][i-1] = 1.0;
		for (j = 1; j <= l; j++)
			a[j-1][i-1] = a[i-1][j-1] = 0.0;
	}
}

/**  Tridiagonal QL algorithm -- Implicit  **********************/
void tqli(double* d, double* e, int n, double** z){
	int		m, l, iter, i, k;
	double	s, r, p, g, f, dd, c, b;
	for (i = 2; i <= n; i++)
		e[i-1-1] = e[i-1];
	e[n-1] = 0.0;
	for (l = 1; l <= n; l++){
		iter = 0;
		do{
			for (m = l; m <= n-1; m++){
				dd = fabs(d[m-1]) + fabs(d[m+1-1]);
				if (fabs(e[m-1]) + dd == dd) break;
			}
			if (m != l){
				if (iter++ == 30){
					//erhand("No convergence in TLQI.");
					std::cerr<<"No convergence in TLQI.\n";
					assert(false);
				};
				g = (d[l+1-1] - d[l-1]) / (2.0 * e[l-1]);
				r = sqrt((g * g) + 1.0);
				g = d[m-1] - d[l-1] + e[l-1] / (g + SIGN(r, g));
				s = c = 1.0;
				p = 0.0;
				for (i = m-1; i >= l; i--){
					f = s * e[i-1];
					b = c * e[i-1];
					if (fabs(f) >= fabs(g)){
						c = g / f;
						r = sqrt((c * c) + 1.0);
						e[i+1-1] = f * r;
						c *= (s = 1.0/r);
					}
					else{
						s = f / g;
						r = sqrt((s * s) + 1.0);
						e[i+1-1] = g * r;
						s *= (c = 1.0/r);
					}
					g = d[i+1-1] - p;
					r = (d[i-1] - g) * s + 2.0 * c * b;
					p = s * r;
					d[i+1-1] = g + p;
					g = c * r - b;
					for (k = 1; k <= n; k++){
						f = z[k-1][i+1-1];
						z[k-1][i+1-1] = s * z[k-1][i-1] + c * f;
						z[k-1][i-1] = c * z[k-1][i-1] - s * f;
					}
				}
				d[l-1] = d[l-1] - p;
				e[l-1] = g;
				e[m-1] = 0.0;
			}
		}  while (m != l);
	}
}
