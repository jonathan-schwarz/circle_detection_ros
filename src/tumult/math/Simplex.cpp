#include "tumult/math/Simplex.h"
#include <cmath>
#include <cassert>

/**
	FROM
	
	Numerical recipes in C : the art of scientific computing / William H. Press ... [et al.]. – 2nd ed.
	ISBN 0-521-43108-5
*/
void simp1(double **a, int mm, int ll[], int nll, int iabf, int *kp,double *bmax){
	int k;
	double test;

	if (nll <= 0)
		*bmax=0.0;
	else {
		*kp=ll[1];
		*bmax=a[mm+1][*kp+1];
		for (k=2;k<=nll;k++) {
			if (iabf == 0)
				test=a[mm+1][ll[k]+1]-(*bmax);
			else
				test=fabs(a[mm+1][ll[k]+1])-fabs(*bmax);
			if (test > 0.0) {
				*bmax=a[mm+1][ll[k]+1];
				*kp=ll[k];
			}
		}
	}
}

void simp2(double **a, int m, int n, double eps, int *ip, int kp){
	int k,i;
	double qp,q0,q,q1;
	
	*ip=0;
	for (i=1;i<=m;i++)
		if (a[i+1][kp+1] < -eps) break;
	if (i>m) return;
	q1 = -a[i+1][1]/a[i+1][kp+1];
	*ip=i;
	for (i=*ip+1;i<=m;i++) {
		if (a[i+1][kp+1] < -eps) {
			q = -a[i+1][1]/a[i+1][kp+1];
			if (q < q1) {
				*ip=i;
				q1=q;
			} else if (q == q1) {
				for (k=1;k<=n;k++) {
					qp = -a[*ip+1][k+1]/a[*ip+1][kp+1];
					q0 = -a[i+1][k+1]/a[i+1][kp+1];
					if (q0 != qp) break;
				}
				if (q0 < qp) *ip=i;
			}
		}
	}
}

void simp3(double **a, int i1, int k1, int ip, int kp){
	int kk,ii;
	double piv;
	piv=1.0/a[ip+1][kp+1];
	for (ii=1;ii<=i1+1;ii++)
		if (ii-1 != ip) {
			a[ii][kp+1] *= piv;
			for (kk=1;kk<=k1+1;kk++)
				if (kk-1 != kp)
					a[ii][kk] -= a[ip+1][kk]*a[ii][kp+1];
		}
	for (kk=1;kk<=k1+1;kk++)
		if (kk-1 != kp) a[ip+1][kk] *= -piv;
	a[ip+1][kp+1]=piv;
}

/**
	double	a[m+3][n+2]:	tableau
	int		m:				number constraints
	int		n:				total number vars
	int		m1:				number constraints less
	int		m2:				number constraints greater
	int		m3:				number constraints equal
	
	
	int		l1[n+2];		internal use, external alloc
	int		l3[m+1];		internal use, external alloc
	
	double	eps:			terminating criterion
	
	On output, the tableau a is indexed by two returned arrays of integers. iposv[j]
	contains, for j= 1 . . . M , the number i whose original variable x i is now represented
	by row j+1 of a. These are thus the left-hand variables in the solution. (The first row
	of a is of course the z-row.) A value i > N indicates that the variable is a y i rather
	than an xi , xN +j -> yj . Likewise, izrov[j] contains, for j= 1 . . . N , the number i
	whose original variable x i is now a right-hand variable, represented by column j+1
	of a. These variables are all zero in the solution. The meaning of i > N is the same
	as above, except that i > N + m 1 + m2 denotes an artificial or slack variable which
	was used only internally and should now be entirely ignored.
	
	int		izrov[n+1]
	int		iposv[m+1]
	
*/
Simplex::Error simplx(double **a, int m, int n, int m1, int m2, int m3, int* l1, int* l3, double eps, int izrov[], int iposv[]){
	int i,ip,is,k,kh,kp,nl1;
	double q1,bmax;
	size_t iter=0;
	size_t iter_max=100;
	
	if (m != (m1+m2+m3)) return Simplex::BAD_INPUT_CONSTRAINT_COUNTS;
	nl1=n;
	for (k=1;k<=n;k++) l1[k]=izrov[k]=k;
	for (i=1;i<=m;i++) {
		if (a[i+1][1] < 0.0) return Simplex::BAD_INPUT_TABLEAU;
		iposv[i]=n+i;
	}
	if (m2+m3) {
		for (i=1;i<=m2;i++) l3[i]=1;
		for (k=1;k<=(n+1);k++) {
			q1=0.0;
			for (i=m1+1;i<=m;i++) q1 += a[i+1][k];
			a[m+2][k] = -q1;
		}
		for (;;) {
			if((iter++) >iter_max){
				std::cerr<<"SIMPLEX::ITER>ITER_MAX in START\n";
				return Simplex::UNFEASIBLE;
			};
			simp1(a,m+1,l1,nl1,0,&kp,&bmax);
			if (bmax <= eps && a[m+2][1] < -eps)
				return Simplex::UNFEASIBLE;
			else if (bmax <= eps && a[m+2][1] <= eps) {
				for (ip=m1+m2+1;ip<=m;ip++) {
					if (iposv[ip] == (ip+n)) {
						simp1(a,ip,l1,nl1,1,&kp,&bmax);
						if (bmax > eps)
							goto one;
					}
				}
				for (i=m1+1;i<=m1+m2;i++)
					if (l3[i-m1] == 1)
						for (k=1;k<=n+1;k++)
							a[i+1][k] = -a[i+1][k];
				break;
			}
			simp2(a,m,n,eps,&ip,kp);
			if (ip == 0)
				return Simplex::UNFEASIBLE;
	one:	simp3(a,m+1,n,ip,kp);
			if (iposv[ip] >= (n+m1+m2+1)) {
				for (k=1;k<=nl1;k++)
					if (l1[k] == kp) break;
				--nl1;
				for (is=k;is<=nl1;is++) l1[is]=l1[is+1];
			} else {
				kh=iposv[ip]-m1-n;
				if (kh >= 1 && l3[kh]) {
					l3[kh]=0;
					++a[m+2][kp+1];
					for (i=1;i<=m+2;i++)
						a[i][kp+1] = -a[i][kp+1];
				}
			}
			is=izrov[kp];
			izrov[kp]=iposv[ip];
			iposv[ip]=is;
		}
	}
	for (;;) {
		if((iter++) >iter_max){
			std::cerr<<"SIMPLEX::ITER>ITER_MAX in MAXIMIZE\n";
			return Simplex::UNFEASIBLE;
		};
		simp1(a,0,l1,nl1,0,&kp,&bmax);
		if (bmax <= eps)
			return Simplex::SUCCESS;
		simp2(a,m,n,eps,&ip,kp);
		if (ip == 0)
			return Simplex::UNBOUNDED;
		simp3(a,m,n,ip,kp);
		is=izrov[kp];
		izrov[kp]=iposv[ip];
		iposv[ip]=is;
	}
}

/************************************************************************************

		wrapper implementation

************************************************************************************/

/************************************************************************************
	descriptor
************************************************************************************/
Simplex::Descriptor::~Descriptor(){
	if(size){
		delete[]	coeffs;
	};
};

Simplex::Descriptor::Descriptor(size_t number_coeffs)
	:	type(BAD_DESCRIPTOR),
		size(number_coeffs),
		coeffs(0),
		value(0.0)
{
	if(size){
		coeffs=new double[size];
	}
	for(size_t i = 0; i<size; ++i){
		coeffs[i]=0.0;
	};
};

Simplex::Descriptor::Descriptor(Type _type, double _value, double*	_coeffs, size_t number_coeffs)
	:	size(number_coeffs),
		coeffs(0)
{
	if(size){
		coeffs = new double[size];
	}
	if(_type==OBJECTIVE){
		_value=0.0;
	};
	if(_value<0){
		value = -_value;
		for(size_t i = 0; i<size; ++i){
			coeffs[i]=-_coeffs[i];
		};
		if(		_type==CONSTRAINT_LESS){	type=CONSTRAINT_GREATER;		}
		else if(_type==CONSTRAINT_GREATER){	type=CONSTRAINT_LESS;			}
		else{								type=_type;			}
	}
	else{
		value = _value;
		for(size_t i = 0; i<size; ++i){
			coeffs[i]=_coeffs[i];
		};
		type=_type;
	};
};
Simplex::Descriptor::Descriptor(const Descriptor& o)
	:	type(o.type),
		size(o.size),
		coeffs(0),
		value(o.value)
{
	if(size){
		coeffs=new double[size];
	}
	for(size_t i = 0; i<size; ++i){
		coeffs[i]=o.coeffs[i];
	};
};
Simplex::Descriptor& Simplex::Descriptor::operator=(const Descriptor& o){
	if(&o==this){
		return *this;
	};
	value=o.value;
	if(o.size!=size){
		if(size){
			delete[]	coeffs;
		}
		coeffs=0;
		size=o.size;
		if(size){
			coeffs=new double[size];
		}
	}
	for(size_t i = 0; i<size; ++i){
		coeffs[i]=o.coeffs[i];
	};
	type=o.type;
	return *this;
}

std::ostream& operator<<(std::ostream& o, const Simplex::Descriptor& d){
	switch(d.type){
	case Simplex::Descriptor::SOLUTION:
		o<<"SOLUTION: z=\t"<<d.value<<"\t;\tx=(\t";
		for(size_t i = 0; i < d.size; ++i){
			o<<d.coeffs[i]<<"\t";
		}
		o<<")";
		break;
	case Simplex::Descriptor::OBJECTIVE:
		o<<"OBJECTIVE: maximise z=\t"<<"(\t";
		for(size_t i = 0; i < d.size; ++i){
			o<<d.coeffs[i]<<"\t";
		}
		o<<") x";
		break;
	case Simplex::Descriptor::CONSTRAINT_LESS:
		o<<"CONSTRAINT:\t\t(\t";
		for(size_t i = 0; i < d.size; ++i){
			o<<d.coeffs[i]<<"\t";
		}
		o<<") x <=\t"<<d.value;
		break;
	case Simplex::Descriptor::CONSTRAINT_GREATER:
		o<<"CONSTRAINT:\t\t(\t";
		for(size_t i = 0; i < d.size; ++i){
			o<<d.coeffs[i]<<"\t";
		}
		o<<") x >=\t"<<d.value;
		break;
	case Simplex::Descriptor::CONSTRAINT_EQUAL:
		o<<"CONSTRAINT:\t\t(\t";
		for(size_t i = 0; i < d.size; ++i){
			o<<d.coeffs[i]<<"\t";
		}
		o<<") x ==\t"<<d.value;
		break;
	default:
		o<<"BAD_DESCRIPTOR";
	}
	o<<"\n";
	return o;
};

/************************************************************************************
	simplex
************************************************************************************/
void Simplex::print_error(std::ostream& o,Error e){
	switch(e){
		case	SUCCESS						:	o<<"SUCCESS";						return;
		case	UNBOUNDED					:	o<<"UNBOUNDED";						return;
		case	UNFEASIBLE					:	o<<"UNFEASIBLE";					return;
		case	BAD_INPUT_CONSTRAINT_COUNTS	:	o<<"BAD_INPUT_CONSTRAINT_COUNTS";	return;
		case	BAD_INPUT_TABLEAU			:	o<<"BAD_INPUT_TABLEAU";				return;
	};
};

void	Simplex::clear(){
	for(size_t i = 0; i < 3; ++i){
		constraints[i].clear();
	};
};

void	Simplex::addDescriptor(const Descriptor& d){
	if(d.type==Descriptor::OBJECTIVE){
		objective=d;
		return;
	}
	constraints[d.type].push_back(d);
};

Simplex::Error Simplex::solve(Descriptor& solution,double eps){
	size_t	m1 = constraints[0].size();
	size_t	m2 = constraints[1].size();
	size_t	m3 = constraints[2].size();
	size_t	M  = m1+m2+m3;
	size_t	n  = objective.size;	//cnt_primary_variables
	size_t	N  = n+m1+m2;			//cnt_varibles
	Matrix<double>	a(M+3,N+2);
	Matrix<int>		l1(1,N+2);
	Matrix<int>		l3(1,M+1);
	Matrix<int>		izrov(1,N+1);
	Matrix<int>		iposv(1,M+1);
	
	//	set_objective
	size_t	row=1;
	a.at(row,1)=0.0;
	for(size_t j = 2; j < n+2; ++j){
		a.at(row,j)=objective.coeffs[j-2];
	};
	++row;
	//set_constraints
	for(size_t k = 0; k < 3; ++k){
		for(size_t i = 0; i < constraints[k].size(); ++i){
			a.at(row,1)=constraints[k][i].value;
			for(size_t j = 2; j < n+2; ++j){
				a.at(row,j)=-constraints[k][i].coeffs[j-2];
			};
			++row;
		};
	};
	
	//prepare_slacks
	for(size_t i = 0; i < M+1; ++i){
		for(size_t j = n+1; j<N+1; ++j){
			a.at(i+1,j+1)=0.0;
		}
	}
	for(size_t i = 1; i < m1+1; ++i){
		a.at(i+1,n+i+1)=-1.0;
	}
	for(size_t i = m1+1; i < m1+m2+1; ++i){
		a.at(i+1,n+i+1)=1.0;
	}
	
	//solve
	Error e = simplx(a.raw(), M, N, m1, m2, m3, l1.raw()[0], l3.raw()[0], eps, izrov.raw()[0], iposv.raw()[0]);
	
	//translate result
	solution=Descriptor(n);
	solution.type=Descriptor::SOLUTION;
	if(e==SUCCESS){
		for(size_t j = 1; j < M+1; ++j)
			if(iposv.at(0,j)>n);
			else
				solution.coeffs[iposv.at(0,j)-1]=a.at(j+1,1);
		solution.value= a.at(1,1);
	};
	return e;
};

