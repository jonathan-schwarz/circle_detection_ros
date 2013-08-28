#pragma once

#include <iostream>
#include <assert.h>
#include "circle_detection/tumult/tools/Vertex.h"

 
template<typename TYPE>
class Matrix
{
private:
	TYPE**		m_data;
	unsigned int	m_columns;
	unsigned int	m_rows;

private:
	void close();
	void copy(const Matrix<TYPE>& other);
	void init();
		
public:
	Matrix(unsigned int rows = 0, unsigned int columns = 0);
	Matrix(const Matrix<TYPE>& other);
	Matrix<TYPE>&	operator=(const Matrix<TYPE>& other);
	~Matrix();
	unsigned int columns()const;
	unsigned int rows()const;
	
	TYPE**	raw(){return m_data;};
	
	TYPE&	at(	unsigned int i,	//	row_index
			unsigned int j	//	column_index
	);
	TYPE	at(	unsigned int i,
			unsigned int j
	)const;
	
	TYPE&	operator()(	unsigned int i,	//	row_index
				unsigned int j	//	column_index
	){return at(i,j);}
	TYPE	operator()(	unsigned int i,
				unsigned int j
	)const{return at(i,j);};
	
	Matrix<TYPE>	transpose()const;
	Matrix<TYPE>	operator*(const Matrix<TYPE>& other)const;
	
	Matrix<TYPE>&	operator*=(TYPE elementModifier(TYPE));
	Matrix<TYPE>	operator*(TYPE elementModifier(TYPE))const;
	Matrix<TYPE>&	operator+=(const Matrix<TYPE>& other);
	Matrix<TYPE>&	operator-=(const Matrix<TYPE>& other);
	Matrix<TYPE>	operator+(const Matrix<TYPE>& other)const;
	Matrix<TYPE>	operator-(const Matrix<TYPE>& other)const;
	
	Matrix<TYPE>&	operator*=(TYPE scalar);
	Matrix<TYPE>	operator*(TYPE scalar)const;
	
	Matrix<TYPE>&	operator/=(TYPE scalar);
	Matrix<TYPE>	operator/(TYPE scalar)const;

	void	print(std::ostream&	out)const;
	
	
	void print_ptr(){
		std::cout<<m_data<<"\t\t";
		for(size_t i = 0; i < m_rows; ++i){
			std::cout<<m_data[i]<<"\t";
		};
		std::cout<<"\n";
	};

};

template<typename TYPE>
inline std::ostream& operator<<(std::ostream& out, const Matrix<TYPE>& matrix)
{
	matrix.print(out);
	return out;
};

template<typename TYPE>
inline void Matrix<TYPE>::close()
{
	if(!m_data){
		m_rows =0;
		m_columns=0;
		return;
	}
	if((m_rows)&&(m_columns))
	{
		for(unsigned int i = 0; i < m_rows; ++i)
		{
			if(m_data[i]){
				delete[] m_data[i];
				m_data[i]=0;
			}
		};
		delete[]	m_data;
	};
	m_data = 0;
	m_rows =0;
	m_columns=0;
};

template<typename TYPE>
inline void Matrix<TYPE>::copy(const Matrix<TYPE>& other)
{
	if(!(other.m_data&&other.columns()&&other.rows())){
		return;
	};
	if(		(	m_columns!=other.columns()	)
		||	(	m_rows!=other.rows()		)
		||	(	m_data==0					)
	  )
	{
		close();
		m_columns	=other.columns();
		m_rows		=other.rows();
		init();
	};
	assert(m_data);
	for(unsigned int i = 0; i < m_rows; ++i)
	{
		assert(m_data[i]);
		for(unsigned int j = 0; j < m_columns; ++j)
		{
			at(i,j)=other.at(i,j);
		};
	};
};
template<typename TYPE>
inline void Matrix<TYPE>::init()
{
	assert(m_data==0);
	if((m_columns==0)||(m_rows==0))
	{
		m_data = 0;
		m_columns=0;
		m_rows=0;
		return;
	};
	
	m_data = new TYPE*[m_rows];
	assert(m_data);
	for(unsigned int i = 0; i < m_rows; ++i)
	{
		m_data[i] = new TYPE[m_columns];
		assert(m_data[i]);
		for(unsigned int j = 0; j < m_columns; ++j)
		{
			m_data[i][j] = ((TYPE)0);
		};
	};
	
};

template<typename TYPE>
inline Matrix<TYPE>::Matrix(unsigned int rows, unsigned int columns)
	:	m_data(0),
		m_columns(columns),
		m_rows(rows)
{
	init();
};
	
template<typename TYPE>
inline Matrix<TYPE>::Matrix(const Matrix<TYPE>& other)
	:	m_data(0),
		m_columns(other.columns()),
		m_rows(other.rows())
{
	init();
	copy(other);
};
	
template<typename TYPE>
inline Matrix<TYPE>&	Matrix<TYPE>::operator=(const Matrix<TYPE>& other)
{
	copy(other);
	return *this;
};
	
template<typename TYPE>
inline Matrix<TYPE>::~Matrix()
{
	close();
};
	
template<typename TYPE>
inline unsigned int Matrix<TYPE>::columns()const
{
	return m_columns;
};
template<typename TYPE>
inline unsigned int Matrix<TYPE>::rows()const
{
	return m_rows;
};
	
template<typename TYPE>
inline TYPE&	Matrix<TYPE>::at(	unsigned int i,unsigned int j)
{
	assert(i<m_rows);
	assert(j<m_columns);
	return m_data[i][j];
};
				  
template<typename TYPE>
inline TYPE	Matrix<TYPE>::at(	unsigned int i,unsigned int j)const
{
	assert(i<m_rows);
	assert(j<m_columns);
	return m_data[i][j];
};
	
template<typename TYPE>
inline Matrix<TYPE>	Matrix<TYPE>::transpose()const
{
	Matrix<TYPE>	C(columns(),rows());
	for(unsigned int i = 0; i < C.rows(); ++i)
	{
		for(unsigned int j = 0; j < C.columns(); ++j)
		{
			C.at(i,j)	=	at(j,i);
		};
	};
	return C;
};
	
template<typename TYPE>
inline Matrix<TYPE>	Matrix<TYPE>::operator*(const Matrix<TYPE>& other)const
{
	const Matrix<TYPE>&	A(*this);
	const Matrix<TYPE>& B(other);
	assert(A.columns()==B.rows());
	Matrix<TYPE>		C(A.rows(),B.columns());

	//	AB=C
	//	c_{i,j} = \sum\limits_{k=1}^{m}{a_{i,k}*b_{k,j}}
	for(unsigned int i = 0; i < C.rows(); ++i)
	{
		for(unsigned int j = 0; j < C.columns(); ++j)
		{
			C.at(i,j) = 0.0;
			for(unsigned int k = 0; k < A.columns(); ++k)
			{
				C.at(i,j)	+= (	A.at(i,k)*B.at(k,j)	);
			};
		};
	};
	return C;
};

template<typename TYPE>
inline Matrix<TYPE>&	Matrix<TYPE>::operator*=(TYPE elementModifier(TYPE))
{
	for(unsigned int i = 0; i < rows(); ++i)
	{
		for(unsigned int j = 0; j < columns(); ++j)
		{
			at(i,j)	=	elementModifier(	at(i,j)	);
		};
	};
	return *this;
};

template<typename TYPE>
inline Matrix<TYPE>	Matrix<TYPE>::operator*(TYPE elementModifier(TYPE))const
{
	Matrix<TYPE>	C(*this);
	C*=elementModifier;
	return C;
};

template<typename TYPE>
inline Matrix<TYPE>&	Matrix<TYPE>::operator*=(TYPE scalar)
{
	for(unsigned int i = 0; i < rows(); ++i)
	{
		for(unsigned int j = 0; j < columns(); ++j)
		{
			at(i,j)	*=	scalar;
		};
	};
	return *this;
};

template<typename TYPE>
inline Matrix<TYPE>	Matrix<TYPE>::operator*(TYPE scalar)const
{
	Matrix<TYPE>	C(*this);
	C*=scalar;
	return C;
};

template<typename TYPE>
inline Matrix<TYPE>&	Matrix<TYPE>::operator/=(TYPE scalar)
{
	return (this->operator*=(1.0/scalar));
};

template<typename TYPE>
inline Matrix<TYPE>	Matrix<TYPE>::operator/(TYPE scalar)const
{
	Matrix<TYPE>	C(*this);
	C/=scalar;
	return C;
};

template<typename TYPE>
inline Matrix<TYPE>&	Matrix<TYPE>::operator+=(const Matrix<TYPE>& other)
{
	assert(rows()==other.rows());
	assert(columns()==other.columns());
	for(unsigned int i = 0; i < rows(); ++i)
	{
		for(unsigned int j = 0; j < columns(); ++j)
		{
			at(i,j)	+=	other.at(i,j);
		};
	};
	return *this;
};
template<typename TYPE>
inline Matrix<TYPE>	Matrix<TYPE>::operator+(const Matrix<TYPE>& other)const
{
	Matrix C(*this);
	C+=other;
	return C;
};

template<typename TYPE>
inline Matrix<TYPE>&	Matrix<TYPE>::operator-=(const Matrix<TYPE>& other)
{
	assert(rows()==other.rows());
	assert(columns()==other.columns());
	for(unsigned int i = 0; i < rows(); ++i)
	{
		for(unsigned int j = 0; j < columns(); ++j)
		{
			at(i,j)	-=	other.at(i,j);
		};
	};
	return *this;
};
template<typename TYPE>
inline Matrix<TYPE>	Matrix<TYPE>::operator-(const Matrix<TYPE>& other)const
{
	Matrix C(*this);
	C-=other;
	return C;
};

template<typename TYPE>
inline void	Matrix<TYPE>::print(std::ostream&	out)const
{
	for(unsigned int i = 0; i < rows(); ++i)
	{
		out	<<	"|\t";
		for(unsigned int j = 0; j < columns(); ++j)
		{
			out	<<	at(i,j)	<<	"\t";
		};
		out	<<	"|\n";
	};
};

