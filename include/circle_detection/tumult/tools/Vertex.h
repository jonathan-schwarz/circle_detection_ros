#ifndef	___VERTEX__H
#define	___VERTEX__H

#include <math.h>
#include <iostream>
#include <assert.h>

#include <iomanip>

template<typename SCALAR_TYPE, size_t DIMENSION>
class Vertex
{
public:
	Vertex();
	Vertex(const SCALAR_TYPE&	init				);
	Vertex(const Vertex<SCALAR_TYPE,DIMENSION>&	other);
	Vertex(SCALAR_TYPE data[DIMENSION]);
	virtual ~Vertex();
		
	SCALAR_TYPE	get(size_t dimension) const;
	
	SCALAR_TYPE	norm()const;
	SCALAR_TYPE	distance(const Vertex<SCALAR_TYPE,DIMENSION>& other)const;
		
	Vertex<SCALAR_TYPE,DIMENSION>	operator+(const Vertex<SCALAR_TYPE,DIMENSION>& other)const;
	Vertex<SCALAR_TYPE,DIMENSION>	operator-(const Vertex<SCALAR_TYPE,DIMENSION>& other)const;
	SCALAR_TYPE						operator*(const Vertex<SCALAR_TYPE,DIMENSION>& other)const;
	Vertex<SCALAR_TYPE,DIMENSION>	operator*(const SCALAR_TYPE& scalar)const;
	Vertex<SCALAR_TYPE,DIMENSION>	operator/(const SCALAR_TYPE& scalar)const;
	
	Vertex<SCALAR_TYPE,DIMENSION>&	operator+=(const Vertex<SCALAR_TYPE,DIMENSION>& other);
	Vertex<SCALAR_TYPE,DIMENSION>&	operator-=(const Vertex<SCALAR_TYPE,DIMENSION>& other);
	Vertex<SCALAR_TYPE,DIMENSION>&	operator*=(const SCALAR_TYPE& scalar);
	Vertex<SCALAR_TYPE,DIMENSION>&	operator/=(const SCALAR_TYPE& scalar);
		
	Vertex<SCALAR_TYPE,DIMENSION>&	operator=(const Vertex<SCALAR_TYPE,DIMENSION>& other);

	SCALAR_TYPE&					operator[](size_t index);
	const SCALAR_TYPE&				operator[](size_t index)const;
	size_t							getDimension()const{	return DIMENSION;	};
	
			
	void print(std::ostream&	out)const;
	
	Vertex<SCALAR_TYPE,DIMENSION>	reScale(const SCALAR_TYPE& newLength)const;


private:
	SCALAR_TYPE	m_data[DIMENSION];

protected:
	SCALAR_TYPE*	getData()	{	return	m_data;	};
};



template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>::Vertex()
{};
template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>::Vertex(const SCALAR_TYPE&	init)
{
	for(size_t i = 0; i < DIMENSION; i++ )
	{
		m_data[i] = init;
	}
};

template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>::Vertex(const Vertex<SCALAR_TYPE,DIMENSION>&	other)
{
	for(size_t i = 0; i < DIMENSION; i++ )
	{
		m_data[i] = other.get(i);
	}
};

template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>::Vertex(SCALAR_TYPE data[DIMENSION])
{
	for(size_t i = 0; i < DIMENSION; i++ )
	{
		m_data[i] = data[i];
	}
};

template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>::~Vertex()
{
//		std::cout	<<	"destroy Vertex:\t"	<<	*this	<<	"\tfinished\n";
};
	
template <typename SCALAR_TYPE, size_t DIMENSION>
SCALAR_TYPE	Vertex<SCALAR_TYPE,DIMENSION>::get(size_t dimension) const
{
	assert(dimension < DIMENSION);
	return m_data[dimension];
};

template <typename SCALAR_TYPE, size_t DIMENSION>
SCALAR_TYPE&	Vertex<SCALAR_TYPE,DIMENSION>::operator[](size_t index)
{
	//assert(index>=0);
	assert(index<DIMENSION);
	return	m_data[index];
};
template <typename SCALAR_TYPE, size_t DIMENSION>
const SCALAR_TYPE&	Vertex<SCALAR_TYPE,DIMENSION>::operator[](size_t index)const
{
	//assert(index>=0);
	assert(index<DIMENSION);
	return	m_data[index];
};

template <typename SCALAR_TYPE, size_t DIMENSION>
SCALAR_TYPE	Vertex<SCALAR_TYPE,DIMENSION>::norm()const
{
	const Vertex<SCALAR_TYPE,DIMENSION>&	me = *this;
	return sqrt( me * me );
};
template <typename SCALAR_TYPE, size_t DIMENSION>
SCALAR_TYPE	Vertex<SCALAR_TYPE,DIMENSION>::distance(const Vertex<SCALAR_TYPE,DIMENSION>& other)const
{
	return Vertex<SCALAR_TYPE,DIMENSION>(*this-other).norm();
};
	
template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>	Vertex<SCALAR_TYPE,DIMENSION>::operator+(const Vertex<SCALAR_TYPE,DIMENSION>& other)const
{
	Vertex<SCALAR_TYPE,DIMENSION>	temp(*this);
	temp+=other;
	return	temp;
};

template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>	Vertex<SCALAR_TYPE,DIMENSION>::operator-(const Vertex<SCALAR_TYPE,DIMENSION>& other)const
{
	Vertex<SCALAR_TYPE,DIMENSION>	temp(*this);
	temp-=other;
	return	temp;
};

template <typename SCALAR_TYPE, size_t DIMENSION>
SCALAR_TYPE	Vertex<SCALAR_TYPE,DIMENSION>::operator*(const Vertex<SCALAR_TYPE,DIMENSION>& other)const
{
	SCALAR_TYPE swp(0);
	for(size_t i = 0; i < DIMENSION; i++)
	{
		swp += get(i) * other.get(i);
	}
	return swp;
};

template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>	Vertex<SCALAR_TYPE,DIMENSION>::operator*(const SCALAR_TYPE& scalar)const
{
	Vertex<SCALAR_TYPE,DIMENSION>	temp(*this);
	temp*=scalar;
	return	temp;
};

template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>	Vertex<SCALAR_TYPE,DIMENSION>::operator/(const SCALAR_TYPE& scalar)const
{
	Vertex<SCALAR_TYPE,DIMENSION>	temp(*this);
	temp/=scalar;
	return	temp;
};

template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>&	Vertex<SCALAR_TYPE,DIMENSION>::operator+=(const Vertex<SCALAR_TYPE,DIMENSION>& other)
{
	for(size_t i = 0; i < DIMENSION; i++)
	{
		m_data[i] += other.get(i);
	}
	return *this;
};
	
template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>&	Vertex<SCALAR_TYPE,DIMENSION>::operator-=(const Vertex<SCALAR_TYPE,DIMENSION>& other)
{
	for(size_t i = 0; i < DIMENSION; i++)
	{
		m_data[i] -= other.get(i);
	}
	return *this;
};

template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>&	Vertex<SCALAR_TYPE,DIMENSION>::operator*=(const SCALAR_TYPE& scalar)
{
	for(size_t i = 0; i < DIMENSION; i++)
	{
		m_data[i] *= scalar;
	}
	return *this;
};

template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>&	Vertex<SCALAR_TYPE,DIMENSION>::operator/=(const SCALAR_TYPE& scalar)
{
	for(size_t i = 0; i < DIMENSION; i++)
	{
		m_data[i] /= scalar;
	}
	return *this;
};

	
template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>&	Vertex<SCALAR_TYPE,DIMENSION>::operator=(const Vertex<SCALAR_TYPE,DIMENSION>& other)
{
	if(&other != this)
	{
		for(size_t i = 0; i < DIMENSION; i++ )
		{
			m_data[i] = other.get(i);
		}
	}
	return *this;
};

template <typename SCALAR_TYPE, size_t DIMENSION>
Vertex<SCALAR_TYPE,DIMENSION>	Vertex<SCALAR_TYPE,DIMENSION>::reScale(const SCALAR_TYPE& newLength)const
{
	Vertex<SCALAR_TYPE,DIMENSION> ret = *this/this->norm();
	ret = ret*newLength;
	return ret;
};


template <typename SCALAR_TYPE, size_t DIMENSION>
void Vertex<SCALAR_TYPE,DIMENSION>::print(std::ostream&	out)const
{
	for(size_t i = 0; i < DIMENSION; i++)
	{
		out << std::setw(8)	<<	std::setprecision(4)<< this->get(i);
		if( i < DIMENSION - 1 )
		{
			out  << "\t";
		}
	}
};


//output...
template <typename SCALAR_TYPE, size_t DIMENSION>
std::ostream& operator<<(std::ostream& out, const Vertex<SCALAR_TYPE,DIMENSION>& vertex)
{
	vertex.print(out);
	return out;
};


#endif	//___VERTEX__H 
