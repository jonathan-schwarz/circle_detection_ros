#pragma once
#include <iostream>
#include <vector>
#include "tumult/math/Matrix.h"


struct Simplex
{
	enum Error{
		SUCCESS,
		UNBOUNDED,
		UNFEASIBLE,
		BAD_INPUT_CONSTRAINT_COUNTS,
		BAD_INPUT_TABLEAU
	};
	
	struct Descriptor{
		enum	Type{
			CONSTRAINT_LESS		=0,
			CONSTRAINT_GREATER	=1,
			CONSTRAINT_EQUAL	=2,
			OBJECTIVE			=3,
			SOLUTION			=4,
			BAD_DESCRIPTOR		=5
		};
		
		Type	type;
		size_t	size;
		double*	coeffs;
		double	value;
		
		~Descriptor();
		Descriptor(size_t number_coeffs=0);
		Descriptor(Type _type, double _value, double*	_coeffs, size_t number_coeffs);
		Descriptor(const Descriptor& o);
		Descriptor& operator=(const Descriptor& o);
	};
	
	Descriptor				objective;
	std::vector<Descriptor>	constraints[3];
	
	Error solve(Descriptor& solution,double eps=1e-6);
	
	void	addDescriptor(const Descriptor& d);
	void	clear();
	void	print_error(std::ostream& o,Error e);
	
	static Descriptor non_non_negative_to_non_negative_descriptor(const Descriptor& desc){
		Descriptor r(2*desc.size);
		r.type = desc.type;
		r.value= desc.value;
		for(size_t i = 0; i < desc.size; ++i){
			r.coeffs[2*i  ] =  desc.coeffs[i];
			r.coeffs[2*i+1] = -desc.coeffs[i];
		};
		return r;
	};
	static Descriptor non_negative_to_non_non_negative_descriptor(const Descriptor& desc){
		Descriptor r(desc.size/2);
		r.type = desc.type;
		r.value= desc.value;
		for(size_t i = 0; i < r.size; ++i){
			r.coeffs[i] =  desc.coeffs[2*i] - desc.coeffs[2*i+1];
		};
		return r;
	}
};

std::ostream& operator<<(std::ostream& o, const Simplex::Descriptor& d);
inline std::ostream& operator<<(std::ostream& o, const Simplex& s){
	o<<s.objective;
	for(size_t j = 0; j < 3; ++j){
		for(size_t i = 0; i < s.constraints[j].size(); ++i){
			o<<s.constraints[j][i];
		}
	};
	return o;
};
