#pragma once
#include <cmath>
#include <iostream>

struct Point{
	double x;
	double y;
	
	Point(const Point& o)
		: x(o.x)
		, y(o.y)
	{}
	
	Point(double _x=0.0, double _y=0.0)
		: x(_x)
		, y(_y)
	{}
	
	double dist(const Point& p)const{
		double dx = p.x-x;
		double dy = p.y-y;
		return sqrt(dx*dx+dy*dy);
	}
};

inline std::ostream& operator<<(std::ostream& o, const Point& p)
{
//    o << "( " << p.x << ", " << p.y << " ) ";
    o << p.x << "\t" << p.y;

    return o;
}
