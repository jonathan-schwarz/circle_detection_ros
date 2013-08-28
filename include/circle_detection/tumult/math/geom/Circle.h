#pragma once
#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include <cstdlib>
#include "circle_detection/tumult/math/geom/Point.h"
#include "circle_detection/tumult/math/DoubleMatrix.h"

struct Circle{
	Point	center;
	double	radius;
	
	Circle(const Circle& o)
		: center(o.center)
		, radius(o.radius)
	{}
	Circle(double _radius=1.0)
		:	radius(_radius)
	{}
	Circle(const Point& _center, double _radius=1.0)
		: center(_center)
		, radius(_radius)
	{}
	
	/**
		nach http://www.arndt-bruenner.de/mathe/scripts/kreis3p.htm
	*/
	static bool from3(Circle& sink, const Point& a, const Point& b, const Point& c){
		// let A*x = y
		Matrix<double> A(3,3);
		A.at(0,0) = 1.0;	A.at(0,1) = -a.x;	A.at(0,2) = -a.y;
		A.at(1,0) = 1.0;	A.at(1,1) = -b.x;	A.at(1,2) = -b.y;
		A.at(2,0) = 1.0;	A.at(2,1) = -c.x;	A.at(2,2) = -c.y;
		Matrix<double> y(3,1);
		y.at(0,0) = -( a.x * a.x + a.y * a.y );
		y.at(1,0) = -( b.x * b.x + b.y * b.y );
		y.at(2,0) = -( c.x * c.x + c.y * c.y );

		try{
			DoubleMatrix::gaussj(A,y);	//	now y is x 
		}catch(...){

			return false;
		}

		double xx = 0.5*y.at(1,0);;
		double yy = 0.5*y.at(2,0);;
		double rr = xx*xx + yy*yy - y.at(0,0);
		if(rr<0.0){
			return false;
		}
		sink.center.x = xx;
		sink.center.y = yy;
		sink.radius = sqrt( rr );
		return true;
	}

	double dist(const Point& p)const{
		return fabs( center.dist(p) - radius );
	}
	
	/**
		Ausgleichskreis aus n punkten.
		  - Ueberpruefung, ob genug punkte da sind fehlt...
		nach http://people.inf.ethz.ch/arbenz/MatlabKurs/node79.html
	*/
	static bool fromN(Circle& sink, const std::vector<Point>& points){
		size_t m = points.size();
		if(m<3){
			return false;
		}
		//aktuelle loesung mit
		//	x[0] -> x des mittelpunktes
		//	x[1] -> y des mittelpunktes
		//	x[2] -> radius
		Matrix<double> x(3,1);
		//fehler an stelle x
		Matrix<double> f(m,1);
		//ableitung des fehlers an stelle x (jacobi-matrix)
		Matrix<double> J(m,3);
		//pseudoinverse der jacobi-matrix
		Matrix<double> Ji;
		
		//f(x+h) \approx f(x)+J(x)h \approx 0
		//f(x)+J(x)h = 0
		//J(x)h = -f(x)
		
		//startloesung
		x.at(0,0) = 0.0;
		x.at(1,0) = 0.0;
		x.at(2,0) = 1.0;
		double thresh = 1e-3;
		double norm = 1e300;
		double delta_norm;
		do{
			for(size_t i = 0; i < m; ++i){
				double dx = x.at(0,0)-points[i].x;
				double dy = x.at(1,0)-points[i].y;
				double den = sqrt( dx*dx+dy*dy);
				f.at(i,0) = den-x.at(2,0);
				J.at(i,0) = dx/den;
				J.at(i,1) = dy/den;
				J.at(i,2) = -1.0;
			}
			DoubleMatrix::pseudo_inverse(Ji,J);
			Matrix<double> h = Ji*(f*(-1.0));
			std::cout<<x<<"\n";
			x+=h;
			double new_norm = (f.transpose()*f).at(0,0);
			delta_norm = fabs(norm-new_norm);
			norm = new_norm;
		}while(delta_norm > thresh);
		
		sink.center.x=x.at(0,0);
		sink.center.y=x.at(1,0);
		sink.radius  =x.at(2,0);
		
		return true;
	}
};
