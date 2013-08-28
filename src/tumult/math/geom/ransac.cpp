#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include <cstdlib>
#include "circle_detection/tumult/math/geom/ransac.h"
#include "circle_detection/tumult/math/geom/Point.h"
#include "circle_detection/tumult/math/geom/Circle.h"

/**
	kreis aus 3 zufaellig gewaehlten punkten
	  - nicht sehr elegant geloest...
*/
Circle generate_random_model(const std::vector<Point>& points){
	assert(points.size()>2);

	Circle r;
	while(
		!Circle::from3(
			r,
			points[rand()%points.size()],
			points[rand()%points.size()],
			points[rand()%points.size()]
		)
	);
	return r;
}
/**
	anzahl der unterstuetzer zu einem kreis bestimmen
*/
size_t check_model(const Circle& model, const std::vector<Point>& points, double eps){
	size_t cnt = 0;
	for(size_t i = 0; i < points.size(); ++i){
		if(model.dist(points[i])<eps){
			++cnt;
		}
	}
	return cnt;
}
/**
	unterstuetzer eines kreises bestimmen
*/
void get_supporters(std::vector<Point>& supporters, const Circle& model, const std::vector<Point>& points, double eps){
	size_t cnt = 0;
	for(size_t i = 0; i < points.size(); ++i){
		if(model.dist(points[i])<eps){
			supporters.push_back(points[i]);
		}
	}
}
/**
	unterstuetzer eines kreises aus den ausgangsdaten loeschen
*/
void remove_supporters(const Circle& model, std::vector<Point>& points, double eps){
	size_t cnt = 0;
	std::vector<Point>::iterator i = points.begin();
	while( i != points.end() ){
		if(model.dist(*i)<eps){
			i = points.erase(i);
		}else{
			++i;
		}
	}
}

/**
	bestimmt den bestunterstuetzten kreis in den ausgangsdaten,
	ausgangsdaten werden nicht veraendert
		- points:		daten
		- eps:			maximale entfernung eines punktes zum kreis
		- iter:			wieviele versuche
		- supporters:	anzahl der unterstuetzer des gefundenen kreises
*/
Circle ransac_1(const std::vector<Point>& points, double eps, size_t iter, size_t& supporters){
	Circle m = generate_random_model(points);
	supporters = check_model(m,points,eps);
	for(size_t i = 1; i < iter; ++i){
		Circle m_tmp = generate_random_model(points);
		size_t c_tmp = check_model(m_tmp,points,eps);
		if(c_tmp>supporters){
			supporters = c_tmp;
			m = m_tmp;
		}
	}
	return m;
}
/**
	sucht wiederholt nach bestunterstuetzten kreisen (solange unterstuetzeranzahl ausreichend)
	  - aus den unterstuetzern des aktuellen kreises wird durch ausgleichsrechung der besteingepasste kreis bestimmt
	  - alle unterstuetzer des ausgleichskreises werden aus den daten entfernt
	  - der ausgleichskreis wird an das ergebnis angehangen
*/
void ransac(std::vector<Circle>& sink, std::vector<Point>& points, double eps, size_t nbModels, size_t min_supporters, size_t maxCircles){
	size_t supporters;
    size_t i=maxCircles;
	while(i--){
		Circle c = ransac_1(points,eps,nbModels,supporters);
        std::cout << c.center.x << "\t" << c.center.y << "\t" << c.radius << "\n";
		if(supporters>=min_supporters){
			std::vector<Point> supporters;
//			get_supporters(supporters, c, points, eps);
//			if(supporters.size()<3){
//				return;
//			}
//			Circle cc;
//			Circle::fromN(cc,supporters);
//			sink.push_back(cc);
//			remove_supporters(cc, points, eps);
            sink.push_back(c);
			remove_supporters(c, points, eps);
		}else{
			return;
		}
	}
}
