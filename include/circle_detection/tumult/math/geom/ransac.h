#pragma once
#include <vector>
#include "circle_detection/tumult/math/geom/Point.h"
#include "circle_detection/tumult/math/geom/Circle.h"

/**
	sucht wiederholt nach bestunterstuetzten kreisen (solange unterstuetzeranzahl ausreichend)
	  - aus den unterstuetzern des aktuellen kreises wird durch ausgleichsrechung der besteingepasste kreis bestimmt
	  - alle unterstuetzer des ausgleichskreises werden aus den daten entfernt
	  - der ausgleichskreis wird an das ergebnis angehangen
*/
void ransac(std::vector<Circle>& sink, std::vector<Point>& points, double eps, size_t nbModels, size_t min_supporters, size_t maxCircles);
