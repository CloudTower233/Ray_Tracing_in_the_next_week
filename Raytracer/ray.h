#pragma once
#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class ray
{
public:
	ray() {}
	ray(const point3 origin, const vec3 direction,double t = 0.0) : orig(origin), dir(direction),_time(t) {}

	point3 origin() const { return orig; }
	vec3 direction() const { return dir; }
	double time() const { return _time; }
	point3 at(double t) const { return orig + t * dir; }
	
public:
	point3 orig;
	vec3 dir;
	double _time;
};

#endif