#pragma once
#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

class vec3
{
public:
	vec3() : e{ 0, 0, 0 } {};		// 无参构造函数
	vec3(double e0, double e1, double e2) : e{ e0, e1, e2 } {};   // 有参构造函数

	double x() const { return e[0]; }
	double y() const { return e[1]; }
	double z() const { return e[2]; }

	// 运算符重载
	vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
	double operator[](int index) const { return e[index]; }
	double& operator[](int index) { return e[index]; }

	vec3& operator+=(vec3 t)
	{
		e[0] += t.x(); e[1] += t.y(); e[2] += t.z();
		return *this;
	}

	vec3& operator*=(double t)
	{
		e[0] *= t; e[1] *= t; e[2] *= t;
		return *this;
	}

	vec3& operator/=(double t)
	{
		e[0] /= t; e[1] /= t; e[2] /= t;
		return *this;
	}

	double length_square() const
	{
		return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
	}

	double length() const
	{
		return sqrt(length_square());
	}

	inline static vec3 random() {
		return vec3(random_double(), random_double(), random_double());
	}

	inline static vec3 random(double min, double max) {
		return vec3(random_double(min, max), random_double(min, max), random_double(min, max));
	}

	bool near_zero() const {
		const auto s = 1e-8;
		return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
	}

	
	
public:
	double e[3];
};

using point3 = vec3;		// point 3D
using color = vec3;         // RGB

inline std::ostream& operator<< (std::ostream& out, vec3& v)
{
	return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vec3 operator+(const vec3& v, const vec3& t)
{
	return vec3(v.e[0] + t.e[0], v.e[1] + t.e[1], v.e[2] + t.e[2]);
}

inline vec3 operator-(const vec3& v, const vec3& t)
{
	return vec3(v.e[0] - t.e[0], v.e[1] - t.e[1], v.e[2] - t.e[2]);
}

inline vec3 operator*(const vec3& v, const vec3& t)
{
	return vec3(v.e[0] * t.e[0], v.e[1] * t.e[1], v.e[2] * t.e[2]);
}

inline vec3 operator*(double t, const vec3& v)
{
	return vec3(v.e[0] * t, v.e[1] * t, v.e[2] * t);
}

inline vec3 operator*(const vec3& v, double t)
{
	return t * v;
}

inline vec3 operator/(const vec3& v, double t)
{
	return 1 / t * v;
}

inline double dot(const vec3& v, const vec3& u)		// 点乘
{
	return v.e[0] * u.e[0] + v.e[1] * u.e[1] + v.e[2] * u.e[2];
}

inline vec3 cross(const vec3& v, const vec3& u)		// 叉乘（按照叉乘公式写）
{
	return vec3(
		v.e[1] * u.e[2] - v.e[2] * u.e[1],
		v.e[2] * u.e[0] - v.e[0] * u.e[2],
		v.e[0] * u.e[1] - v.e[1] * u.e[0]
	);
}

inline vec3 unit_vector(const vec3& v)				// 求单位向量
{
	return v / v.length();
}

inline vec3 random_in_unit_sphere() {
	while (true) {
		auto p = vec3::random(-1, 1);
		if (p.length_square() >= 1) continue;
		//std::cerr << "\n"<<p[0]<<" "<<p[1]<<" "<<p[2]<<"\n";
		return p;
	}
}

inline vec3 random_unit_vector() {
	return unit_vector(random_in_unit_sphere());
}

inline vec3 random_on_hemisphere(const vec3& normal) {
	vec3 on_unit_sphere = random_unit_vector();
	if (dot(on_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
		return on_unit_sphere;
	else
		return -on_unit_sphere;
}

inline vec3 random_in_unit_disk() {
	while (true) {
		auto p = vec3(random_double(-1, 1), random_double(-1, 1), 0);
		if (p.length_square() >= 1) continue;
		return p;
	}
}

inline vec3 reflect(const vec3& v, const vec3& n) {
		return v - 2 * dot(v, n) * n;
}

inline vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
	auto cos_theta = fmin(dot(-uv, n), 1.0);
	vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
	vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_square())) * n;
	return r_out_perp + r_out_parallel;
}




#endif
