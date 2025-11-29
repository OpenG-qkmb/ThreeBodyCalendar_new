#include "3dv.h"
#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>

bool _3dv::iszero(double val) const
{
	return (abs(val) < SMALL_NUM);
}

_3dv& _3dv::operator=(const _3dv& v)
{
	x = v.x; y = v.y; z = v.z;
	return *this;
}
_3dv _3dv::operator+(const _3dv& v) const
{
	return _3dv(x + v.x, y + v.y, z + v.z);
}
_3dv& _3dv::operator+=(const _3dv& v)
{
	x += v.x; y += v.y; z += v.z;
	return *this;
}
_3dv _3dv::operator-(const _3dv& v) const
{
	return _3dv(x - v.x, y - v.y, z - v.z);
}
_3dv& _3dv::operator-=(const _3dv& v)
{
	x -= v.x; y -= v.y; z -= v.z;
	return *this;
}
_3dv _3dv::operator*(double scalar) const // 数乘
{
	return _3dv(x * scalar, y * scalar, z * scalar);
}
_3dv& _3dv::operator*=(double scalar)
{
	x *= scalar; y *= scalar; z *= scalar;
	return *this;
}
_3dv operator*(double scalar, const _3dv& v) // 数乘（实数左乘）
{
	return v * scalar;
}
_3dv _3dv::operator/(double scalar) const
{
	if (iszero(scalar))
		return _3dv(0.0, 0.0, 0.0);
	return _3dv(x / scalar, y / scalar, z / scalar);
}
_3dv& _3dv::operator/=(double scalar)
{
	if (iszero(scalar))
	{
		x = y = z = 0.0;
		return *this;
	}
	x /= scalar; y /= scalar; z /= scalar;
	return *this;
}
bool _3dv::operator==(const _3dv& v) const
{
	return (x == v.x && y == v.y && z == v.z);
}
bool _3dv::operator!=(const _3dv& v) const
{
	return !(*this == v);
}
std::ostream& operator<<(std::ostream& os, const _3dv& v) // 输出流
{
	os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
	return os;
}

double _3dv::dot(const _3dv& v) const // 点积
{
	return x * v.x + y * v.y + z * v.z;
}
_3dv _3dv::cross(const _3dv& v) const // 叉积
{
	return _3dv(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
}
double _3dv::mag_2() const // 模长之平方
{
	return x * x + y * y + z * z;
}
double _3dv::mag() const // 模长
{
	return sqrt(x * x + y * y + z * z);
}
_3dv _3dv::_e() const // 求单位向量
{
	double len = mag();
	if (iszero(len)) return _3dv(0.0, 0.0, 0.0);
	return _3dv(x / len, y / len, z / len);
}
double _3dv::distance_2(const _3dv& v) const // 求两点间距离之平方
{
	return (x - v.x) * (x - v.x) + (y - v.y) * (y - v.y) + (z - v.z) * (z - v.z);
}
double _3dv::distance(const _3dv& v) const // 求两点间距离
{
	return sqrt((x - v.x) * (x - v.x) + (y - v.y) * (y - v.y) + (z - v.z) * (z - v.z));
}

std::string _3dv::to_str() const
{
	return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";
}

bool _3dv::is_zero() const
{
	return iszero(x) && iszero(y) && iszero(z);
}

double _3dv::angle_with(const _3dv& v) const // 求夹角（弧度制）
{
	double dot_prod = this->dot(v);
	double mags = this->mag() * v.mag();
	if (iszero(mags))
		return 0.0;
	double cos_theta = dot_prod / mags;
	if (cos_theta > 1.0)
		cos_theta = 1.0;
	if (cos_theta < -1.0)
		cos_theta = -1.0;
	return acos(cos_theta);
}
_3dv _3dv::proj_onto(const _3dv& v) const // 向量投影
{
	double v_mag_2 = v.mag_2();
	if (iszero(v_mag_2))
		return _3dv(0.0, 0.0, 0.0);
	double scalar = this->dot(v) / v_mag_2;
	return v * scalar;
}
_3dv _3dv::reflect(const _3dv& normal) const // 向量反射（给出反射面之法向量）
{
	_3dv e = normal._e();
	double dot_prod = this->dot(e);
	return *this - e * (2.0 * dot_prod);
}