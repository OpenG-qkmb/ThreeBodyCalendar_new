#pragma once

#ifndef _3DV_H_
#define _3DV_H_

#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>

constexpr auto SMALL_NUM = 1e-7;

class _3dv // 3D Vector
{
private:
	bool iszero(double val) const;
public:
	double x, y, z;
	
	// 构造函数

	_3dv() : x(0.0), y(0.0), z(0.0) {}
	_3dv(double x, double y, double z) : x(x), y(y), z(z) {}
	_3dv(double mag) : x(mag), y(mag), z(mag) {}
	_3dv(const _3dv& v) : x(v.x), y(v.y), z(v.z) {}
	~_3dv() = default;

	// 运算符重载

	_3dv& operator=(const _3dv& v);
	_3dv operator+(const _3dv& v) const;
	_3dv& operator+=(const _3dv& v);
	_3dv operator-(const _3dv& v) const;
	_3dv& operator-=(const _3dv& v);
	_3dv operator*(double scalar) const; // 数乘
	_3dv& operator*=(double scalar);
	friend _3dv operator*(double scalar, const _3dv& v); // 数乘（实数左乘）
	_3dv operator/(double scalar) const;
	_3dv& operator/=(double scalar);
	bool operator==(const _3dv& v) const;
	bool operator!=(const _3dv& v) const;
	friend std::ostream& operator<<(std::ostream& os, const _3dv& v); // 输出流
	
	// 常见运算

	double dot(const _3dv& v) const; // 点积
	_3dv cross(const _3dv& v) const; // 叉积
	double mag_2() const; // 模长之平方
	double mag() const; // 模长
	_3dv _e() const; // 求单位向量
	double distance_2(const _3dv& v) const; // 求两点间距离之平方
	double distance(const _3dv& v) const; // 求两点间距离

	// 展示向量

	std::string to_str() const;

	// 检查零向量

	bool is_zero() const;

	// 向量求夹角、投影和反射

	double cos_angle_with(const _3dv& v) const; // 求夹角（cos）
	double angle_with(const _3dv& v) const; // 求夹角（弧度制）
	_3dv proj_onto(const _3dv& v) const; // 向量投影
	_3dv reflect(const _3dv& normal) const; // 向量反射（给出反射面之法向量）
};

#endif