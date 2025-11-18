#pragma once

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include <iostream>
#include "3dv.h"
#include <vector>
#include <string>

enum _type{STAR, PLANET}; // 恒星 行星

// 常用常量

constexpr double G = 6.67430e-11;        // 引力常数 (m³/kg/s²)
constexpr double SOLAR_MASS = 1.989e30;  // 太阳质量 (kg)
constexpr double EARTH_MASS = 5.972e24;  // 地球质量 (kg)
constexpr double AU = 1.496e11;          // 天文单位 (m)
constexpr double DAY_SECONDS = 86400.0;  // 一天的秒数
constexpr double STD_YEAR = 365.25;      // 一年的天数

class _obj // 天体
{
public:
	_type type;
	double m; // 质量（视为质点）
	_3dv pos, v, a; // 位置 速度 加速度
	std::string id; // 名称（不建议重复，后续用户交互时设置限制）

	// 构造函数

	_obj() : type(STAR), m(SOLAR_MASS), pos(_3dv()), v(_3dv()), a(_3dv()), id("sun") {} // 不要使用这种方法初始化，否则面临重名风险
	_obj(_type type, double m, const _3dv& pos, const _3dv& v, const _3dv& a, const std::string& st) : type(type), m(m), pos(pos), v(v), a(a), id(st) {}
	_obj(const _obj& o) : type(o.type), m(o.m), pos(o.pos), v(o.v), a(o.a), id(o.id) {}
	~_obj() = default;

	// 运算符重载

	bool operator==(const _obj& o) const // == 判断id是否相等
	{
		return (id == o.id);
	}
	bool operator!=(const _obj& o) const
	{
		return !(*this == o);
	}
	friend std::ostream& operator<<(std::ostream& os, const _obj& o) // 输出流
	{
		os	<< "Object \"" << o.id << "\":"
			<< std::endl << "Type=" << (o.type == STAR ? "STAR" : "PLANET")
			<< std::endl << "Mass=" << o.m
			<< std::endl << "Position=" << o.pos
			<< std::endl << "Velocity=" << o.v
			<< std::endl << "Acceleration=" << o.a << std::endl;
		return os;
	}

};

class _state
{
public:
	std::vector<_obj> objs;
	double time;
	bool available = true; // 状态是否有效（天体未相撞）

	// 构造函数

	_state() : objs(), time(0.0), available(true) {}
	/*_state(size_t size = 0)
	*{
	*	objs.resize(size);
	*	time = 0.0;
	*	available = true;
	*}*/
	~_state() = default;

	// 运算符重载

	friend std::ostream& operator<<(std::ostream& os, const _state& s) // 输出流
	{
		os << "State at time " << s.time << " s:" << std::endl;
		for (const _obj& o : s.objs)
		{
			os << "  " << o << std::endl;
		}
		return os;
	}

	// 有关操作

	size_t size() const
	{
		return objs.size();
	}

	void add(const _obj& o)
	{
		objs.push_back(o);
	}

	void set_a() // 设置加速度
	{
		for (_obj& i : objs)
		{
			i.a = _3dv(0);
			for (_obj& j : objs)
			{
				if (i == j)
					continue;
				_3dv r_ij = j.pos - i.pos;
				double dist_2 = r_ij.mag_2();
				if (r_ij.is_zero())
				{
					// 防止除以零，后续作状态检查时直接宣布天体相撞爆炸，结束模拟
					available = false;
					continue; // 合体了，不再计算引力
				}
				i.a += (G * j.m / dist_2) * r_ij._e();
			}
		}
		return;
	}

	// man, what can i say!!!
};

#endif