#pragma once

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include <iostream>
#include "3dv.h"
#include <vector>
#include <string>
#include <cmath>

enum _type{STAR, PLANET}; // 恒星 行星

// 常用常量

constexpr double G = 6.67430e-11;        // 引力常数 (m³/kg/s²)
constexpr double SOLAR_MASS = 1.989e30;  // 太阳质量 (kg)
constexpr double EARTH_MASS = 5.972e24;  // 地球质量 (kg)
constexpr double AU = 1.496e11;          // 天文单位 (m)
constexpr double DAY_SECONDS = 86400.0;  // 一天的秒数
constexpr double STD_YEAR = 365.25;      // 一年的天数
constexpr double V_SUN = 22000.0;        // 太阳系绕银河系中心公转速度 (m/s)

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

	_obj operator+(const _obj& o) const // 合并天体
	{
		_obj res;
		res.m = m + o.m;
		res.pos = (m * pos + o.m * o.pos) / res.m; // 质心位置
		res.v = (m * v + o.m * o.v) / res.m; // 动量守恒
		res.a = _3dv(0); // 加速度作废
		res.id = id + " + " + o.id;
		res.type = (type == STAR || o.type == STAR) ? STAR : PLANET; // 有恒星则为恒星
		return res;
	}
	friend _obj operator+=(const _obj& lhs, const _obj& rhs)
	{
		return lhs + rhs;
	}
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
		os	<< "[Object] \"" << o.id << "\": "
			<< std::endl << "Type = " << (o.type == STAR ? "[STAR]" : "[PLANET]")
			<< std::endl << "Mass = " << o.m
			<< std::endl << "Position = " << o.pos
			<< std::endl << "Velocity = " << o.v
			<< std::endl << "Acceleration = " << o.a << std::endl;
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

	// 有关操作

	size_t size() const
	{
		return objs.size();
	}

	bool isexist(std::string uid) const
	{
		for (const _obj& o : objs)
		{
			if (o.id == uid)
				return true;
		}
		return false;
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

	void merge() // 合并相撞天体
	{
		if (available)
			return;
		while (true)
		{
			bool tag = false;
			for (size_t i = 0; i < objs.size(); ++i)
			{
				for (size_t j = i + 1; j < objs.size(); ++j)
				{
					if (objs[i].pos.distance_2(objs[j].pos) < SMALL_NUM) // 距离过近
					{
						_obj new_obj = objs[i] + objs[j];
						objs.erase(objs.begin() + j);
						objs.erase(objs.begin() + i);
						objs.push_back(new_obj);
						tag = true;
						break;
					}
				}
				if (tag)
					break;
			}
			if (!tag)
				break;
		}
		available = true;
		set_a();
		return;
	}

	void merge_compulsory()
	{
		available = false;
		merge();
		return;
	}

	bool orbit_me(std::string sunid, _obj& o, double r)
	{
		_obj sun = _obj();
		bool flag = false;
		for (const _obj& obj : objs)
		{
			if (obj.id == sunid)
			{
				sun = obj;
				flag = true;
				break;
			}
		}
		if (!flag || sun.type == PLANET)
			return false;
		double v = sqrt(G * sun.m / r);
		o.pos = sun.pos + _3dv(r, 0, 0);
		o.v = sun.v + _3dv(0, v, 0);
		o.a = _3dv();
		add(o);
		return true;
	}

	double get_energy() // 获取系统总能量
	{
		merge_compulsory(); // 强制合并，防止除以零
		double ek = 0.0, ep = 0.0;
		for (size_t i = 0; i < objs.size(); ++i)
		{
			ek += 0.5 * objs[i].m * objs[i].v.mag_2();
			for (size_t j = i + 1; j < objs.size(); ++j)
			{
				double dist = objs[i].pos.distance(objs[j].pos);
				if (dist < SMALL_NUM)
					continue; // 防止除以零
				ep -= G * objs[i].m * objs[j].m / dist;
			}
		}
		return ek + ep;
	}

	// 运算符重载

	friend std::ostream& operator<<(std::ostream& os, _state& s) // 输出流
	{
		os << "State at time " << s.time << " s:" << std::endl;
		for (const _obj& o : s.objs)
		{
			os << "  " << o << std::endl;
		}
		os << "Energy: " << s.get_energy() << std::endl;
		return os;
	}

	// man, what can i say!!!
} state;

#endif