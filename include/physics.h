#pragma once

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include <iostream>
#include "3dv.h"
#include <vector>
#include <string>
#include <functional>
#include <cmath>

enum _type{STAR, PLANET}; // 恒星 行星

// 常用常量


namespace oriphy
{
	constexpr double G = 6.67430e-11;        // 引力常数 (m³/kg/s²)
	constexpr double SOLAR_MASS = 1.989e30;  // 太阳质量 (kg)
	constexpr double EARTH_MASS = 5.965e24;  // 地球质量 (kg)
	constexpr double AU = 1.496e11;          // 天文单位 (m)
	constexpr double DAY_SECONDS = 86400.;  // 一天的秒数
	constexpr double STD_YEAR = 365.25;      // 一年的天数
	constexpr double V_SUN = 22000.;        // 太阳系绕银河系中心公转速度 (m/s)
}

// 事实是：使用上面这样带原始单位的物理量作模拟，性能较差且不易于交互，故作如下处理：

namespace phy
{
	constexpr double AU = 1.;
	constexpr double EARTH_MASS = 1.;
	constexpr double HOUR = 1.;
	constexpr double DAY = 24.; // 小时
	constexpr double YEAR = 365.25 * DAY;
	constexpr double SOLAR_MASS = oriphy::SOLAR_MASS / oriphy::EARTH_MASS; // = 333445.0964 (EARTH_MASS)
	constexpr double G = oriphy::G * (3600. * 3600.) * (oriphy::EARTH_MASS) / (oriphy::AU * oriphy::AU * oriphy::AU); // (AU³/EARTH_MASS/h²)
	constexpr double V_SUN = oriphy::V_SUN * (3600.) / (oriphy::AU); // (AU/h)
	constexpr double CRASH = 0.09; // 撞击阈值
}

// 这样做的另一个好处：天体撞击合并不会出现因数值过大未判断而穿过对方的情况

class _obj // 天体
{
public:
	_type type;
	double m; // 质量（视为质点）
	_3dv pos, v, a; // 位置 速度 加速度
	std::string id; // 名称（不重复，用户交互设置限制）

	// 构造函数

	_obj() : type(STAR), m(phy::SOLAR_MASS), pos(_3dv()), v(_3dv()), a(_3dv()), id("sun_default_NULL_OBJ") {} // 不要使用这种方法初始化，否则面临重名风险
	_obj(_type type, double m, const _3dv& pos, const _3dv& v, const _3dv& a, const std::string& st) : type(type), m(m), pos(pos), v(v), a(a), id(st) {}
	_obj(const _obj& o) : type(o.type), m(o.m), pos(o.pos), v(o.v), a(o.a), id(o.id) {}
	~_obj() = default;

	// 运算符重载

	_obj operator+(const _obj& o) const; // 合并天体
	friend _obj operator+=(const _obj& lhs, const _obj& rhs);
	bool operator==(const _obj& o) const; // == 判断id是否相等
	bool operator!=(const _obj& o) const;
	friend std::ostream& operator<<(std::ostream& os, const _obj& o); // 输出流

};

class _state
{
public:
	std::vector<_obj> objs;
	_obj NULL_OBJ; // 空天体
	double time;
	bool available = true, obj_still_exist = true; // 状态是否有效（天体未相撞）
	std::reference_wrapper<_obj> analyse_obj = std::ref(NULL_OBJ);

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

	size_t size() const;
	bool isexist(std::string uid) const;
	void add(const _obj& o);

	void set_a(); // 设置加速度
	void merge(); // 合并相撞天体
	void merge_compulsory();

	bool orbit_me(std::string sunid, _obj& o, double r);

	double get_energy(); // 获取系统总能量

	_obj& get_analyse_obj();

	// 运算符重载

	friend std::ostream& operator<<(std::ostream& os, _state& s); // 输出流

};

#endif