#include "physics.h"
#include <iostream>
#include "3dv.h"
#include <vector>
#include <string>
#include <cmath>
#include <functional>
#include <algorithm>

// _obj

_obj _obj::operator+(const _obj& o) const // 合并天体
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
_obj operator+=(const _obj& lhs, const _obj& rhs)
{
	return lhs + rhs;
}
bool _obj::operator==(const _obj& o) const // == 判断id是否相等
{
	return (id == o.id);
}
bool _obj::operator!=(const _obj& o) const
{
	return !(*this == o);
}
std::ostream& operator<<(std::ostream& os, const _obj& o) // 输出流
{
	os << "[Object] \"" << o.id << "\": "
		<< std::endl << "Type = " << (o.type == STAR ? "[STAR]" : "[PLANET]")
		<< std::endl << "Mass = " << o.m
		<< std::endl << "Position = " << o.pos
		<< std::endl << "Velocity = " << o.v
		<< std::endl << "Acceleration = " << o.a << std::endl;
	return os;
}


// _state

size_t _state::size() const
{
	return objs.size();
}
bool _state::isexist(std::string id) const
{
	for (const _obj& o : objs)
	{
		if (o.id == id)
			return true;
	}
	return false;
}
void _state::add(const _obj& o)
{
	objs.push_back(o);
}

void _state::set_a() // 设置加速度
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
			if (r_ij.mag() < phy::CRASH)
			{
				available = false;
				continue; // 合体，不再计算引力
			}
			i.a += (phy::G * j.m / dist_2) * r_ij._e();
		}
	}
	return;
}

void _state::merge() // 合并相撞天体
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
				if (objs[i].pos.distance(objs[j].pos) < phy::CRASH) // 距离过近
				{
					_obj new_obj = objs[i] + objs[j];
					bool analyse_obj_tag = true;
					try
					{
						_obj rubbish = analyse_obj.get();
					}
					catch (...)
					{
						obj_still_exist = false;
						analyse_obj_tag = false;
					}
					if (analyse_obj_tag && analyse_obj.get() == NULL_OBJ)
					{
						obj_still_exist = false;
						analyse_obj_tag = false;
					}
					if (analyse_obj_tag && (analyse_obj.get() == objs[i] || analyse_obj.get() == objs[j]))
					{
						obj_still_exist = false;
					}
					objs.erase(objs.begin() + j);
					objs.erase(objs.begin() + i);
					objs.push_back(new_obj);
					if (analyse_obj_tag)
					{
						analyse_obj = std::ref(objs.back());
					}
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

void _state::merge_compulsory()
{
	available = false;
	merge();
	return;
}

bool _state::orbit_me(std::string sunid, _obj& o, double r)
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
	double v = sqrt(phy::G * sun.m / r);
	o.pos = sun.pos + _3dv(r, 0, 0);
	o.v = sun.v + _3dv(0, v, 0);
	o.a = _3dv();
	add(o);
	return true;
}

double _state::get_energy() // 获取系统总能量
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
				continue; // 防止除以零，不过这个跟合并无关
			ep -= phy::G * objs[i].m * objs[j].m / dist;
		}
	}
	return ek + ep;
}

_obj& _state::get_analyse_obj()
{
	try
	{
		_obj o = analyse_obj.get();
	}
	catch (...)
	{
		return NULL_OBJ;
	}
	return analyse_obj.get();
}

// 运算符重载

std::ostream& operator<<(std::ostream& os, _state& s) // 输出流
{
	os << "Object under analysis: " << s.get_analyse_obj().id << std::endl
		<< "State at time " << s.time << " h:" << std::endl << std::endl;
	for (const _obj& o : s.objs)
	{
		os << "  " << o << std::endl;
	}
	os << "Energy: " << s.get_energy() << std::endl;
	return os;
}


// 判断轨道运行状态的有关基础方法

double _state::get_orbit_energy(const _obj& o) // 相对于给定的analyse_obj
{
	_3dv r(0), v(0);
	_obj& anlys = get_analyse_obj();
	if (anlys == NULL_OBJ)
		return 0.;
	r = anlys.pos - o.pos;
	v = anlys.v - o.v;
	return 0.5 * anlys.m * v.mag_2() - phy::G * anlys.m * o.m / r.mag();
}
double _state::get_orbit_energy(const std::string& id) // 重载
{
	auto it = std::find(objs.begin(), objs.end(), _obj(id));
	if (it == objs.end())
		return 0.;
	return get_orbit_energy(*it);
}

_3dv _state::get_angmom(const _obj& o) // 角动量
{
	_obj& anlys = get_analyse_obj();
	if (anlys == NULL_OBJ)
		return _3dv(0);
	_3dv r = anlys.pos - o.pos;
	_3dv v = anlys.v - o.v;
	return r.cross(anlys.m * v);
}
_3dv _state::get_angmom(const std::string& id)
{
	auto it = std::find(objs.begin(), objs.end(), _obj(id));
	if (it == objs.end())
		return _3dv(0);
	return get_angmom(*it);
}

double _state::get_semi_a(const _obj& o) // 半长轴
{
	_obj& anlys = get_analyse_obj();
	if (anlys == NULL_OBJ)
		return 0.;
	return -phy::G * anlys.m * o.m / (2.0 * get_orbit_energy(o));
}
double _state::get_semi_a(const std::string& id)
{
	auto it = std::find(objs.begin(), objs.end(), _obj(id));
	if (it == objs.end())
		return 0.;
	return get_semi_a(*it);
}

_3dv _state::get_eccent(const _obj& o) // 离心率矢量
{
	_obj& anlys = get_analyse_obj();
	if (anlys == NULL_OBJ)
		return _3dv(0);
	_3dv r = anlys.pos - o.pos, v = anlys.v - o.v;
	_3dv e = v.cross(get_angmom(o)) / (phy::G * anlys.m * o.m) - r._e();
	return e;
}
_3dv _state::get_eccent(const std::string& id)
{
	auto it = std::find(objs.begin(), objs.end(), _obj(id));
	if (it == objs.end())
		return _3dv(0);
	return get_eccent(*it);
}

double _state::get_T(const _obj& o) // 周期
{
	_obj& anlys = get_analyse_obj();
	if (anlys == NULL_OBJ)
		return 0.;
	double a = get_semi_a(o);
	if (a <= 0.)
		return 0.;
	return 2.0 * phy::PI * std::sqrt(a * a * a / (phy::G * o.m));
}
double _state::get_T(const std::string& id)
{
	auto it = std::find(objs.begin(), objs.end(), _obj(id));
	if (it == objs.end())
		return 0.;
	return get_T(*it);
}

double _state::get_hill_radius(const _obj& thisobj, const _obj& thatobj) // 希尔半径
{
	if (get_analyse_obj() == NULL_OBJ)
		return 0.;
	double a = get_semi_a(thisobj);
	if (a <= 0.)
		return 0.;
	return a * (1 - get_eccent(thisobj).mag()) * std::cbrt(thisobj.m / (3.0 * thatobj.m));
}
double _state::get_hill_radius(const std::string& id1, const std::string& id2)
{
	auto it1 = std::find(objs.begin(), objs.end(), _obj(id1));
	if (it1 == objs.end())
		return 0.;
	auto it2 = std::find(objs.begin(), objs.end(), _obj(id2));
	if (it2 == objs.end())
		return 0.;
	return get_hill_radius(*it1, *it2);
}