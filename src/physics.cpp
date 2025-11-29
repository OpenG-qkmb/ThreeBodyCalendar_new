#include "physics.h"
#include <iostream>
#include "3dv.h"
#include <vector>
#include <string>
#include <cmath>

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
bool _state::isexist(std::string uid) const
{
	for (const _obj& o : objs)
	{
		if (o.id == uid)
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

// 运算符重载

std::ostream& operator<<(std::ostream& os, _state& s) // 输出流
{
	os << "State at time " << s.time << " h:" << std::endl;
	for (const _obj& o : s.objs)
	{
		os << "  " << o << std::endl;
	}
	os << "Energy: " << s.get_energy() << std::endl;
	return os;
}