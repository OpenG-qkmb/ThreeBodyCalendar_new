#include "physics.h"
#include <iostream>
#include "3dv.h"
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <functional>
#include <algorithm>
#include <omp.h>

#define TEST_SUNRANK

// 其实，三体问题也可以有稳定的周期解。
// 最有代表性的就是欧拉-拉格朗日族、8字族和布鲁克-赫农族这三族(family)解。
// https://mp.weixin.qq.com/s?__biz=MjM5NDA1Njg2MA==&mid=2652050750&idx=1&sn=f42979abb4d1a6ae71b7727b3063b009&chksm=bd6a2acd8a1da3db65943e7419b3d40b49430f16ff06e04a2368a33ec8d5aec00f84eee4a665&scene=27

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

bool _state::operator==(const _state& s)
{
	return state_ver == s.state_ver;
}
bool _state::operator!=(const _state& s)
{
	return !(*this == s);
}

size_t _state::size() const
{
	return objs.size();
}
bool _state::isexist(std::string id) const
{
	auto it = std::find(objs.begin(), objs.end(), _obj(id));
	if (it == objs.end())
		return false;
	return true;
}
_obj& _state::get_obj_by_id(const std::string& id)
{
	auto it = std::find(objs.begin(), objs.end(), _obj(id));
	if (it == objs.end())
		return NULL_OBJ;
	return *it;
}
void _state::add(const _obj& o)
{
	objs.push_back(o);
}

void _state::set_a() // 设置加速度
{
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < objs.size(); ++i)
	{
#pragma omp critical
		{
			objs[i].a = _3dv(0);
#pragma omp	parallel for reduction(+ : objs[i].a)
			for (int j = 0; j < objs.size(); ++j)
			{
				_3dv add = _3dv(0);
				if (i != j)
				{
					_3dv r_ij = objs[j].pos - objs[i].pos;
					double dist_2 = r_ij.mag_2();
					if (r_ij.mag() < phy::CRASH)
					{
						available = false;
						// 合体，不再计算引力
					}
					else
					{
						add = (phy::G * objs[j].m / dist_2) * r_ij._e();
					}
				}
				objs[i].a += add;
			}
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
#pragma omp parallel for schedule(dynamic)
		for (size_t i = 0; i < objs.size(); ++i)
		{
			for (size_t j = i + 1; j < objs.size(); ++j)
			{
#pragma omp critical
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
	auto it = std::find(objs.begin(), objs.end(), _obj(sunid));
	if (it == objs.end())
		return false;
	sun = *it;
	if (sun.type == PLANET)
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
	merge(); // 合并，防止除以零
	double ek = 0.0, ep = 0.0;
#pragma omp parallel for reduction (+ : ek)
	for (size_t i = 0; i < objs.size(); ++i)
	{
		double add0 = 0.5 * objs[i].m * objs[i].v.mag_2();
		ek += add0;
#pragma omp parallel for reduction (- : ep)
		for (size_t j = i + 1; j < objs.size(); ++j)
		{
			double dist = objs[i].pos.distance(objs[j].pos);
			double add = 0;
			if (dist >= SMALL_NUM) // 防止除以零，不过这个跟合并无关
				add = phy::G * objs[i].m * objs[j].m / dist;
			ep -= add;
		}
	}
	return ek + ep;
}

_obj& _state::get_analyse_obj()
{
	try
	{
		_obj& o = analyse_obj.get();
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
	// os << "Energy: " << s.get_energy() << std::endl; // 能量守恒已验证
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

//double _state::rank_sun(const _obj& o, std::map<std::string, bool>& bool_list) // 为捕获的强度打分
//{
//	if (get_analyse_obj() == NULL_OBJ)
//		return -8000.;
//	_obj anlys = get_analyse_obj();
//
//	if (o == anlys)
//		return -4000.;
//
//	double totrank = 0.;
//	double orbit_energy, cos_angle;
//
//	bool_list.clear();
//
//#ifdef TEST_SUNRANK
//	std::cout << std::endl << "[Debug] Object " << o.id << ":" << std::endl;
//#endif
//
//	// 以下权重均有待验证
//
//	// 轨道能量（0.2）
//	orbit_energy = get_orbit_energy(o);
//	if (orbit_energy >= 0.)
//	{
//		// totrank -= 2500.;
//		totrank -= orbit_energy * 1e7 * 0.75;
//#ifdef TEST_SUNRANK
//		std::cout << "[Debug] Rank of orbit energy: " << -orbit_energy * 1e7 * 0.75 << std::endl;
//#endif
//	}
//	else
//	{
//		totrank -= orbit_energy * 1e7 * 0.75; // 经验之谈：地球绕太阳此项数值大约为2.5e-7，下同
//#ifdef TEST_SUNRANK
//		std::cout << "[Debug] Rank of orbit energy: " << -orbit_energy * 1e7 * 0.75 << std::endl;
//#endif
//	}
//
//	// 角动量与位矢应几近垂直（0.3）
//	cos_angle = (anlys.pos - o.pos).cos_angle_with(get_angmom(o)); // 1e-25
//	totrank -= cos_angle * 1e25 * 3; // 不垂直则扣分
//#ifdef TEST_SUNRANK
//	std::cout << "[Debug] Rank of cos_angle: " << -cos_angle * 1e25 * 3 << std::endl;
//#endif
//
//	// 希尔球（0.5）
//
//	for (const _obj& otherobj : objs)
//	{
//		if (otherobj == o || otherobj == anlys)
//			continue;
//		double hill_radius = get_hill_radius(o, otherobj);
//		bool hill_ball_all = true; // 对所有其他天体，是否均在该恒星希尔球内
//
//		if (hill_radius <= 0.)
//		{
//#ifndef TEST_SUNRANK
//			bool_list.insert(std::pair<std::string, bool>(o.id, false));
//			totrank -= 1;
//			continue;
//#else
//			std::cout << "[Debug] Hill sphere: Main=" << o.id << "; Other=" << otherobj.id << "; Radius<0: " << hill_radius << std::endl;
//#endif
//		}
//
//		if (anlys.pos.distance(o.pos) > hill_radius) // 在本恒星希尔球外
//		{
//			hill_ball_all = false;
//			totrank -= 0.5 / ((size() > 2) ? (static_cast<double>(size() - 2)) : 1.);
//#ifdef TEST_SUNRANK
//			std::cout << "[Debug] Hill sphere: Main=" << o.id << "; Other=" << otherobj.id << "; deltaR=" << -0.5 / ((size() > 2) ? (static_cast<double>(size() - 2)) : 1.) << std::endl;
//#endif
//		}
//		else
//		{
//#ifdef TEST_SUNRANK
//			std::cout << "[Debug] Hill sphere: Main=" << o.id << "; Other=" << otherobj.id << "; ===IN===" << std::endl;
//#endif
//		}
//		if (hill_ball_all)
//			bool_list.insert(std::pair<std::string, bool>(o.id, true));
//		else
//			bool_list.insert(std::pair<std::string, bool>(o.id, false));
//	}
//#ifdef TEST_SUNRANK
//	std::cout << "[Debug] Total rank for object " << o.id << ": " << totrank << std::endl;
//#endif
//	return totrank;
//}
//
//_obj& _state::get_sun()
//{
//	if (get_analyse_obj() == NULL_OBJ)
//		return NULL_OBJ;
//	_obj anlys = get_analyse_obj();
//	
//	std::map<std::string, bool> is_in_hill_ball;
//	std::map<std::string, double>ranks;
//	for (const _obj& o : objs)
//	{
//		if (o == anlys)
//			continue;
//		ranks.insert(std::pair<std::string, double>(o.id, rank_sun(o, is_in_hill_ball)));
//	}
//	return NO_SUN_OBJ;
//}