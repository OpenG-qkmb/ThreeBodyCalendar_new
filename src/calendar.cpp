#include "3dv.h"
#include "calendar.h"
#include "physics.h"
#include "integrator.h"
#include "test.h"
#include <functional>
#include <map>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <graphics.h>
#include <ege.h>
#include <omp.h>


_ranks _ranks::operator+(const _ranks& r) const
{
	return _ranks{ v + r.v, hill + r.hill, ecc + r.ecc };
}
_ranks& _ranks::operator+=(const _ranks& r)
{
	v += r.v;
	hill += r.hill;
	ecc += r.ecc;
	return *this;
}
_ranks _ranks::operator-(const _ranks& r) const
{
	return _ranks{ v - r.v, hill - r.hill, ecc - r.ecc };
}
_ranks& _ranks::operator-=(const _ranks& r)
{
	v -= r.v;
	hill -= r.hill;
	ecc -= r.ecc;
	return *this;
}
_ranks _ranks::operator*(const _ranks& r) const
{
	return _ranks{ v * r.v, hill * r.hill, ecc * r.ecc };
}
_ranks& _ranks::operator*=(const _ranks& r)
{
	v *= r.v;
	hill *= r.hill;
	ecc *= r.ecc;
	return *this;
}
_ranks& _ranks::operator*=(double val)
{
	v *= val;
	hill *= val;
	ecc *= val;
	return *this;
}
_ranks _ranks::operator*(double val) const
{
	return _ranks{ v * val, hill * val, ecc * val };
}
_ranks operator*(double val, const _ranks& r)
{
	return r * val;
}
std::ostream& operator<<(std::ostream& os, const _ranks& rk)
{
	os << "V: " << rk.v << ", HILL: " << rk.hill << ", ECC: " << rk.ecc;
	return os;
}
double _ranks::sum() const
{
	return v + hill + ecc;
}
double _ranks::sum_halfecc() const
{
	return v + hill + ecc * 0.5;
}

_state& _calendar::get_state()
{
	try
	{
		_state& s = this_state.get();
	}
	catch (...)
	{
		return NULL_STATE;
	}
	return this_state.get();
}
bool _calendar::is_state_set()
{
	return get_state() != NULL_STATE;
}
void _calendar::set_state(_state& s)
{
	state_true = this_state = std::ref(s);
	current_sun = std::ref(s.NO_SUN_OBJ);
	return;
}
_obj& _calendar::get_current_sun()
{
	try
	{
		_obj& o = current_sun.get();
	}
	catch (...)
	{
		return get_state().NULL_OBJ;
	}
	return current_sun.get();
}

std::string _calendar::era_name()
{
	std::string res;
	switch (current_era)
	{
	case STABLE: res = "Stable"; break;
	case CHAOTIC: res = "Chaotic"; break;
	case UNCERTAIN:
	default: res = "Not yet determined"; break;
	}
	return res;
}


void _calendar::sample_energy(_state& s, const _obj& sun)
{
	if (s == NULL_STATE)
		return;
	if (sun == s.NULL_OBJ || sun.type != STAR)
		return;
	double energy = s.get_orbit_energy(sun);
	_cnt_lf& _ave = ave_energy[sun.id];
	_cnt_lf& _var = var_energy[sun.id];
	++_ave.first;
	_ave.second += energy;
	++_var.first;
	_var.second += (energy - _ave.second) * (energy - _ave.second);
	return;
}
double _calendar::rank_energy(_state& s, const _obj& sun)
{
	_cnt_lf &_var = var_energy[sun.id], &_ave = ave_energy[sun.id];
	if (_var.first == 0 || _ave.first == 0)
		return 0.;
	double stdvar = std::sqrt(_var.second / _var.first), ave = _ave.second / _ave.first, sgn = -1; // 直接改用标准差
	if (ave < 0.)
	{
		ave = -ave;
		sgn = 1;
	}
	if (ave < SMALL_NUM)
		ave += SMALL_NUM; // 防止除以零
	double frac = stdvar * 10. + ave;
	if (std::abs(frac) < SMALL_NUM * 100)
		frac += SMALL_NUM * 100;
	return ave * sgn / frac; // 能量为负说明行星捕获较强，故添加符号
}

void _calendar::sample_angmom(_state& s, const _obj& sun)
{
	if (s == NULL_STATE)
		return;
	if (sun == s.NULL_OBJ || sun.type != STAR)
		return;
	_3dv& prev = prev_ang_mom[sun.id];
	if (prev.is_zero())
	{
		prev = s.get_angmom(sun);
		dtheta_ranks[sun.id] = 0;
		return;
	}
	_3dv&& angmom = s.get_angmom(sun);
	double dtheta = angmom.cos_angle_with(prev);
	prev = angmom;
	if (dtheta > 0.9)
		dtheta_ranks[sun.id] += 0.1;
	else if (dtheta > 0.7)
		dtheta_ranks[sun.id] += 0.05;
	return;
}
double _calendar::rank_angmom(_state& s, const _obj& sun)
{
	return dtheta_ranks[sun.id];
}

double _calendar::rank_f_status(_state& s, const _obj& sun)
{
	if (s == NULL_STATE)
		return 0;
	if (sun == s.NULL_OBJ || sun.type != STAR)
		return 0;
	set_f_mags(s);
	double fm = (f_mags > SMALL_NUM) ? f_mags : SMALL_NUM;
	return s.force_with(sun).mag() / fm;
}

void _calendar::sample_all(_state& s)
{
	if (s == NULL_STATE)
		return;
	ranked = false;
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < s.objs.size(); ++i)
	{
		if (s.objs[i].type != STAR || s.objs[i] == s.NULL_OBJ)
			continue;
#pragma omp critical
		{
			sample_energy(s, s.objs[i]);
			sample_angmom(s, s.objs[i]);
		}
	}
	return;
}

void _calendar::set_f_mags(_state& s)
{
	f_mags = 0;
	if (s == NULL_STATE)
		return;
#pragma omp parallel for reduction(+ : f_mags)
	for (int i = 0; i < s.objs.size(); ++i)
	{
		if (s.objs[i].type != STAR || s.objs[i] == s.NULL_OBJ)
			continue;
		f_mags += s.force_with(s.objs[i]).mag();
	}
	return;
}

void _calendar::rank_all(_state& s)
{
	if (s == NULL_STATE)
		return;
	s.set_a();
	ranklist.clear();
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < s.objs.size(); ++i)
	{
		if (s.objs[i].type != STAR || s.objs[i] == s.NULL_OBJ)
			continue;
#pragma omp critical
		{
			double thisrank = 0;
			thisrank += rank_energy(s, s.objs[i]) * 36; // 此数字为实验调整得到，或有更优值，下同
			// thisrank += rank_angmom(s.objs[i]) * 0.003; // 不计入评分：改为（稳定？）判断
			thisrank += rank_f_status(s, s.objs[i]) * 0.5;
			ranklist[s.objs[i].id] = thisrank;
		}
	}
	if (ranklist.empty())
	{
		current_sun = std::ref(s.NO_SUN_OBJ);
		sun_cnt = 0;
		leading_rate = 0;
	}
	else
	{
		auto best_sun = ranklist.begin();
		double mrank = -1e9, sec_rank = -2e9, minr = 1e9;
		bool sec_init = false;
		sun_cnt = 0; // 每次重新计算，以防天体合并
		for (auto it = ranklist.begin(); it != ranklist.end(); ++it)
		{
			if (it->second > mrank)
			{
				mrank = it->second;
				best_sun = it;
			}
			else if (it->second > sec_rank)
			{
				sec_rank = it->second;
				sec_init = true;
			}
			if (it->second < minr)
				minr = it->second;
			++sun_cnt;
		}
		if (!sec_init)
			sec_rank = minr;
		current_sun = std::ref(s.get_obj_by_id(best_sun->first));
		if ((mrank - minr) > SMALL_NUM)
			leading_rate = (mrank - sec_rank) / (mrank - minr);
		else if (sun_cnt == 1)
			leading_rate = 1;
		else
			leading_rate = 0;
	}
	var_energy.clear();
	ave_energy.clear();
	dtheta_ranks.clear(); // 不要清除前面记录的上一步角动量
	ranked = true;
	return;
}

// 从实践的角度看：行星被恒星稳定捕获时几乎可完全根据能量判断；
// 反之，所有轨道能量均非负时，行星“乱飘”，此时几乎无有效判据，直接令NO_SUN_OBJ即可


void _calendar::helper_printrank(std::ostream& os)
{
	if (ranklist_ranks.empty())
		return;
	os << "Among the " << sun_cnt << " sun(s):" << std::endl;
	for (auto it = ranklist_ranks.begin(); it != ranklist_ranks.end(); ++it)
	{
		os << "Rank of star [" << it->first << "] = " << it->second << std::endl;
	}
	os << "Current sun: " << get_current_sun().id/*  << " (rate = " << leading_rate << ")"*/ << std::endl;
	return;
}


// 以上rank的函数全部暂弃用（hepler_printrank除外）


double _calendar::rank_v(_state& s, _obj& me, const _obj& sun)
{
	double dist = (me.pos - sun.pos).mag(), v = (me.v - sun.v).mag();
	double v_esc = std::sqrt(2 * phy::G * sun.m / dist);
	if (v > v_esc)
		return 0;
	return (1 - (v * v / (v_esc * v_esc))); // 1 - (v / v_esc) ^ 2
}
double _calendar::rank_hradius(_state& s, _obj& me, const _obj& sun)
{
	double dist = (me.pos - sun.pos).mag(), min_hillr = INFINITY, max_hillr = -1e9;
	int sun_cnt = 1;
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < s.objs.size(); ++i)
	{
		if (s.objs[i] == sun || s.objs[i].type != STAR || s.objs[i] == s.NULL_OBJ)
			continue;
		double hill_r = s.get_hill_radius(sun, s.objs[i]);
		if (hill_r > max_hillr)
			max_hillr = hill_r;
		if (hill_r < min_hillr)
			min_hillr = hill_r;
		++sun_cnt;
	}
	if (dist <= (min_hillr + SMALL_NUM))
		return 1;
	if (dist >= (max_hillr - SMALL_NUM))
		return 0;
	return ((dist - min_hillr) / (max_hillr - min_hillr));
}
double _calendar::rank_e(_state& s, _obj& me, const _obj& sun)
{
	double ecc = s.get_eccent(sun).mag();
	if (ecc >= 1 - SMALL_NUM)
		return 0;
	return (1 - ecc) * (1 - ecc) * (1 + 2 * ecc);
}

_ranks _calendar::rank_this(const _obj& sun)
{
	_state& s = get_state();
	if (s == NULL_STATE)
		return _ranks{ 0, 0, 0 };
	if (sun == s.NULL_OBJ || sun == s.NO_SUN_OBJ)
		return _ranks{ 0, 0, 0 };
	_obj& me = s.get_analyse_obj();
	if (me == s.NULL_OBJ)
		return _ranks{ 0, 0, 0 };
	double rank = 0;
	return _ranks{ rank_v(s, me, sun), rank_hradius(s, me, sun), rank_e(s, me, sun)};
}

void _calendar::rank_all_new()
{
	_state& s = get_state();
	if (s == NULL_STATE)
		return;
	ranklist_ranks.clear();
	sun_cnt = 0;
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < s.objs.size(); ++i)
	{
		if (s.objs[i].type != STAR || s.objs[i] == s.NULL_OBJ)
			continue;
#pragma omp critical
		{
			ranklist_ranks[s.objs[i].id] = rank_this(s.objs[i]);
			++sun_cnt;
		}
	}
	if (ranklist_ranks.empty())
	{
		current_sun = std::ref(s.NO_SUN_OBJ);
	}
	else
	{
		auto best_sun = ranklist_ranks.begin();
		double mrank = -1e9;
		for (auto it = ranklist_ranks.begin(); it != ranklist_ranks.end(); ++it)
		{
			if ((it->second.hill + it->second.v + it->second.ecc * 0.5) > mrank)
			{
				mrank = it->second.hill + it->second.v + it->second.ecc * 0.5;
				best_sun = it;
			}
		}
		if (mrank < 0.05)
		{
			current_sun = std::ref(s.NO_SUN_OBJ);
			ranklist_ranks[s.NO_SUN_OBJ.id] = _ranks{ 0, 0, 0 };
		}
		else
		{
			current_sun = std::ref(s.get_obj_by_id(best_sun->first));
		}
	}
	ranked = true;
	return;
}
void _calendar::reset_ranked()
{
	ranked = false;
	return;
}
void _calendar::pushrank()
{
	if (!ranked)
		rank_all_new();
	_obj& sun = get_current_sun();
	_era_info tmp;
	tmp.time = get_state().time;
	tmp.sun_id = sun.id;
	tmp.ranks = ranklist_ranks[sun.id];
	periodlist.push_back(tmp);
	ranked = false;
	return;
}

_era_type _calendar::get_era()
{
	if (periodlist.empty())
	{
		current_era = UNCERTAIN;
		return UNCERTAIN;
	}
	
	_ranks ave = _ranks{ 0, 0, 0 }, var = _ranks{ 0, 0, 0 }; // 均值 方差
	bool alter = false;
	std::string prev_id;

	prev_id = periodlist[0].sun_id;
	int n = periodlist.size();

#pragma omp parallel for reduction(+ : ave)
	for (int i = 0; i < n; ++i)
	{
		ave += periodlist[i].ranks;
		if (prev_id != periodlist[i].sun_id)
		{
			alter = true;
			break;
		}
		prev_id = periodlist[i].sun_id;
	}

	if (alter)
	{
		current_era = CHAOTIC;
		return CHAOTIC;
	}

	double scal = 1. / static_cast<double>(n);
	ave *= scal;

#pragma omp parallel for reduction(+ : var)
	for (int i = 0; i < n; ++i)
	{
		_ranks rk = ave - periodlist[i].ranks;
		var += rk * rk;
	}

	var *= scal;

	// 唉，经验数值，其实从经验上根本没什么好用的经验
	// 恒纪元其实不需要这么复杂的判断，乱纪元判断了也没必要
	if (ave.sum_halfecc() < 0.45 || ave.ecc < 5e-3)
	{
		current_era = CHAOTIC;
		return CHAOTIC;
	}
	if (var.sum_halfecc() > 1.25 || var.v > 0.25 || var.hill > 0.25 || var.ecc > 0.25)
	{
		current_era = CHAOTIC;
		return CHAOTIC;
	}
	current_era = STABLE;
	ranklist_ranks.clear();
	periodlist.clear();
	return STABLE;
}

void _calendar::get_stdyear(_state& s, const std::string& method)
{
	double time0 = 0;
	state_true = this_state;
	state_getyear = std::ref(s);
	this_state = state_getyear;
	periodlist.clear();
	std_year = phy::YEAR;
	while (s.time < gyear::timelen)
	{
		if (s.time - time0 >= gyear::piece)
		{
			_era_type current_type = get_era();
			_obj& cur_sun = get_current_sun();
			if (current_type == STABLE && cur_sun != s.NO_SUN_OBJ && cur_sun != s.NULL_OBJ)
			{
				std_year = s.get_T(cur_sun);
				break;
			}
		}
		integrate_dt(s, gyear::dt, method);
		rank_all_new();
		pushrank();
	}
	this_state = state_true;
	state_getyear = std::ref(NULL_STATE);
	current_era = UNCERTAIN;
	periodlist.clear();
	ranklist_ranks.clear();
	reset_ranked();
	return;
}

void _calendar::printstate(std::ostream& os)
{
	_state& s = get_state();
	_obj& cur_sun = get_current_sun();
	os << s << std::endl;
	os << "Current sun: " << cur_sun.id << std::endl;
	os << "Year " << static_cast<int>(s.time / std_year) << ": " << era_name() << std::endl;
	os << _LINE_LONG << std::endl;
	return;
}