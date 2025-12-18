#include "3dv.h"
#include "calendar.h"
#include "physics.h"
#include <functional>
#include <map>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <omp.h>


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
	this_state = std::ref(s);
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


void _calendar::sample_energy(const _obj& sun)
{
	_state& s = get_state();
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
	_var.second += energy * energy;
	return;
}
double _calendar::rank_energy(const _obj& sun)
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

void _calendar::sample_angmom(const _obj& sun)
{
	_state& s = get_state();
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
	_3dv& angmom = s.get_angmom(sun);
	double dtheta = angmom.cos_angle_with(prev);
	prev = angmom;
	if (dtheta > 0.9)
		dtheta_ranks[sun.id] += 0.1;
	else if (dtheta > 0.7)
		dtheta_ranks[sun.id] += 0.05;
	return;
}
double _calendar::rank_angmom(const _obj& sun)
{
	return dtheta_ranks[sun.id];
}

double _calendar::rank_f_status(const _obj& sun)
{
	_state& s = get_state();
	if (s == NULL_STATE)
		return 0;
	if (sun == s.NULL_OBJ || sun.type != STAR)
		return 0;
	_obj& o = s.get_analyse_obj();
	double dist2 = o.pos.distance(sun.pos);
	if (o.pos.distance(sun.pos) < phy::CRASH)
		dist2 += SMALL_NUM;
	return ((phy::G * sun.m / dist2) * (o.pos - sun.pos)._e())._e().dot(o.a._e());
}

void _calendar::sample_all()
{
	_state& s = get_state();
	if (s == NULL_STATE)
		return;
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < s.objs.size(); ++i)
	{
		if (s.objs[i].type != STAR || s.objs[i] == s.NULL_OBJ)
			continue;
#pragma omp critical
		{
			sample_energy(s.objs[i]);
			sample_angmom(s.objs[i]);
		}
	}
	return;
}

void _calendar::rank_all()
{
	_state& s = get_state();
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
			thisrank += rank_energy(s.objs[i]) * 36; // 此数字为实验调整得到，或有更优值，下同
			thisrank += rank_angmom(s.objs[i]) * 0.003;
			thisrank += rank_f_status(s.objs[i]) * 0.3;
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
	return;
}



void _calendar::helper_printrank(std::ostream& os)
{
	if (ranklist.empty())
		return;
	os << "Among the " << sun_cnt << " sun(s):" << std::endl;
	for (auto it = ranklist.begin(); it != ranklist.end(); ++it)
	{
		os << "Rank of star [" << it->first << "] = " << it->second << std::endl;
	}
	os << "Current sun: " << get_current_sun().id  << " (rate = " << leading_rate << ")" << std::endl;
	return;
}