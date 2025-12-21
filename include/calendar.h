#pragma once

#ifndef SHOW_CONSOLE
#define SHOW_CONSOLE
#endif

#ifndef _CALENDAR_H_
#define _CALENDAR_H_

#include "3dv.h"
#include "physics.h"
#include <functional>
#include <map>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

enum _era_type {STABLE, SEMI, CHAOTIC, UNCERTAIN};
// 恒纪元 过渡纪元 乱纪元 无法确定 (现天体均为质点，暂不考虑温度等因素)

using _era_pair = std::pair<int, _era_type>;
using _cnt_lf = std::pair<int, double>;

struct _era_info
{
	double time, e_orbit, leading_rate; // 结束时刻
	std::string sun_id;
};

class _calendar
{
private:
	_state NULL_STATE = _state(-1);
	_obj NO_STATE_NULL_OBJ = _obj("NO_STATE_NULL_OBJ");
	std::reference_wrapper<_state> this_state = std::ref(NULL_STATE);
	std::reference_wrapper<_obj> current_sun = std::ref(NO_STATE_NULL_OBJ);
	double leading_rate = 0; // 取值范围：[0, 1] (= (max - second_max) / (max - min))
	int sun_cnt = 0;

	bool ranked = false;
	double f_mags = 0;

	// samples（暂弃用）
	std::map<std::string, _cnt_lf> ave_energy; // 算术平均（实际存总和，取值时再除）
	std::map<std::string, _cnt_lf> var_energy; // 方差（存储同上）

	std::map<std::string, _3dv> prev_ang_mom;
	std::map<std::string, double> dtheta_ranks;


	// 尝试：采用瞬时判定方法，以下方法返回归一化评分
	double rank_v(_state& s, _obj& me, const _obj& sun);
	double rank_hradius(_state& s, _obj& me, const _obj& sun);


public:

	_calendar() : this_state(std::ref(NULL_STATE)), current_sun(std::ref(NO_STATE_NULL_OBJ)), leading_rate(0), sun_cnt(0), ranked(false), f_mags(0) {}
	_calendar(_state& s) : this_state(std::ref(s)), current_sun(std::ref(s.NO_SUN_OBJ)), leading_rate(0), sun_cnt(0), ranked(false), f_mags(0) {}
	~_calendar() = default;

	// 现放于public便于即时验证
	std::unordered_map<std::string, double> ranklist;
	std::vector<_era_info> info_list;
	std::vector<_era_pair> era_list;
	
	_state& get_state();
	_obj& get_current_sun();
	bool is_state_set();
	void set_state(_state& s);

	// 判别绕行的依据（暂弃用）

	void sample_energy(const _obj& sun); // 取样
	double rank_energy(const _obj& sun); // 多次取样后计算某天体能量之得分
	
	void sample_angmom(const _obj& sun);
	double rank_angmom(const _obj& sun);

	double rank_f_status(const _obj& sun);

	void set_f_mags();

	void sample_all();
	void rank_all(); // 评分的同时设置最佳太阳


	// 尝试：采用瞬时判定方法
	double rank_this(const _obj& sun);
	void rank_all_new();

	
	_era_type get_era();

	void helper_printrank(std::ostream& os); // 仅调试使用

};

#endif