#pragma once

#ifndef SHOW_CONSOLE
#define SHOW_CONSOLE
#endif

#ifndef _CALENDAR_H_
#define _CALENDAR_H_

#include "3dv.h"
#include "physics.h"
#include "integrator.h"
#include "test.h"
#include <functional>
#include <map>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

enum _era_type {STABLE/*, SEMI*/, CHAOTIC, UNCERTAIN};
// 恒纪元 过渡纪元 乱纪元 无法确定 (现天体均为质点，暂不考虑温度等因素)
// 决定：弃用SEMI（过渡纪元）但凡不稳就是乱

using _era_pair = std::pair<int, _era_type>;
using _cnt_lf = std::pair<int, double>;

struct _ranks
{
	double v, hill, ecc;

	_ranks operator+(const _ranks& r) const;
	_ranks& operator+=(const _ranks& r);
	_ranks operator-(const _ranks& r) const;
	_ranks& operator-=(const _ranks& r);
	_ranks operator*(const _ranks& r) const;
	_ranks& operator*=(const _ranks& r);
	_ranks& operator*=(double val);
	_ranks operator*(double val) const;
	friend _ranks operator*(double val, const _ranks& r);
	friend std::ostream& operator<<(std::ostream& os, const _ranks& rk);

	double sum() const;
	double sum_halfecc() const;
};

struct _era_info
{
	double time = 0; // 结束时刻
	_ranks ranks = _ranks{ 0, 0, 0 };
	std::string sun_id;
};

namespace gyear
{
	constexpr double dt = 1.;
	constexpr double timelen = 1500;
	constexpr double piece = 500;
}

class _calendar
{
private:
	_state NULL_STATE = _state(-1);
	_obj NO_STATE_NULL_OBJ = _obj("NO_STATE_NULL_OBJ");
	std::reference_wrapper<_state> this_state = std::ref(NULL_STATE);
	std::reference_wrapper<_state> state_getyear = std::ref(NULL_STATE);
	std::reference_wrapper<_state> state_true = std::ref(NULL_STATE);
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
	double rank_e(_state& s, _obj& me, const _obj& sun);



	void helper_printrank(std::ostream& os); // 仅调试使用

public:

	_calendar() : state_true(std::ref(NULL_STATE)), this_state(std::ref(NULL_STATE)), current_sun(std::ref(NO_STATE_NULL_OBJ)), leading_rate(0), sun_cnt(0), ranked(false), f_mags(0) {}
	_calendar(_state& s) : state_true(std::ref(s)), this_state(std::ref(s)), current_sun(std::ref(s.NO_SUN_OBJ)), leading_rate(0), sun_cnt(0), ranked(false), f_mags(0) {}
	~_calendar() = default;

	// 现放于public便于即时验证
	std::unordered_map<std::string, double> ranklist; // 弃用
	std::unordered_map<std::string, _ranks> ranklist_ranks; // 现用
	std::vector<_era_info> periodlist; // 现用 存储一个时间段（例如一年）内的“太阳成绩”
	std::vector<_era_pair> era_list; // 弃用
	
	_state& get_state();
	_obj& get_current_sun();
	bool is_state_set();
	void set_state(_state& s);

	double std_year = phy::YEAR;
	_era_type current_era = UNCERTAIN;
	std::string era_name();

	// 判别绕行的依据（暂弃用）

	void sample_energy(_state& s, const _obj& sun); // 取样
	double rank_energy(_state& s, const _obj& sun); // 多次取样后计算某天体能量之得分
	
	void sample_angmom(_state& s, const _obj& sun);
	double rank_angmom(_state& s, const _obj& sun);

	double rank_f_status(_state& s, const _obj& sun);

	void set_f_mags(_state& s);

	void sample_all(_state& s);
	void rank_all(_state& s); // 评分的同时设置最佳太阳


	// 尝试：采用瞬时判定方法
	_ranks rank_this(const _obj& sun);
	void rank_all_new();
	void pushrank();
	void reset_ranked();

	
	_era_type get_era();

	void get_stdyear(_state& s, const std::string& method);

	void printstate(std::ostream& os);

};

#endif