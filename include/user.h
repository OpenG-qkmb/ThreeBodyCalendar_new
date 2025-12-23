#pragma once

#ifndef SHOW_CONSOLE
#define SHOW_CONSOLE
#endif

#ifndef _USER_H_
#define _USER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <cstdlib>
#include <ctime>
#include <random>
#include <chrono>
#include <functional>
#include <graphics.h>
#include <ege.h>
#include "3dv.h"
#include "physics.h"
#include "calendar.h"


// command format:
// 
// initialize -manual
// add <type (= star/sun)> <id> <mass> <x> <y> <z> <vx> <vy> <vz> (add sun, mass: 1 = SOLAR_MASS)
// add <type (= planet/earth)> <id> <mass> <x> <y> <z> <vx> <vy> <vz> (add planet - 1, mass: 1 = EARTH_MASS)
// add <type (= planet/earth)> <id> <mass> orbit <id_sun> <r> (add planet - 2)
// end initialize
// 
// initialize -rand
// rand <id> <mass> (type could only be STAR)
// (planet cannot be added randomly, use "add ... orbit ...")
// end initialize
// 
// initialize -stable
// import <filename> (file format should be identical to the "initialize -manual" section)
// <end initialize> (don't type this line when importing from file)
// 
// step <sample_dt (= 1h in simulation)> <steps (between every output, = 1h)>
// method <method (= euler/verlet/rk4)>
// print2screen (default) / print2file <filename>
// timelen <total_time (= 1y)>
// timelen unlimited
// analyse <id>
// display [<height (= 1024)> <width (= 768)>]
// 
// start/begin

typedef std::queue<std::string>  _Q;

class _user
{
private:
	
	const char* _INVALID = "[Error] Invalid input.";

	bool rand_init = false, is_planet_added = false;
	std::random_device rd;
	std::mt19937 gen;

	inline int mymin(const int& a, const int& b) const;
	void _clear(_Q& q);
	void set_lower(std::string& st);

	_Q _get(std::istream& s_in = std::cin);
	inline _3dv rand_v(double mark);
	_3dv rand_v_new(double mark);

	// 可视化
	double zoom = 7.3143; // 与rand的范围相应，768 / 105
	_3dv getpos_ori(const _3dv& pos);
	inline double get_z(double z);
	bool zoomed = false;
	double std_rcircle = 3;

	const COLORS STAR_COLOR = RED, PLANET_COLOR = WHITE;

	long double gettime();

public:

	bool screen_print = true;
	std::string filename = "";
	std::ofstream fout;

	double dt = 1., sample_dt = 1., steps = 24.;
	std::string method = "verlet", analyse_id = "";
	double timelen = phy::YEAR;
	bool unlimited = false, finished = false;



	bool display = false;
	int width = 1024, height = 768;
	

	void initialize(_state& state, const char& mode);
	bool setsteps(_Q& q);
	bool setmethod(_Q& q);
	bool setp(_Q& q);
	bool set_tlen(_Q& q);
	bool set_analyse_id(_state& state, _Q& q);
	bool set_display(_Q& q);
	void read_cmd(_state& state/*, bool put_prompt = true*/);

	_user()
	{
		// srand(static_cast<unsigned int>(time(NULL)));
		gen = std::mt19937(rd());
		rand_init = true;
	}
	~_user()
	{
		if (fout.is_open())
			fout.close();
	}


	// 可视化
	void show(_state& s, _calendar& cal);
};


#endif