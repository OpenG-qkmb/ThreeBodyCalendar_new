#define SHOW_CONSOLE

#include <iostream>
#include <fstream>
#include <graphics.h>
#include <ege.h>
// #include <ctime>
#include "3dv.h"
#include "physics.h"
#include "integrator.h"
#include "user.h"
#include "calendar.h"
#include "test.h"

_state state(1);
_user user;
_calendar cal(state);

int main()
{
	// int time0 = time(NULL);
	double time0 = 0., time_sample = 0.;
	user.read_cmd(state);

	if (!cal.is_state_set())
		cal.set_state(state);

	std::ifstream fin;
	std::ofstream fout;

	if (!user.finished)
	{
		fin.open("input.txt");
		if (!fin.is_open())
		{
			fout.open("error.log");
			if (!fout.is_open())
			{
				std::cerr << "[Error] Cannot open error log file." << std::endl;
				return 0;
			}
			fout << "[Error] Simulation not started by user command." << std::endl;
			fout.close();
			return 0;
		}
		std::cin.rdbuf(fin.rdbuf());
		fout.open("std_output.log");
		if (!fout.is_open())
		{
			std::cerr << "[Error] Cannot open output log file." << std::endl;
			return -1;
		}
		std::cout.rdbuf(fout.rdbuf());
		std::cerr.rdbuf(fout.rdbuf());
		user.read_cmd(state);
		if (!user.finished)
		{
			std::cerr << "[Error] Simulation not started by user command." << std::endl;
			fin.close();
			fout.close();
			return 0;
		}
	}
	
	// 以上勿改：可视化后如不显示控制台，则控制台IO重定向至文件

	// test

	if (user.display)
	{
		initgraph(user.width, user.height);
		setcaption("ThreeBodyCalendar");
		setbkcolor(BLACK);
		cleardevice();
		xyprintf(0, 0, "Simulator");
	}

	{
		_state teststate = state;
		cal.get_stdyear(teststate, user.method);
	}

	if ((!user.display) && (!user.screen_print))
	{
		std::cout << "[Info] Please wait patiently..." << std::endl;
	}

	while (user.finished && (state.get_analyse_obj() != state.NULL_OBJ) && (state.time < user.timelen || user.unlimited))
	{
		if (state.time - time0 >= cal.std_year)
		{
			(void)cal.get_era();
			if (user.screen_print)
			{
				cal.printstate(std::cout);
			}
			if (user.fout.is_open())
			{
				cal.printstate(user.fout);
			}
			time0 = state.time;
		}
		integrate_dt(state, user.dt, user.method);
		cal.rank_all_new();
		cal.pushrank();
		if (user.display)
		{
			user.show(state, cal);
		}
		// cal.reset_ranked();  // pushrank()自带reset
	}
	if (user.display)
		closegraph();
	if (fin.is_open())
		fin.close();
	if (fout.is_open())
		fout.close();
	return 0;
}