#define SHOW_CONSOLE
#include <iostream>
#include <fstream>
// #include <graphics.h>
// #include <ege.h>
// #include <ctime>
#include "3dv.h"
#include "physics.h"
#include "integrator.h"
#include "user.h"
#include "test.h"

_state state;
_user user;

int main()
{
	// int time0 = time(NULL);
	double time0 = 0.;
	user.read_cmd(state);

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

	while (user.finished && (state.get_analyse_obj() != state.NULL_OBJ) && (state.time < user.timelen || user.unlimited))
	{
		if (/*time(NULL)*/state.time - time0 >= user.steps)
		{
			if (user.screen_print)
			{
				std::cout << state << std::endl;
				for (_obj& o : state.objs)
				{
					if (o.type == PLANET || o == state.get_analyse_obj() || state.get_orbit_energy(o) > 0)
						continue;
					std::cout << "Orbiting object: " << o.id << std::endl;
					std::cout << "  Orbital energy: " << state.get_orbit_energy(o) << std::endl
						<< "  Angular momentum: " << state.get_angmom(o) << std::endl
						<< "  Semi-major axis: " << state.get_semi_a(o) << std::endl
						<< "  Eccentricity: " << state.get_eccent(o).mag() << std::endl
						<< "  Period: " << state.get_T(o) << std::endl;
				}
				std::cout << _LINE_LONG << std::endl;
			}
			if (user.fout.is_open())
			{
				user.fout << state << std::endl;
				for (_obj& o : state.objs)
				{
					if (o.type == PLANET || o == state.get_analyse_obj() || state.get_orbit_energy(o) > 0)
						continue;
					user.fout << "Orbiting object: " << o.id << std::endl;
					user.fout << "  Orbital energy: " << state.get_orbit_energy(o) << std::endl
						<< "  Angular momentum: " << state.get_angmom(o) << std::endl
						<< "  Semi-major axis: " << state.get_semi_a(o) << std::endl
						<< "  Eccentricity: " << state.get_eccent(o).mag() << std::endl
						<< "  Period: " << state.get_T(o) << std::endl;
				}
				user.fout << _LINE_LONG << std::endl;
			}
			time0 = /*time(NULL)*/state.time;
		}
		integrate_dt(state, user.dt, user.method);
	}

	/*initgraph(1024, 768);
	* setcaption("ThreeBodyCalendar");
	* setbkcolor(BLACK);
	* cleardevice();
	* xyprintf(0, 0, "ThreeBodyCalendar [Author: OpenG-qkmb]");
	* while (!(kbhit() && getch() == 27))
	* {
	* 	
	* }
	* closegraph();*/
	if (fin.is_open())
		fin.close();
	if (fout.is_open())
		fout.close();
	return 0;
}