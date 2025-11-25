#include <iostream>
// #include <ctime>
#include "3dv.h"
#include "physics.h"
#include "integrator.h"
#include "user.h"
#include "test.h"

int main()
{
	// int time0 = time(NULL);
	double time0 = 0.;
	user.read_cmd(state);
	
	// test

	// state_goto(state, user.timelen, user.dt, user.method);

	while (state.time < user.timelen || user.unlimited)
	{
		if (/*time(NULL)*/state.time - time0 >= user.steps)
		{
			if (user.screen_print)
			{
				std::cout << state << std::endl;
				std::cout << _LINE_LONG << std::endl;
			}
			if (user.fout.is_open())
			{
				user.fout << state << std::endl;
				user.fout << _LINE_LONG << std::endl;
			}
			time0 = /*time(NULL)*/state.time;
		}
		integrate_dt(state, user.dt, user.method);
	}
	return 0;
}