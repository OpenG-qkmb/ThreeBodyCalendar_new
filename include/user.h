#pragma once

#ifndef _USER_H_
#define _USER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <cstdlib>
#include <ctime>
#include "3dv.h"
#include "physics.h"
#include "integrator.h"


// command format:
// 
// initialize -manual
// add <type (= star/sun)> <id> <mass> <x> <y> <z> <vx> <vy> <vz> (add sun)
// add <type (= planet/earth)> <id> <mass> <x> <y> <z> <vx> <vy> <vz> (add planet - 1)
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
// auto (default)
// step <dt (= 3600s in simulation)> <steps (between every output, = 3s in real world)>
// method <method (= euler/verlet/rk4)>
// print2screen (default) / print2file <filename> (default: output_<time>.txt)
// timelen <total_time (= 1y)>
// timelen unlimited
// 
// start/begin

//class _prompt
//{
//public:
//	const std::map<int, std::string> PROMPT_LIST = { // int p_id, std::string prompt
//		{-1, "Invalid input"},
//		{0, "[Error] "},
//		{1, "Initialize and start to simulate. See README.md for help."},
//		{2, "Initialization completed."},
//		{3, "wtf"}
//	};
//};

typedef std::queue<std::string>  _Q;

class _user
{
private:
	
	const char* _INVALID = "[Error] Invalid input.";

	inline int mymin(const int& a, const int& b) const
	{
		return (a < b) ? a : b;
	}

	void _clear(_Q& q)
	{
		while (!q.empty())
			q.pop();
		return;
	}
	
	void set_lower(std::string& st)
	{
		for (int i = 0; i < st.size(); ++i)
		{
			if (st[i] >= 'A' && st[i] <= 'Z')
				st[i] += 'a' - 'A';
		}
		return;
	}

	_Q _get(std::istream& s_in = std::cin)
	{
		std::string cmd, single;
	read:
		if (!std::getline(s_in, cmd))
		{
			if (s_in.eof())
				return _Q();
			if (s_in.fail())
			{
				std::cerr << "[Error] Input read operation failed." << std::endl;
				return _Q();
			}
			std::cerr << "[Error] Unexpected error." << std::endl;
			return _Q();
		}
		if (cmd.empty())
			goto read;
		_Q q = {};
		set_lower(cmd);
		// _clear(q);
		std::stringstream ss(cmd);
		while (ss >> single)
		{
			q.push(single);
		}
		return q;
	}

	inline _3dv rand_v(double mark)
	{
		return _3dv(
			static_cast<double>(rand()) * mark / RAND_MAX,
			static_cast<double>(rand()) * mark / RAND_MAX,
			static_cast<double>(rand()) * mark / RAND_MAX
		);
	}

public:

	bool screen_print = true;
	std::string filename = "";
	std::ofstream fout;

	double dt = 3600, steps = 3;
	std::string method = "verlet";
	double timelen = STD_YEAR * DAY_SECONDS;
	bool unlimited = false;

	
	void initialize(_state& state, const char& mode)
	{
		_Q q;
		std::ifstream fin;
		std::string filename;
		const bool s = (mode == 's');
		if (s)
		{
		read_filename:
			do
			{
				q = _get();
				if (q.front().front() != 'i' || q.size() < 2)
				{
					std::cerr << _INVALID << std::endl;
				}
				else
				{
					q.pop();
					break;
				}
			} while (true);
			while (!q.empty())
			{
				filename += q.front() + " ";
				q.pop();
			}
			
			try
			{
				filename.pop_back();
				fin.open(filename);
			}
			catch (...)
			{
				std::cerr << "[Error] Cannot open the file." << std::endl;
				filename.clear();
				goto read_filename;
			}
			if(!fin.is_open())
			{
				std::cerr << "[Error] No such file." << std::endl;
				filename.clear();
				goto read_filename;
			}
		}
		do
		{
			bool isrand = false;
			_obj sample;
			_clear(q);
			if (s)
				q = _get(fin);
			else
				q = _get();
			if (q.front().front() == 'e' || q.empty())
				break;
			if (q.front().front() == 'r')
			{
				isrand = true;
				if (q.size() < 3)
				{
					std::cerr << _INVALID << std::endl;
					continue;
				}
				q.pop();
			}

			if (isrand)
				goto skip_id;

			if (q.front().front() != 'a' || q.size() < 7)
			{
				std::cerr << _INVALID << std::endl;
				continue;
			}

			q.pop();
			
			// select type

			switch (q.front().front())
			{
			case 'p':
			case 'e': sample.type = PLANET; break;
			default: sample.type = STAR;
			}
			q.pop();

			// set id

		skip_id:

			sample.id = q.front();
			if (state.isexist(sample.id))
			{
				std::cerr << "[Error] ID must be unique." << std::endl;
				continue;
			}
			q.pop();

			// mass

			// sample.m = std::stod(q.front());
			try
			{
				sample.m = std::stod(q.front());
			}
			catch (...)
			{
				std::cerr << _INVALID << std::endl;
				continue;
			}
			q.pop();

			if (isrand)
			{
				sample.pos = rand_v(AU * 100);
				sample.v = rand_v(V_SUN);
				sample.a = _3dv();
				state.add(sample);
				continue;
			}

			// 自此起分叉 orbit单独处理 否则设置三维向量

			if (q.front().front() == 'o')
			{
				q.pop();
				std::string sunid = q.front();
				double r = AU;
				if (!state.isexist(sunid))
				{
					std::cerr << "[Error] Cannot find the object." << std::endl;
					continue;
				}
				q.pop();
				try
				{
					r = std::stod(q.front());
				}
				catch (...)
				{
					std::cerr << _INVALID << std::endl;
					continue;
				}
				q.pop();
				if (!state.orbit_me(sunid, sample, r))
				{
					std::cerr << "[Error] Planet should orbit around stars." << std::endl;
					continue;
				}
			}
			else
			{
				if (q.size() < 6)
				{
					std::cerr << _INVALID << std::endl;
					continue;
				}
				double xyz[6] = {0}, tmp = .0; // 这里的6是无可辩驳的6，无需单独设作常量
				bool flag = true;
				for (int i = 0; i < 6; ++i)
				{
					try
					{
						tmp = std::stod(q.front());
					}
					catch (...)
					{
						std::cerr << _INVALID << std::endl;
						flag = false;
						continue;
					}
					xyz[i] = tmp;
					q.pop();
				}
				if (!flag)
					continue;
				sample.pos = _3dv(xyz[0], xyz[1], xyz[2]);
				sample.v = _3dv(xyz[3], xyz[4], xyz[5]);
				sample.a = _3dv();
				state.add(sample);
			}
		} while (true);
		if (fin.is_open())
			fin.close();
		state.merge_compulsory();
		return;
	}

	bool setsteps(_Q& q)
	{
		if (q.size() < 2)
			return false;
		try
		{
			dt = std::stod(q.front());
		}
		catch (...)
		{
			return false;
		}
		q.pop();
		try
		{
			steps = std::stod(q.front());
		}
		catch (...)
		{
			return false;
		}
		q.pop();
		return true;
	}

	bool setmethod(_Q& q)
	{
		if (q.empty())
			return false;
		switch (q.front().front())
		{
		case 'e':
		case 'v':
		case 'r': break;
		default: return false;
		}
		method = q.front();
		q.pop();
		return true;
	}

	bool setp(_Q& q)
	{
		if (q.empty())
			return false;
		if (q.front().back() == 'e')
		{
			q.pop();
			if (q.empty())
				return false;
			filename.clear();
			while (!q.empty())
			{
				filename += q.front() + " ";
				q.pop();
			}
			try
			{
				filename.pop_back();
				fout.open(filename);
			}
			catch (...)
			{
				return false;
			}
			screen_print = false;
		}
		_clear(q);
		return true;
	}

	bool set_tlen(_Q& q)
	{
		if (q.empty())
			return false;
		if (q.front().front() == 'u')
		{
			unlimited = true;
			return true;
		}
		try
		{
			timelen = std::stod(q.front());
		}
		catch (...)
		{
			return false;
		}
		q.pop();
		return true;
	}

	void read_cmd(_state& state/*, bool put_prompt = true*/)
	{
		_Q q;
		
		srand(time(NULL));

		// initialize
		std::cout << "[Tip] Initialize before simulation. See README.md for help." << std::endl;
	start_initialize:
		do
		{
			_clear(q);
			q = _get();
			if (q.front().front() != 'i' || q.size() < 2)
			{
				std::cerr << "[Error] No initialization mode selected." << std::endl;
			}
			else
			{
				q.pop();
				break;
			}
		} while (true);
		switch (q.front().front())
		{
		case 'm':
		case 'r':
		case 's': initialize(state, q.front().front()); break;
		default: std::cerr << _INVALID << std::endl; goto start_initialize;
		}

		std::cout << "[Tip] Confirm settings before simulation. See README.md for help." << std::endl;

		// 防止用户硬要输入超过预期数量的 "end initialize"
		do
		{
			_clear(q);
			q = _get();
		} while (q.front().front() == 'e');


		// settings
		do
		{
			bool success = false, jump = false;
			if (!q.empty() && q.front().back() == 't')
				break;
			switch (q.front().front())
			{
			case 'b':
			case 'e': jump = true; // don't break
			case 'a': success = true; break;
			case 's': q.pop(); success = setsteps(q); break;
			case 'm': q.pop(); success = setmethod(q); break;
			case 'p': /*don't pop*/ success = setp(q); break;
			case 't': q.pop(); success = set_tlen(q); break;
			}
			if (!success)
				std::cerr << _INVALID << std::endl;
			if (jump)
				break;
			q = _get();
			if (q.front() == "start")
				break;
		} while (q.front().back() != 't' && q.front().front() != 'b' && q.front().front() != 'e');
		return;
	}

	~_user()
	{
		if (fout.is_open())
			fout.close();
	}
} user;


#endif