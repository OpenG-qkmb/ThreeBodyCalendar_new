#include "user.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <cstdlib>
#include <ctime>
#include <random>
#include <functional>
#include "3dv.h"
#include "physics.h"

inline int _user::mymin(const int& a, const int& b) const
{
	return (a < b) ? a : b;
}

void _user::_clear(_Q& q)
{
	while (!q.empty())
		q.pop();
	return;
}

void _user::set_lower(std::string& st)
{
	for (int i = 0; i < st.size(); ++i)
	{
		if (st[i] >= 'A' && st[i] <= 'Z')
			st[i] += 'a' - 'A';
	}
	return;
}

_Q _user::_get(std::istream& s_in)
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

inline _3dv _user::rand_v(double mark)
{
	if (!rand_init)
	{
		srand(static_cast<unsigned int>(time(NULL)));
		rand_init = true;
	}
	return _3dv(
		static_cast<double>(rand()) * mark / RAND_MAX,
		static_cast<double>(rand()) * mark / RAND_MAX,
		static_cast<double>(rand()) * mark / RAND_MAX
	);
}

_3dv _user::rand_v_new(double mark)
{
	if (!rand_init)
	{
		gen = std::mt19937(rd());
		rand_init = true;
	}
	std::uniform_real_distribution<double> dist(-mark, mark);
	return _3dv(
		dist(gen),
		dist(gen),
		dist(gen)
	);
}


void _user::initialize(_state& state, const char& mode)
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
		if (!fin.is_open())
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
		case 'e': sample.type = PLANET; is_planet_added = true; break;
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
		if (sample.type == STAR)
			sample.m *= phy::SOLAR_MASS;
		/*else
		*	sample.m *= phy::EARTH_MASS;*/

		if (isrand)
		{
			sample.pos = rand_v_new(phy::AU * 50); // 姑且先这么考虑
			sample.v = rand_v_new(phy::V_SUN * 10);
			sample.a = _3dv();
			state.add(sample);
			continue;
		}

		// 自此起分叉 orbit单独处理 否则设置三维向量

		if (q.front().front() == 'o')
		{
			q.pop();
			std::string sunid = q.front();
			double r = phy::AU;
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
			double xyz[6] = { 0 }, tmp = .0; // 这里的6是无可辩驳的6，无需单独设作常量
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

bool _user::setsteps(_Q& q)
{
	if (q.size() < 2)
		return false;
	try
	{
		sample_dt = std::stod(q.front());
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

bool _user::setmethod(_Q& q)
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

bool _user::setp(_Q& q)
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

bool _user::set_tlen(_Q& q)
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

bool _user::set_analyse_id(_state& state, _Q& q)
{
	if (q.empty())
		return false;
	while (!q.empty())
	{
		analyse_id += q.front() + " ";
		q.pop();
	}
	analyse_id.pop_back();
	for (int i = 0; i < state.objs.size(); ++i) // NO ITERATOR
	{
		if (state.objs[i].id == analyse_id && state.objs[i].type == PLANET)
		{
			state.analyse_obj = std::ref(state.objs[i]);
			/*std::cerr << "[Info] Added!" << std::endl << state.analyse_obj.get() << std::endl;*/
			return true;
		}
	}
	analyse_id.clear();
	return false;
}

void _user::read_cmd(_state& state/*, bool put_prompt = true*/)
{
	_Q q;
	bool console = true;
	char mode = '\0';
	// srand(static_cast<unsigned int>(time(NULL)));

	// initialize
	std::cout << "[Tip] Initialize before simulation. See README.md for help." << std::endl;
start_initialize:
	do
	{
		if (mode != '\0')
			break;
		_clear(q);
		q = _get();
		if (q.empty())
		{
			std::cerr << "[Error] Console unavailable." << std::endl;
			console = false;
			break;
		}
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
	if (!console)
		return;
	if (mode == '\0')
		mode = q.front().front();
	switch (mode)
	{
	case 'm':
	case 'r':
	case 's': initialize(state, mode); break;
	default: std::cerr << _INVALID << std::endl; goto start_initialize;
	}
	if (!is_planet_added)
	{
		std::cerr << "[Error] No PLANET has yet been added to the system." << std::endl;
		goto start_initialize;
	}

	std::cout << "[Tip] Confirm settings before simulation. See README.md for help." << std::endl;

	// 防止用户硬要输入超过预期数量的 "end initialize"
	do
	{
		_clear(q);
		q = _get();
	} while (q.front().front() == 'e');

	bool analyse_set = false;

settings:

	// settings
	do
	{
		bool success = false, jump = false;
		if (!q.empty() && q.front().back() == 't')
			break;
		switch (q.front().front())
		{
		case 'b':
		case 'e': jump = true; success = true; break;
		case 'a': q.pop(); success = analyse_set = set_analyse_id(state, q); break; // already changed to "analyse_id"
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
	} while (q.front().back() != 't' && q.front().front() != 'b' && q.front().front() != 'e');
	if (!analyse_set)
	{      
		if (!is_planet_added)
		{
			std::cerr << "[Error] No PLANET has yet been added for analysis. Restart the program to address this issue." << std::endl;
			q = _get();
			goto settings;
		}
		for (int i = 0; i < state.objs.size(); ++i) // NO ITERATOR
		{
			if (state.objs[i].type == PLANET)
			{
				state.analyse_obj = std::ref(state.objs[i]);
				analyse_set = true;
				break;
			}
		}
	}
	dt = ::std::max(1., timelen * 1e-6); // 经验之谈：不这样，那么用时巨大。积分十万步用时约数秒
	sample_dt = ::std::max(dt, ::std::max(static_cast<int>(sample_dt / dt), 1) * dt); // 整倍数化，较为必要
	steps = ::std::max(sample_dt, ::std::max(static_cast<int>(steps / sample_dt), 1) * sample_dt);
	finished = true;
	return;
}