#include "user.h"
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
#include <format>
#include "3dv.h"
#include "physics.h"
#include "calendar.h"
#include <omp.h>

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

//inline _3dv _user::rand_v(double mark) // 弃用
//{
//	if (!rand_init)
//	{
//		gen = std::mt19937(rd());
//		rand_init = true;
//	}
//	return _3dv(
//		static_cast<double>(rand()) * mark / RAND_MAX,
//		static_cast<double>(rand()) * mark / RAND_MAX,
//		static_cast<double>(rand()) * mark / RAND_MAX
//	);
//}

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
			try
			{
				q = _get();
			}
			catch (...)
			{
				return;
			}
			if (q.empty())
				return;
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
		try
		{
			if (s)
				q = _get(fin);
			else
				q = _get();
		}
		catch (...)
		{
			return;
		}
		if (q.empty())
			return;
		if (q.empty() || q.front().front() == 'e')
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
			sample.pos = rand_v_new(phy::AU * 5); // 姑且先这么考虑
			sample.v = rand_v_new(phy::V_SUN * 3);
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
	unlimited = false;
	timelen = phy::YEAR * 100;
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

bool _user::set_display(_Q& q)
{
	if (q.empty())
	{
		display = true;
		return true;
	}
	display = false;
	if (q.size() < 2)
		return false;
	int w, h;
	try
	{
		w = std::stoi(q.front());
		q.pop();
		h = std::stoi(q.front());
		q.pop();
	}
	catch (...)
	{
		return false;
	}
	display = true;
	width = w;
	height = h;
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
		if (mode == 's' || mode == 'm' || mode == 'r')
			break;
		try
		{
			_clear(q);
			q = _get();
		}
		catch (...)
		{
			break;
		}
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
	if (mode != 's' && mode !='m' && mode !='r')
		mode = q.front().front();
	switch (mode)
	{
	case 'm':
	case 'r':
	case 's': initialize(state, mode); break;
	default: std::cerr << _INVALID << std::endl; goto start_initialize; // wtf
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
		try
		{
			_clear(q);
			q = _get();
		}
		catch (...)
		{
			break;
		}
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
		case 'd': q.pop(); success = set_display(q); break;
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
	if (display)
		dt = 1.; // display下积分步长可细一些
	else
		dt = (std::max)(1., timelen * 2e-6); // 至多积分五十万步以控制输出时长
	sample_dt = (std::max)(dt, (std::max)(static_cast<int>(sample_dt / dt), 1) * dt); // 整倍数化，较为必要 // 这个设置项已被废弃，没有设置的意义了
	steps = (std::max)(sample_dt, (std::max)(static_cast<int>(steps / sample_dt), 1) * sample_dt); // 这个设置项也已被废弃
	timelen += dt * 5; // 确保最后一年的数据输出
	finished = true;
	return;
}


long double _user::gettime()
{
	auto duration = std::chrono::high_resolution_clock::now().time_since_epoch();
	auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
	return static_cast<long double>(ns);
}

inline double _user::get_z(double z)
{
	return 2. / (1. + exp(-0.5 * z));
}

_3dv _user::getpos_ori(const _3dv& pos)
{
	_3dv res;
	res.x = static_cast<double>(width / 2) + pos.x / zoom;
	res.y = static_cast<double>(height / 2) - pos.y / zoom;
	res.z = get_z(pos.z);
	return res;
}

void _user::show(_state& s, _calendar& cal)
{
	static int zoom_times = 0, times = 0;
	static bool is_shrink = false;
	//static long double time = 0;
	if (!display)
		return;
	long double current_time = gettime();
	/*if (current_time - time >= 3e7)
	{
		time = gettime();
		return;
	}*/
	if (times < 5)
	{
		++times;
		return;
	}
	times = 0;
	if (!zoomed)
	{
		zoom = static_cast<double>((std::min)(height, width)) / 180.;
		std_rcircle = 3.; // 放弃根据屏幕大小调整
		zoomed = true;
	}
	double xmax = -INFINITY, ymax = -INFINITY;
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < s.objs.size(); ++i)
	{
		if (xmax < std::abs(s.objs[i].pos.x))
			xmax = std::abs(s.objs[i].pos.x);
		if (ymax < std::abs(s.objs[i].pos.y))
			ymax = std::abs(s.objs[i].pos.y);
	}
	if ((xmax * 2 / zoom) > width || (ymax * 2 / zoom) > height)
	{
		++zoom_times;
		zoom *= 1 +  1 / std::sqrt(static_cast<double>(zoom_times));
		is_shrink = true;
	}
	while ((!is_shrink) && ((xmax * 5 / zoom) < width && (ymax * 5 / zoom) < width))
		zoom /= 2;
	if (!is_shrink)
		is_shrink = true;
	cleardevice();
	setcolor(WHITE);
	std::string id = cal.get_current_sun().id;
	std::string title = std::format("Current sun: {} (rate = {:f}, {:f}, {:f}) \n Time: {:f} \n YEAR_LEN = {:f} \n Year {:d}: {}", id, cal.ranklist_ranks[id].v, cal.ranklist_ranks[id].hill, cal.ranklist_ranks[id].ecc, s.time, cal.std_year, static_cast<int>(s.time / cal.std_year), cal.era_name());
	xyprintf(0, 0, title.c_str());
	/*std::string test = std::format("a = {:f}; ecc = {:f}; energy = {:f}", s.get_semi_a(id), s.get_eccent(id).mag(), s.get_orbit_energy(id));
	xyprintf(0, height / 2, test.c_str());*/
	for (int i = 0; i < s.objs.size(); ++i)
	{
		_3dv screen = getpos_ori(s.objs[i].pos);
		int radius = static_cast<int>(std::ceil(screen.z * std_rcircle)); // !!!
		COLORS col = (s.objs[i].type == STAR) ? STAR_COLOR : PLANET_COLOR;
		setfillcolor(col);
		setcolor(col);
		fillellipse(static_cast<int>(std::round(screen.x)), static_cast<int>(std::round(screen.y)), radius, radius);
	}
	delay_ms(16);
	//time = gettime();
	return;
}