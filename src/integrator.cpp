#include "integrator.h"
#include "3dv.h"
#include "physics.h"
#include <vector>
#include <queue>
#include <map>
#include <string>
#include <string_view>
#include <cstdlib>
#include <omp.h>

void euler_step(_state& state, double dt) // 欧拉法前进dt
{
	state.set_a();
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < state.objs.size(); ++i)
	{
#pragma omp critical
		{
			state.objs[i].v += state.objs[i].a * dt;
			state.objs[i].pos += state.objs[i].v * dt;
		}
	}
	state.time += dt;
	state.set_a();
	return;
}

void verlet_step(_state& state, double dt) // 速度Verlet法前进dt
{

	state.set_a();
	std::queue<_3dv> a_old;
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < state.objs.size(); ++i)
	{
#pragma omp critical
		{
			a_old.push(state.objs[i].a);
			state.objs[i].pos += state.objs[i].v * dt + state.objs[i].a * (0.5 * dt * dt);
		}
	}
	state.set_a();
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < state.objs.size(); ++i)
	{
#pragma omp critical
		{
			state.objs[i].v += (a_old.front() + state.objs[i].a) * (0.5 * dt);
			a_old.pop();
		}
	}
	state.time += dt;
	state.set_a();
	return;
}
// https://zhuanlan.zhihu.com/p/1941592778577519357（文末）

void rk4_step(_state& state, double dt) // RK4法前进dt（若需，考虑优化？）
{
	size_t n = state.size();
	state.set_a();
	_state state_copy = state;
	std::vector<_3dv> k1_v(n), k1_a(n);
	std::vector<_3dv> k2_v(n), k2_a(n);
	std::vector<_3dv> k3_v(n), k3_a(n);
	std::vector<_3dv> k4_v(n), k4_a(n);
	// k1 = state (original)
	//state_copy.set_a();（先set_a再赋值）
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < n; ++i)
	{
		k1_v[i] = state.objs[i].v;
		k1_a[i] = state.objs[i].a;
	}
	// k2
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < n; ++i)
	{
#pragma omp critical
		{
			state_copy.objs[i].pos += k1_v[i] * (0.5 * dt);
			state_copy.objs[i].v += k1_a[i] * (0.5 * dt);
		}
	}
	state_copy.set_a();
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < n; ++i)
	{
		k2_v[i] = state_copy.objs[i].v;
		k2_a[i] = state_copy.objs[i].a;
	}
	// k3
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < n; ++i)
	{
#pragma omp critical
		{
			state_copy.objs[i].pos += k2_v[i] * (0.5 * dt);
			state_copy.objs[i].v += k2_a[i] * (0.5 * dt);
		}
	}
	state_copy.set_a();
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < n; ++i)
	{
		k3_v[i] = state_copy.objs[i].v;
		k3_a[i] = state_copy.objs[i].a;
	}
	// k4
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < n; ++i)
	{
#pragma omp critical
		{
			state_copy.objs[i].pos += k3_v[i] * dt;
			state_copy.objs[i].v += k3_a[i] * dt;
		}
	}
	state_copy.set_a();
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < n; ++i)
	{
		k4_v[i] = state_copy.objs[i].v;
		k4_a[i] = state_copy.objs[i].a;
	}
	// 结果
#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < n; ++i)
	{
#pragma omp critical
		{
			state.objs[i].pos += (k1_v[i] + 2.0 * k2_v[i] + 2.0 * k3_v[i] + k4_v[i]) * (dt / 6.0);
			state.objs[i].v += (k1_a[i] + 2.0 * k2_a[i] + 2.0 * k3_a[i] + k4_a[i]) * (dt / 6.0);
		}
	}
	state.time += dt;
	state.set_a();
	return;
}
// 翻百科

void integrate_dt(_state& state, double dt, std::string_view method) // 前进dt
{

	if (method.empty())
	{
		DEFAULT_INTEGRATOR(state, dt);
		state.merge();
		return;
	}
	char first_char = method[0];
	if (first_char >= 'A' && first_char <= 'Z')
		first_char += 'a' - 'A';
	auto it = FUNC_MAP.find(first_char);
	if (it != FUNC_MAP.end())
		it->second(state, dt);
	else
		DEFAULT_INTEGRATOR(state, dt);
	return;
}

void state_goto(_state& state, double t_tar, double dt, std::string_view method) // 积分前进至t_tar
{
	if (t_tar <= state.time)
		return;
	double t_rem = t_tar - state.time;
	while (abs(t_rem) >= SMALL_NUM)
	{
		double step_dt = (t_rem < dt) ? t_rem : dt;
		integrate_dt(state, step_dt, method);
		state.merge();
		t_rem -= step_dt;
	}
	state.set_a();
	state.merge();
	return;
}