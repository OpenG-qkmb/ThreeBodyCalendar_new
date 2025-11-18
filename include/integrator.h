#pragma once

#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "3dv.h"
#include "physics.h"
#include <vector>
#include <queue>
#include <map>
#include <string_view>
#include <cstdlib>

using _func_ptr = void(*)(_state&, double);

// 积分方法（单步）

void euler_step(_state& state, double dt) // 欧拉法前进dt
{
	state.set_a();
	for (_obj& o : state.objs)
	{
		o.v += o.a * dt;
		o.pos += o.v * dt;
	}
	state.time += dt;
	return;
}

void verlet_step(_state& state, double dt) // 显式Verlet法前进dt，无需欧拉法初始化
{
	state.set_a();
	std::queue<_3dv> a_old;
	for (_obj& o : state.objs)
	{
		a_old.push(o.a);
		o.pos += o.v * dt + o.a * (0.5 * dt * dt);
	}
	state.set_a();
	for (_obj& o : state.objs)
	{
		o.v += (a_old.front() + o.a) * (0.5 * dt);
		a_old.pop();
	}
	state.time += dt;
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
	for (size_t i = 0; i < n; ++i)
	{
		k1_v[i] = state.objs[i].v;
		k1_a[i] = state.objs[i].a;
	}
	// k2
	for (size_t i = 0; i < n; ++i)
	{
		state_copy.objs[i].pos += k1_v[i] * (0.5 * dt);
		state_copy.objs[i].v += k1_a[i] * (0.5 * dt);
	}
	state_copy.set_a();
	for (size_t i = 0; i < n; ++i)
	{
		k2_v[i] = state_copy.objs[i].v;
		k2_a[i] = state_copy.objs[i].a;
	}
	// k3
	for (size_t i = 0; i < n; ++i)
	{
		state_copy.objs[i].pos += k2_v[i] * (0.5 * dt);
		state_copy.objs[i].v += k2_a[i] * (0.5 * dt);
	}
	state_copy.set_a();
	for (size_t i = 0; i < n; ++i)
	{
		k3_v[i] = state_copy.objs[i].v;
		k3_a[i] = state_copy.objs[i].a;
	}
	// k4
	for (size_t i = 0; i < n; ++i)
	{
		state_copy.objs[i].pos += k3_v[i] * dt;
		state_copy.objs[i].v += k3_a[i] * dt;
	}
	state_copy.set_a();
	for (size_t i = 0; i < n; ++i)
	{
		k4_v[i] = state_copy.objs[i].v;
		k4_a[i] = state_copy.objs[i].a;
	}
	// 得出结果
	for (size_t i = 0; i < n; ++i)
	{
		state.objs[i].pos += (k1_v[i] + 2.0 * k2_v[i] + 2.0 * k3_v[i] + k4_v[i]) * (dt / 6.0);
		state.objs[i].v += (k1_a[i] + 2.0 * k2_a[i] + 2.0 * k3_a[i] + k4_a[i]) * (dt / 6.0);
	}
	state.time += dt;
	return;
}
// 翻百科能找到这个积分方法的公式


// 后续应手动确定默认积分方法
const std::map<char, _func_ptr> FUNC_MAP = {
	{'e', euler_step},
	{'v', verlet_step},
	{'r', rk4_step}
}; // 便于用户手动切换积分方法

constexpr _func_ptr DEFAULT_INTEGRATOR = verlet_step; // 默认积分方法

void integrate_dt(_state& state, double dt, std::string_view method = "verlet") // 前进dt，method指定积分方法
{
	
	if (method.empty())
	{
		DEFAULT_INTEGRATOR(state, dt);
		return;
	}
	auto it = FUNC_MAP.find(method[0]);
	if (it != FUNC_MAP.end())
		it->second(state, dt);
	else
		DEFAULT_INTEGRATOR(state, dt);
	return;
}

void state_goto(_state& state, double t_tar, double dt, std::string_view method = "verlet") // 积分前进至t_tar
{
	if (t_tar <= state.time)
		return;
	double t_rem = t_tar - state.time;
	while (abs(t_rem) >= SMALL_NUM)
	{
		double step_dt = (t_rem < dt) ? t_rem : dt;
		integrate_dt(state, step_dt, method);
		t_rem -= step_dt;
	}
	state.set_a();
	return;
}

#endif