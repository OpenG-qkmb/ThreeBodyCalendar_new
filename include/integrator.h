#pragma once

#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "3dv.h"
#include "physics.h"
#include <vector>
#include <queue>
#include <map>
#include <string>
#include <string_view>
#include <cstdlib>

using _func_ptr = void(*)(_state&, double);

// 积分方法（单步）

void euler_step(_state& state, double dt); // 欧拉法前进dt
void verlet_step(_state& state, double dt); // 显式Verlet法前进dt，无需欧拉法初始化
void rk4_step(_state& state, double dt); // RK4法前进dt（若需，考虑优化？）


const std::map<char, _func_ptr> FUNC_MAP = {
	{'e', euler_step},
	{'v', verlet_step},
	{'r', rk4_step}
};

constexpr _func_ptr DEFAULT_INTEGRATOR = verlet_step; // 默认积分方法

void integrate_dt(_state& state, double dt, std::string_view method = "verlet"); // 前进dt
void state_goto(_state& state, double t_tar, double dt, std::string_view method = "verlet"); // 积分前进至t_tar

#endif