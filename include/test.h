#pragma once

#ifndef _TEST_H_
#define _TEST_H_

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <chrono>
#include <corecrt.h>
#include "3dv.h"
#include "physics.h"
#include "integrator.h"

constexpr const char* _LINE_LONG = "=======================";
constexpr const char* _LINE_MID = "=================";
constexpr const char* _LINE = "===========";

long double gettime()
{
	auto duration = std::chrono::high_resolution_clock::now().time_since_epoch();
	auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
	return static_cast<long double>(ns);
}

void test_3dv()
{
	_3dv v1, v2;
	std::cin >> v1.x >> v1.y >> v1.z;
	std::cin >> v2.x >> v2.y >> v2.z;
	std::cout << "v1 + v2 = " << (v1 + v2) << std::endl;
	std::cout << "v1 - v2 = " << (v1 - v2) << std::endl;
	std::cout << "v1 * 5.0 = " << (v1 * 5.0) << std::endl;
	std::cout << "5.0 * v1 = " << (5.0 * v1) << std::endl;
	std::cout << "v1 / 5.0 = " << (v1 / 5.0) << std::endl;
	std::cout << "v1 · v2 = " << v1.dot(v2) << std::endl;
	std::cout << "v1 × v2 = " << v1.cross(v2) << std::endl;
	std::cout << "|v1| = " << v1.mag() << " = sqrt(" << v1.mag_2() << ')' << std::endl;
	std::cout << "|v1 - v2| = " << v1.distance(v2) << std::endl;
	std::cout << "^v1 = " << v1._e() << std::endl;
	std::cout << "<v1, v2> = " << v1.angle_with(v2) << " rad" << std::endl;
	std::cout << "proj_v2(v1) = " << v1.proj_onto(v2) << std::endl;
	std::cout << "refl_v2(v1) = " << v1.reflect(v2) << std::endl;
	return;
}

void test_intrg()
{
	std::cout << "Test: Two-body problem (Earth orbiting the Sun)]" << std::endl;
	_state solar_sys, sys[3];
	long double t = 0.0, tm[3]{};
	double v_earth = sqrt(G * SOLAR_MASS / AU); // 计算地球公转速度大小
	solar_sys.add(_obj(STAR, SOLAR_MASS, _3dv(0), _3dv(0), _3dv(0), "sun"));
	solar_sys.add(_obj(PLANET, EARTH_MASS, _3dv(AU, 0, 0), _3dv(0, v_earth, 0), _3dv(0), "earth"));
	solar_sys.set_a();
	sys[0] = sys[1] = sys[2] = solar_sys;
	t = gettime();
	state_goto(sys[0], DAY_SECONDS * 30.0, 3600.0, "euler");
	tm[0] = gettime() - t;
	std::cout << _LINE << std::endl << "Euler:" << std::endl << sys[0]
		<< "Time (ns): " << tm[0] << std::endl << _LINE_MID << std::endl << std::endl;
	tm[0] = gettime();
	state_goto(sys[1], DAY_SECONDS * 30.0, 3600.0, "verlet");
	tm[1] = gettime() - tm[0];
	std::cout << _LINE << std::endl << "Verlet:" << std::endl << sys[1]
		<< "Time (ns): " << tm[1] << std::endl << _LINE_MID << std::endl << std::endl;
	tm[1] = gettime();
	state_goto(sys[2], DAY_SECONDS * 30.0, 3600.0, "rk4");
	tm[2] = gettime() - tm[1];
	std::cout << _LINE << std::endl << "RK4:" << std::endl << sys[2]
		<< "Time (ns): " << tm[2] << std::endl << _LINE_MID << std::endl << std::endl;
	return;
}

inline void printlnn()
{
	std::cout << std::endl << _LINE_LONG << _LINE_LONG << std::endl << std::endl;
	system("pause");
	return;
}

void test_all()
{
	std::cout << "[Test 3D vector: input two vectors]" << std::endl;
	test_3dv();
	printlnn();
	std::cout << "[Test integrator]" << std::endl;
	test_intrg();
	printlnn();
	return;
}

#endif