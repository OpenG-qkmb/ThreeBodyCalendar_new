#pragma once

#ifndef _TEST_H_
#define _TEST_H_

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <string>
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

void test_intrg_two_body()
{
	std::cout << "[Test]: Two-body problem (Earth orbiting the Sun)]" << std::endl;
	constexpr int N = 3;
	_state solar_sys, sys;
	long double t0, t1, e0 = 0, e1 = 0;
	const double v_earth = sqrt(G * SOLAR_MASS / AU); // 计算地球公转速度大小
	const std::string method[N] = {"Euler", "Verlet", "RK4"};
	solar_sys.add(_obj(STAR, SOLAR_MASS, _3dv(0), _3dv(0), _3dv(0), "sun"));
	solar_sys.add(_obj(PLANET, EARTH_MASS, _3dv(AU, 0, 0), _3dv(0, v_earth, 0), _3dv(0), "earth"));
	solar_sys.set_a();
	e0 = solar_sys.get_energy();
	for (int i = 0; i < N; ++i)
	{
		sys = solar_sys;
		t0 = gettime();
		state_goto(sys, DAY_SECONDS * 30.0, 3600.0, method[i]);
		t1 = gettime() - t0;
		e1 = sys.get_energy();
		std::cout << _LINE << std::endl << method[i] << ":" << std::endl << sys
			<< "Time (ns): " << t1 << std::endl
			<< "Energy0: " << e0 << "   " << "Energy1: " << e0 << std::endl 
			<< _LINE_MID << std::endl << std::endl;
	}
	return;
}

void test_intrg_merge() // FAILED
{
	std::cout << std::endl << "[Test]: Merge (Earth dropping into the Sun)]" << std::endl;
	_state test_sys;
	double v_earth = sqrt(G * SOLAR_MASS / AU);
	test_sys.add(_obj(STAR, SOLAR_MASS, _3dv(0), _3dv(0), _3dv(0), "sun"));
	test_sys.add(_obj(PLANET, EARTH_MASS, _3dv(AU, 0, 0), _3dv(0, 0, 0), _3dv(0), "earth"));
	test_sys.set_a();
	std::cout << _LINE << std::endl << "Before:" << std::endl << test_sys << std::endl;
	state_goto(test_sys, DAY_SECONDS * STD_YEAR , 60, "verlet");
	test_sys.merge();
	std::cout << _LINE << std::endl << "After:" << std::endl << test_sys << std::endl;
	return;
}

void test_intrg()
{
	test_intrg_two_body();
	// std::cout << _LINE_MID << std::endl;
	// test_intrg_merge();
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
	std::cout << "[Test 3D vector] Input two vectors:" << std::endl;
	test_3dv();
	printlnn();
	std::cout << "[Test integrator]" << std::endl;
	test_intrg();
	printlnn();
	return;
}

#endif