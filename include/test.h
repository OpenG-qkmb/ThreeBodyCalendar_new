#pragma once

#ifndef _TEST_H_
#define _TEST_H_

#include <iostream>
#include "3dv.h"

constexpr char _LINE[] = "=======================";

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
	std::cout << "v1 ¡¤ v2 = " << v1.dot(v2) << std::endl;
	std::cout << "v1 ¡Á v2 = " << v1.cross(v2) << std::endl;
	std::cout << "|v1| = " << v1.mag() << " = sqrt(" << v1.mag_2() << ')' << std::endl;
	std::cout << "|v1 - v2| = " << v1.distance(v2) << std::endl;
	std::cout << "^v1 = " << v1._e() << std::endl;
	std::cout << "<v1, v2> = " << v1.angle_with(v2) << " rad" << std::endl;
	std::cout << "proj_v2(v1) = " << v1.proj_onto(v2) << std::endl;
	std::cout << "refl_v2(v1) = " << v1.reflect(v2) << std::endl;
	return;
}

inline void println()
{
	std::cout << std::endl << _LINE << std::endl << std::endl;
	return;
}

void test_all()
{
	std::cout << "Test 3D vector: input two vectors" << std::endl;
	test_3dv();
	println();
	return;
}

#endif