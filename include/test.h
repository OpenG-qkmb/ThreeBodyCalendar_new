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

long double gettime();
void test_3dv();
void test_intrg_two_body_comp();
void test_intrg_two_body_timeline();
void test_intrg_merge();
void test_intrg();
inline void printlnn();
void test_all();

#endif