/**
 * @file stdfunc.h
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief Head file of some commonly used functions
 * @version 1.0
 * @date 2019-08-25
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#ifndef STDFUNC_H
#define STDFUNC_H
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;
#include <math.h>
// include head file for argument parse of main function
// ==================Perhaps only works in Mac or Linux system===============

#include "textcolor.h"  //defined text color
#ifdef _WIN32
// define something for Windows (32-bit and 64-bit, this part is common)
// cout << "32λwinϵͳ" << endl;
#endif
#ifdef _WIN64
// cout << "64λwinϵͳ" << endl;
#include "windows.h"
#else
	#include <sys/ioctl.h>
	
	#include "unistd.h"     //for linux and macos
#endif
// ==========================================================================

// commonly used math macro
#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MOD			/* Knuth-style modulo function (remainder after floored division) */
#define MOD(x, y) (x - y * floor((double)(x)/(double)(y)))
#endif


/**
 * @brief Determine whether a string is a number
 * 
 * @param str 
 * @return true 
 * @return false 
 */
bool isNum(string str);


/**
 * @brief Cross of two 2D-vectors, return \f$ \vec{pa} \times \vec{pb} \f$
 * 
 * @param p start point
 * @param a endpoint1
 * @param b endpoint2
 * @return double 
 */
double TwoCross(double p[2],double a[2],double b[2]);
/**
 * @brief Cross of two 2D-vectors, return \f$ \vec{v1} \times \vec{v2} \f$
 * 
 * @param v1 
 * @param v2 
 * @return double 
 */
double TwoCross(double v1[2],double v2[2]);
/**
 * @brief Dot of two 2D-vectors, return \f$ \vec{pa} \cdot \vec{pb} \f$
 * 
 * @param a 
 * @param b 
 * @return double 
 */
double TwoDot(double a[2],double b[2]);
#endif