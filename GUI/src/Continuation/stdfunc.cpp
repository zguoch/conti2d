/**
 * @file stdfunc.cpp
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief Source file of some commonly used functions
 * @version 1.0
 * @date 2019-08-25
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "stdfunc.h"

bool isNum(string str)
{
    stringstream sin(str);
    double d;
    char c;
    if(!(sin >> d))
        return false;
    if (sin >> c)
        return false;
    return true;
}

/**
 * @brief The cross of two vectors on a 2D plane. \f$ \vec{pa} \times \vec{pb} \f$ .
 * The direction will be parallel with z axis. positive means along z direction, negative means -z direction.
 * 
 * @param p 
 * @param a 
 * @param b 
 * @return double 
 */
double TwoCross(double p[2],double a[2],double b[2])
{
    return (a[0]-p[0])*(b[1]-p[1])-(a[1]-p[1])*(b[0]-p[0]);
}
/**
 * @brief calculate cross of two vectors v1 x v2
 * 
 * @param v1 
 * @param v2 
 * @return double 
 */
double TwoCross(double v1[2],double v2[2])
{
    return (v1[0]*v2[1]-v1[1]*v2[0]);
}
double TwoDot(double a[2],double b[2])
{
    return (a[0]*b[0]+a[1]*b[1]);
}