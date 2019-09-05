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
#include "unistd.h"     //for linux and macos
#include "textcolor.h" //defined text color
#include <sys/ioctl.h> 
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

//static function
static void Text_Axis()
{
    cout<<"   /y"<<endl;
    cout<<"  /"<<endl;
    cout<<" /"<<endl;
    cout<<"/0____________x"<<endl;
    cout<<"|"<<endl;
    cout<<"|"<<endl;
    cout<<"|"<<endl;
    cout<<"|z(-)"<<endl;
}
static void StartText()
{
    // cout << "***************************************************\n";
    // cout << "*                 program conti2d                 *\n";
    // cout << "*                 ~~~~~~~ ~~~~~~~                 *\n";
    // cout << "*                                                 *\n";
    // cout << "* Analytical continuation of potential field data *\n";
    // cout << "*  - Upward Continuation from plane to plane      *\n";
    // cout << "*  - Upward Continuation from plane to surface    *\n";
    // cout << "*  - Downward Continuation from plane to plane    *\n";
    // cout << "*  - Downward Continuation from surface to plane  *\n";
    // cout << "*                                                 *\n";
    // cout << "* See head file for input/output details.         *\n";
    // cout << "* (c) Zhikui Guo, CUG, Aug 2016, Hangzhou         *\n";
    // cout << "*                                                 *\n";
    // cout << "***************************************************\n";       
    string version="2.0";
    string author="Zhikui Guo (zhikuiguo@live.cn)";
    string locus="SIO, Hangzhou";
    unsigned int wordWidth=10;
    string now_str="Aug 8, 2016";
    string namepro="conti2d";
    cout<<"===================== conti2d ======================"<<endl;
    cout<<"Analytical continuation of potential field data"<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Author "<<GREEN<<author<<COLOR_DEFALUT<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Locus "<<GREEN<<locus<<COLOR_DEFALUT<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Date "<<GREEN<<now_str<<COLOR_DEFALUT<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Version "<<GREEN<<version<<COLOR_DEFALUT<<endl;
    Text_Axis();                                                                                                                                                                                     
}

static void helpINFO()
{
    string version="2.0";
    string author="Zhikui Guo (zhikuiguo@live.cn)";
    string locus="SIO, Hangzhou";
    unsigned int wordWidth=10;
    string now_str="Aug 8, 2016";
    string namepro="conti2d";
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Synopsis:"<<endl;
    cout<<"conti2d [-I<input file>] [-O<output file> [-C<value|file>]]"<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Options:"<<endl;
    // option *
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<RED<<"  -I "
    <<COLOR_DEFALUT<<": Input file name. The file contains at least 4 columns, "
    <<BLUE<<"x y z field "<<COLOR_DEFALUT<<", separated by space or tab"<<endl;
    // option *
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<RED<<"  -O "
    <<COLOR_DEFALUT<<": Output file name. File format can be identified by the file extentional name, "
    <<BLUE<<".txt, .msh"<<COLOR_DEFALUT<<" ."<<endl;
    // option *
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<RED<<"  -C "
    <<COLOR_DEFALUT<<": Continuation option. The supported option are: "<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<BLUE<<"    -C"<<BLUE<<"value"<<COLOR_DEFALUT<<" : [float] elevation of continuation plane."<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<BLUE<<"    -C"<<BLUE<<"file"<<COLOR_DEFALUT<<" : continuation points file(ASCII) includes three columns, "<<BLUE<<" x y z "<<COLOR_DEFALUT<<"separated by space or tab"<<endl;
    // option *
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<RED<<"  -D "
    <<COLOR_DEFALUT<<": Discritization of the scatter observation points. Write (1) triangle mesh (.msh file) connect each points; (2) vertex centered polygon geometry (.geo file)"
    <<COLOR_DEFALUT<<"."<<endl;
    
    cout<<"==========================================================================================================="<<endl<<endl;
    exit(0);
}
static void StartText_artASCII()
{
    cout<<GREEN<<" ________          ________          ________           _________        ___           _______          ________     \n"
    <<"|\\   ____\\        |\\   __  \\        |\\   ___  \\        |\\___   ___\\     |\\  \\         /  ___  \\        |\\   ___ \\    \n"
    <<"\\ \\  \\___|        \\ \\  \\|\\  \\       \\ \\  \\\\ \\  \\       \\|___ \\  \\_|     \\ \\  \\       /__/|_/  /|       \\ \\  \\_|\\ \\   \n"
    <<" \\ \\  \\            \\ \\  \\\\\\  \\       \\ \\  \\\\ \\  \\           \\ \\  \\       \\ \\  \\      |__|//  / /        \\ \\  \\ \\\\ \\  \n"
    <<"  \\ \\  \\____        \\ \\  \\\\\\  \\       \\ \\  \\\\ \\  \\           \\ \\  \\       \\ \\  \\         /  /_/__        \\ \\  \\_\\\\ \\ \n"
    <<"   \\ \\_______\\       \\ \\_______\\       \\ \\__\\\\ \\__\\           \\ \\__\\       \\ \\__\\       |\\________\\       \\ \\_______\\\n"
    <<"    \\|_______|        \\|_______|        \\|__| \\|__|            \\|__|        \\|__|        \\|_______|        \\|_______|\n"
    <<COLOR_DEFALUT<<endl;    
}
static void HelloConti()
{
    struct winsize w;
    ioctl(0, TIOCGWINSZ, &w);
    if(w.ws_col>119)
    {
        StartText_artASCII();
    }else
    {
        StartText();
    }
}
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