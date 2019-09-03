/**
 * @file stdfunc.h
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief Definition of some basic functions.
 * @version 1.0
 * @date 2019-09-03
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
#include <string.h>
using namespace std;
#include <math.h>
#include "MultiProgressBar.h"
#include <netcdf.h>
#ifdef _WIN32
#include "windows.h"
#endif


#define ERRCODE 2
#define ERR(e) {cout<<"Error:"<<nc_strerror(e)<<endl; exit(ERRCODE);}

#define PI 3.141592653

#define DWC_CGLS 1
#define DWC_TIKHONOV 2
#define DWC_INTEGRALITERATION 3
#define DWC_LANDWEBER 4

#define ERROR_COUT "["<<"\033[31mError"<<"\033[0m] "
#define WARN_COUT "["<<"\033[33mWarning"<<"\033[0m] "
#define ERROR_COUT "["<<"\033[31mError"<<"\033[0m] "
#define PURPLE "\033[35m"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define YELLOW "\033[33m"
#define BLUE "\033[34m"
#define COLOR_DEFALUT "\033[0m"
#include <unistd.h>
#define MOVEUP(x) printf("\033[%dA", (x))
// clean screen
#define CLEAR() printf("\033[2J") 
// move up cursor
#define MOVEUP(x) printf("\033[%dA", (x)) 
// move down cursor
#define MOVEDOWN(x) printf("\033[%dB", (x)) 
// move left cursor
#define MOVELEFT(y) printf("\033[%dD", (y)) 
// move right cursor
#define MOVERIGHT(y) printf("\033[%dC",(y)) 
// locate cursor
#define MOVETO(x,y) printf("\033[%d;%dH", (x), (y)) 
// reset cursor
#define RESET_CURSOR() printf("\033[H") 
// hide cursor
#define HIDE_CURSOR() printf("\033[?25l") 
// show cursor
#define SHOW_CURSOR() printf("\033[?25h") 
#define HIGHT_LIGHT() printf("\033[7m")
#define UN_HIGHT_LIGHT() printf("\033[27m")
////////////////////////////////////

/**
 * @brief Head of Golden Software Surfer format 6
 * 
 */
struct GrdHead
{
    int cols, rows;
    double bounds[6];
};

/**
 * @brief Read a ASCII Surfer format grid (6)
 * 
 * @param filename 
 * @param grdhead 
 * @param extNum 
 * @return double* 
 */
double* ReadGrd(string filename,GrdHead& grdhead,int extNum);

/**
 * @brief Save a grid data as ASCII surfer grid file.
 * 
 * @param filename 
 * @param grdhead 
 * @param data 
 * @param extNum 
 * @param savexxyz 
 * @param isInfo 
 * @return true 
 * @return false 
 */
bool SaveGrd(string filename, GrdHead grdhead,double* data,int extNum, bool savexxyz=false,bool isInfo=true);



/**
 * @brief Save grid data to netCDF format. see https://www.unidata.ucar.edu/software/netcdf/ for details.
 * 
 * @param filename 
 * @param grdhead 
 * @param data 
 * @param extNum 
 * @param isInfo 
 * @return true 
 * @return false 
 */
bool SaveGrd2netCDF(string filename, GrdHead grdhead,double* data,int extNum,bool isInfo=true);


// Calculate H2 norm of a vector: ||x||2
double Norm2(double* x,const int num);

/**
 * @brief Get file name from a full path, e.g. figures/figure_uwc_p2p.ps  -> figure_uwc_p2p
 * 
 * @param filepath 
 * @return string 
 */
string Path_GetBaseName(string filepath);

/**
 * @brief Get file extension name from a full file path. e.g. figures/figure_uwc_p2p.ps -> ps
 * 
 * @param filepath 
 * @return string 
 */
string Path_GetExtName(string filepath);

/**
 * @brief Get file path from a full path. e.g. figures/figure_uwc_p2p.ps -> figures/figure_uwc_p2p
 * 
 * @param filepath 
 * @return string 
 */
string Path_GetPath(string filepath);

/**
 * @brief Get file name from a full path. e.g. figures/figure_uwc_p2p.ps -> figures
 * 
 * @param filepath 
 * @return string 
 */
string Path_GetFileName(string filepath);

/**
 * @brief Get 2nd norm of the gradient of a grid data
 * 
 * @param result 
 * @param grdhead 
 * @return double 
 */
double Norm2_Gradient(double* result,GrdHead grdhead);

/**
 * @brief Save results of at each iteration step as vtk file
 * 
 * @param outputfile 
 * @param grdhead 
 * @param data 
 * @return int 
 */
int SaveGrd2VTK(string outputfile,GrdHead grdhead,double* data,double z=0);

/**
 * @brief Save a grid field data on a topography as vtk format (3d)
 * 
 * @param outputfile 
 * @param grdhead 
 * @param data N=mxn array
 * @param topo N=mxn array,  the dimension must be the save  as data
 * @return int 
 */
int SaveGrd2VTK_topo(string outputfile,GrdHead grdhead,double* data,double* topo);

bool SaveGrd2xyz(string filename, GrdHead grdhead, double* data,int extNum,bool isInfo=true);

static void StartText()
{
    //30: black  31:red  32:green  33:yellow  34:blue  35:purple  36:darkgreen
    cout<<"\033[33m";       //print text in yellow color
    cout << "***************************************************\n";
    cout << "*                 program conti2d                 *\n";
    cout << "*                 ~~~~~~~ ~~~~~~~                 *\n";
    cout << "*                                                 *\n";
    cout << "* Analytical continuation of potential field data *\n";
    cout << "*  - Upward Continuation from plane to plane      *\n";
    cout << "*  - Upward Continuation from plane to surface    *\n";
    cout << "*  - Downward Continuation from plane to plane    *\n";
    cout << "*  - Downward Continuation from surface to plane  *\n";
    cout << "*                                                 *\n";
    cout << "* See head file for input/output details.         *\n";
    cout << "* (c) Zhikui Guo, CUG, Aug 2016, Hangzhou         *\n";
    cout << "*                                                 *\n";
    cout << "***************************************************\n";
    cout << "\n\n";
    cout<<"\033[0m";
                                                                                                                                                                                                                      
}
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
static void helpINFO()
{
    string version="2.0";
    string author="Zhikui Guo";
    string locus="SIO, Hangzhou";
    unsigned int wordWidth=20;
    // time_t now=time(0);
    // char* now_str=ctime(&now);
    string now_str="Aug 8, 2016";

    //30: black  31:red  32:green  33:yellow  34:blue  35:purple  36:darkgreen
    cout<<"===================== conti2d ======================"<<endl;
    cout<<"Analytical continuation of potential field data"<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Author \033[032m"<<author<<"\033[0m"<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Locus \033[032m"<<locus<<"\033[0m"<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Date \033[032m"<<now_str<<"\033[0m"<<endl;
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Version \033[032m"<<version<<"\033[0m"<<endl;
    Text_Axis();
    cout<<setw(wordWidth)<<setiosflags(ios::left)<<"Synopsis:\033[031m"
        <<"conti2d "
        <<GREEN<<"input.grd "
        <<PURPLE<<"-G"<<GREEN<<"output.grd "
        <<PURPLE<<"-E"<<GREEN<<"extBoundNum "
        <<PURPLE<<"-H"<<GREEN<<"level_outputdata|.grd "
        <<PURPLE<<"-T"<<GREEN<<"level_inputdata|.grd "
        <<PURPLE<<"-D"<<COLOR_DEFALUT<<"["
        <<GREEN<<"+T"<<COLOR_DEFALUT<<"Tik_Par"<<PURPLE<<"|"
        <<GREEN<<"+I"<<COLOR_DEFALUT<<"It_Par"<<PURPLE<<"|"
        <<GREEN<<"+L"<<COLOR_DEFALUT<<"Landweber_Par"<<PURPLE<<"|"
        <<GREEN<<"+C"<<COLOR_DEFALUT<<"CGLS_Par"<<PURPLE<<""
        <<COLOR_DEFALUT<<"] "
        <<PURPLE<<"-t"<<GREEN<<"threads "
        <<PURPLE<<"-f"<<GREEN<<"[FrequencyDomain] "
        <<COLOR_DEFALUT<<endl;
    cout<<"===================================================="<<endl<<endl;
    // exit(1);
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
//error info
static void OutputErrorinfo(string info)
{
    cout<<"\033[31m";       //print error infor in red color
    cout<<info<<endl;
    cout<<"\033[0m";
    exit(0);
}
//warning info
static void OutputWarninginfo(string info)
{
    cout<<"\033[34m";       //print warning info as blue color
    cout<<info<<endl;
    cout<<"\033[0m";
}

#endif