//
//  Conti2D.h
//  Conti2D
//
//  Created by Zhikui Guo on 2017/2/1.

#ifndef Conti3D_h
#define Conti3D_h

#ifdef _WIN32
//define something for Windows (32-bit and 64-bit, this part is common)
//cout << "32位win系统" << endl;
#endif
#ifdef _WIN64
//cout << "64位win系统" << endl;
#include "windows.h"
#else
//define something for Windows (32-bit only)
#endif

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string.h>
using namespace std;
#include <math.h>
#include "omp.h"
#include "FFTN.h"
#include "MultiProgressBar.h"
#include "netcdf.h"
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
//#include <unistd.h>

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
 * @brief head information of a grid data, e.g. Goldernsoftward Surfer Grid.
 * 
 */

struct GrdHead
{
    int cols, rows;
    double bounds[6];
};
/**
 * @brief Read a Surfer Grid data from file
 * 
 * @param filename 
 * @param grdhead 
 * @param extNum 
 * @return double* return data array
 */
double* ReadGrd(string filename,GrdHead& grdhead,int extNum);
/**
 * @brief Save a Surfer Grid data to file
 * 
 * @param filename 
 * @param grdhead 
 * @param data 
 * @param extNum 
 * @param savexxyz true: save xyz format simultaneously
 * @param isInfo 
 * @return true 
 * @return false 
 */
bool SaveGrd(string filename, GrdHead grdhead,double* data,int extNum, bool savexxyz=false,bool isInfo=true);


bool SaveGrd2netCDF(string filename, GrdHead grdhead,double* data,int extNum,bool isInfo=true);


/**
 * @brief Get the Kernal Matrix object
 * 
 * @param grdhead 
 * @param G 
 * @param rph 
 * @param num_thread 
 */ 

void GetKernalMatrix(GrdHead grdhead, double* G, const double rph, int num_thread=4);
/**
 * @brief Get the kernel matrix using the new developed formula
 * 
 * @param grdhead 
 * @param G a 1-d array to store the kernel matrix
 * @param rph height of continuation
 * @param num_thread 
 */
void GetKernalMatrix_new(GrdHead grdhead, double* G, const double rph, int num_thread=4);
/**
 * @brief Get the value at i row and the j column of kernel matrix from the first row of the matrix only.
 * 
 * @param i 
 * @param j 
 * @param firstRow The first row of the kernel matrix
 * @param grdhead 
 * @return double 
 */
double GetGij(const int i, const int j, double* firstRow, const GrdHead grdhead);

/**
 * @brief Upward continue from plane to plane
 * 
 * @param grdhead 
 * @param rph 
 * @param kernel 
 * @param num_thread 
 * @return int 
 */
int Getkernel_p2p(GrdHead grdhead, double rph, double** kernel, int num_thread);
int Getkernel_p2p_new(GrdHead grdhead, double rph, double* kernel_firstRow, int num_thread);
int Getkernel_p2p_new(GrdHead grdhead, double rph, double** kernel, int num_thread);
void UWC_p2p(string inputfilename,string outputfilename,
    double height1,double height2,int extNum, int num_thread, bool isProgress,bool useOldKernel,
    string filename_exact);
void UWC_p2p_f(string inputfilename,string outputfilename,double height1,double height2,
                int extNum,string filename_exact);
int UWC_p2p_f(double* inputdata,double* result,GrdHead grdhead,double rph);

//2. UWC: uneven surface to plane
/**
 * @brief Upward continue from surface to plane
 * 
 * @param grdhead 
 * @param terrain1 2-D array of terrain of the surface
 * @param height2 Height of upward continuation
 * @param kernel 
 * @param num_thread 
 * @return int 
 */
int Getkernel_u2p(GrdHead grdhead, double* terrain1, double height2, double** kernel, int num_thread);
int Getkernel_u2p_new(GrdHead grdhead, double* terrain1, double height2,
                    double* nx,double* ny,double* nz, 
                    double** kernel, int num_thread);
int Getkernel_p2s_new(GrdHead grdhead, double h1,double* topo2, double** kernel, int num_thread);
// void UWC_u2p(string inputfilename,string outputfilename,
    // string terrainfile1, double height2, int num_thered, bool isProgress);
void GetPmnij(double* Pmnij,int rows,int cols,double dx,double dy,double rph,double xm,double ym);
void UWC(double* datain, double* dataout, GrdHead grdhead,double** G);
void UWC_Gij(double* b, double* G,double* x, GrdHead grdhead, int num_thread=1);
void UWC_Gji(double* b, double* G,double* x, GrdHead grdhead, int num_thread=1);
void UWC_Gij(double* b,double** G,double* x, int modelnum,int num_thread=1);
//compute: b=(GT*G+lamubda*I)*x == b=M*x
void UWC_G_CGLS_Tik(double* b,double** G,double* x, int modelnum,double lambda2,int num_thread=1);
void UWC_Gji(double* b,double** G,double* x, int modelnum,int num_thread=1);
// 3. DWC: Plane to plane
/**
 * @brief Downward continue from plane to plane
 * 
 * @param inputfilename Input file name of grid data 
 * @param outputfilename Output file name of the downward continuation result
 * @param height1 
 * @param filename_topo2 
 */
void DWC_u2p(string inputfilename,string outputfilename,double height1,string filename_topo2);
void DWC_p2p_f(string inputfilename,string outputfilename,double height1,
    double height2,int extNum,double TRP, int num_thread, bool isProgress,bool useOldKernel);
void DWC_p2p(string inputfilename,string outputfilename,double height1,
    double height2,int extNum,double DWC_parameter,int DWC_method, int num_thread, bool isProgress,bool useOldKernel,
    string filename_exact);
/**
 * @brief Upward continuation from plane to surface
 * 
 * @param inputfilename 
 * @param outputfilename 
 * @param height1 
 * @param topoFile 
 * @param extNum 
 * @param num_thread 
 * @param isProgress 
 * @param useOldKernel 
 * @param filename_exact 
 */
void UWC_p2s(string inputfilename,string outputfilename,
    double height1,string topoFile,int extNum, int num_thread, bool isProgress,
    bool useOldKernel,string filename_exact);
/**
 * @brief downward continuation from surface to plane
 * 
 * @param inputfilename 
 * @param outputfilename 
 * @param topo1 
 * @param height2 
 * @param extNum 
 * @param DWC_parameter 
 * @param DWC_method 
 * @param num_thread 
 * @param isProgress 
 * @param useOldKernel 
 * @param filename_exact 
 */
void DWC_s2p(string inputfilename,string outputfilename,string topo1,
    double height2,int extNum,double DWC_parameter,int DWC_method, int num_thread, bool isProgress,bool useOldKernel,
    string filename_exact);

void DWC_Tikhonov_old(double* G_firstRow,double* dataout,double* indata,double TRP,int kmax,double daierta,GrdHead grdhead,int num_thread);
void DWC_p2p_CGLS(double* G_firstRow,double* x, double* b,GrdHead grdhead,int extNum,double delta,int num_thread);
void DWC_s2p_CGLS(double** G,double* x, double* b,GrdHead grdhead,int extNum,double delta,int num_thread);
void DWC_s2p_Tikhonov(double** G,double* x, double* b,GrdHead grdhead,int extNum,double lambda,int num_thread);

void DWC_s2p_ItegrationIter(double** G,double* x, double* b,
    GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number,double* ExactSolution=NULL);
void DWC_p2p_ItegrationIter(double* G,double* x, double* b,
    GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number,double* ExactSolution=NULL);

void DWC_p2p_LandweberIter(double* G,double* x, double* b,
    GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number,double* ExactSolution=NULL);
void DWC_s2p_LandweberIter(double** G,double* x, double* b,
    GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number,double* ExactSolution=NULL);

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
//保存中间结果为时间序列的vtk格式
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
