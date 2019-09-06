/**
 * @file Conti2D.h 
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief Definition of sub-functions of continuation of potential field data.
 * @version 1.0
 * @date 2019-09-03
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#ifndef Conti3D_h
#define Conti3D_h

#include "stdfunc.h"
#include "omp.h"
#include "FFTN.h"


//Upward continuation from plane to plane
void UWC_p2p(string inputfilename,string outputfilename,double height1,double height2,int extNum, int num_thread, bool isProgress,string filename_exact);

//Upward continuation from plane to plane in frequency domain
void UWC_p2p_f(string inputfilename,string outputfilename,double height1,double height2,int extNum,string filename_exact);

//Downward continuation from plane to plane
void DWC_p2p(string inputfilename,string outputfilename,double height1,double height2,int extNum,double DWC_parameter,int DWC_method, int num_thread, bool isProgress,string filename_exact);

//Upward continuation from surface to plane
void UWC_p2s(string inputfilename,string outputfilename,double height1,string topoFile,int extNum, int num_thread, bool isProgress,string filename_exact);

//Downward continuation from surface to plane
void DWC_s2p(string inputfilename,string outputfilename,string topo1,double height2,int extNum,double DWC_parameter,int DWC_method, int num_thread, bool isProgress,string filename_exact);

// Get i row j column element of kernel matrix
double GetGij(const int i, const int j, double* firstRow, const GrdHead grdhead);
 
// Get first row of new kernal matrix 
int Getkernel_p2p_new(GrdHead grdhead, double rph, double* kernel_firstRow, int num_thread);

// Get first row of new kernal matrix 
int Getkernel_p2p_new(GrdHead grdhead, double rph, double** kernel, int num_thread);

// Calculate upward continuation in frequency domain
int UWC_p2p_f(double* inputdata,double* result,GrdHead grdhead,double rph);

//K of equation 5
int Getkernel_p2s_new(GrdHead grdhead, double h1,double* topo2, double** kernel, int num_thread);

//Equation 4
void GetPmnij(double* Pmnij,int rows,int cols,double dx,double dy,double rph,double xm,double ym);

// calculate upward continuation: b=Gx
void UWC_Gij(double* b, double* G,double* x, GrdHead grdhead, int num_thread=1);

//calculate upward continuation: b=Gx
void UWC_Gij(double* b,double** G,double* x, int modelnum,int num_thread=1);

//Downward continuation from plane to plane using Landweber iteration method (Hansen, 2002)
void DWC_p2p_LandweberIter(double* G,double* x, double* b,GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number,double* ExactSolution=NULL);

void DWC_s2p_LandweberIter(double** G,double* x, double* b,GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number,double* ExactSolution=NULL);

#endif 
