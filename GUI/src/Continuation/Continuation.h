/**
 * @file cContinuation.h
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief The C++ class of cContinuation is defined in this head file
 * @version 1.0
 * @date 2019-08-26
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#ifndef CONTINUATION_H
#define CONTINUATION_H
#include "stdfunc.h"
#include "GModel.h"
 #include "omp.h"
// =========triangle===========
#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */
#include "triangle.h"
#include "trimesh.h"
//File format definition
#define FMT_UNKNOWN 0
#define FMT_GRID_SURF 11
#define FMT_NETCDF 12
#define FMT_GMSH 13
#define FMT_VTK 14
#define FMT_XYZ 15

#define UWC_P2P 0
#define UWC_P2S 1
#define DWC_P2P 2
#define DWC_S2P 3

#define DWC_CGLS 1
#define DWC_TIKHONOV 2
#define DWC_INTEGRALITERATION 3
#define DWC_LANDWEBER 4

// constants definition
#define PI 3.141592653

struct GrdHead
{
    int cols, rows;
    double bounds[6];
};
// // define actions for conti main program
// #define ACTION_CONTINUATION 1  //continuation of potential field data, if -C option is set
// #define ACTION_DISCRITIZATION 2 //discritization of observation scatter points, if -D is set
class cContinuation
{
private:
    bool m_isRegular; //whether the input data is on regular grid
    string m_extName_Input; 
    int m_typeFileFormat;
    int m_InfoLevel;
public:
    cContinuation();
    cContinuation(string );
    ~cContinuation();
    int ReadFile();
    int WriteFile();
    string getExtName(string ); /**< e.g., input xx/xxx/file.txt, return txt */
    string getFileName(string );/**< e.g., input xx/xxx/file.txt, return xx/xxx/file */
    string getPath(string );/**< e.g., input xx/xxx/file.txt, return xx/xxx */
    string getBaseName(string );/**< e.g., input xx/xxx/file.txt, return xx/xxx/file */
    int ReadXYZ(string FileName, cTriMesh& trimesh);
    double* ReadGrd(string filename, GrdHead& grdhead,int extNum=0);
    bool SaveGrd(string filename, GrdHead grdhead,double* data,int extNum, bool savexxyz=false,bool isInfo=true);
    int SaveGrd2VTK(string outputfile,GrdHead grdhead,double* data,double z=0);
    bool SaveGrd2xyz(string filename, GrdHead grdhead, double* data,int extNum,bool isInfo=true);
    int SaveGrd2VTK_topo(string outputfile,GrdHead grdhead,double* data,double* topo);
    string Grid2Gmsh_plane(string FileName,string Name_data="Original Field",double z0=0);//surfer grid file to gmsh
    string Grid2Gmsh_surface(string FileName,string Name_data="Original Field",string Topography="");//surfer grid file to gmsh 
    int WriteGmsh();
    int CheckOpts();
    int update();
    void UWC_p2p(string inputfilename,string outputfilename,double height1,double height2,int extNum,int num_thread,bool show=true);
    void UWC_p2s(string inputfilename,string outputfilename,double height1,string topoFile, int extNum, int num_thread);
    void DWC_p2p(string inputfilename,string outputfilename,double height1,double height2,int extNum,double DWC_parameter,int DWC_method, int num_thread);
    void DWC_s2p_LandweberIter(double** G,double* x, double* b,GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number);
    void DWC_s2p(string inputfilename,string outputfilename,string topo1,double height2,int extNum,double DWC_parameter,int DWC_method, int num_thread);

    //kernel
    //Equation 4
    void GetPmnij(double* Pmnij,int rows,int cols,double dx,double dy,double rph,double xm,double ym);
    int Getkernel_p2p_new(GrdHead grdhead, double rph, double* kernel_firstRow, int num_thread);
    int Getkernel_p2p_new(GrdHead grdhead, double rph, double** kernel, int num_thread);
    int Getkernel_p2s_new(GrdHead grdhead, double h1,double* topo2, double** kernel, int num_thread);
    // calculate upward continuation: b=Gx
    void UWC_Gij(double* b, double* G,double* x, GrdHead grdhead, int num_thread=1);
    //calculate upward continuation: b=Gx
    void UWC_Gij(double* b,double** G,double* x, int modelnum,int num_thread=1);
    // Get i row j column element of kernel matrix
    double GetGij(const int i, const int j, double* firstRow, const GrdHead grdhead);
    double Norm2(double* x,const int num);
    double Norm2_Gradient(double* result,GrdHead grdhead);
public: /*public data*/
    cTriMesh m_triMESH;
    string m_InputFile;
    string m_OutputFile;
    bool m_Action_Continuation;
    bool m_Action_Discritization;
private:
    int triangle(cTriMesh& trimesh);
    int getTypeofFileFormat(string extName);
    // vector<double>m_x,m_y,m_z,m_field;
    int Info(string info, int level=1);
    inline void ERR(string err){cout<<RED<<"Error: "<<COLOR_DEFALUT<<err<<endl;exit(0);}
    void init_triangulateio(struct triangulateio& tri );
    void free_triangulateio(struct triangulateio& tri );
    int init();
    
};

#endif