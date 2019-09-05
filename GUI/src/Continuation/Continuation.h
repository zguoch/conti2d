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

// constants definition
#define PI 3.141592653

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
    int WriteGmsh();
    int CheckOpts();
    int update();
    int run(); /**< run the main program according to the options */
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
    /**
     * @brief Calculate continuation coefficient, \f$ C_i \f$ , of a vertex centered polygon
     * 
     * \f$ \iint\limits_{{S_i}} {\frac{{U(\alpha ,\beta ,{z_0})}}{{{R^3}}}d\alpha d\beta } = U(\alpha ,\beta ,{z_0})\sum\limits_{j = 1}^n {({A_{j2}} - {A_{j1}} + {B_{j2}} - {B_{j1}})}  = {U_i}{C_i} \f$
     * @param alpha x-like coordinate of observation point \f$ Q_i \f$ on the plane for a upward continuation quation (Figure 2 of Guo and Zhang(2019))
     * @param beta x-like coordinate ...
     * @param polygon vertex (\f$ Q_i \f$) centered polygon, contains nodes, connection of the polygon
     * @return double coefficient of \f$ C_i \f$
     */
    double kernel(int index_Q,double deltaZ,double P[3], const cPolyMesh* polygon);
    /**
     * @brief calculate coefficience of a edge. The observation point is Q, and the start point of the edge vector is V1.
     * 
     * @param Q 
     * @param V1 
     * @param V2
     * @param deltaZ 
     * @return double 
     */
    double kernel_edge(double deltaZ,const double V1[2],const double V2[2],const double P[3]);
};

#endif