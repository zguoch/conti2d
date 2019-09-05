/**
 * @file cContinuation.cpp
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief The cContinuation class contains all the algorithm of potential field data cContinuation, including upward and downword, from a arbitrary surface to another
 * @version 1.0
 * @date 2019-08-26
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "continuation.h"

cContinuation::cContinuation(string FileName):
m_typeFileFormat(FMT_UNKNOWN)
{
    init();
    m_InputFile=FileName;
    update();
}
cContinuation::cContinuation():
m_typeFileFormat(FMT_UNKNOWN)
{
    init();
}
int cContinuation::init()
{
    m_InfoLevel=2;
    m_Action_Continuation=false;
    m_Action_Discritization=false;
    m_InputFile="";
    m_OutputFile="";

    return 0;
}
int cContinuation::update()
{
    m_extName_Input=getExtName(m_InputFile);
    m_typeFileFormat=getTypeofFileFormat(m_extName_Input);
    if(m_typeFileFormat==FMT_UNKNOWN)
    {
        cout<<ERROR_COUT<<"Unknown file format: "<<m_InputFile<<endl;
        exit(0);
    }
    // 1. read file
    Info("Reading file ...");
    ReadFile(); //m_TriMesh will be updated after reading input file
    // 2. calculate vertex centered polygon mesh, make every observation point at the center of each polygon
    // Info("Generating vertex centered polygon mesh ..");
    // m_triMESH.GenerateVertexCenteredPolygon();
    // Info("Generate vertex centered polygon mesh done!");
    return 0;
}
cContinuation::~cContinuation()
{
}

int cContinuation::ReadFile()
{
    switch (m_typeFileFormat)
    {
    case FMT_XYZ:
        {
            ReadXYZ(m_InputFile,m_triMESH);
            //if input file is in xyz format, the trianglar meshing is required
            m_triMESH.numberofelements=triangle(m_triMESH);
            Info("Discritization of observation points is done.");
        }
        
        break;
    default:
        break;
    }
    return 0;
}

string cContinuation::getFileName(string filepath)
{
    int pos_ext=filepath.find_last_of('.');
    string filename=filepath.substr(0,pos_ext);
    return filename;
}
string cContinuation::getPath(string filepath)
{
    int pos_ext=filepath.find_last_of('/');
    string path=filepath.substr(0,pos_ext);
    return path;
}
string cContinuation::getBaseName(string filepath)
{
    int pos_ext=filepath.find_last_of('.');
    string pathname=filepath.substr(0,pos_ext);
    int pos_basename=pathname.find_last_of('/');
    string basename=pathname.substr(pos_basename+1,pathname.length());
    return basename;
}
string cContinuation::getExtName(string filepath)
{
    int pos_ext=filepath.find_last_of('.');
    string extname=filepath.substr(pos_ext+1,filepath.length());
    return extname;
}

/**
 * @brief Delaunay triangulation of Cartesian data. 
 * It find how the points should be connected to give the most equilateral triangulation possible 
 * 
 * @return vector<int> 
 */
int cContinuation::triangle(cTriMesh& trimesh)
{
    uint64_t n=trimesh.numberofpoints;
    struct triangulateio In, Out, vorOut;
    init_triangulateio(In);
    init_triangulateio(Out);
    init_triangulateio(vorOut);
    In.numberofpoints = (int)n;
    In.pointlist = (REAL *) malloc(In.numberofpoints * 2 * sizeof(REAL));
    
	/* Copy x,y points to In structure array */
    int i=0,j=0;
    
	for (i = j = 0; i < n; i++) {
		In.pointlist[j++] = trimesh.m_nodes[i].x;
		In.pointlist[j++] = trimesh.m_nodes[i].y;
	}
    
    string cmd="zIQBDn";//command n means calculate neighboring triangles
    triangulate ((char*)cmd.data(), &In, &Out, &vorOut);
    int nEdges=Out.numberoftriangles*3;
    trimesh.m_triElements=new TRIANGLE[Out.numberoftriangles];
    for(int i=0;i<Out.numberoftriangles;i++)
    {
        trimesh.m_triElements[i].isBoundary=false;
        for(int j=0;j<3;j++)
        {
            trimesh.m_triElements[i].iNodes[j]=Out.trianglelist[i*3+j];
            trimesh.m_triElements[i].iNeighbour[j]=Out.neighborlist[i*3+j];
            if(Out.neighborlist[i*3+j]==-1)
            {
                trimesh.m_triElements[i].isBoundary=true;
                //determine which vertex is the boundary vertex
                if(j==2)
                {
                    trimesh.m_nodes[Out.trianglelist[i*3]].isBoundary=true;// magic!!! -1 flat means the next vertex as boundary point
                }else
                {
                    trimesh.m_nodes[Out.trianglelist[i*3+j+1]].isBoundary=true;
                }
            }
        }
    }
    //get number of boundary points
    trimesh.numberofpoints_bd=0;
    for(int i=0;i<trimesh.numberofpoints;i++)
    {
        if(trimesh.m_nodes[i].isBoundary)trimesh.numberofpoints_bd++;
    }
    //get number of boundary elements
    trimesh.numberofelements_bd=0;
    for(int i=0;i<trimesh.numberofelements;i++)
    {
        if(trimesh.m_triElements[i].isBoundary)trimesh.numberofelements_bd++;
    }
    // -----------------------
    free_triangulateio(In);
    free_triangulateio(Out);
    free_triangulateio(vorOut);
    return Out.numberoftriangles;
}

int cContinuation::getTypeofFileFormat(string extName)
{
    if(extName=="vtk")return FMT_VTK;
    if(extName=="xyz")return FMT_XYZ;
    if(extName=="grd")return FMT_GRID_SURF;
    if(extName=="msh")return FMT_GMSH;
    if(extName=="txt")return FMT_XYZ;
    if(extName=="dat")return FMT_XYZ;
    return FMT_UNKNOWN;
}
/**
 * @brief 
 * 
 * @param FileName 
 * @param x 
 * @param y 
 * @param z 
 * @param field 
 * @return int 
 */
int cContinuation::ReadXYZ(string FileName, cTriMesh& trimesh)
{
    ifstream fin(FileName);
    if(!fin)
    {
        cout<<ERROR_COUT<<"Open file failed: "<<FileName<<endl;
        exit(0);
    }
    vector<double>x,y,z,field;
    double tmp_x,tmp_y,tmp_z,tmp_f;
    while (!fin.eof())
    {
    fin>>tmp_x>>tmp_y>>tmp_z>>tmp_f;
    x.push_back(tmp_x);
    y.push_back(tmp_y);
    z.push_back(tmp_z);
    field.push_back(tmp_f);
    }
    // if the last row is empty, the x,y,z would be duplicate, need to remove
    for(int i=x.size()-1;i>0;i--)
    {
        if((x[i]==x[i-1]) && (y[i]==y[i-1]) && (z[i]==z[i-1]))
        {
            x.pop_back();
            y.pop_back();
            z.pop_back();
            Info("Remove one duplicate point at last line");
        }else
        {
            break;
        }
    }
    fin.close();
    // copy data to trimesh
    int numberofpoints=x.size();
    trimesh.numberofpoints=numberofpoints;
    trimesh.m_nodes=new NODES[numberofpoints];
    trimesh.field=new double[numberofpoints];
    for(int i=0;i<numberofpoints;i++)
    {
        trimesh.m_nodes[i].x=x[i];
        trimesh.m_nodes[i].y=y[i];
        trimesh.m_nodes[i].z=z[i];
        trimesh.m_nodes[i].isBoundary=false;
        trimesh.field[i]=field[i];
    }
    

    Info("xyz file is loaded.");
    return 0;
}

int cContinuation::Info(string info, int level)
{
    if(level==0)return 0;
    if(level<=m_InfoLevel)
    {
        cout<<PURPLE<<"Info: "<<COLOR_DEFALUT<<info<<endl;
    }
    return 0;
}
/**
 * @brief Initialize elements (variable and dynamic arrays) of a triangulateio variable
 * 
 * @param tri 
 */
void cContinuation::init_triangulateio(struct triangulateio& tri)
{
    tri.pointlist=(REAL *) NULL;                                               /* In / out */
    tri.pointattributelist=(REAL *) NULL;                                      /* In / out */
    tri.pointmarkerlist=0;                                          /* In / out */
    tri.numberofpoints=0;                                            /* In / out */
    tri.numberofpointattributes=0;                                   /* In / out */

    tri.trianglelist=(int*)NULL;                                             /* In / out */
    tri.triangleattributelist=(REAL *) NULL;                                   /* In / out */
    tri.trianglearealist=(REAL *) NULL;                                         /* In only */
    tri.neighborlist=(int*)NULL;                                             /* Out only */
    tri.numberoftriangles=0;                                         /* In / out */
    tri.numberofcorners=0;                                           /* In / out */
    tri.numberoftriangleattributes=0;                                /* In / out */

    tri.segmentlist=(int*)NULL;                                              /* In / out */
    tri.segmentmarkerlist=(int*)NULL;                                      /* In / out */
    tri.numberofsegments=0;                                         /* In / out */

    tri.holelist=(REAL *)NULL;                        /* In / pointer to array copied out */
    tri.numberofholes=0;                                      /* In / copied out */

    tri.regionlist=(REAL *)NULL;                     /* In / pointer to array copied out */
    tri.numberofregions=0;                                   /* In / copied out */

    tri.edgelist=(int*)NULL;                                                 /* Out only */
    tri.edgemarkerlist=(int*)NULL;            /* Not used with Voronoi diagram; out only */
    tri.normlist=(REAL *)NULL;               /* Used only with Voronoi diagram; out only */
    tri.numberofedges=0;   
}
/**
 * @brief Free memory of dynamic arrays of a triangulateio variable
 * 
 * @param tri 
 */
void cContinuation::free_triangulateio(struct triangulateio& tri)
{
    free(tri.pointlist);
    free(tri.pointattributelist);
    free(tri.trianglelist);
    free(tri.triangleattributelist);
    free(tri.trianglearealist);
    free(tri.neighborlist);
    free(tri.segmentlist);
    free(tri.segmentmarkerlist);
    free(tri.holelist);
    free(tri.regionlist);
    free(tri.edgelist);
    free(tri.edgemarkerlist);
    free(tri.normlist);
}

/**
 * @brief Check the input arguments and print some error, warning or information
 * 
 * @return int 
 */
int cContinuation::CheckOpts()
{
    if(m_InputFile=="")ERR("You have to give the input file by -I option.");
    if(m_OutputFile=="")ERR("You have to give the output file by -O option.");
    return 0;
}

/**
 * @brief Write mesh and field data to gmsh file
 * 
 * @return int 
 */
int cContinuation::WriteGmsh()
{
    ofstream fout(m_OutputFile);
    fout<<"$MeshFormat"<<endl;
    fout<<"2.2 0 8"<<endl;
    fout<<"$EndMeshFormat"<<endl;
    fout<<"$Nodes"<<endl;
    fout<<m_triMESH.numberofpoints<<endl;
    for(int i=0;i<m_triMESH.numberofpoints;i++)
    {
      fout<<i+1<<" "<<m_triMESH.m_nodes[i].x<<" "
                    <<m_triMESH.m_nodes[i].y<<" "
                    <<m_triMESH.m_nodes[i].z<<endl;
    }
    fout<<"$EndNodes"<<endl;
    fout<<"$Elements"<<endl;
    fout<<m_triMESH.numberofelements<<endl;
    for(int i=0;i<m_triMESH.numberofelements;i++)
    {
      fout<<i+1<<" 2 2 1 1 ";
      for(int j=0;j<3;j++)fout<<m_triMESH.m_triElements[i].iNodes[j]+1<<" ";
      fout<<endl;
    }
    fout<<"$EndElements"<<endl;
    // field data
    fout<<"$NodeData"<<endl;
    fout<<1<<endl;
    fout<<"\""<<m_triMESH.fieldname<<"\""<<endl;
    fout<<1<<endl;
    fout<<0<<endl;
    fout<<3<<endl;
    fout<<0<<endl;
    fout<<1<<endl;
    fout<<m_triMESH.numberofpoints<<endl;
    for(int i=0;i<m_triMESH.numberofpoints;i++)
    {
        fout<<i+1<<" "<<m_triMESH.field[i]<<endl;
    }
    fout<<"$EndElementData"<<endl;
    fout.close();
    Info("Write mesh and field data to .msh file successfully: "+m_OutputFile);
    return 0;
}

int cContinuation::WriteFile()
{
    string extName=getExtName(m_OutputFile);
    int type_fmt=getTypeofFileFormat(extName);
    switch (type_fmt)
    {
    case FMT_GMSH:
        WriteGmsh();    
        break;
    default:
        break;
    }
    return 0;
}
