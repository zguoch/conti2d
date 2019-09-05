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
#include "Continuation.h"

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
    case FMT_GRID_SURF:
        {

        }
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

string cContinuation::Grid2Gmsh(string FileName,string Name_data)
{
    string gmshfile=getFileName(FileName)+".gmsh";
    Msg::Info("Reading surfer grid file and display on GUI");
    ifstream fin;
    fin.open(FileName, ios::in);
    if (!fin)
    {
        Msg::Error("Open file failed: %s",FileName.c_str());
        return "";
    }
    string dsaa;
    int nx,ny;
    double bounds[6];
    //read head
    fin >> dsaa;
    if(dsaa=="DSAA")
    {
        Msg::Info("The grd file in ASCII");
    }else{
        Msg::Error("Input grd file is Binary. Please convert it to ascii: https://convert.goldensoftware.com/Application/Conversion");
        return "";
    }
    fin >> nx;
    fin >> ny;

    for (int i = 0; i < 6; i++)
    {
        fin >> bounds[i];
    }
    double dx=(bounds[1]-bounds[0])/(nx-1);
    double dy=(bounds[3]-bounds[2])/(ny-1);
    //read data
    double* data = new double[nx*ny];
    for (int i = 0; i < ny; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            fin >> data[j+nx*i];
        }
    }
    fin.close();
    //transform to gmsh
    ofstream fout(gmshfile);
    if (!fin)
    {
        Msg::Error("Open file failed: %s",FileName.c_str());
        return "";
    }else
    {
        Msg::Info("Save grd file to gmsh: %s",gmshfile.c_str());
    }
    fout<<"$MeshFormat"<<endl;
    fout<<"2.2 0 8"<<endl;
    fout<<"$EndMeshFormat"<<endl;
    fout<<"$Nodes"<<endl;
    fout<<nx*ny<<endl;
    for(int i=0;i<ny;i++)
    {
        for(int j=0;j<nx;j++)
        {
            fout<<j+i*nx+1<<" "
                <<bounds[0]+dx*j<<" "
                <<bounds[2]+dy*i<<" "
                <<0<<endl;
        }
    }
    fout<<"$EndNodes"<<endl;
    fout<<"$Elements"<<endl;
    int nx_el=nx-1,ny_el=ny-1;
    int nel=nx_el*ny_el;
    fout<<nel<<endl;
    for(int i=0;i<ny_el;i++)
    {
        for(int j=0;j<nx_el;j++)
        {
            int index_LL=j+i*nx+1;
            int index_TL=j+(i+1)*nx+1;
            fout<<j+i*nx_el+1<<" "
                <<3<<" 2 1 1 "
                <<index_LL<<" "
                <<index_LL+1<<" "
                <<index_TL+1<<" "
                <<index_TL<<endl;
        }
    }
    fout<<"$EndElements"<<endl;
    //write data
    fout<<"$NodeData"<<endl;
    fout<<1<<endl;
    fout<<"\""<<Name_data<<"\""<<endl;
    fout<<"1 \n0 \n3 \n0 \n1"<<endl;
    fout<<nx*ny<<endl;
    for(int i=0;i<ny*nx;i++)
    {
        // for(int j=0;j<nx;j++)
        {
            fout<<i+1<<" "<<data[i]<<endl;
        }
    }
    fout<<"$EndNodeData"<<endl;
    fout.close();
    delete[]  data;
    return gmshfile;
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

double* cContinuation::ReadGrd(string filename, GrdHead& grdhead,int extNum)
{
    ifstream fin;
    fin.open(filename, ios::in);
    if (!fin)
    {
        Msg::Error("Open file failed: %s",filename.c_str());
        return NULL;
    }
    string dsaa;
    //read head
    fin >> dsaa;
    if(dsaa=="DSAA")
    {
        Msg::Info("The grd file in ASCII");
    }else{
        Msg::Error("Input grd file is Binary. Please convert it to ascii: https://convert.goldensoftware.com/Application/Conversion");
        return NULL;
    }
    fin >> grdhead.cols;
    fin >> grdhead.rows;
    for (int i = 0; i < 6; i++)
    {
        fin >> grdhead.bounds[i];
    }
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    grdhead.bounds[0]-=extNum*dx;
    grdhead.bounds[1]+=extNum*dx;
    grdhead.bounds[2]-=extNum*dy;
    grdhead.bounds[3]+=extNum*dy;
    grdhead.cols+=2*extNum;
    grdhead.rows+=2*extNum;
    //read data
    double* data = new double[grdhead.rows*grdhead.cols];
    for (int i = extNum; i < grdhead.rows-extNum; i++)
    {
        for (int j = extNum; j < grdhead.cols-extNum; j++)
        {
            fin >> data[j+grdhead.cols*i];
        }
    }
    //extend as zero
    int index0=extNum;
    for(int i=0;i<extNum;i++)
    {
        for(int j=extNum;j<grdhead.rows-extNum;j++)
        {
            data[j+grdhead.cols*i]=data[index0+grdhead.cols*i];
        }
    }
    index0=grdhead.cols-extNum-1;
    for(int i=grdhead.cols-extNum;i<grdhead.cols;i++)
    {
        for(int j=extNum;j<grdhead.rows-extNum;j++)
        {
            data[j+grdhead.cols*i]=data[index0+grdhead.cols*i];
        }
    }
    //extend y direction (row direction)
    index0=extNum;
    for(int i=0;i<grdhead.cols;i++)
    {
        for(int j=0;j<extNum;j++)
        {
            data[j+grdhead.cols*i]=data[j+grdhead.cols*index0];
        }
    }
    index0=grdhead.rows-extNum-1;
    for(int i=0;i<grdhead.cols;i++)
    {
        for(int j=grdhead.rows-extNum;j<grdhead.rows;j++)
        {
            data[j+grdhead.cols*i]=data[j+grdhead.cols*index0];
        }
    }
    fin.close();
    return data;
}

string cContinuation::UWC_p2p(string inputfilename,string outputfilename,double height1,double height2,int extNum,int num_thread)
{
    cout << "***************************************************\n";
    cout << " Upward continuation from plane to plane:"<<height1<<"->"<<height2<<endl;
    cout << "***************************************************\n";
    //1. read grd data
    GrdHead grdhead;
    double* indata = NULL;
    if (!(indata = ReadGrd(inputfilename, grdhead,extNum)))return "";
    int modelnum = grdhead.rows*grdhead.cols;
    //2. height of continuation
    double rph=fabs(height1-height2);

    //3. kernel mat.
    double* G_firstRow=new double[modelnum];
    cout<<"calculating kernal matrix"<<endl;
    Getkernel_p2p_new(grdhead,rph,G_firstRow,num_thread);
    
    //4. outdata
    double *outdata=new double[modelnum];
    //compute UWC
    cout<<"calculating uwc: "<<num_thread<<" threads are used"<<endl;
    UWC_Gij(outdata,G_firstRow, indata,grdhead,num_thread);
    cout << "Finished\n";
    // 4. write result
    // string basename_infile=Path_GetBaseName(inputfilename);
    string ext_outputfile=getExtName(outputfilename);
    if(ext_outputfile=="grd" || ext_outputfile=="GRD")
    {
        // cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        if (!SaveGrd(outputfilename, grdhead, outdata,extNum))return "";
    }
    else if(ext_outputfile=="vtk" || ext_outputfile=="VTK")
    {
        // cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2VTK(outputfilename,grdhead,outdata,height2);
        // string filename_outfile=Path_GetFileName(outputfilename);
        // SaveGrd2VTK(filename_outfile+"_origin.vtk",grdhead,indata,height1);//save the origin data as well
    }else if(ext_outputfile=="xyz" || ext_outputfile=="txt" || ext_outputfile=="dat")
    {
        // cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2xyz(outputfilename, grdhead, outdata,extNum);
    }
    // else if(ext_outputfile=="nc")
    // {
    //     cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
    //     SaveGrd2netCDF(outputfilename, grdhead, outdata,extNum);
    // }
    else
    {
        cout<<RED<<"The output file format is not supported: "<<ext_outputfile<<COLOR_DEFALUT<<endl;
    }
    
    //delete data array
    delete[] indata;
    delete[] outdata;
    delete[] G_firstRow;

    return "";
}

/**
 * @brief Ger complete kernel matrix for upward continuation from plane to plane.
 * 
 * @param grdhead 
 * @param rph (\f$ \Delta z \f$ in equation 4) Height or vertical distance for the upward continuation. 
 * @param kernel 
 * @param num_thread 
 * @return int 
 */
int cContinuation::Getkernel_p2p_new(GrdHead grdhead, double rph, double** kernel, int num_thread)
{
    cout<<"new kernel of plane to surface: "<<num_thread<<" threads are used"<<endl;
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    int datanum=grdhead.rows*grdhead.cols;
    //omp parallel
     omp_set_num_threads(num_thread);
    // ProgressBar bar0(grdhead.rows);
    // MultiProgressBar multibar(grdhead.rows,COLOR_BAR_BLUE);
    #pragma omp parallel for shared(grdhead,dx,dy,kernel) 
    for(int irow=0; irow<grdhead.rows; irow++)
    {
        double yup=irow*dy;
        for(int icol=0; icol<grdhead.cols; icol++)
        {
            double* vector_row_K=new double[datanum];
            //=============================calculate Pmnij=======================
            int index_col_k=icol+irow*grdhead.cols;
            double xup=icol*dx;
            // double rph=topo2[index_col_k]-h1;
            GetPmnij(vector_row_K,grdhead.rows,grdhead.cols,dx,dy,rph,xup,yup);
            for(int k=0;k<datanum;k++)
            kernel[index_col_k][k]=vector_row_K[k];
            //====================================================================
            delete[] vector_row_K;
        }
        //#pragma omp critical
        // multibar.Update();
    }
	cout<<"\n";
    return 0;
}

/**
 * @brief Ger the first row of kernel matrix for upward continuation from plane to plane.
 * 
 * @param grdhead 
 * @param rph \f$ \Delta z \f$
 * @param kernel_firstRow 
 * @param num_thread 
 * @return int 
 */
int cContinuation::Getkernel_p2p_new(GrdHead grdhead, double rph, double* kernel_firstRow, int num_thread)
{
    cout<<"calculating first row of new kernel"<<endl;
    // double xup=0,yup=0;
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    int datanum=grdhead.rows*grdhead.cols;
    
    GetPmnij(kernel_firstRow,grdhead.rows,grdhead.cols,dx,dy,rph,0,0);
    
    return 1;
}

/**
 * @brief Get the kernel vector for a calculation point \f$ (x_m, y_m) \f$. 
 * 
 * See equation 4 in the manuscript.
 * 
 * \f$ \begin{gathered}
    {P_{m,n:i,j}} = \Delta z({x_m},{y_n})\int_{{\alpha _i} - \Delta x/2}^{{\alpha _i} + \Delta x/2} {\int_{{\beta _j} - \Delta y/2}^{{\beta _j} + \Delta y/2} {\frac{{d\alpha d\beta }}{{{R^3}}}} } {\text{  }} \\ 
    \;\;\;\; = \arctan \left( {\frac{{(\alpha  - {x_m})(\beta  - {y_n})}}{{\Delta z({x_m},{y_n})\sqrt {{{(\alpha  - {x_m})}^2} + {{(\beta  - {y_m})}^2} + \Delta z{{({x_m},{y_n})}^2}} }}} \right)|_{{\alpha _i} - \Delta x/2}^{{\alpha _i} + \Delta x/2}|_{{\beta _j} - \Delta y/2}^{{\beta _j} + \Delta y/2} \\ 
    \end{gathered} \f$
 * 
 * @param Pmnij 
 * @param rows Points number in y direction.
 * @param cols Points number in x direction.
 * @param dx \f$ \Delta x \f$ in equation (4)
 * @param dy \f$ \Delta y \f$ in equation (4)
 * @param rph \f$ \Delta z \f$ in equation (4)
 * @param xm \f$ x_m \f$ in equation (4)
 * @param ym \f$ x_m \f$ in equation (4)
 */
void cContinuation::GetPmnij(double* Pmnij,int rows,int cols,double dx,double dy,double rph,double xm,double ym)
{
    double xup=xm;
    double yup=ym;
    double dx_2=dx/2.0;
    double dy_2=dy/2.0;
    double rph_2=rph*rph;
    int index=0;
    double p11,p12,p21,p22;
    double x1,x2,y1,y2,x0,y0;
    double x1x1,x2x2,y1y1,y2y2;
    double ydown,xdown;
    for (int i=0;i<rows;i++)
    {
        ydown=i*dy;
        for(int j=0;j<cols;j++)
        {
            xdown=j*dx;
            //------
            x0=xup-xdown;
            y0=yup-ydown;
            x1=x0-dx_2;
            x2=x0+dx_2;
            y1=y0-dy_2;
            y2=y0+dy_2;
            
            x1x1=x1*x1;
            y1y1=y1*y1;
            x2x2=x2*x2;
            y2y2=y2*y2;
            //--------
            p22=atan2(x2*y2,rph*sqrt(rph_2+x2x2+y2y2));
            p11=atan2(x1*y1,rph*sqrt(rph_2+x1x1+y1y1));
            p12=atan2(x1*y2,rph*sqrt(rph_2+x1x1+y2y2));
            p21=atan2(x2*y1,rph*sqrt(rph_2+x2x2+y1y1));
            Pmnij[index]=(p22+p11-p21-p12)/(2*PI);
            index++;
        }
    }
}

/**
 * @brief Calculate upward continuation as \f$ b=\mathbf{G}x \f$, see equation (5). 
 * This function is used several times in Landweber iteration steps.
 * 
 * @param b 
 * @param G 
 * @param x 
 * @param modelnum 
 * @param num_thread 
 */
void cContinuation::UWC_Gij(double* b,double** G,double* x, int modelnum,int num_thread)
{
     omp_set_num_threads(num_thread);
    #pragma omp parallel for shared(b,x,modelnum)
    for (int i = 0; i < modelnum; i++)
    {
        b[i] = 0;
        for (int j = 0; j < modelnum; j++)
        {
            b[i] += G[i][j]*x[j];
        }
    }
}

/**
 * @brief Calculate upward continuation \f$ b=\mathbf{G}x \f$ just using the first row of the kernel matrix.
 * Only valid for case of plane to plane.
 * 
 * @param b 
 * @param G 
 * @param x 
 * @param grdhead 
 * @param num_thread 
 */
void cContinuation::UWC_Gij(double* b, double* G,double* x, GrdHead grdhead, int num_thread)
{
    int modelnum = grdhead.rows*grdhead.cols;

     omp_set_num_threads(num_thread);
    #pragma omp parallel for 				
    for (int i = 0; i < modelnum; i++)
    {
        b[i] = 0;
        for (int j = 0; j < modelnum; j++)
        {
            b[i] += GetGij(i,j,G,grdhead)*x[j];
        }
    }
}

/**
 * @brief Get the i row and j column element  \f$ K_{ij} \f$ from the first row of the kernel matrix.
 * 
 * Corresponding to equation (7) in the manuscript.
 * 
 * \f$ \begin{array}{*{20}{c}}
    {k\left( {m,n:i,j} \right) = k\left( {1,1:\left| {{i_0}} \right|,\left| {{j_0}} \right|} \right) \cdot sign({i_0}) \cdot sign({j_0})} \\ 
    {{i_0} = (i - m),\;\;\;\;{j_0} = (j - n)} 
    \end{array},\;sign(x) = \left\{ {\begin{array}{*{20}{c}}
    {1,\;if\;x > 0} \\ 
    {0,\;if\;x = 0} \\ 
    { - 1,\;if\;x < 0} 
    \end{array}} \right.\f$
 * 
 * @param i 
 * @param j 
 * @param firstRow 
 * @param grdhead 
 * @return double 
 */
double cContinuation::GetGij(const int i, const int j, double* firstRow, const GrdHead grdhead)
{
    double Gij;
    int row_block, colum_block, count_block;
    int modelnum = grdhead.rows*grdhead.cols;
    int index_row=0;
    int index_col=0;  
    row_block = i / grdhead.cols;						
    colum_block = j / grdhead.cols;					
    index_row=i-row_block*grdhead.cols;                 
    index_col=j-colum_block*grdhead.cols;               
    int index_block=abs(row_block-colum_block);         
    int col0=abs(index_row-index_col);                 
    Gij=firstRow[index_block*grdhead.cols+col0];
    return Gij;
}

bool cContinuation::SaveGrd(string filename, GrdHead grdhead, double* data,int extNum, bool savexxyz,bool isInfo)
{
    if(isInfo)cout<<filename<<GREEN<<" saved"<<COLOR_DEFALUT<<endl;
    ofstream fout(filename);
    if (!fout)
    {
        cout<<RED << "Open file false: "<<COLOR_DEFALUT << filename << "\n";
        return false;
    }
    //compute zmin and zmax
    double zmin=data[0];
    double zmax=zmin;
    int num_data=grdhead.rows*grdhead.cols;
    int index0=0;
    for (int i=extNum; i<grdhead.rows-extNum; i++) 
    {
        for(int j=extNum;j<grdhead.cols-extNum;j++)
        {
            if (zmin>data[index0]) {
            zmin=data[index0];
            }
            if (zmax<data[index0]) {
                zmax=data[index0];
            }
        }
    }
    grdhead.bounds[4]=zmin;
    grdhead.bounds[5]=zmax;
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    //write head
    fout << "DSAA\n";
    fout << grdhead.cols-2*extNum<<"\t";
    fout << grdhead.rows-2*extNum<<"\n";
    fout << grdhead.bounds[0]+extNum*dx<<"\t";
    fout << grdhead.bounds[1]-extNum*dx<<"\n";
    fout << grdhead.bounds[2]+extNum*dy<<"\t";
    fout << grdhead.bounds[3]-extNum*dy<<"\n";
    fout << grdhead.bounds[4]<<"\t";
    fout << grdhead.bounds[5]<<"\n";
    //write data
    for (int i = extNum; i < grdhead.rows-extNum; i++)
    {
        for (int j = extNum; j < grdhead.cols-extNum; j++)
        {
            fout << data[j+grdhead.cols*i]<<" ";
        }
        fout << "\n";
    }
    fout.close();
    
    return true;
}

 int cContinuation::SaveGrd2VTK(string outputfile,GrdHead grdhead,double* data,double z)
{
    //vtk format
    ofstream fout(outputfile);
    if(!fout)
    {
        // cout<<"["<<RED<<"error"<<COLOR_DEFALUT<<"]open file failed: "<<outputfile+".vtk"<<endl;
        exit(0);
    }
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    int num_pt=grdhead.cols*grdhead.rows;
    fout<<"# vtk DataFile Version 2.0"<<endl;
    fout<<"VTK from vtk2gm"<<endl;
    fout<<"ASCII"<<endl;
    fout<<"DATASET STRUCTURED_GRID"<<endl;
    fout<<"DIMENSIONS "<<grdhead.cols<<" "<<grdhead.rows<<" 1"<<endl;
    fout<<"POINTS "<<num_pt<<" float"<<endl;
    int index=0;
    for(int iy=0;iy<grdhead.rows;iy++)
    {
        for(int ix=0;ix<grdhead.cols;ix++)
        {
            fout<<grdhead.bounds[0]+ix*dx<<" "
                <<grdhead.bounds[2]+iy*dy<<" "
                <<z<<" ";
            index++;
        }
    }fout<<endl;
    fout<<"POINT_DATA "<<num_pt<<endl;
    fout<<"SCALARS PotentialField float"<<endl;
    fout<<"LOOKUP_TABLE default"<<endl;
    for(int i=0;i<num_pt;i++)
    {
        fout<<data[i]<<" ";
    }fout<<endl;
    fout.close();
    // cout<<"Output VTK file of result finished"<<endl;
    return 0;
}


bool cContinuation::SaveGrd2xyz(string filename, GrdHead grdhead, double* data,int extNum,bool isInfo)
{
    //save xyz format as well
    ofstream fxyz(filename);
    if (!fxyz)
    {
        // cout<<RED << "Open file false: "<<COLOR_DEFALUT << filename << "\n";
        return false;
    }else
    {
        if(isInfo)cout<<filename<<GREEN<<" saved"<<COLOR_DEFALUT<<endl;
    }
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    for(int i=extNum;i<grdhead.rows-extNum;i++)
    {
        for(int j=extNum;j<grdhead.cols-extNum;j++)
        {
            fxyz<<grdhead.bounds[0]+j*dx<<"\t"
                <<grdhead.bounds[2]+i*dy<<"\t"
                <<data[j+grdhead.cols*i]
                <<endl;
        }
    }
    fxyz.close();
    return true;
}
