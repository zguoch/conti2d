/**
 * @file Conti2D.cpp
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief Implementation of the sub-functions of continuation of potential field in spatial domain.
 * @version 1.0
 * @date 2019-09-03
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#include "Conti2D.h"

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
double GetGij(const int i, const int j, double* firstRow, const GrdHead grdhead)
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
void GetPmnij(double* Pmnij,int rows,int cols,double dx,double dy,double rph,double xm,double ym)
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
 * @brief Get complete kernel matrix for upward continuation from plane to surface
 * 
 * @param grdhead 
 * @param h1 Elevation of the observation field.
 * @param topo2 Topography of the continuation points(or calculation points)
 * @param kernel Save kernel matrix to a 2D array
 * @param num_thread How many threads used to calculate
 * @return int 
 */
int Getkernel_p2s_new(GrdHead grdhead, double h1,double* topo2, double** kernel, int num_thread)
{
    cout<<"new kernel of plane to surface: "<<num_thread<<" threads are used"<<endl;
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    int datanum=grdhead.rows*grdhead.cols;
    //omp parallel
    omp_set_num_threads(num_thread);
    // ProgressBar bar0(grdhead.rows);
    MultiProgressBar multibar(grdhead.rows,COLOR_BAR_BLUE);
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
            double rph=topo2[index_col_k]-h1;
            GetPmnij(vector_row_K,grdhead.rows,grdhead.cols,dx,dy,rph,xup,yup);
            for(int k=0;k<datanum;k++)
            kernel[index_col_k][k]=vector_row_K[k];
            //====================================================================
            delete[] vector_row_K;
        }
        #pragma omp critical
        multibar.Update();
    }cout<<"\n";
    
    return 1;
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
int Getkernel_p2p_new(GrdHead grdhead, double rph, double** kernel, int num_thread)
{
    cout<<"new kernel of plane to surface: "<<num_thread<<" threads are used"<<endl;
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    int datanum=grdhead.rows*grdhead.cols;
    //omp parallel
    omp_set_num_threads(num_thread);
    // ProgressBar bar0(grdhead.rows);
    MultiProgressBar multibar(grdhead.rows,COLOR_BAR_BLUE);
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
        #pragma omp critical
        multibar.Update();
    }cout<<"\n";
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
int Getkernel_p2p_new(GrdHead grdhead, double rph, double* kernel_firstRow, int num_thread)
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
 * @brief Upward continuation: plane to surface(space domain)
 * 
 * @param inputfilename 
 * @param outputfilename 
 * @param height1 Elevation of the observation field data, which read from file named inputfilename
 * @param topoFile Topography file contains x,y,z of calculation points.
 * @param extNum 
 * @param num_thread 
 * @param isProgress 
 * @param filename_exact 
 */
void UWC_p2s(string inputfilename,string outputfilename,double height1,string topoFile, int extNum, int num_thread, bool isProgress,string filename_exact)
{
    cout << "***************************************************\n";
    cout << " Upward continuation from plane to surface:"<<height1<<"-> Topography"<<endl;
    cout << "***************************************************\n";
    //1. read field grd data
    GrdHead grdhead;
    double* indata = NULL;
    if (!(indata = ReadGrd(inputfilename, grdhead,extNum)))return ;
    int modelnum = grdhead.rows*grdhead.cols;
    //2. read topography grd data
    GrdHead grdhead_topo;
    double* topo = NULL;
    if (!(topo = ReadGrd(topoFile, grdhead_topo,extNum)))return ;
    if((grdhead_topo.cols != grdhead.cols) && (grdhead_topo.rows != grdhead.rows))
    {
        cout<<"["<<RED<<"Error"<<COLOR_DEFALUT<<"]: Input field and upward continued topography data have different dimensions"
            <<endl;
        exit(0);
    }
    //3. kernel mat.
    double** G=new double*[modelnum];
    for(int i=0; i<modelnum; i++)
    {
        G[i]=new double[modelnum];
    }
    cout<<"calculating kernal matrix"<<endl;
    Getkernel_p2s_new(grdhead,height1,topo,G,num_thread);
    // //4. outdata
    double *outdata=new double[modelnum];
    //compute UWC
    cout<<"calculating uwc from plane to surface"<<endl;
    UWC_Gij(outdata,G,indata,modelnum,num_thread);

    cout << "Finished\n";
    
    //4. write result
    string ext_outputfile=Path_GetExtName(outputfilename);
    if(ext_outputfile=="grd" || ext_outputfile=="GRD")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        if (!SaveGrd(outputfilename, grdhead, outdata,extNum))return ;
    }else if(ext_outputfile=="vtk" || ext_outputfile=="VTK")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2VTK(outputfilename,grdhead,outdata);
    }else if(ext_outputfile=="xyz" || ext_outputfile=="txt" || ext_outputfile=="dat")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2xyz(outputfilename, grdhead, outdata,extNum);
    }else if(ext_outputfile=="nc")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2netCDF(outputfilename, grdhead, outdata,extNum);
    }
    else
    {
        cout<<RED<<"The output file format is not supported: "<<ext_outputfile<<COLOR_DEFALUT<<endl;
    }
    //5. compute error if exact solution is given
    if(filename_exact!="")
    {
        double* exactsolution = NULL;
        if (!(exactsolution = ReadGrd(filename_exact, grdhead,extNum)))return ;
        for(int i=0;i<modelnum;i++)
        {
            outdata[i]=outdata[i]-exactsolution[i];
        }
        if (!SaveGrd(outputfilename+"_error.grd", grdhead, outdata,extNum))return ;
        delete[] exactsolution;
    }
    
    //delete data array
    delete[] indata;
    delete[] outdata;
    delete[] topo;
    for(int i=0; i<modelnum; i++)
    {
        delete[] G[i];
    }delete[] G;
}

/**
 * @brief Upward continuation: plane to plane(space domain)
 * 
 * @param inputfilename 
 * @param outputfilename 
 * @param height1 
 * @param height2 
 * @param extNum 
 * @param num_thread 
 * @param isProgress 
 * @param filename_exact 
 */
void UWC_p2p(string inputfilename,string outputfilename,double height1,double height2,int extNum,int num_thread, bool isProgress,string filename_exact)
{
    cout << "***************************************************\n";
    cout << " Upward continuation from plane to plane:"<<height1<<"->"<<height2<<endl;
    cout << "***************************************************\n";
    //1. read grd data
    GrdHead grdhead;
    double* indata = NULL;
    if (!(indata = ReadGrd(inputfilename, grdhead,extNum)))return ;
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
    //4. write result
    // string basename_infile=Path_GetBaseName(inputfilename);
    string ext_outputfile=Path_GetExtName(outputfilename);
    if(ext_outputfile=="grd" || ext_outputfile=="GRD")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        if (!SaveGrd(outputfilename, grdhead, outdata,extNum))return ;
    }else if(ext_outputfile=="vtk" || ext_outputfile=="VTK")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2VTK(outputfilename,grdhead,outdata,height2);
        string filename_outfile=Path_GetFileName(outputfilename);
        SaveGrd2VTK(filename_outfile+"_origin.vtk",grdhead,indata,height1);//save the origin data as well
    }else if(ext_outputfile=="xyz" || ext_outputfile=="txt" || ext_outputfile=="dat")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2xyz(outputfilename, grdhead, outdata,extNum);
    }else if(ext_outputfile=="nc")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2netCDF(outputfilename, grdhead, outdata,extNum);
    }
    else
    {
        cout<<RED<<"The output file format is not supported: "<<ext_outputfile<<COLOR_DEFALUT<<endl;
    }
    
    //5. compute error if exact solution is given
    if(filename_exact!="")
    {
        double* exactsolution = NULL;
        if (!(exactsolution = ReadGrd(filename_exact, grdhead,extNum)))return ;
        for(int i=0;i<modelnum;i++)
        {
            outdata[i]=outdata[i]-exactsolution[i];
        }
        if (!SaveGrd(outputfilename+"_error.grd", grdhead, outdata,extNum))return ;
        delete[] exactsolution;
    }
    
    //delete data array
    delete[] indata;
    delete[] outdata;
    delete[] G_firstRow;
}
 
/**
 * @brief Upward continuation: plane to plane(frequency domain)
 * 
 * @param inputfilename 
 * @param outputfilename 
 * @param height1 
 * @param height2 
 * @param extNum 
 * @param filename_exact 
 */
void UWC_p2p_f(string inputfilename,string outputfilename,double height1,double height2,int extNum,string filename_exact)
{
    cout << "*********************************************************\n";
    cout << "Upward continuation from plane to plane: frequency domain\n";
    cout << "*********************************************************\n";
    //1. read grd data
    GrdHead grdhead;
    double* indata = NULL;
    if (!(indata = ReadGrd(inputfilename, grdhead,extNum)))return ;
    int modelnum = grdhead.rows*grdhead.cols;

    //2. height of continuation
    double rph=fabs(height1-height2);

    //3. call uwc function in frequency domain
    double *outdata=new double[modelnum];
    UWC_p2p_f(indata,outdata,grdhead,rph);

    //4. write result
    // if (!SaveGrd(outputfilename, grdhead, outdata,extNum))return ;
    string ext_outputfile=Path_GetExtName(outputfilename);
    if(ext_outputfile=="grd" || ext_outputfile=="GRD")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        if (!SaveGrd(outputfilename, grdhead, outdata,extNum))return ;
    }else if(ext_outputfile=="vtk" || ext_outputfile=="VTK")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2VTK(outputfilename,grdhead,outdata,height2);
        string filename_outfile=Path_GetFileName(outputfilename);
        SaveGrd2VTK(filename_outfile+"_origin.vtk",grdhead,indata,height1);//save the origin data as well
    }else if(ext_outputfile=="xyz" || ext_outputfile=="txt" || ext_outputfile=="dat")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2xyz(outputfilename, grdhead, outdata,extNum);
    }else if(ext_outputfile=="nc")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2netCDF(outputfilename, grdhead, outdata,extNum);
    }
    else
    {
        cout<<RED<<"The output file format is not supported: "<<ext_outputfile<<COLOR_DEFALUT<<endl;
    }
    //5. compute error if exact solution is given
    if(filename_exact!="")
    {
        double* exactsolution = NULL;
        if (!(exactsolution = ReadGrd(filename_exact, grdhead,extNum)))return ;
        for(int i=0;i<modelnum;i++)
        {
            outdata[i]=outdata[i]-exactsolution[i];
        }
        if (!SaveGrd(outputfilename+"_error.grd", grdhead, outdata,extNum))return ;
        delete[] exactsolution;
    }
    //delete data array
    delete[] indata;
    delete[] outdata;
    
}

/**
 * @brief Calculate upward continuation: plane to plane(frequency domain)
 * 
 * @param inputdata 
 * @param result 
 * @param grdhead 
 * @param rph 
 * @return int 
 */
int UWC_p2p_f(double* inputdata,double* result,GrdHead grdhead,double rph)
{
    double tp = 2 * PI;
	int number_x=grdhead.cols;
    int number_y=grdhead.rows;
	int N = number_x*number_y;
	double* originaldata = inputdata;
	fftw_complex* out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	FFT2d(number_y, number_x, originaldata, out_c);
	double *bounds = grdhead.bounds;
	double dy = fabs(bounds[3] - bounds[2]) / (number_y - 1), dx = fabs(bounds[1] - bounds[0]) / (number_x - 1);
	double en = dy*(number_y), em = dx*(number_x);
	double u, v;
	double factor = -tp*rph;
	double factor_conti;
	int number_x_half = (int)(number_x / 2);
	int number_y_half = (int)(number_y / 2);
	int index_ij;

	for (int i = 0; i <= number_y_half; i++)
	{
		v = i / en;
		for (int j = 0; j <= number_x_half; j++)
		{
			u = j / em;
			index_ij = i*number_x + j;
			factor_conti = exp(factor*sqrt(u*u + v*v));
			out_c[index_ij][0] = out_c[index_ij][0] * factor_conti;
			out_c[index_ij][1] = out_c[index_ij][1] * factor_conti;
		}

		for (int j = number_x_half + 1; j < number_x; j++)
		{
			u = (j - number_x) / em;
			index_ij = i*number_x + j;
			factor_conti = exp(factor*sqrt(u*u + v*v));
			out_c[index_ij][0] = out_c[index_ij][0] * factor_conti;
			out_c[index_ij][1] = out_c[index_ij][1] * factor_conti;
		}
	}
	for (int i = number_y_half + 1; i <number_y; i++)
	{
		v = (i - number_y) / en;
		for (int j = 0; j <= number_x_half; j++)
		{
			u = j / em;
			index_ij = i*number_x + j;
			factor_conti = exp(factor*sqrt(u*u + v*v));
			out_c[index_ij][0] = out_c[index_ij][0] * factor_conti;
			out_c[index_ij][1] = out_c[index_ij][1] * factor_conti;
		}

		for (int j = number_x_half + 1; j < number_x; j++)
		{
			u = (j - number_x) / em;
			index_ij = i*number_x + j;
			factor_conti = exp(factor*sqrt(u*u + v*v));
			out_c[index_ij][0] = out_c[index_ij][0] * factor_conti;
			out_c[index_ij][1] = out_c[index_ij][1] * factor_conti;
		}
	}
	IFFT2d(number_y, number_x, out_c, result);
    return 1;
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
void UWC_Gij(double* b,double** G,double* x, int modelnum,int num_thread)
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
void UWC_Gij(double* b, double* G,double* x, GrdHead grdhead, int num_thread)
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
 * @brief Downward continuation: uneven surface to plane(space domain)
 * 
 * @param inputfilename Grid file of the observation field.
 * @param outputfilename 
 * @param topoFile Topography grid file of the observation field. 
 * @param height2 Elevation of the continuation plane.
 * @param extNum 
 * @param DWC_parameter 
 * @param DWC_method 
 * @param num_thread 
 * @param isProgress 
 * @param filename_exact 
 */
void DWC_s2p(string inputfilename,string outputfilename,string topoFile,
    double height2,int extNum,double DWC_parameter,int DWC_method, int num_thread, bool isProgress,
    string filename_exact)
{
    cout << "***************************************************\n";
    cout << " Downward continuation from surface to plane: topography"<<"->"<<height2<<endl;
    cout << "***************************************************\n";
    // double TRP = 0.01;
    int kmax = 500;
    const double daierta = 0.00001;										//Iterative termination condition
    cout<<"The maximum iteration number: "<<kmax<<", Iterative termination condition: "<<daierta<<endl;
    //1. read grd data
    GrdHead grdhead;
    double* indata = NULL;
    if (!(indata = ReadGrd(inputfilename, grdhead,extNum)))return ;
    int modelnum = grdhead.rows*grdhead.cols;
    double* exactsolution = NULL;
    if(filename_exact!="")
    {
        if (!(exactsolution = ReadGrd(filename_exact, grdhead,extNum)))return ;
    }
    //2. read topography grd data
    GrdHead grdhead_topo;
    double* topo = NULL;
    if (!(topo = ReadGrd(topoFile, grdhead_topo,extNum)))return ;
    if((grdhead_topo.cols != grdhead.cols) && (grdhead_topo.rows != grdhead.rows))
    {
        cout<<"["<<RED<<"Error"<<COLOR_DEFALUT<<"]: Input field and upward continued topography data have different dimensions"
            <<endl;
        exit(0);
    }
    //2. height of continuation
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    // cout << "*****************************************************\n";
    // cout << "       Downward continuation from plane to plane(space)     \n";
    cout <<"Data size: "<<grdhead.rows<<" X "<<grdhead.cols<<endl;
    cout <<"Size of kernel mat.: "<<grdhead.rows*grdhead.cols<<" X "<<grdhead.rows*grdhead.cols<<endl;
    cout <<"Memory of kernel mat.: "<<sizeof(double)*pow(grdhead.rows*grdhead.cols,2)/pow(1024.0,2)<< "Mb"<<endl;
    cout<<"Point spacing dx: "<<dx<<"  dy: "<<dy<<endl;
    cout <<"Downward continue: ";
    cout<<""<<"\033[31m"<<((grdhead_topo.bounds[5]-height2)/dx)<<"\033[0m"<<" point spacing";
    cout<<" to "<<"\033[31m"<<((grdhead_topo.bounds[4]-height2)/dx)<<"\033[0m"<<" point spacing"<<endl;
    cout << "*****************************************************\n";
    
    //3. kernel mat.
    double** G=new double*[modelnum];
    for(int i=0; i<modelnum; i++)
    {
        G[i]=new double[modelnum];
    }
    cout<<"calculating kernal matrix"<<endl;
    Getkernel_p2s_new(grdhead,height2,topo,G,num_thread);
    //4. outdata
    double *dataout=new double[modelnum];
    //call downward continuation subfunction
    switch(DWC_method)
    {
        case DWC_LANDWEBER:
            DWC_s2p_LandweberIter(G,dataout,indata,grdhead,extNum,num_thread,outputfilename,DWC_parameter,exactsolution);
            break;
        default:
            cout<<RED<<"The method of downward continuation is not available yet."<<COLOR_DEFALUT<<endl;
            exit(0);
    }

    //4. write result
    // if (!SaveGrd(outputfilename, grdhead, dataout,extNum))return ;
    string ext_outputfile=Path_GetExtName(outputfilename);
    if(ext_outputfile=="grd" || ext_outputfile=="GRD")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        if (!SaveGrd(outputfilename, grdhead, dataout,extNum))return ;
    }else if(ext_outputfile=="vtk" || ext_outputfile=="VTK")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2VTK(outputfilename,grdhead,dataout,height2);
        string filename_outfile=Path_GetFileName(outputfilename);
        SaveGrd2VTK_topo(filename_outfile+"_origin.vtk",grdhead,indata,topo);//save the origin data as well
    }else if(ext_outputfile=="xyz" || ext_outputfile=="txt" || ext_outputfile=="dat")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2xyz(outputfilename, grdhead, dataout,extNum);
    }else if(ext_outputfile=="nc")
    {
        cout<<GREEN<<"Output file format is : "<<ext_outputfile<<COLOR_DEFALUT<<endl;
        SaveGrd2netCDF(outputfilename, grdhead, dataout,extNum);
    }
    else
    {
        cout<<RED<<"The output file format is not supported: "<<ext_outputfile<<COLOR_DEFALUT<<endl;
    }

    if(filename_exact!="")
    {
        cout<<"calculate error according to exact solution file"<<endl;
        for(int i=0;i<modelnum;i++)
        {
            dataout[i]=dataout[i]-exactsolution[i];
        }
        if (!SaveGrd(outputfilename+"_error.grd", grdhead, dataout,extNum))return ;
        delete[] exactsolution;
    }
    //delete array
    delete[] indata;
    delete[] dataout;
    delete[] topo;
    // if(!exactsolution)delete[] exactsolution;
    for (int i = 0; i < modelnum; i++)
    {
        delete[] G[i];
    }delete[] G;
}

/**
 * @brief Downward continuation: plane to plane(space domain)
 * 
 * @param inputfilename 
 * @param outputfilename 
 * @param height1 Elevation of the observation field.
 * @param height2 Elevation of the continuation field.
 * @param extNum 
 * @param DWC_parameter 
 * @param DWC_method 
 * @param num_thread 
 * @param isProgress 
 * @param filename_exact 
 */
void DWC_p2p(string inputfilename,string outputfilename,double height1,double height2,int extNum,double DWC_parameter,int DWC_method,int num_thread, bool isProgress,string filename_exact)
{
    cout << "***************************************************\n";
    cout << " Downward continuation from plane to plane:"<<height1<<"->"<<height2<<endl;
    cout << "***************************************************\n";
    // double TRP = 0.01;
    int kmax = 500;

    //1. read grd data
    GrdHead grdhead;
    double* indata = NULL;
    if (!(indata = ReadGrd(inputfilename, grdhead,extNum)))return ;
    int modelnum = grdhead.rows*grdhead.cols;

    double* exactsolution = NULL;
    if(filename_exact!="")
    {
        cout<<"Read exact solution file"<<endl;
        if (!(exactsolution = ReadGrd(filename_exact, grdhead,extNum)))return ;
    }
    //2. height of continuation
    double rph=fabs(height1-height2);
    //print general information
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    // cout << "*****************************************************\n";
    // cout << "       Downward continuation from plane to plane(space)     \n";
    cout<<num_thread<<" threads"<<endl;
    cout <<"Size of data: "<<grdhead.rows<<" X "<<grdhead.cols<<endl;
    cout <<"Size of kernel mat.: "<<grdhead.rows*grdhead.cols<<" X "<<grdhead.rows*grdhead.cols<<endl;
    cout <<"Only save the first row of kernel mat: "<<sizeof(double)*pow(grdhead.rows*grdhead.cols,1)/pow(1024.0,2)<< "Mb"<<endl;
    cout <<"Downward continuation:  ";
    cout <<"\033[31m";       
    cout <<(rph/dx);
    cout <<"\033[0m";
    cout <<"  point spacing\n";
    cout << "*****************************************************\n";

    double** G=new double*[modelnum];
    for(int i=0; i<modelnum; i++)
    {
        G[i]=new double[modelnum];
    }

    Getkernel_p2p_new(grdhead,rph,G,num_thread);
    //4. outdata
    double *dataout=new double[modelnum];

    switch(DWC_method)
    {
        case DWC_LANDWEBER:
            DWC_s2p_LandweberIter(G,dataout,indata,grdhead,extNum,num_thread,outputfilename,DWC_parameter,exactsolution);
            break;
        default:
            cout<<RED<<"The downward continuation method is not available yet"<<COLOR_DEFALUT<<endl;
            exit(0);
    }
    //5. write result
    if (!SaveGrd(outputfilename, grdhead, dataout,extNum))return ;
    //5. compute error if exact solution is given
    if(filename_exact!="")
    {
        cout<<"calculate error according to exact solution file"<<endl;
        for(int i=0;i<modelnum;i++)
        {
            dataout[i]=dataout[i]-exactsolution[i];
        }
        if (!SaveGrd(outputfilename+"_error.grd", grdhead, dataout,extNum))return ;
        delete[] exactsolution;
    }
    //delete 
    delete[] indata;
    delete[] dataout;
    if(!exactsolution)delete[] exactsolution;
}

/**
 * @brief Calculate downward continuation from plane to plane using Landweber iteration method. 
 * 
 * Solve linear equations of \f$ \mathbf{G}x=b \f$ for \f$ b \f$.
 * 
 * @param G First row of the kernel matrix.
 * @param x Result
 * @param b Observation field
 * @param grdhead 
 * @param extNum 
 * @param num_thread 
 * @param outputfile 
 * @param iter_number 
 * @param ExactSolution Read and compare with exactsolution if it is known. e.g. synethetic test.
 */
void DWC_p2p_LandweberIter(double* G,double* x, double* b,GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number,double* ExactSolution)
{
    int kmax=100;
    if(iter_number>0)kmax=int(iter_number);
    cout<<"Landweber iteration method: kmax="<<kmax<<endl;
    double REmin=1E-8;
    int modelnum=grdhead.cols*grdhead.rows;
    for(int i=0;i<modelnum;i++)x[i]=b[i];
    double* residual=new double[modelnum];
    for(int i=0;i<modelnum;i++)residual[i]=0;
    double* x_uwc=new double[modelnum];
    double* tempArray=new double[modelnum];

    string path_tempResult=Path_GetFileName(outputfile)+"_Landweber";
    string command="mkdir "+path_tempResult;
    if(system(command.c_str()))
    {
        command="rm "+path_tempResult+"/*.grd";
        if(!system(command.c_str()))
        {
            cout<<GREEN<<"Clean .grd file in temporary folder: "<<COLOR_DEFALUT<<path_tempResult<<endl;
        }
    }else{
        cout<<GREEN<<"Create directory to save temporary resut: "<<COLOR_DEFALUT<<path_tempResult<<endl;
    }
    string filename_log=path_tempResult+"/LandweberInteration_log.txt";
    ofstream fout_log(filename_log);
    if(!fout_log)
    {
        cout<<RED<<"[DWC_p2p_LandweberIter]Open file failed: "<<filename_log<<COLOR_DEFALUT<<endl;
        exit(0);
    }
    if(ExactSolution)
    {
        fout_log<<"Iter\tNorm2(Ax-b)/Norm2(b)\tNorm2(x)\tNorm2(grad_x)\tNorm2(Exact-x)"<<endl;
    }
    else{
        fout_log<<"Iter\tNorm2(Ax-b)/Norm2(b)\tNorm2(x)\tNorm2(grad_x)"<<endl;
    }
    //iteration
    vector<double>bar_left;bar_left.push_back(0);bar_left.push_back(0);
    vector<double>bar_right;bar_right.push_back(double(kmax));bar_right.push_back(log10(REmin));
    vector<double>bar_pos; bar_pos.push_back(bar_left[0]);bar_pos.push_back(bar_left[1]);
    vector<string>bar_title; bar_title.push_back("Iteration");bar_title.push_back("RE(log10):");
    MultiProgressBar multiBar(bar_left,bar_right,bar_title);
    for(int k=0;k<kmax;k++)
    {
        for(int i=0;i<modelnum;i++)x[i]=x[i]+residual[i];//xn=x_n1+residual
        UWC_Gij(x_uwc,G,x,grdhead,num_thread);//x_uwc=A*xn
        for(int i=0;i<modelnum;i++)residual[i]=b[i]-x_uwc[i];//b-A*xn
        UWC_Gij(x_uwc,G,residual,grdhead,num_thread);//residual=A*(b-A*xn)
        for(int i=0;i<modelnum;i++)residual[i]=x_uwc[i];

        double norm2_residual=Norm2(residual,modelnum);
        double norm2_x=Norm2(x,modelnum);
        double norm2_gradient=Norm2_Gradient(x,grdhead);
        double RE=norm2_residual/Norm2(b,modelnum);
        fout_log<<k<<"\t"<<RE<<"\t"<<norm2_x<<"\t"<<norm2_gradient;
        if(ExactSolution)
        {
            for(int i=0;i<modelnum;i++)tempArray[i]=ExactSolution[i]-x[i];
            fout_log<<"\t"<<Norm2(tempArray,modelnum)<<endl;
        }else{
            fout_log<<endl;
        }
        //5. write temporary result
        if (!SaveGrd(path_tempResult+"/"+std::to_string(k)+".grd", grdhead, x,extNum,false,false))return ;
        SaveGrd2VTK(path_tempResult+"/result_"+std::to_string(k)+".vtk",grdhead,x);
        //update progressbar
        bar_pos[0]=k+1;
        bar_pos[1]=log10(RE);
        multiBar.Update(bar_pos);
    }
    fout_log.close();
    delete[] residual;
    delete[] x_uwc;
    delete[] tempArray;
}

/**
 * @brief Calculate downward continuation of potential field from surface to a plane.
 * 
 * @param G The complete kernel matrix.
 * @param x Result
 * @param b Observation field data.
 * @param grdhead 
 * @param extNum 
 * @param num_thread 
 * @param outputfile 
 * @param iter_number 
 * @param ExactSolution 
 */
void DWC_s2p_LandweberIter(double** G,double* x, double* b,
    GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number,
    double* ExactSolution)
{
    int kmax=100;
    if(iter_number>0)kmax=int(iter_number);
    cout<<"Landweber iteration method: kmax="<<kmax<<endl;
    double REmin=1E-8;
    int modelnum=grdhead.cols*grdhead.rows;
    for(int i=0;i<modelnum;i++)x[i]=b[i];
    double* residual=new double[modelnum];
    for(int i=0;i<modelnum;i++)residual[i]=0;
    double* x_uwc=new double[modelnum];
    double* tempArray=new double[modelnum];

    string path_tempResult=Path_GetFileName(outputfile)+"_Landweber";
    string command="mkdir "+path_tempResult;
    if(system(command.c_str()))
    {
        command="rm "+path_tempResult+"/*.grd";
        if(!system(command.c_str()))
        {
            cout<<GREEN<<"Clean .grd file in temporary folder: "<<COLOR_DEFALUT<<path_tempResult<<endl;
        }
    }else{
        cout<<GREEN<<"Create directory to save temporary resut: "<<COLOR_DEFALUT<<path_tempResult<<endl;
    }

    string filename_log=path_tempResult+"/LandweberIteration_log.txt";
    ofstream fout_log(filename_log);
    if(!fout_log)
    {
        cout<<RED<<"[DWC_p2p_ItegrationIter]Open file failed: "<<filename_log<<COLOR_DEFALUT<<endl;
        exit(0);
    }
    if(ExactSolution)
    {
        fout_log<<"Iter\tNorm2(Ax-b)/Norm2(b)\tNorm2(x)\tNorm2(grad_x)\tNorm2(Exact-x)"<<endl;
    }
    else{
        fout_log<<"Iter\tNorm2(Ax-b)/Norm2(b)\tNorm2(x)\tNorm2(grad_x)"<<endl;
    }
    //iteration
    vector<double>bar_left;bar_left.push_back(0);bar_left.push_back(0);
    vector<double>bar_right;bar_right.push_back(double(kmax));bar_right.push_back(log10(REmin));
    vector<double>bar_pos; bar_pos.push_back(bar_left[0]);bar_pos.push_back(bar_left[1]);
    vector<string>bar_title; bar_title.push_back("Iteration");bar_title.push_back("RE(log10):");
    MultiProgressBar multiBar(bar_left,bar_right,bar_title);
    for(int k=0;k<kmax;k++)
    {
        for(int i=0;i<modelnum;i++)x[i]=x[i]+residual[i];//xn=xn1+residual
        UWC_Gij(x_uwc,G,x,modelnum,num_thread);//x_uwc=A*xn
        for(int i=0;i<modelnum;i++)residual[i]=b[i]-x_uwc[i];//residual=b-A*xn
        //5. write temporary result
        if (!SaveGrd(path_tempResult+"/"+std::to_string(k)+".grd", grdhead, x,extNum,false,false))return ;
        SaveGrd2VTK(path_tempResult+"/result_"+std::to_string(k)+".vtk",grdhead,x,0,residual);
        
        UWC_Gij(x_uwc,G,residual,modelnum,num_thread);//x_uwc=A*residual
        for(int i=0;i<modelnum;i++)residual[i]=x_uwc[i];
        
        double norm2_residual=Norm2(residual,modelnum);
        double norm2_x=Norm2(x,modelnum);
        double norm2_gradient=Norm2_Gradient(x,grdhead);
        double RE=norm2_residual/Norm2(b,modelnum);
        fout_log<<k<<"\t"<<RE<<"\t"<<norm2_x<<"\t"<<norm2_gradient;
        if(ExactSolution)
        {
            for(int i=0;i<modelnum;i++)tempArray[i]=ExactSolution[i]-x[i];
            fout_log<<"\t"<<Norm2(tempArray,modelnum)<<endl;
        }else{
            fout_log<<endl;
        }
        
        if(RE<REmin)break;
        //update progressbar
        bar_pos[0]=k+1;
        bar_pos[1]=log10(RE);
        multiBar.Update(bar_pos);
    }
    fout_log.close();
    delete[] residual;
    delete[] x_uwc;
    delete[] tempArray;
}
