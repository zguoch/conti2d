
#include "Conti2D.h"

string Path_GetFileName(string filepath)
{
    int pos_ext=filepath.find_last_of('.');
    string filename=filepath.substr(0,pos_ext);
    return filename;
}
string Path_GetPath(string filepath)
{
    int pos_ext=filepath.find_last_of('/');
    string path=filepath.substr(0,pos_ext);
    return path;
}
string Path_GetExtName(string filepath)
{
    int pos_ext=filepath.find_last_of('.');
    string extname=filepath.substr(pos_ext+1,filepath.length());
    return extname;
}
string Path_GetBaseName(string filepath)
{
    int pos_ext=filepath.find_last_of('.');
    string pathname=filepath.substr(0,pos_ext);
    int pos_basename=pathname.find_last_of('/');
    string basename=pathname.substr(pos_basename+1,pathname.length());
    return basename;
}
/*===================================
 Get first row of kernal mat
 ====================================*/
void GetKernalMatrix(GrdHead grdhead, double* G, const double rph, int num_thread )
{
    int index_row=0;
    int index_col=0;
    double xdown=0,ydown=0,xup=0,yup=0;
    double r;
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    double xishu=dx*dy*rph/(2*PI);
    //omp parallel
    omp_set_num_threads(num_thread);
    #pragma omp parallel for 
    for(int irow=0; irow<1; irow++)
    {
        index_row=irow*grdhead.cols;
        yup=irow*dy;
        for(int icol=0; icol<1; icol++)
        {
            xup=icol*dx;
            for(int i=0; i<grdhead.rows; i++)
            {
                index_col=i*grdhead.cols;
                ydown=i*dy;
                for(int j=0; j<grdhead.cols; j++)
                {
                    xdown=j*dx;
                    r=pow((xup-xdown),2.0)+pow((yup-ydown),2.0)+rph*rph;
                    // kernel[index_row][index_col]=xishu/pow(r,1.5);
                    G[index_col]=xishu/pow(r,1.5);
                    index_col++;
                }
            }
            index_row++;
        }
    }

}

/*===================================
 Get first row of kernal mat
 ====================================*/
 void GetKernalMatrix_new(GrdHead grdhead, double* G, const double rph, int num_thread)
 {
    int index_row=0;
    int index_col=0;
    double xdown=0,ydown=0,xup=0,yup=0;
    double r;
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    double dx_2=dx/2.0, dy_2=dy/2.0;
    double xishu=dx*dy*rph/(2*PI);
    double x1=0,x2=0,y1=0,y2=0;
    double x1x1=0,x2x2=0,y1y1=0,y2y2=0;//平方
    double x0=0,y0=0;
    double p22=0,p11=0,p12=0,p21=0;
    double rph_2=rph*rph;
    //omp parallel
    omp_set_num_threads(num_thread);
    #pragma omp parallel for 
    for(int irow=0; irow<1/*grdhead.rows*/; irow++)
    {
        index_row=irow*grdhead.cols;
        yup=irow*dy;
        for(int icol=0; icol<1/*grdhead.cols*/; icol++)
        {
            xup=icol*dx;
            for(int i=0; i<grdhead.rows; i++)
            {
                index_col=i*grdhead.cols;
                ydown=i*dy;
                for(int j=0; j<grdhead.cols; j++)
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
                    //-----------
                    // kernel[index_row][index_col]=(p22+p11-p21-p12)/(2*PI);
                    G[index_col]=(p22+p11-p21-p12)/(2*PI);
                    index_col++;
                }
            }
            index_row++;
        }
    }
 }
 
 /*===================================
  Get component Gij from the first row
  ====================================*/
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
    // cout<<"row: "<<i<<"*"<<j<<endl;
    // cout<<"block: "<<row_block<<"*"<<colum_block<<endl;
    // cout<<"coordinate: "<<index_row<<"*"<<index_col<<endl;
    int index_block=abs(row_block-colum_block);         
    int col0=abs(index_row-index_col);                 
    Gij=firstRow[index_block*grdhead.cols+col0];
     return Gij;
 }
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
int Getkernel_p2p(GrdHead grdhead, double rph, double** kernel, int num_thread)
{
    cout<<"old kernel"<<endl;
    int index_row=0;
    int index_col=0;
    double xdown=0,ydown=0,xup=0,yup=0;
    double r;
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    double xishu=dx*dy*rph/(2*PI);
    //omp parallel
    MultiProgressBar bar0(grdhead.rows);
    omp_set_num_threads(num_thread);
    #pragma omp parallel for 
    for(int irow=0; irow<grdhead.rows; irow++)
    {
        index_row=irow*grdhead.cols;
        yup=irow*dy;
        for(int icol=0; icol<grdhead.cols; icol++)
        {
            xup=icol*dx;
            for(int i=0; i<grdhead.rows; i++)
            {
                index_col=i*grdhead.cols;
                ydown=i*dy;
                for(int j=0; j<grdhead.cols; j++)
                {
                    xdown=j*dx;
                    r=pow((xup-xdown),2.0)+pow((yup-ydown),2.0)+rph*rph;
                    kernel[index_row][index_col]=xishu/pow(r,1.5);
                    index_col++;
                }
            }
            index_row++;
        }
        bar0.Update();
    }cout<<"\n";
    
    return 1;
}
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
/*===============================================================================
 Upward continuation: plane to surface(space domain)
 ================================================================================*/
void UWC_p2s(string inputfilename,string outputfilename,
    double height1,string topoFile, int extNum,
    int num_thread, bool isProgress,bool useOldKernel,string filename_exact)
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
    if(useOldKernel)
    {
        cout<<"Old kernel p2s is not added now"<<endl;
        exit(0);
        // Getkernel_p2s(grdhead,rph,G,num_thread);
    }else
    {
        Getkernel_p2s_new(grdhead,height1,topo,G,num_thread);
    }
    // //4. outdata
    double *outdata=new double[modelnum];
    //compute UWC
    cout<<"calculating uwc from plane to surface"<<endl;
    UWC_Gij(outdata,G,indata,modelnum,num_thread);

    cout << "Finished\n";
    
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
/*===============================================================================
 Upward continuation: plane to plane(space domain)
 ================================================================================*/
void UWC_p2p(string inputfilename,string outputfilename,
    double height1,double height2,int extNum,
     int num_thread, bool isProgress,bool useOldKernel,string filename_exact)
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
    if(useOldKernel)
    {
        cout<<"The old kernel function is  depressed "<<endl;
        exit(0);
        // Getkernel_p2p(grdhead,rph,G,num_thread);
    }else
    {
        Getkernel_p2p_new(grdhead,rph,G_firstRow,num_thread);
    }
    
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
/*===============================================================================
 Upward continuation: plane to plane(frequency domain)
 ================================================================================*/
void UWC_p2p_f(string inputfilename,string outputfilename,double height1,double height2,
                int extNum,string filename_exact)
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
int UWC_p2p_f(double* inputdata,double* result,GrdHead grdhead,double rph)
{
    double tp = 2 * PI;
	int number_x=grdhead.cols;
    int number_y=grdhead.rows;
	int N = number_x*number_y;
	double* originaldata = inputdata;
	fftw_complex* out_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)* N);
	FFT2d(number_y, number_x, originaldata, out_c);		
	double dy;
	dy = fabs(grdhead.bounds[3] - grdhead.bounds[2]) / (number_y - 1);
	double dx = fabs(grdhead.bounds[1] - grdhead.bounds[0]) / (number_x - 1);
	double en = dy * (number_y);
	double em = dx * (number_x);
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

int Getkernel_u2p(GrdHead grdhead, double* terrain1, double height2, double** kernel, int num_thread)
{
    int index_row=0;
    int index_col=0;
    double xdown=0,ydown=0,xup=0,yup=0;
    double r;
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    double xishu=dx*dy/(2*PI);
    double height_bottom=0;
    double rph=0;
    //omp parallel
    omp_set_num_threads(num_thread);
    #pragma omp parallel for 
    for(int irow=0; irow<grdhead.rows; irow++)
    {
        index_row=irow*grdhead.cols;
        yup=irow*dy;
        for(int icol=0; icol<grdhead.cols; icol++)
        {
            xup=icol*dx;
            for(int i=0; i<grdhead.rows; i++)
            {
                index_col=i*grdhead.cols;
                ydown=i*dy;
                for(int j=0; j<grdhead.cols; j++)
                {
                    height_bottom=terrain1[j+grdhead.cols*i];
                    rph=(height2-height_bottom);
                    xdown=j*dx;
                    r=pow((xup-xdown),2.0)+pow((yup-ydown),2.0)+rph*rph;
                    kernel[index_row][index_col]=xishu*rph/pow(r,1.5);
                    index_col++;
                }
            }
            index_row++;
        }
    }
    return 1;
}

//refer (Syberg, 1971; H. Granser, 1983; Bhattacharyya and Chan, 1977)
// maybe wrong, not used in the program
int Getkernel_u2p_new(GrdHead grdhead, double* terrain1, double height2,
                    double* nx,double* ny,double* nz, 
                    double** kernel, int num_thread)
{
    int index_row=0;
    int index_col=0;
    double xdown=0,ydown=0,xup=0,yup=0;
    double r;
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    double dx_2=dx/2.0, dy_2=dy/2.0;
    double x1=0,x2=0,y1=0,y2=0;
    double x1x1=0,x2x2=0,y1y1=0,y2y2=0;//
    double fabs_x1,fabs_x2,fabs_y1,fabs_y2;//
    double x0=0,y0=0;
    double p22=0,p11=0,p12=0,p21=0;
    double rph_2=0;
    double rph=0;
    double height_bottom=0;
    double P1=0,P2=0,P3=0;
    int index_nxyz=0;
    //omp parallel
    omp_set_num_threads(num_thread);
    #pragma omp parallel for 
    for(int irow=0; irow<grdhead.rows; irow++)
    {
        index_row=irow*grdhead.cols;
        yup=irow*dy;
        for(int icol=0; icol<grdhead.cols; icol++)
        {
            xup=icol*dx;
            for(int i=0; i<grdhead.rows; i++)
            {
                index_col=i*grdhead.cols;
                ydown=i*dy;
                for(int j=0; j<grdhead.cols; j++)
                {
                    height_bottom=terrain1[j+grdhead.cols*i];
                    rph=(height2-height_bottom);
                    rph_2=rph*rph;
                    xdown=j*dx;
                    index_nxyz=j+grdhead.cols*i;
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
                    //P1--------
                    p22=atan2(x2*y2,rph*sqrt(rph_2+x2x2+y2y2));
                    p11=atan2(x1*y1,rph*sqrt(rph_2+x1x1+y1y1));
                    p12=atan2(x1*y2,rph*sqrt(rph_2+x1x1+y2y2));
                    p21=atan2(x2*y1,rph*sqrt(rph_2+x2x2+y1y1));
                    P1=(p22+p11-p21-p12)*nz[index_nxyz];//
                    //-----------
                    fabs_y1=fabs(y1);
                    fabs_y2=fabs(y2);
                    //P2--------
                    p22=atan2(sqrt(x2x2+y2y2+rph_2),fabs_y2)/fabs_y2;
                    p11=atan2(sqrt(x1x1+y1y1+rph_2),fabs_y1)/fabs_y1;
                    p12=atan2(sqrt(x1x1+y2y2+rph_2),fabs_y2)/fabs_y2;
                    p21=atan2(sqrt(x2x2+y1y1+rph_2),fabs_y1)/fabs_y1;
                    P2=(p22+p11-p21-p12)*nx[index_nxyz];
                    //-----------
                    fabs_x1=fabs(x1);
                    fabs_x2=fabs(x2);
                    //P3--------
                    p22=atan2(sqrt(x2x2+y2y2+rph_2),fabs_x2)/fabs_x2;
                    p11=atan2(sqrt(x1x1+y1y1+rph_2),fabs_x1)/fabs_x1;
                    p12=atan2(sqrt(x1x1+y2y2+rph_2),fabs_x1)/fabs_x1;
                    p21=atan2(sqrt(x2x2+y1y1+rph_2),fabs_x2)/fabs_x2;
                    P3=(p22+p11-p21-p12)*ny[index_nxyz];
                    //-----------
                    kernel[index_row][index_col]=(P1+P2+P3)/(2*PI);
                    index_col++;
                }
            }
            index_row++;
        }
    }
    
    return 1;
}

/*===============================================================================
 Read surfer grd text file
 Delete data array,you can see function: bool DeleteGrdData(GRD_surfer grddata)
 ================================================================================*/
double* ReadGrd(string filename, GrdHead& grdhead,int extNum)
{
    ifstream fin;
    fin.open(filename, ios::in);
    if (!fin)
    {
        cout<<RED << "Open file false: "<<COLOR_DEFALUT << filename << "\n";
        return NULL;
    }
    string dsaa;
    //read head
    fin >> dsaa;
    if(dsaa=="DSAA")
    {
        cout<<"ASCII"<<endl;
    }else{
        cout<<BLUE<<"Input grd file is Binary\n"<<COLOR_DEFALUT
            <<RED<<"Please convert it to ascii: \n"<<COLOR_DEFALUT
            <<GREEN<<"https://convert.goldensoftware.com/Application/Conversion"
            <<COLOR_DEFALUT<<endl;
        exit(0);
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

/*===================================
 Save Grd data as surfer grd text file
 ====================================*/
bool SaveGrd(string filename, GrdHead grdhead, double* data,int extNum, bool savexxyz,bool isInfo)
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
bool SaveGrd2xyz(string filename, GrdHead grdhead, double* data,int extNum,bool isInfo)
{
    //save xyz format as well
    ofstream fxyz(filename);
    if (!fxyz)
    {
        cout<<RED << "Open file false: "<<COLOR_DEFALUT << filename << "\n";
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
/*===================================
 Save Grd data as netCDF format
 ====================================*/
bool SaveGrd2netCDF(string filename, GrdHead grdhead,double* data,int extNum,bool isInfo)
{
     double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
     double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
     const int NDIMS=2;
     int ncid, x_dimid, y_dimid, varid,x_varid,y_varid;
     int dimids[NDIMS];
     // nx and ny
     int NX=grdhead.rows-2*extNum;
     int NY=grdhead.cols-2*extNum;
     // data
     double *x=new double[NX];
     double *y=new double[NY];
     double** data_out=new double* [NX];
     for(int i=0;i<NX;i++)
     {
         data_out[i]=new double[NY];
     }
     // attribute, e.g. unit
     char field_units[]="unit of potential field";
     char name_x[]="x";
     char name_Y[]="Y";
     int retval;/* Error handling. */
     for(int i=0;i<NX;i++)
     {
         x[i]=grdhead.bounds[0]+i*dx;
     }
     for(int i=0;i<NY;i++)
     {
         y[i]=grdhead.bounds[2]+i*dy;
     }
     //write data
     for (int i = extNum; i < grdhead.rows-extNum; i++)
     {
         for (int j = extNum; j < grdhead.cols-extNum; j++)
         {
             data_out[i-extNum][j-extNum]=data[j+grdhead.cols*i];
         }
     }
     /* Create the  nc file. */
     if ((retval = nc_create(filename.c_str(), NC_CLOBBER, &ncid)))
       ERR(retval);
     /* Define the dimensions. NetCDF will hand back an ID for each. */
    if ((retval = nc_def_dim(ncid, name_x, NX, &x_dimid)))
       ERR(retval);
    if ((retval = nc_def_dim(ncid, name_Y, NY, &y_dimid)))
       ERR(retval);
     /* Define coordinate netCDF variables.*/
    if ((retval = nc_def_var(ncid, name_x, NC_DOUBLE, 1, &x_dimid,&x_varid)))
       ERR(retval);
    if ((retval = nc_def_var(ncid, name_Y, NC_DOUBLE, 1, &y_dimid,&y_varid)))
       ERR(retval);
     /* The dimids array is used to pass the IDs of the dimensions of the variable. */
    dimids[0] = x_dimid;
    dimids[1] = y_dimid;
    /* Define the variable.  */
    if ((retval = nc_def_var(ncid, "Result by conti2d", NC_DOUBLE, NDIMS,
 			    dimids, &varid)))
       ERR(retval);
     /* End define mode. This tells netCDF we are done defining metadata. */
    if ((retval = nc_enddef(ncid)))
       ERR(retval);
     if ((retval = nc_put_var_double(ncid, x_varid, &x[0])))
       ERR(retval);
    if ((retval = nc_put_var_double(ncid, y_varid, &y[0])))
       ERR(retval);
    if ((retval = nc_put_var_double(ncid, varid, &data_out[0][0])))
       ERR(retval);
     /* Close the file. */
    if ((retval = nc_close(ncid)))
       ERR(retval);

     // delete pointer
     delete[] x,y;
     for(int i=0;i<NX;i++)delete[] data_out[i];
     delete[] data_out;
    return true;
}


void UWC_Gji(double* b,double** G,double* x, int modelnum,int num_thread)
{
    omp_set_num_threads(num_thread);
    #pragma omp parallel for shared(b,x,modelnum)
    for (int i = 0; i < modelnum; i++)
    {
        b[i] = 0;
        for (int j = 0; j < modelnum; j++)
        {
            b[i] += G[j][i]*x[j];
        }
    }
    // cout<<"\n";
}

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

void UWC_Gji(double* b, double* G,double* x, GrdHead grdhead, int num_thread)
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
void UWC(double* datain, double* dataout, GrdHead grdhead,double** G)
{
    int modelnum = grdhead.rows*grdhead.cols;

    #pragma omp parallel for 				
    for (int i = 0; i < modelnum; i++)
    {
        dataout[i] = 0;
        for (int j = 0; j < modelnum; j++)
        {
            dataout[i] += G[i][j]*datain[j];
        }

    }

}
/*===============================================================================
 Downward continuation: uneven surface to plane(frequency domain)
 ================================================================================*/
 void DWC_p2p_f(string inputfilename,string outputfilename,double height1,double height2,
int extNum, double TRP,int num_thread, bool isProgress,bool useOldKernel)
 {

    int kmax = 500;
    const double daierta = 0.00001;										//Iterative termination condition
    
    //1. read grd data
    GrdHead grdhead;
    double* indata = NULL;
    if (!(indata = ReadGrd(inputfilename, grdhead,extNum)))return ;
    int modelnum = grdhead.rows*grdhead.cols;

    //2. height of continuation
    double rph=fabs(height1-height2);
    //print general information
    double dx=(grdhead.bounds[1]-grdhead.bounds[0])/(grdhead.cols-1);
    double dy=(grdhead.bounds[3]-grdhead.bounds[2])/(grdhead.rows-1);
    cout << "*****************************************************\n";
    cout << "       Downward continuation from plane to plane(frequency)     \n";
    cout <<"        Data: "<<grdhead.rows<<" X "<<grdhead.cols<<endl;
    cout <<"        Size of kernel matrix: "<<grdhead.rows*grdhead.cols<<" X "<<grdhead.rows*grdhead.cols<<endl;
    cout <<"        Memory: "<<sizeof(double)*pow(grdhead.rows*grdhead.cols,1)/pow(1024.0,2)<< "Mb"<<endl;
    cout <<"        Downward continuation  ";
    cout <<"\033[31m";       //red
    cout <<(rph/dx);
    cout <<"\033[0m";
    cout <<"  point spacing\n";
    cout << "*****************************************************\n";
    
    
    double* G = new double[modelnum];	

    //get the first row of the matrix
    GetKernalMatrix_new(grdhead, G, rph);
    // outdata
    double *dataout=new double[modelnum];
   
    double* gk = new double[modelnum], *dk = new double[modelnum], *yk1 = new double[modelnum], *sk1 = new double[modelnum];//yk1表示yk-1
    double* Adk = new double[modelnum], *gk1 = new double[modelnum], *mk1 = new double[modelnum];//gk1表示gk-1
    double* B = new double[modelnum], *ATAdk = new double[modelnum];
    double thetak, betak1, aerfak;
    UWC_p2p_f(indata,B,grdhead,rph);

    for (int i = 0; i<modelnum; i++)
    {
        mk1[i] = 0;								 
        gk1[i] = -B[i];						
        dk[i] = -gk1[i];			
    }
    //====================calculate Adk=============
    UWC_p2p_f(dk, Adk, grdhead,rph);
    UWC_p2p_f(Adk, ATAdk,grdhead,rph);
    for (int i = 0; i < modelnum; i++)
    {
        ATAdk[i] += TRP*dk[i];													
    }
    //=========================================
    double fenzi = 0, fenmu = 0;
    for (int i = 0; i<modelnum; i++)
    {
        fenzi += gk1[i] * dk[i];
        fenmu += dk[i] * ATAdk[i];
    }
    aerfak = -fenzi / fenmu;

    for (int i = 0; i<modelnum; i++)
    {
        dataout[i] = mk1[i] + aerfak*dk[i];
        gk[i] = gk1[i] + aerfak*ATAdk[i];
    }
    //1.5 loop and iteration
    double* error_vector = new double[kmax];
    for (int k = 0; k<kmax; k++)
    {
        //1.1 yk-1,sk-1
        for (int i = 0; i<modelnum; i++)
        {
            yk1[i] = gk[i] - gk1[i];
            sk1[i] = dataout[i] - mk1[i];
        }
        //1.2 thetak
        double fenzi = 0, fenmu = 0, fenzi2 = 0, fenmu2 = 0;
        for (int i = 0; i<modelnum; i++)
        {
            fenzi += pow(yk1[i], 2.0);
            fenzi2 += gk[i] * sk1[i];
            fenmu += gk[i] * yk1[i];
            fenmu2 += sk1[i] * yk1[i];
        }
        thetak = 1.0 - fenzi*fenzi2 / (fenmu2*fenmu);
        if (thetak <= 0.25)
        {
            thetak = 1;
        }
        //1.3 betak-1
        double fenzi3 = 0;
        fenzi = 0, fenmu = 0; fenzi2 = 0; fenmu2 = 0;
        for (int i = 0; i<modelnum; i++)
        {
            fenzi += gk[i] * yk1[i];
            fenmu += sk1[i] * yk1[i];
            
            fenzi2 += pow(yk1[i], 2.0);
            fenzi3 += gk[i] * sk1[i];
            fenmu2 += sk1[i] * yk1[i];
        }
        betak1 = fenzi / fenmu - fenzi2*fenzi3 / pow(fenmu2, 2.0);//printf("betak-1 %lf\t",betak1);
        //1.4 dk
        for (int i = 0; i<modelnum; i++)
        {
            dk[i] = -thetak*gk[i] + betak1*sk1[i];
        }

        UWC_p2p_f(dk, Adk, grdhead,rph);
        UWC_p2p_f(Adk, ATAdk,grdhead,rph);
        for (int i = 0; i < modelnum; i++)
        {
            ATAdk[i] += TRP*dk[i];
        }
        //==============================================================
        
        fenzi = 0; fenmu = 0;
        for (int i = 0; i<modelnum; i++)
        {
            fenzi += gk[i] * dk[i];
            fenmu += dk[i] * ATAdk[i];
        }
        aerfak = -fenzi / fenmu;
        //1.6 mk+1
        for (int i = 0; i<modelnum; i++)											//把mk给mk-1,gk给gk-1
        {
            mk1[i] = dataout[i];
            gk1[i] = gk[i];
        }
        for (int i = 0; i<modelnum; i++)
        {
            dataout[i] = mk1[i] + aerfak*dk[i];
        }
        //1.7 gk
        double error_MHS = 0;
        for (int i = 0; i<modelnum; i++)
        {
            gk[i] = gk1[i] + aerfak*ATAdk[i];
            error_MHS += pow(gk[i], 2.0);
        }
        printf("The %dth iteration  ,error:", k);
        cout<<"\033[32m"; 
        printf("%lf\t", error_MHS);
        cout<<"\033[0m";
        printf("%lf\t%lf\n", aerfak, thetak);
        error_vector[k] = error_MHS;
        if (error_MHS<daierta)
        {
            break;
        }
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
        SaveGrd2VTK(filename_outfile+"_origin.vtk",grdhead,indata,height1);//save the origin data as well
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
    //delete array
    delete[] gk;
    delete[] dk;
    delete[] yk1;
    delete[] sk1;
    delete[] Adk;
    delete[] gk1;
    delete[] mk1;
    delete[] B;
    delete[] ATAdk;
    delete[] indata;
    delete[] dataout;
    delete[] G;
 }
 /*===============================================================================
 Downward continuation: uneven surface to plane(space domain)
 ================================================================================*/
void DWC_s2p(string inputfilename,string outputfilename,string topoFile,
    double height2,int extNum,double DWC_parameter,int DWC_method, int num_thread, bool isProgress,bool useOldKernel,
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
    if(useOldKernel)
    {
        cout<<"Old kernel p2s is depressed"<<endl;
        exit(0);
    }else
    {
        Getkernel_p2s_new(grdhead,height2,topo,G,num_thread);
    }
    //4. outdata
    double *dataout=new double[modelnum];
    //call downward continuation subfunction
    switch(DWC_method)
    {
        case DWC_CGLS:
            DWC_s2p_CGLS(G,dataout,indata,grdhead,extNum,DWC_parameter,num_thread);
            break;
        case DWC_INTEGRALITERATION:
            DWC_s2p_ItegrationIter(G,dataout,indata,grdhead,extNum,num_thread,outputfilename,DWC_parameter,exactsolution);
            break;
        case DWC_LANDWEBER:
            DWC_s2p_LandweberIter(G,dataout,indata,grdhead,extNum,num_thread,outputfilename,DWC_parameter,exactsolution);
            break;
        case DWC_TIKHONOV:
            cout<<RED<<"The Tikhonov method is not implemented yet "<<COLOR_DEFALUT<<endl;
            // DWC_s2p_Tikhonov(G,dataout,indata,grdhead,extNum,DWC_parameter,num_thread);
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

void DWC_s2p_ItegrationIter(double** G,double* x, double* b,
    GrdHead grdhead,int extNum,int num_thread,string outputfile,double iter_number,
    double* ExactSolution)
{
    int kmax=100;
    if(iter_number>0)kmax=int(iter_number);
    cout<<"Integral iteration method: kmax="<<kmax<<endl;
    double REmin=1E-8;
    int modelnum=grdhead.cols*grdhead.rows;
    for(int i=0;i<modelnum;i++)x[i]=b[i];
    double* residual=new double[modelnum];
    for(int i=0;i<modelnum;i++)residual[i]=0;
    double* x_uwc=new double[modelnum];
    double* tempArray=new double[modelnum];

    string path_tempResult=Path_GetFileName(outputfile)+"_Integral";
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

    string filename_log=path_tempResult+"/IntegralInteration_log.txt";
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
        for(int i=0;i<modelnum;i++)x[i]=x[i]+residual[i];
        UWC_Gij(x_uwc,G,x,modelnum,num_thread);
        for(int i=0;i<modelnum;i++)residual[i]=b[i]-x_uwc[i];
        
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
        
        //5. write temporary result,save as vtk
        // if (!SaveGrd(path_tempResult+"/"+std::to_string(k)+".grd", grdhead, x,extNum,false,false))return ;
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
//CGLS
//P. C. Hansen, Deconvolution and regularization with Toeplitz matrices, 2002, 323-378
//iteration form: 
void DWC_p2p_CGLS(double* G_firstRow,double* x, double* b,GrdHead grdhead,int extNum,
double delta,int num_thread)
{
    
    int kmax=500;
    if(delta<=0)delta=0.001;
    cout<<"CGLS method: kmax="<<kmax<<" delta="<<delta<<endl;
    //////////////////////////////////////
    int modelnum=grdhead.cols*grdhead.rows;
    double* rk1=new double[modelnum];//k-1
    double* dk1=new double[modelnum];
    double* xk1=new double[modelnum];
    
    double* rk=new double[modelnum];//k
    double* dk=new double[modelnum]; 
    // double* xk=new double[modelnum];xk=x
    double* temp_Ab=new double[modelnum];//temprary variable to store: Matrix * vector
    double alphak1=0,alphak2=0,alphak=0;
    double betak1=0,betak2=0,betak=0;
    //1. first guess of x
    for(int i=0;i<modelnum;i++)
    {
        xk1[i]=0;
        rk1[i]=b[i];//r0=b-A*x0, but x0=0,so r0=b
    }
    UWC_Gij(dk1,G_firstRow,rk1,grdhead,num_thread);//d0=AT * r0
    //iteration
    vector<double>bar_left;bar_left.push_back(0);bar_left.push_back(log10(1E6));
    vector<double>bar_right;bar_right.push_back(double(kmax));bar_right.push_back(log10(delta));
    vector<double>bar_pos; bar_pos.push_back(bar_left[0]);bar_pos.push_back(bar_left[1]);
    vector<string>bar_title; bar_title.push_back("Iteration");bar_title.push_back("Convergence:log");
    MultiProgressBar multiBar(bar_left,bar_right,bar_title);

    for(int k=0;k<kmax;k++)
    {
        //compute alpha_k
        UWC_Gji(temp_Ab,G_firstRow,rk1,grdhead,num_thread);//AT * r_k-1
        alphak1=Norm2(temp_Ab,modelnum);
        UWC_Gij(temp_Ab,G_firstRow,dk1,grdhead,num_thread);//A * d_k-1
        alphak2=Norm2(temp_Ab,modelnum);
        alphak=alphak1/alphak2;
        alphak=alphak*alphak;//square

        //compute x_k
        for(int i=0;i<modelnum;i++)x[i]=xk1[i]+alphak*dk1[i]; //x_k = x_k-1 + alpha_k * d_k-1
        
        //compute r_k
        UWC_Gij(temp_Ab,G_firstRow,dk1,grdhead,num_thread);//A * d_k-1
        for(int i=0;i<modelnum;i++)rk[i]=rk1[i]-alphak*temp_Ab[i];//r_k = r_k-1 - alpha_k * (A*d_k-1)

        //compute beta_k
        UWC_Gji(temp_Ab,G_firstRow,rk,grdhead,num_thread);
        betak1=Norm2(temp_Ab,modelnum);
        UWC_Gji(temp_Ab,G_firstRow,rk1,grdhead,num_thread);
        betak2=Norm2(temp_Ab,modelnum);
        betak=betak1/betak2;
        betak=betak*betak;

        //compute d_k
        UWC_Gji(temp_Ab,G_firstRow,rk,grdhead,num_thread);//AT * r_k
        for(int i=0;i<modelnum;i++)dk[i]=temp_Ab[i]+betak*dk1[i];
        double residual=Norm2(temp_Ab,modelnum);
        //update xk1,dk1,rk1
        for(int i=0;i<modelnum;i++)
        {
            xk1[i]=x[i];
            dk1[i]=dk[i];
            rk1[i]=rk[i];
        }
        //5. write temporary result
        // if (!SaveGrd("dwc_temp_"+std::to_string(k)+".grd", grdhead, x,extNum,false,false))return ;
        SaveGrd2VTK("dwc_temp_"+std::to_string(k)+".vtk",grdhead,x);
        //update progressbar
        bar_pos[0]=k+1;
        bar_pos[1]=log10(residual);
        multiBar.Update(bar_pos);
        //check convergence
        if(residual<delta || k>=kmax)break;
    }
    cout<<endl;
    
    //delete
    delete[] rk1,dk1,xk1,rk,dk,temp_Ab;
}
//CGLS:
//P. C. Hansen, Deconvolution and regularization with Toeplitz matrices, 2002, 323-378
//iteration form: 
void DWC_s2p_CGLS(double** G,double* x, double* b,GrdHead grdhead,int extNum,double delta,int num_thread)
{
    if(delta<=0)delta=0.001;
    int kmax=500;
    cout<<"CGLS method: kmax="<<kmax<<" delta:"<<delta<<endl;
    //////////////////////////////////////
    int modelnum=grdhead.cols*grdhead.rows;
    double* rk1=new double[modelnum];//k-1
    double* dk1=new double[modelnum];
    double* xk1=new double[modelnum];
    
    double* rk=new double[modelnum];//k
    double* dk=new double[modelnum]; 
    // double* xk=new double[modelnum];xk=x
    double* temp_Ab=new double[modelnum];//temprary variable to store: Matrix * vector
    double alphak1=0,alphak2=0,alphak=0;
    double betak1=0,betak2=0,betak=0;
    //1. first guess of x
    for(int i=0;i<modelnum;i++)
    {
        xk1[i]=0;
        rk1[i]=b[i];//r0=b-A*x0, but x0=0,so r0=b
    }
    UWC_Gij(dk1,G,rk1,modelnum,num_thread);//d0=AT * r0
    //iteration
    vector<double>bar_left;bar_left.push_back(0);bar_left.push_back(log10(1E6));
    vector<double>bar_right;bar_right.push_back(double(kmax));bar_right.push_back(log10(delta));
    vector<double>bar_pos; bar_pos.push_back(bar_left[0]);bar_pos.push_back(bar_left[1]);
    vector<string>bar_title; bar_title.push_back("Iteration");bar_title.push_back("Convergence:log");
    MultiProgressBar multiBar(bar_left,bar_right,bar_title);

    for(int k=0;k<kmax;k++)
    {
        //compute alpha_k
        UWC_Gji(temp_Ab,G,rk1,modelnum,num_thread);//AT * r_k-1
        alphak1=Norm2(temp_Ab,modelnum);
        UWC_Gij(temp_Ab,G,dk1,modelnum,num_thread);//A * d_k-1
        alphak2=Norm2(temp_Ab,modelnum);
        alphak=alphak1/alphak2;
        alphak=alphak*alphak;//square

        //compute x_k
        for(int i=0;i<modelnum;i++)x[i]=xk1[i]+alphak*dk1[i]; //x_k = x_k-1 + alpha_k * d_k-1
        
        //compute r_k
        UWC_Gij(temp_Ab,G,dk1,modelnum,num_thread);//A * d_k-1
        for(int i=0;i<modelnum;i++)rk[i]=rk1[i]-alphak*temp_Ab[i];//r_k = r_k-1 - alpha_k * (A*d_k-1)

        //compute beta_k
        UWC_Gji(temp_Ab,G,rk,modelnum,num_thread);
        betak1=Norm2(temp_Ab,modelnum);
        UWC_Gji(temp_Ab,G,rk1,modelnum,num_thread);
        betak2=Norm2(temp_Ab,modelnum);
        betak=betak1/betak2;
        betak=betak*betak;

        //compute d_k
        UWC_Gji(temp_Ab,G,rk,modelnum,num_thread);//AT * r_k
        for(int i=0;i<modelnum;i++)dk[i]=temp_Ab[i]+betak*dk1[i];
        double residual=Norm2(temp_Ab,modelnum);
        //update xk1,dk1,rk1
        for(int i=0;i<modelnum;i++)
        {
            xk1[i]=x[i];
            dk1[i]=dk[i];
            rk1[i]=rk[i];
        }
        //5. write temporary result
        // if (!SaveGrd("dwc_temp_"+std::to_string(k)+".grd", grdhead, x,extNum,false,false))return ;
        SaveGrd2VTK("dwc_temp_"+std::to_string(k)+".vtk",grdhead,x);
        //update progressbar
        bar_pos[0]=k+1;
        bar_pos[1]=log10(Norm2(rk,modelnum));
        multiBar.Update(bar_pos);
        //check convergence
        if(residual<delta || k>=kmax)break;
    }
    cout<<endl;
    
    //delete
    delete[] rk1,dk1,xk1,rk,dk,temp_Ab;
}
void UWC_G_CGLS_Tik(double* b,double** G,double* x, int modelnum,double lambda2,int num_thread)
{
    double* temp=new double[modelnum];
    UWC_Gij(temp,G,x,modelnum,num_thread);//G*x
    UWC_Gji(b,G,temp,modelnum,num_thread);//b=GT*G*x
    for(int i=0;i<modelnum;i++)b[i]=b[i]+lambda2*x[i];
    //delete
    delete[] temp;
}
//CGLS+Tikhonov regularization
//(AT*A+lambda^2*I)*x = AT*b
void DWC_s2p_Tikhonov(double** G,double* x, double* b0,GrdHead grdhead,int extNum,double lambda,int num_thread)
{
    if(lambda<=0)lambda=0.001;
    cout<<"Tikhonov+CG method: "<<lambda<<endl;
    int kmax=500;
    double delta=0.01;
    double lambda2=lambda*lambda;
    //////////////////////////////////////
    int modelnum=grdhead.cols*grdhead.rows;
    double* b=new double[modelnum];
    UWC_Gji(b,G,b0,modelnum,num_thread);//b=GT*b0
    double* rk1=new double[modelnum];//k-1
    double* dk1=new double[modelnum];
    double* xk1=new double[modelnum];
    
    double* rk=new double[modelnum];//k
    double* dk=new double[modelnum]; 
    // double* xk=new double[modelnum];xk=x
    double* temp_Ab=new double[modelnum];//temprary variable to store: Matrix * vector
    double alphak1=0,alphak2=0,alphak=0;
    double betak1=0,betak2=0,betak=0;
    //1. first guess of x
    for(int i=0;i<modelnum;i++)
    {
        xk1[i]=0;
        rk1[i]=b[i];//r0=b-A*x0, but x0=0,so r0=b
    }
    UWC_G_CGLS_Tik(dk1,G,rk1,modelnum,lambda2,num_thread);
    //iteration
    vector<double>bar_left;bar_left.push_back(0);bar_left.push_back(log10(1E6));
    vector<double>bar_right;bar_right.push_back(double(kmax));bar_right.push_back(log10(delta));
    vector<double>bar_pos; bar_pos.push_back(bar_left[0]);bar_pos.push_back(bar_left[1]);
    vector<string>bar_title; bar_title.push_back("Iteration");bar_title.push_back("Convergence:log");
    MultiProgressBar multiBar(bar_left,bar_right,bar_title);

    for(int k=0;k<kmax;k++)
    {
        //compute alpha_k
        UWC_G_CGLS_Tik(temp_Ab,G,rk1,modelnum,lambda2,num_thread);//(GT*G + lambda*I) * r_k-1
        alphak1=Norm2(temp_Ab,modelnum);
        UWC_G_CGLS_Tik(temp_Ab,G,dk1,modelnum,lambda2,num_thread);//(GT*G + lambda*I) * d_k-1
        alphak2=Norm2(temp_Ab,modelnum);
        alphak=alphak1/alphak2;
        alphak=alphak*alphak;//square

        //compute x_k
        for(int i=0;i<modelnum;i++)x[i]=xk1[i]+alphak*dk1[i]; //x_k = x_k-1 + alpha_k * d_k-1
        
        //compute r_k
        UWC_G_CGLS_Tik(temp_Ab,G,dk1,modelnum,lambda2,num_thread);//(GT*G + lambda*I) * d_k-1
        for(int i=0;i<modelnum;i++)rk[i]=rk1[i]-alphak*temp_Ab[i];//r_k = r_k-1 - alpha_k * (A*d_k-1)

        //compute beta_k
        UWC_G_CGLS_Tik(temp_Ab,G,rk,modelnum,lambda2,num_thread);
        betak1=Norm2(temp_Ab,modelnum);
        UWC_G_CGLS_Tik(temp_Ab,G,rk1,modelnum,lambda2,num_thread);
        betak2=Norm2(temp_Ab,modelnum);
        betak=betak1/betak2;
        betak=betak*betak;

        //compute d_k
        UWC_G_CGLS_Tik(temp_Ab,G,rk,modelnum,lambda2,num_thread);//(GT*G + lambda*I) * r_k
        for(int i=0;i<modelnum;i++)dk[i]=temp_Ab[i]+betak*dk1[i];
        double residual=Norm2(temp_Ab,modelnum);
        //update xk1,dk1,rk1
        for(int i=0;i<modelnum;i++)
        {
            xk1[i]=x[i];
            dk1[i]=dk[i];
            rk1[i]=rk[i];
        }
        //5. write temporary result
        // if (!SaveGrd("dwc_temp_"+std::to_string(k)+".grd", grdhead, x,extNum,false,false))return ;
        SaveGrd2VTK("dwc_temp_"+std::to_string(k)+".vtk",grdhead,x);
        //update progressbar
        bar_pos[0]=k+1;
        bar_pos[1]=log10(Norm2(rk,modelnum));
        multiBar.Update(bar_pos);
        //check convergence
        if(residual<delta || k>=kmax)break;
    }
    cout<<endl;
    
    //delete
    delete[] rk1,dk1,xk1,rk,dk,temp_Ab,b;
}
 /*===============================================================================
 Downward continuation: plane to plane(space domain)
 ================================================================================*/
 void DWC_p2p(string inputfilename,string outputfilename,
 double height1,double height2,int extNum,double DWC_parameter,int DWC_method,
 int num_thread, bool isProgress,bool useOldKernel,string filename_exact)
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
        case DWC_CGLS:
            // DWC_p2p_CGLS(G_firstRow,dataout,indata,grdhead,extNum,DWC_parameter,num_thread);
            DWC_s2p_CGLS(G,dataout,indata,grdhead,extNum,DWC_parameter,num_thread);
            
            break;
        case DWC_INTEGRALITERATION:
            // DWC_p2p_ItegrationIter(G_firstRow,dataout,indata,grdhead,extNum,num_thread,outputfilename,DWC_parameter,exactsolution);
            DWC_s2p_ItegrationIter(G,dataout,indata,grdhead,extNum,num_thread,outputfilename,DWC_parameter,exactsolution);
            break;
        case DWC_LANDWEBER:
            // DWC_p2p_LandweberIter(G_firstRow,dataout,indata,grdhead,extNum,num_thread,outputfilename,DWC_parameter,exactsolution);
            DWC_s2p_LandweberIter(G,dataout,indata,grdhead,extNum,num_thread,outputfilename,DWC_parameter,exactsolution);
            break;
        case DWC_TIKHONOV:
            // DWC_s2p_Tikhonov(G,dataout,indata,grdhead,DWC_parameter,num_thread);
            cout<<RED<<"The Tikhonov method of downward continuation is not implemented yet"<<endl;
            exit(0);
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
    // delete[] G_firstRow;
    delete[] indata;
    delete[] dataout;
    if(!exactsolution)delete[] exactsolution;
 }
void DWC_p2p_ItegrationIter(double* G,double* x, double* b,
    GrdHead grdhead,int extNum,int num_thread,string outputfile,
    double iter_number,double* ExactSolution)
{
    int kmax=100;
    if(iter_number>0)kmax=int(iter_number);
    cout<<"Integral iteration method: kmax="<<kmax<<endl;
    double REmin=1E-8;
    int modelnum=grdhead.cols*grdhead.rows;
    for(int i=0;i<modelnum;i++)x[i]=b[i];
    double* residual=new double[modelnum];
    for(int i=0;i<modelnum;i++)residual[i]=0;
    double* x_uwc=new double[modelnum];
    double* tempArray=new double[modelnum];

    string path_tempResult=Path_GetFileName(outputfile)+"_Integral";
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
    string filename_log=path_tempResult+"/IntegralInteration_log.txt";
    ofstream fout_log(filename_log);
    if(!fout_log)
    {
        cout<<RED<<"[DWC_s2p_ItegrationIter]Open file failed: "<<filename_log<<COLOR_DEFALUT<<endl;
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
        for(int i=0;i<modelnum;i++)x[i]=x[i]+residual[i];
        UWC_Gij(x_uwc,G,x,grdhead,num_thread);
        for(int i=0;i<modelnum;i++)residual[i]=b[i]-x_uwc[i];
        
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
        // if (!SaveGrd(path_tempResult+"/"+std::to_string(k)+".grd", grdhead, x,extNum,false,false))return ;
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
void DWC_p2p_LandweberIter(double* G,double* x, double* b,
    GrdHead grdhead,int extNum,int num_thread,string outputfile,
    double iter_number,double* ExactSolution)
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
        cout<<RED<<"[DWC_s2p_ItegrationIter]Open file failed: "<<filename_log<<COLOR_DEFALUT<<endl;
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
        // if (!SaveGrd(path_tempResult+"/"+std::to_string(k)+".grd", grdhead, x,extNum,false,false))return ;
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
        
        //5. write temporary result
        // if (!SaveGrd(path_tempResult+"/"+std::to_string(k)+".grd", grdhead, x,extNum,false,false))return ;
        SaveGrd2VTK(path_tempResult+"/result_"+std::to_string(k)+".vtk",grdhead,x);
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
 void DWC_Tikhonov_old(double* G_firstRow,double* dataout,double* indata,double TRP,int kmax,double daierta,GrdHead grdhead,int num_thread)
 {
     int modelnum=grdhead.rows*grdhead.cols;

    double* gk = new double[modelnum], *dk = new double[modelnum], *yk1 = new double[modelnum], *sk1 = new double[modelnum];//yk1表示yk-1
    double* Adk = new double[modelnum], *gk1 = new double[modelnum], *mk1 = new double[modelnum];//gk1表示gk-1
    double* B = new double[modelnum], *ATAdk = new double[modelnum];
    double thetak, betak1, aerfak;

    UWC_Gij(B, G_firstRow, indata,grdhead,num_thread);

    for (int i = 0; i<modelnum; i++)
    {
        mk1[i] = 0;								 
        gk1[i] = -B[i];						
        dk[i] = -gk1[i];							
    }
    UWC_Gij(Adk, G_firstRow,dk, grdhead, num_thread);
    UWC_Gij(ATAdk,G_firstRow,Adk, grdhead,num_thread);

    for (int i = 0; i < modelnum; i++)
    {
        ATAdk[i] += TRP*dk[i];													
    }
    //==========================================
    double fenzi = 0, fenmu = 0;
    for (int i = 0; i<modelnum; i++)
    {
        fenzi += gk1[i] * dk[i];
        fenmu += dk[i] * ATAdk[i];
    }
    aerfak = -fenzi / fenmu;

    for (int i = 0; i<modelnum; i++)
    {
        dataout[i] = mk1[i] + aerfak*dk[i];
        gk[i] = gk1[i] + aerfak*ATAdk[i];
    }

    double* error_vector = new double[kmax];
    for (int k = 0; k<kmax; k++)
    {
        for (int i = 0; i<modelnum; i++)
        {
            yk1[i] = gk[i] - gk1[i];
            sk1[i] = dataout[i] - mk1[i];
        }
        double fenzi = 0, fenmu = 0, fenzi2 = 0, fenmu2 = 0;
        for (int i = 0; i<modelnum; i++)
        {
            fenzi += pow(yk1[i], 2.0);
            fenzi2 += gk[i] * sk1[i];
            fenmu += gk[i] * yk1[i];
            fenmu2 += sk1[i] * yk1[i];
        }
        thetak = 1.0 - fenzi*fenzi2 / (fenmu2*fenmu);
        if (thetak <= 0.25)
        {
            thetak = 1;
        }
        double fenzi3 = 0;
        fenzi = 0, fenmu = 0; fenzi2 = 0; fenmu2 = 0;
        for (int i = 0; i<modelnum; i++)
        {
            fenzi += gk[i] * yk1[i];
            fenmu += sk1[i] * yk1[i];
            
            fenzi2 += pow(yk1[i], 2.0);
            fenzi3 += gk[i] * sk1[i];
            fenmu2 += sk1[i] * yk1[i];
        }
        betak1 = fenzi / fenmu - fenzi2*fenzi3 / pow(fenmu2, 2.0);//printf("betak-1 %lf\t",betak1);
        for (int i = 0; i<modelnum; i++)
        {
            dk[i] = -thetak*gk[i] + betak1*sk1[i];
        }

        UWC_Gij(Adk, G_firstRow,dk,  grdhead,num_thread);
        UWC_Gij(ATAdk,G_firstRow,Adk,  grdhead,num_thread);

        for (int i = 0; i < modelnum; i++)
        {
            ATAdk[i] += TRP*dk[i];
        }
        
        fenzi = 0; fenmu = 0;
        for (int i = 0; i<modelnum; i++)
        {
            fenzi += gk[i] * dk[i];
            fenmu += dk[i] * ATAdk[i];
        }
        aerfak = -fenzi / fenmu;

        for (int i = 0; i<modelnum; i++)		
        {
            mk1[i] = dataout[i];
            gk1[i] = gk[i];
        }
        for (int i = 0; i<modelnum; i++)
        {
            dataout[i] = mk1[i] + aerfak*dk[i];
        }

        double error_MHS = 0;
        for (int i = 0; i<modelnum; i++)
        {
            gk[i] = gk1[i] + aerfak*ATAdk[i];
            error_MHS += pow(gk[i], 2.0);
        }
        printf("The %dth iteration  ,error:", k);
        cout<<"\033[32m";    
        printf("%lf\t", error_MHS);
        cout<<"\033[0m";
        printf("%lf\t%lf\n", aerfak, thetak);
        error_vector[k] = error_MHS;
        if (error_MHS<daierta)
        {
            break;
        }
    }
    
    //delete array
    delete[] gk;
    delete[] dk;
    delete[] yk1;
    delete[] sk1;
    delete[] Adk;
    delete[] gk1;
    delete[] mk1;
    delete[] B;
    delete[] ATAdk;

 }

 double Norm2(double* x,const int num)
 {
     double sum=0;
     for(int i=0;i<num;i++)
     {
         sum=sum+x[i]*x[i];
     }
     return sqrt(sum);
 }
 double Norm2_Gradient(double* result,GrdHead grdhead)
 {
    double fx=0,fy=0,sum=0;
    //compute fx
    int index_ij=0;
    for(int i=0;i<grdhead.rows-1;i++)
    {
        for(int j=0;j<grdhead.cols-1;j++)
        {
            index_ij=j+i*grdhead.cols;
            fx=result[index_ij+1]-result[index_ij];
            fy=result[index_ij+grdhead.cols]-result[index_ij];
            sum=fx*fx+fy*fy;
        }
    }

    return sqrt(sum);
 }

 int SaveGrd2VTK(string outputfile,GrdHead grdhead,double* data,double z)
{
    //vtk format
    ofstream fout(outputfile);
    if(!fout)
    {
        cout<<endl;
        cout<<"["<<RED<<"error"<<COLOR_DEFALUT<<"]open file failed: "<<outputfile+".vtk"<<endl;
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

int SaveGrd2VTK_topo(string outputfile,GrdHead grdhead,double* data,double* topo)
{
    //vtk format
    ofstream fout(outputfile);
    if(!fout)
    {
        cout<<endl;
        cout<<"["<<RED<<"error"<<COLOR_DEFALUT<<"]open file failed: "<<outputfile+".vtk"<<endl;
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
                <<topo[iy*grdhead.cols+ix]<<" ";
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