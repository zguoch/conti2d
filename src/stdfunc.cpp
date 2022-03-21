/**
 * @file stdfunc.cpp
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief Implementation of some basic functions. For example read and write file.
 * @version 1.0
 * @date 2019-09-03
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#include "stdfunc.h"

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

 /**
  * @brief Save Grd data as netCDF format. see https://www.unidata.ucar.edu/software/netcdf/ for details.
  * 
  * @param filename 
  * @param grdhead 
  * @param data 
  * @param extNum 
  * @param isInfo 
  * @return true 
  * @return false 
  */
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
    delete[] x;
    delete[] y;
    for(int i=0;i<NX;i++)delete[] data_out[i];
    delete[] data_out;
    return true;
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

 int SaveGrd2VTK(string outputfile,GrdHead grdhead,double* data,double z,double* err)
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
    // write residual if exist
    if(err)
    {
        fout<<"SCALARS Residual_UWC float"<<endl;
        fout<<"LOOKUP_TABLE default"<<endl;
        for(int i=0;i<num_pt;i++)
        {
            fout<<err[i]<<" ";
        }fout<<endl;
    }
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