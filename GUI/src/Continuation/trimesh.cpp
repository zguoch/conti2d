/**
 * @file trimesh.cpp 
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief 
 * @version 1.0
 * @date 2019-08-27
 * 
 * @copyright Copyright (c) 2019 by Zhikui Guo and Yi Zhang
 * 
 */
#include "trimesh.h"

/**
 * @brief Construct a new c Tri Mesh::c Tri Mesh object
 * 
 */
cTriMesh::cTriMesh()
{
    m_nodes=NULL;
    field=NULL;
    field_conti=NULL;
    m_triElements=NULL;
    fieldname="Original Field";
    numberofelements=0;
    numberofpoints=0;
    numberofpoints_bd=0;
    numberofelements_bd=0;
}
cTriMesh::~cTriMesh()
{
    if(m_nodes)delete[] m_nodes;
    if(field)delete[] field;
    if(m_triElements)delete[] m_triElements;
}

int cTriMesh::GenerateVertexCenteredPolygon()
{
    m_vertexCenteredPolygons.nNodes=numberofelements+numberofpoints_bd;
    m_vertexCenteredPolygons.nCells=numberofpoints;
    // 1. 
    m_vertexCenteredPolygons.x=new double[m_vertexCenteredPolygons.nNodes];
    m_vertexCenteredPolygons.y=new double[m_vertexCenteredPolygons.nNodes];
    getCentroids(m_vertexCenteredPolygons.x,m_vertexCenteredPolygons.y);
    // 2. get index of elements wich share the same vertex
    
    // vector<vector<int> > nod2ele;
    m_vertexCenteredPolygons.el2nod.resize(numberofpoints);
    for(int i=0;i<numberofelements;i++)
    {
        for(int j=0;j<3;j++)
        {
            m_vertexCenteredPolygons.el2nod[m_triElements[i].iNodes[j]].push_back(i);
        }
    }
    // 3. reconnect the centroids related to a vertex and extend the boundary points
    int index_pt=0;
    vector<double>extendX,extendY;
    for(int i=0;i<numberofpoints;i++)
    {
        if(!m_nodes[i].isBoundary)
        {
            reConnect(m_vertexCenteredPolygons.el2nod[i],m_nodes[i].x,m_nodes[i].y,m_vertexCenteredPolygons.x,m_vertexCenteredPolygons.y);
        }else
        {
            reConnect_Boundary(m_vertexCenteredPolygons.el2nod[i],m_nodes[i].x,m_nodes[i].y,m_vertexCenteredPolygons.x,m_vertexCenteredPolygons.y);
            // extension of boundary points. 
            double p0[2]={m_nodes[i].x,m_nodes[i].y};
            //(1) Only two centroids, mirror the center of these two centroids to the vertex
            //(2) Three or more, mirror the central index of the centroids  
            if(m_vertexCenteredPolygons.el2nod[i].size()==2) //only two centroids
            {
                double p1[2]={(m_vertexCenteredPolygons.x[m_vertexCenteredPolygons.el2nod[i][0]]+m_vertexCenteredPolygons.x[m_vertexCenteredPolygons.el2nod[i][1]])/2.0, 
                              (m_vertexCenteredPolygons.y[m_vertexCenteredPolygons.el2nod[i][0]]+m_vertexCenteredPolygons.y[m_vertexCenteredPolygons.el2nod[i][0]])/2.0};
                double newx,newy;
                MirrorPoint2(p1,p0,newx,newy);
                extendX.push_back(newx);
                extendY.push_back(newy);
            }else
            {
                int index=ceil(m_vertexCenteredPolygons.el2nod[i].size()/2);
                double p1[2]={m_vertexCenteredPolygons.x[m_vertexCenteredPolygons.el2nod[i][index]], 
                              m_vertexCenteredPolygons.y[m_vertexCenteredPolygons.el2nod[i][index]]};
                double newx,newy;
                MirrorPoint2(p1,p0,newx,newy);
                extendX.push_back(newx);
                extendY.push_back(newy);
            }
            m_vertexCenteredPolygons.el2nod[i].push_back(numberofelements-1+extendX.size());
        }
    }
    // copy the extend point to the end of the centroids 
    for(int i=0;i<numberofpoints_bd;i++)
    {
        m_vertexCenteredPolygons.x[i+numberofelements]=extendX[i];
        m_vertexCenteredPolygons.y[i+numberofelements]=extendY[i];
    }
    
    return 0;
}
/**
 * @brief The same as cTriMesh.reConnect, just do some special process for boundary points. e.g., boundary points extension.
 * 
 * @param sharedTriangles 
 * @param x0 
 * @param y0 
 * @param xx 
 * @param yy 
 */
void cTriMesh::reConnect_Boundary(vector<int>& sharedTriangles,double x0,double y0,const double* xx,const double* yy)
{
    int tmp;
    double p0[2]={x0,y0},p1[2],p2[2];//**< \f$ \vec{p_0p_1} \times \vec{p_0p_2} \f$ must be positive to meet right hand rule */
    // 1. if only two centroids, just determin the connection direction accroding to right hand rule
    if(sharedTriangles.size()==2)
    {
        p1[0]=xx[sharedTriangles[0]]; p1[1]=yy[sharedTriangles[0]];
        p2[0]=xx[sharedTriangles[1]]; p2[1]=yy[sharedTriangles[1]];
        if(TwoCross(p0,p1,p2)<0)
        {
            tmp=sharedTriangles[0];
            sharedTriangles[0]=sharedTriangles[1];
            sharedTriangles[1]=tmp;
        }
        return;
    }
    // 2. greater than 2 centroids, need to sort the connection to a right hand polygon
    //the first on, must be one of the centroid of a boundary triangle
    //the number of boundary triangle must be 2, if not equal to 2, print information and the program need to update!!!!
    vector<int> index_bds;
    for(int i=0;i<sharedTriangles.size();i++)
    {
        if(m_triElements[sharedTriangles[i]].isBoundary==true)index_bds.push_back(i);
    }
    if(index_bds.size()!=2)
    {
        cout<<"Error report: for one vertex on boundary, the number of boundary triangle among shared-triangles are greater than 2, that's impossible!!"<<endl;
        cout<<"The vertex coordinate is ( "<<x0<<", "<<y0<<")"<<endl;
        cout<<"The shared triangle index are ";
        for(int i=0;i<sharedTriangles.size();i++)cout<<sharedTriangles[i]<<" ";
        cout<<endl;
        cout<<"The boundary triangle index are ";
        for(int i=0;i<index_bds.size();i++)cout<<index_bds[i]<<" ";
        cout<<endl;
    }
    p1[0]=xx[sharedTriangles[index_bds[0]]]; p1[1]=yy[sharedTriangles[index_bds[0]]];
    p2[0]=xx[sharedTriangles[index_bds[1]]]; p2[1]=yy[sharedTriangles[index_bds[1]]];
    if(TwoCross(p0,p1,p2)>0)
    {
        tmp=sharedTriangles[0];
        sharedTriangles[0]=sharedTriangles[index_bds[0]];
        sharedTriangles[index_bds[0]]=tmp;
    }else
    {
        tmp=sharedTriangles[0];
        sharedTriangles[0]=sharedTriangles[index_bds[1]];
        sharedTriangles[index_bds[1]]=tmp;
    }
    //the second one
    for(int i=0;i<3;i++)
    {
        for(int j=1;j<sharedTriangles.size();j++)
        {
            // make sure the connection meet the right hand rule, for calculation of equation * in the paper of Guo and Zhang(2019)
            if(m_triElements[sharedTriangles[0]].iNeighbour[i] == sharedTriangles[j])
            {
                p1[0]=xx[sharedTriangles[0]]; p1[1]=yy[sharedTriangles[0]];
                p2[0]=xx[sharedTriangles[j]]; p2[1]=yy[sharedTriangles[j]];
                if(TwoCross(p0,p1,p2)>0)
                {
                    tmp=sharedTriangles[1];
                    sharedTriangles[1]=sharedTriangles[j];
                    sharedTriangles[j]=tmp;
                }
            }
        }
    }
    // the remains
    for(int i=1;i<sharedTriangles.size()-1;i++)
    {
        for(int j=0;j<3;j++)
        {
            for(int jj=i+1;jj<sharedTriangles.size();jj++)
            {
                if( (m_triElements[sharedTriangles[i]].iNeighbour[j] == sharedTriangles[jj]) && 
                    (m_triElements[sharedTriangles[i]].iNeighbour[j] != sharedTriangles[i-1]) )
                    {
                        tmp=sharedTriangles[i+1];
                        sharedTriangles[i+1]=sharedTriangles[jj];
                        sharedTriangles[jj]=tmp;
                    }
            }
        }
    }
}
/**
 * @brief Sort the index of centroids relate to a vertex, make sure the connection of these centroids is a right hand polygon
 * 
 * @param sharedTriangles 
 * @param x0 
 * @param y0 
 * @param xx 
 * @param yy 
 */
void cTriMesh::reConnect(vector<int>& sharedTriangles,double x0,double y0,const double* xx,const double* yy)
{
    int tmp;
    double p0[2]={x0,y0},p1[2],p2[2];//**< \f$ \vec{p_0p_1} \times \vec{p_0p_2} \f$ must be positive to meet right hand rule */
    //the second one
    for(int i=0;i<3;i++)
    {
        for(int j=1;j<sharedTriangles.size();j++)
        {
            // make sure the connection meet the right hand rule, for calculation of equation * in the paper of Guo and Zhang(2019)
            if(m_triElements[sharedTriangles[0]].iNeighbour[i] == sharedTriangles[j])
            {
                p1[0]=xx[sharedTriangles[0]]; p1[1]=yy[sharedTriangles[0]];
                p2[0]=xx[sharedTriangles[j]]; p2[1]=yy[sharedTriangles[j]];
                if(TwoCross(p0,p1,p2)>0)
                {
                    tmp=sharedTriangles[1];
                    sharedTriangles[1]=sharedTriangles[j];
                    sharedTriangles[j]=tmp;
                }
                
            }
        }
    }
    // the remains
    for(int i=1;i<sharedTriangles.size()-1;i++)
    {
        for(int j=0;j<3;j++)
        {
            for(int jj=i+1;jj<sharedTriangles.size();jj++)
            {
                if( (m_triElements[sharedTriangles[i]].iNeighbour[j] == sharedTriangles[jj]) && 
                    (m_triElements[sharedTriangles[i]].iNeighbour[j] != sharedTriangles[i-1]) )
                    {
                        tmp=sharedTriangles[i+1];
                        sharedTriangles[i+1]=sharedTriangles[jj];
                        sharedTriangles[jj]=tmp;
                    }
            }
        }
    }
}
/**
 * @brief Get centroid of each triangles
 * 
 * @param cx 1D dynamic array with size of cTriMesh.numberofelements, 
 * 
 * @param cy 
 */
void cTriMesh::getCentroids(double* cx,double* cy)
{
    for (int i = 0; i < numberofelements; i++)
    {
        cx[i]=0; cy[i]=0;
        for(int j=0;j<3;j++)
        {
            cx[i]+=m_nodes[m_triElements[i].iNodes[j]].x;
            cy[i]+=m_nodes[m_triElements[i].iNodes[j]].y;
        }
        cx[i]=cx[i]/3;
        cy[i]=cy[i]/3;
    }
}

void cTriMesh::MirrorPoint2(double p1[2],double p2[2],double& newX,double& newY)
{
    newX=p2[0]+(p2[0]-p1[0]);
    newY=p2[1]+(p2[1]-p1[1]);
}