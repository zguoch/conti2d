/**
 * @file polygonmesh.cpp
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief 
 * @version 1.0
 * @date 2019-08-27
 * 
 * @copyright Copyright (c) 2019 by Zhikui Guo and Yi Zhang
 * 
 */
#include "polygonmesh.h"

/**
 * @brief Construct a new c Poly Mesh::c Poly Mesh object
 * 
 */
cPolyMesh::cPolyMesh()
{
    nNodes=0;        
    nCells=0;        
    x=NULL, y=NULL,m_cData=NULL;            
    DataName="";        
}
cPolyMesh::~cPolyMesh()
{
    if(x)delete[] x;
    if(y)delete[] y;
    if(m_cData)delete[] m_cData;
}
/**
 * @brief Write polygonmesh to gmsh's .geo file to visualize the geometry. 
 * Because the edge number of the polygonmesh is not fixed,
 * and gmsh doesn't support arbritrary polygon type. Therefore just write point and line of each edge.
 * 
 * @param outfile 
 */
void cPolyMesh::WriteGeometry2Gmsh(string outfile)
{
    ofstream fout(outfile);
    if(!fout)
    {
        cout<<"Open file failed: "<<outfile<<endl;
        exit(0);
    }
    // write points
    for(int i=0;i<nNodes;i++)
    {
        fout<<"Point("<<i<<")={"
        <<x[i]<<","
        <<y[i]<<",0,1};"<<endl;
    }
    // write edge line
    for(int i=0;i<el2nod.size()-1;i++)
    {
        for(int j=0;j<el2nod[i].size()-1;j++)
        {
            fout<<"Line(newl)={"
            <<el2nod[i][j]<<","
            <<el2nod[i][j+1]<<"};"<<endl;
        }
        // make it close
        fout<<"Line(newl)={"
            <<el2nod[i][el2nod[i].size()-1]<<","
            <<el2nod[i][0]<<"};"<<endl;
    }
    fout.close();
}
