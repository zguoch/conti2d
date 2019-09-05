/**
 * @file trimesh.h defines triangle mesh
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief 
 * @version 1.0
 * @date 2019-08-27
 * 
 * @copyright Copyright (c) 2019 by Zhikui Guo and Yi Zhang
 * 
 */
#ifndef TRIMESH_H
#define TRIMESH_H
#include "polygonmesh.h"
#include "stdfunc.h"
#include <fstream>
#include <vector>
using namespace std;

struct NODES
{
    double x,y,z; //**< Coordinages of a node*/
    bool isBoundary; //* Whether a node is on the boundary*/
};

struct TRIANGLE
{
    int iNodes[3]; /**< index of three vertex*/
    int iNeighbour[3]; /**< index of three neighbour, -1 means boundary*/
    bool isBoundary; /**< whethere this triangle is a boundary element*/
};
/**
 * @brief cTriMesh defines mesh structure, coordinate and field data for observation points, field and continuation field. 
 * 
 * \code{.cpp}
  The coordinate system is described as following, 
        /y
       /
      /
    0/__________-x
     |
     |
     |
     |z(-)
  \endcode
 */
class cTriMesh
{
    public:
        cTriMesh();
        ~cTriMesh();
    public:
        int numberofpoints;
        int numberofpoints_bd; /**< boundary points number */
        // double* x;          /**< 1D array with size of cTriMesh.numberofpoints. Coordinate x of observation points */
        // double* y;          /**< 1D array with size of cTriMesh.numberofpoints. Coordinate y of observation points */
        // double* z;          /**< 1D array with size of cTriMesh.numberofpoints. Coordinate z of observation points */
        NODES* m_nodes;
        int numberofelements; /**< Number of elements of the mesh */
        int numberofelements_bd; /**< boundary elements number*/
        TRIANGLE* m_triElements; /**< triangle elements*/
        double* field;      /**< 1D array with size of cTriMesh.numberofpoints. Potential field data on the observation points (cTriMesh.x, cTriMesh.y, cTriMesh.z) */
        string fieldname;   /**< Name of the potential field, default is "Original Field" */
        static const int nv_el=3;          /**< Number of verticex of element of the triangle mesh */
        
        double* field_conti; /**< Continuation field array */
        /**
         * @brief Polygons which centered at each observation points. 
         * Therefore  the m_vertexCenteredPolygons.nNodes must be equal to 
         * cTriMesh.numberofelements, 
         * but edge number of each polygon could not be equal, 
         * for example an observation point is shared with 5 triangles, 
         * this point centered polygon has 5 edges.
         * 
         */
        cPolyMesh m_vertexCenteredPolygons; 
    public:
        /**
         * @brief Generates vertex centered polygon mesh according to coordinates of vertex and connection of each triangles.
         * 
         * @return int 
         */
        int GenerateVertexCenteredPolygon();
    private:
        void getCentroids(double* cx, double* cy);
        /**
         * @brief Get the Connection VertexCenteredPolygon object. 
         * Make sure the points have correct connection for a right hand polygon.
         * 
         * @param sharedTriangles 
         * @param x0 x coordinate of the vertex of a triangle
         * @param y0 y coordinate of the vertex of a triangle
         * @param xx x coordinate array of centroids of triangles
         * @param yy y coordinate array of centroids of triangles
         */
        void reConnect(vector<int>& sharedTriangles,double x0,double y0,const double* xx,const double* yy);
        /**
         * @brief Get the connection of centroids of shared triangles to construct a right hand polygon.
         * But for the vertex on boundary, the extension should be made to locate the vertex in the polygon.
         * 
         * @param sharedTriangles 
         * @param x0 
         * @param y0 
         * @param xx 
         * @param yy 
         */
        void reConnect_Boundary(vector<int>& sharedTriangles,double x0,double y0,const double* xx,const double* yy);
        void MirrorPoint2(double p1[2],double p2[2],double& newX,double& newY);
};
#endif