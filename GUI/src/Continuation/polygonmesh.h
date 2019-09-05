/**
 * @file polygonmesh.h defines class of the cPolyMesh
 * @author Zhikui Guo (zhikuiguo@live.cn)
 * @brief 
 * @version 1.0
 * @date 2019-08-27
 * 
 * @copyright Copyright (c) 2019 by Zhikui Guo and Yi Zhang
 * 
 */

#ifndef POLYGONMESH_H
#define POLYGONMESH_H
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;
/**
 * @brief Polygon mesh structure definition. Only defined on a horizontal plane.
 * 
 */
class cPolyMesh
{
private:
    /* data */
public:
    cPolyMesh(/* args */);
    ~cPolyMesh();
public:
    int nNodes;         /**< Number of nodes of the polyMESH */
    int nCells;         /**< Number of elements(cells) of the polyMESH */
    double* x;          /**< 1D array with size of polyMESH.nNodes. Coordinate x of observation points */
    double* y;          /**< 1D array with size of polyMESH.nNodes. Coordinate y of observation points */
    double* m_cData;      /**< 1D array with size of polyMESH.nCells. Cell data*/
    string DataName;    /**< Name of the cell data */
    // int* nEdges;        /**< 1D array with size of polyMESH.nCells, represents how many edges of a cell*/
    // int* el2nod;        /**< 1D array with size of sum(polyMESH.nEdges). Node index of connection of a cell */
    vector< vector<int> > el2nod;
    void WriteGeometry2Gmsh(string outfile);
};

#endif