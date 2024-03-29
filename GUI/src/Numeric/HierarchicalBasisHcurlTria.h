// Gmsh - Copyright (C) 1997-2019 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/gmsh/issues.
//
// Contributed by Ismail Badia.
// Reference :  "Higher-Order Finite Element  Methods"; Pavel Solin, Karel
// Segeth ,
//                 Ivo Dolezel , Chapman and Hall/CRC; Edition : Har/Cdr (2003).

#ifndef HIERARCHICAL_BASIS_HCURL_TRIA_H
#define HIERARCHICAL_BASIS_HCURL_TRIA_H

#include "HierarchicalBasisHcurl.h"
#include <math.h>
#include <algorithm>

/*
 * MTriangle
 *
 *   v
 *   ^
 *   |
 *   2
 *   |`\
 *   |  `\
 *   |    `\
 *   |      `\
 *   v        `\
 *   0---------->1 --> u
 *
 *
 * Oriented Edges:
 *  e0={v0;v1}    e1={v1;v2}  e2={v2;v0}
 *  pe0,pe1,pe2<=pf
 *
 */
class HierarchicalBasisHcurlTria : public HierarchicalBasisHcurl {
public:
  HierarchicalBasisHcurlTria(int order);
  virtual ~HierarchicalBasisHcurlTria();
  virtual void generateBasis(double const &u, double const &v, double const &w,
                             std::vector<std::vector<double> > &vertexBasis,
                             std::vector<std::vector<double> > &edgeBasis,
                             std::vector<std::vector<double> > &faceBasis,
                             std::vector<std::vector<double> > &bubbleBasis,
                             std::string typeFunction)
  {
    if(typeFunction == "HcurlLegendre") {
      generateHcurlBasis(u, v, w, edgeBasis, faceBasis, bubbleBasis);
    }
    else if("CurlHcurlLegendre" == typeFunction) {
      generateCurlBasis(u, v, w, edgeBasis, faceBasis, bubbleBasis);
    }
    else {
      throw std::string("unknown typeFunction");
    }
  }
  virtual void orientEdge(int const &flagOrientation, int const &edgeNumber,
                          std::vector<std::vector<double> > &edgeBasis);
  virtual void orientFace(double const &u, double const &v, double const &w,
                          int const &flag1, int const &flag2, int const &flag3,
                          int const &faceNumber,
                          std::vector<std::vector<double> > &faceFunctions,
                          std::string typeFunction);

  virtual void getKeysInfo(std::vector<int> &functionTypeInfo,
                           std::vector<int> &orderInfo);

private:
  int _pf; // face function order
  int _pOrderEdge[3]; // Edge functions order (pOrderEdge[0] matches the edge 0
                      // order)
  static double
  _affineCoordinate(int const &j, double const &u,
                    double const &v); // affine coordinate lambdaj j=1..3

  // edgeBasis=[phie0_{0},...phie0_{pe0},phie1_{0},...phie1_{pe1}...; edge-based
  // bubble functions ] faceBasis=[ genuine bubble functions]
  virtual void
  generateHcurlBasis(double const &u, double const &v, double const &w,
                     std::vector<std::vector<double> > &edgeBasis,
                     std::vector<std::vector<double> > &faceBasis,
                     std::vector<std::vector<double> > &bubbleBasis);
  virtual void
  generateCurlBasis(double const &u, double const &v, double const &w,
                    std::vector<std::vector<double> > &edgeBasis,
                    std::vector<std::vector<double> > &faceBasis,
                    std::vector<std::vector<double> > &bubbleBasis);

  static double dotProduct(const std::vector<double> &u,
                           const std::vector<double> &v);
};

#endif
