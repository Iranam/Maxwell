#ifndef PML_H
#define PML_H

struct PMLpoint{
  Vector3d sigma;
  Vector3d E;
  Vector3d H;
};

class PMLayer{
public:
  PMLayer(uint thickness,Real sigma,Vector3i fieldsize,BoundaryCondition bc[6]);//fieldsize is size of updated area, excluding borders
  ~PMLayer();
  uint thickness;
  PMLpoint* points;
  uint offset[5];
  //PM layer surrounds all 6 sides of the bounding box
  //all the necessary data for all sides is stored in array "points"
  //"points" array is divided into 6 parts: one for each side of the cube
  //"offset" stores offsets to these parts
};

#endif
