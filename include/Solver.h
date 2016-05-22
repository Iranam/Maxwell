#ifndef SOLVER_H
#define SOLVER_H

#include<Eigen/Core>
using Eigen::Vector3d;
using Eigen::Vector3i;
#include"Typedefs.h"
#include"Field.h"

enum BoundaryCondition{CONST,CYCLIC,REFLECTIVE,PML};
struct SolverParams{
  SolverParams();
  Vector3i size;
  BoundaryCondition bcond[6];
  Real pml_sigma;
  uint pml_thickness;
};

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

//Solver:
//works on a linear regular grid; tau:=h/c;
//E is updated on steps 1,2,3,4...
//H is updated on steps 1/2,1+1/2,2+1/2,...
//E is defined at (0,0,0), (0,0,1), (0,1,0)...
//H is defined at (-1/2,-1/2,-1/2), (-1/2,-1/2,+1/2)...
//Each field has a 1 width boundary which is not updated normally, it is used for boundary conditions;
//for E that is upper boundaries (size.x,y,z)
//for H that is lower boundaries (0,y,z);

class Solver{
public:
  Solver(SolverParams P);
  ~Solver();
  Vector3i size;
  Field E,H;
  BoundaryCondition bcond[6];
  /*b is condition at
   *b[0] - z=0
   *b[1] - z=max
   *b[2] - y=0
   *b[3] - y=max
   *b[4] - x=0
   *b[5] - x=max
   */
  void setE(Vector3d (*f)(Real,Real,Real));//f is a vector function in range [0,1]^3;
  void setH(Vector3d (*f)(Real,Real,Real));
  void step();//changes both E and H
  void print();
private:
  PMLayer* pml;
  char field_sign;//=1 for E and =-1 for H
	uint mx,my,mz;
	void copyExy();
	void copyExz();
	void copyEyz();
	void copyHxy();
	void copyHxz();
	void copyHyz();
	void reflectExy();
	void reflectExz();
	void reflectEyz();
	void reflectHxy();
	void reflectHxz();
	void reflectHyz();
  void applyPML();
  void applyPMLat(uint x,uint y,uint z,PMLpoint*& p);
	void handleBoundary();
	void handleBoundaryE();
	void handleBoundaryH();
};

#endif // SOLVER_H
