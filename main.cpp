#include <math.h>
#include <iostream>
using namespace std;
#include <Eigen/Core>
using Eigen::Vector3d;
using Eigen::Vector3i;
#include"Field.h"
#include"Solver.h"
#include<mgl2/qt.h>
#include"Visualizer.h"

bool DEBUG=true;

void setTestWaves(Solver* solver){
  Vector3i size=solver->size;
	uint mx=size.x();
	uint my=size.y();
	uint mz=size.z();
  for(uint z=0;z<mz;z++)
    for(uint y=0;y<my;y++)
      for(uint x=0;x<mx;x++){
        //const double b=1./sqrt(2);
        double c=cos(2*M_PI*((double)x/mx));
        solver->E(x,y,z)=Vector3d(0,0,c);
        solver->H(x+1,y+1,z+1)=Vector3d(0,c,0);
      }
}
const Real
  IMPACT_MIN=0.3,
  IMPACT_MAX=0.7,
  IMPACT_WIDTH=IMPACT_MAX-IMPACT_MIN,
  AMPLITUDE=1;

Vector3d wavesE(Real x,Real y, Real z){
  return AMPLITUDE*Vector3d(0,0,cos(2*M_PI*x));
}

Vector3d wavesH(Real x,Real y, Real z){
  return AMPLITUDE*Vector3d(0,cos(2*M_PI*x),0);
}

Vector3d impactE(Real x,Real y,Real z){
  Real d=(x>IMPACT_MIN&&x<IMPACT_MAX)?sin(M_PI*(x-IMPACT_MIN)/IMPACT_WIDTH):0;
  return AMPLITUDE*Vector3d(0,0,d);
}

Vector3d impactH(Real x,Real y,Real z){
  Real d=(x>IMPACT_MIN&&x<IMPACT_MAX)?sin(M_PI*(x-IMPACT_MIN)/IMPACT_WIDTH):0;
  return AMPLITUDE*Vector3d(0,d,0);
}

Vector3d tableE(Real x,Real y,Real z){
  return AMPLITUDE*Vector3d(0,0,(x>IMPACT_MIN&&x<IMPACT_MAX)?1:0);
}

Vector3d tableH(Real x,Real y,Real z){
  return AMPLITUDE*Vector3d(0,(x>IMPACT_MIN&&x<IMPACT_MAX)?1:0,0);
}

Vector3d constant(Real x,Real y, Real z){
  return AMPLITUDE*Vector3d(1,0,0);
}

int main(){
  SolverParams P;
  P.size=Vector3i(16,16,16);
  P.bcond[0]=CYCLIC;
  P.bcond[1]=CYCLIC;
  P.bcond[2]=CYCLIC;
  P.bcond[3]=CYCLIC;
  P.bcond[4]=PML;
  P.bcond[5]=PML;
	P.pml_sigma=0.1;
	P.pml_thickness=5;
  Solver* solver=new Solver(P);
  solver->setE(tableE);
  solver->setH(tableH);
  /*for(uint i=0;i<1;i++){
    solver->step();
    solver->print();
  }*/
  Visualizer vis(solver);
  vis.fieldtype=E;
  vis.framedelay=200;
  vis.makeGif("Waves.gif",32);
  delete solver;
  return 0;
}
