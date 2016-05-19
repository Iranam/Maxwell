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
  return AMPLITUDE*Vector3d(0,0,0);
}

void printFields(Solver* solver){
  Vector3i size=solver->size;
	cout<<"E"<<endl;
	for(int z=0;z<size.z()+1;z++)
    for(int y=0;y<size.y()+1;y++){
      for(int x=0;x<size.x()+1;x++){
        if(solver->E(x,y,z)!=Vector3d::Zero())
          cout<<"("<<x<<":"<<y<<":"<<z<<"):"<<(solver->E(x,y,z)).transpose()<<endl;
      }
		}
	cout<<"H"<<endl;
	for(int z=0;z<size.z()+1;z++)
    for(int y=0;y<size.y()+1;y++){
      for(int x=0;x<size.x()+1;x++){
        if(solver->H(x,y,z)!=Vector3d::Zero()||true)
          cout<<"("<<x<<":"<<y<<":"<<z<<"):"<<(solver->H(x,y,z)).transpose()<<endl;
      }
		}
}


int main(){
  SolverParams P;
  P.size=Vector3i(8,8,8);
  P.bcond[0]=CYCLIC;
  P.bcond[1]=CYCLIC;
  P.bcond[2]=CYCLIC;
  P.bcond[3]=CYCLIC;
  P.bcond[4]=CYCLIC;
  P.bcond[5]=CYCLIC;
	P.pml_sigma=1;
	P.pml_thickness=1;
  Solver* solver=new Solver(P);
  solver->setE(impactE);
  solver->setH(impactH);
  printFields(solver);
  Visualizer vis(solver);
  vis.fieldtype=E;
  vis.framedelay=200;
  vis.makeGif("Waves.gif",32);
  delete solver;
  return 0;
}
