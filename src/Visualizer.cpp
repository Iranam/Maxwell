#include <Eigen/Core>
using Eigen::Vector3d;
using Eigen::Vector3i;
#include"Typedefs.h"
#include"Field.h"
#include"Solver.h"
#include<mgl2/qt.h>
#include"Visualizer.h"
#include<iostream>
using namespace std;

Visualizer::Visualizer(Solver* nsolver):solver{nsolver}{
  fieldtype=E;
  framedelay=100;
}
int Visualizer::Draw(mglGraph *gr){
  mglData v[3];
  uint mx,my,mz;
  {Vector3i size=solver->size;
  mx=size.x();
  my=size.y();
  mz=size.z();}
  for(uint j=0;j<3;j++)v[j].Create(mx,my,mz);
  for(uint z=0;z<mz;z++)
  for(uint y=0;y<my;y++)
  for(uint x=0;x<mx;x++){
    Vector3d field=(fieldtype==E)?solver->E(x,y,z):solver->H(x,y,z);
    uint i0=(z*my+y)*mx+x;
    for(uint j=0;j<3;j++)v[j].a[i0]=field(j);
  }
  v[0].a[0]=0.5;//I don't know, somehow it fixes the mgl's vector scaling
  gr->Rotate(60,130);
  gr->Axis("_xyz");
  gr->Vect(v[0],v[1],v[2]);
  gr->Box();
  return 0;
}
void Visualizer::makeGif(const char* file,const uint Nframes){
  mglGraph gr;
  gr.StartGIF(file,framedelay);
  cout<<"Making gif; frame:"<<endl;
  for(uint n=0;n<Nframes;n++){
    cout<<'\r'<<n<<'/'<<Nframes<<flush;
    if(n!=0)solver->step();
    gr.NewFrame();
    Draw(&gr);
    gr.EndFrame();
  }
  gr.CloseGIF();
}
void Visualizer::Test(const char* file,const uint Nframes){
  mglGraph gr;
  mglData v[3]=mglData(12,12,12);
  v[0].a[0]=1;
  gr.StartGIF(file,framedelay);
  cout<<"Making gif; frame:"<<endl;
  float len=0.01;
  for(uint n=0;n<Nframes;n++){
    cout<<'\r'<<n<<'/'<<Nframes<<flush;
    if(n!=0)solver->step();
    gr.NewFrame();
    v[0].a[12*12*5+12*5+5]=10.+1./n;
    gr.Rotate(60,60);
    gr.Axis("_xyz");
    gr.Vect(v[0],v[1],v[2]);
    gr.Box();
    gr.EndFrame();
  }
  gr.CloseGIF();
}
int Visualizer::run(){
  mglQT window(this,"MathGL examples");
  return window.Run();
}
