#include <Eigen/Core>
using Eigen::Vector3d;
using Eigen::Vector3i;
#include "Typedefs.h"
#include "Field.h"

#include<iostream>

Field::Field(Vector3i nsize):size{nsize}{
  uint N=size.x()*size.y()*size.z();
  data=new Vector3d[N];
  for(uint i=0;i<N;i++)data[i].setZero();
}

Field::~Field(){
  delete[] data;
}

Vector3d& Field::operator()(uint j)const{
  return data[j];
}

Vector3d& Field::operator()(int x,int y,int z)const{
  return data[(z*size.y()+y)*size.x()+x];
}

Vector3d& Field::operator()(Vector3i p)const{
  return data[(p.z()*size.y()+p.y())*size.x()+p.x()];
}

Vector3i Field::getSize()const{return size;}

Vector3d Field::rot(Vector3i pos)const{
  Vector3d v[8];{
    uint x[2],y[2],z[2];
    x[0]=pos.x();
    y[0]=pos.y();
    z[0]=pos.z();
    x[1]=x[0]+1;
    y[1]=y[0]+1;
    z[1]=z[0]+1;
    for(int i=0;i<0b1000;i++){//it's magic
      int a=i&1;
      int b=(i>>1)&1;
      int c=(i>>2)&1;
      v[i]=(*this)(x[a],y[b],z[c]);
    }
  }
  Vector3d ret(0,0,0);
  for(int a=0;a<3;a++){
    for(int d=0;d<2;d++){
      int b=(a+d+1)%3;
      for(int i=0;i<8;i++){
        ret(a)+=v[i](b)*(((i>>(3-a-b))^d)&1?-1:1);
      }
    }
  }
  ret*=sigma/4;
  return ret;
}

//std::ostream& operator<<(std::ostream &os,myclass const &m);


