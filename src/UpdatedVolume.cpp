//TODO
//UpdatedVolume module is WIP

#include <Eigen/Core>
using Eigen::Vector3i;
#include"Typedefs.h"
#include"UpdatedVolume.h"

UpdatedVolume::UpdatedVolume(VolumeType type):volume_type{type}{
  scalar_fields=nullptr;
}

UpdatedVolume::UpdatedVolume(VolumeType type,Vector3i npos,Vector3i nsize):volume_type{type},pos{npos},size{nsize}{
  switch(volume_type){
    case PML_X:{
      scalar_fields=new Real*[2];
      uint N=size.x()*size.y()*size.z();
      scalar_fields[0]=new Real[N];
      scalar_fields[1]=new Real[N];
    }
    default:scalar_fields=nullptr;
  }
}

UpdatedVolume::~UpdatedVolume(){
  switch(volume_type){
    case PML_X:{
      delete[] scalar_fields[1];
      delete[] scalar_fields[0];
      delete[] scalar_fields;
    }
  }
}
