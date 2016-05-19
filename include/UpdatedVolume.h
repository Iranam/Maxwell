#ifndef UPDATEDVOLUME_H
#define UPDATEDVOLUME_H

enum VolumeType{VACUUM,CYCLIC_BOUNDARY,PML_X}

class UpdatedVolume{
  public:
    ~UpdatedVolume();
    UpdatedVolume(VolumeType=VACUUM);
    UpdatedVolume(VolumeType,Vector3i pos,Vector3i size)
    Vector3i pos,size;
    VolumeType volume_type;
    Real** scalar_fields;
    /*
      VACUUM and BOUNDARY volumes do not have any additional scalar fields
      scalar fields such as permittivity and conductivity must be introduced in matter, but are not yet implemented
      PML also introduces several auxilary fields
    */
};

#endif // UPDATEDVOLUME_H
