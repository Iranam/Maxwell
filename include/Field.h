#ifndef FIELD_H
#define FIELD_H

const double sigma=1.;//=c*dt/h

class Field{//Vector field
public:
  Field(Vector3i);//Constructs a field with given size, filled with 0-vectors;
  ~Field();
  Vector3d& operator()(uint j)const;
  Vector3d& operator()(int x,int y,int z)const;
  Vector3d& operator()(Vector3i)const;
  Vector3i getSize()const;
  Vector3d rot(Vector3i)const;//Returns rotor of the vector field at point (x+1/2,y+1/2,z+1/2)
    //uses 8 vectors at positions (x,y,z) to (x+1,y+1,z+1)
    //warning: no range checking is performed, e.g. do not call rot(size);
  Vector3i size;
  Vector3d* data;
};
//std::ostream& operator<<(std::ostream,myclass const&);

#endif // FIELD_H
