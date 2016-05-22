#include"Solver.h"
#include<iostream>
using std::cout;
using std::endl;

extern bool DEBUG;

SolverParams::SolverParams(){
  size=Vector3i(16,16,16);
  bcond[0]=CYCLIC;
  bcond[1]=CYCLIC;
  bcond[2]=CYCLIC;
  bcond[3]=CYCLIC;
  bcond[4]=CYCLIC;
  bcond[5]=CYCLIC;
  pml_sigma=1;
  pml_thickness=3;
};

PMLayer::PMLayer(uint t,Real sigma,Vector3i s,BoundaryCondition bc[6]):thickness{t}{
  uint mx=s.x(),my=s.y(),mz=s.z(),t2=t*2;
  uint N=mx*my*mz-(mx-t2)*(my-t2)*(mz-t2);
  points=new PMLpoint[N];
  offset[0]=t*mx*my;
  offset[1]=offset[0]*2;
  offset[2]=offset[1]+t*mx*(mz-t2);
  offset[3]=offset[2]+t*mx*(mz-t2);
  offset[4]=offset[3]+t*(my-t2)*(mz-t2);
	PMLpoint* p=points;
	for(uint z=0;z<t;z++)
	for(uint y=0;y<my;y++)
	for(uint x=0;x<mx;x++){
    #include"SetPML.inc"
	}
	for(uint z=mz-t;z<mz;z++)
	for(uint y=0;y<my;y++)
	for(uint x=0;x<mx;x++){
    #include"SetPML.inc"
	}
	for(uint z=t;z<mz-t;z++)
	for(uint y=0;y<t;y++)
	for(uint x=0;x<mx;x++){
    #include"SetPML.inc"
	}
	for(uint z=t;z<mz-t;z++)
	for(uint y=my-t;y<my;y++)
	for(uint x=0;x<mx;x++){
    #include"SetPML.inc"
	}
	for(uint z=t;z<mz-t;z++)
	for(uint y=t;y<my-t;y++)
	for(uint x=0;x<t;x++){
    #include"SetPML.inc"
	}
	for(uint z=t;z<mz-t;z++)
	for(uint y=t;y<my-t;y++)
	for(uint x=mx-t;x<mx;x++){
    #include"SetPML.inc"
	}
}

PMLayer::~PMLayer(){
  delete[] points;
}

Solver::Solver(SolverParams P):size{P.size},E{P.size+Vector3i(1,1,1)},H{P.size+Vector3i(1,1,1)}{
	for(uint i=0;i<6;i++)bcond[i]=P.bcond[i];
	handleBoundaryE();
	handleBoundaryH();
	pml=nullptr;
	for(uint i=0;i<6;i++){
    if(bcond[i]==PML){
      pml=new PMLayer(P.pml_thickness,P.pml_sigma,size,bcond);
      break;
    }
  }
}

Solver::~Solver(){
	if(!pml)delete[] pml;
}

void Solver::setE(Vector3d (*f)(Real,Real,Real)){
  for(int z=0;z<size.z();z++)
  for(int y=0;y<size.y();y++)
  for(int x=0;x<size.x();x++)
		E(x,y,z)=f((Real)x/size.x(),(Real)y/size.y(),(Real)z/size.z());
	handleBoundaryE();
}

void Solver::setH(Vector3d (*f)(Real,Real,Real)){
  for(int z=1;z<=size.z();z++)
  for(int y=1;y<=size.y();y++)
  for(int x=1;x<=size.x();x++)
		H(x,y,z)=f((Real)x/size.x(),(Real)y/size.y(),(Real)z/size.z());
	handleBoundaryH();
}

void Solver::copyExy(){
	for(int y=0;y<=size.y();y++)
	for(int x=0;x<=size.x();x++)
		E(x,y,size.z())=E(x,y,0);
}
void Solver::copyExz(){
	for(int z=0;z<size.z();z++)
	for(int x=0;x<=size.x();x++)
		E(x,size.y(),z)=E(x,0,z);
}
void Solver::copyEyz(){
	for(int z=0;z<size.z();z++)
	for(int y=0;y<size.y();y++)
		E(size.x(),y,z)=E(0,y,z);
}
void Solver::copyHxy(){
	for(int y=0;y<=size.y();y++)
	for(int x=0;x<=size.x();x++)
		H(x,y,0)=H(x,y,size.z());
}
void Solver::copyHxz(){
	for(int z=1;z<=size.z();z++)
	for(int x=0;x<=size.x();x++)
		H(x,0,z)=H(x,size.y(),z);
}
void Solver::copyHyz(){
	for(int z=1;z<=size.z();z++)
	for(int y=1;y<=size.y();y++)
		H(0,y,z)=H(size.x(),y,z);
}

void Solver::reflectExy(){
	for(int y=0;y<=size.y();y++)
	for(int x=0;x<=size.x();x++){
		E(x,y,size.z()).x()=0;
		E(x,y,size.z()).y()=0;
		E(x,y,size.z()).z()=E(x,y,size.z()-1).z();
	}
}
void Solver::reflectExz(){
	for(int z=0;z<size.z();z++)
	for(int x=0;x<=size.x();x++){
		E(x,size.y(),z).x()=0;
		E(x,size.y(),z).y()=E(x,size.y()-1,z).y();
		E(x,size.y(),z).z()=0;
	}
}
void Solver::reflectEyz(){
	for(int z=0;z<size.z();z++)
	for(int y=0;y<size.y();y++){
		E(size.x(),y,z).x()=E(size.x()-1,y,z).x();
		E(size.x(),y,z).y()=0;
		E(size.x(),y,z).z()=0;
	}
}
void Solver::reflectHxy(){
	for(int y=0;y<=size.y();y++)
	for(int x=0;x<=size.x();x++){
		H(x,y,0).x()=H(x,y,1).x();
		H(x,y,0).y()=H(x,y,1).y();
		H(x,y,0).z()=0;
	}
}
void Solver::reflectHxz(){
	for(int z=1;z<=size.z();z++)
	for(int x=0;x<=size.x();x++){
		H(x,0,z).x()=H(x,1,z).x();
		H(x,0,z).y()=0;
		H(x,0,z).z()=H(x,1,z).z();
	}
}
void Solver::reflectHyz(){
	for(int z=1;z<=size.z();z++)
	for(int y=1;y<=size.y();y++){
		H(0,y,z).x()=0;
		H(0,y,z).y()=H(1,y,z).y();
		H(0,y,z).z()=H(1,y,z).z();
	}
}

void Solver::applyPMLat(uint x,uint y,uint z,PMLpoint*& p){
	/*Perfectly matched layer is applied in 4 stages:
	 *1.Update E without PML, as in E+=rot(B)
	 *2.Update auxilary PML vector Psi
	 *3.E+=Psi
	 *4.Ex/=1+sy+sz
	 */
	Vector3d *Aref,*psiref;
	Field* Bref;
  if(field_sign==1){
    Aref=&E(x,y,z);
		psiref=&(p->E);
    Bref=&H;
  }else{
    Aref=&H(x+1,y+1,z+1);
		psiref=&(p->H);
    Bref=&E;
  }
	Vector3d& A=*Aref;
	Vector3d& psi=*psiref;
	Field& B=*Bref;
	Vector3d
    B000=B(x+0,y+0,z+0),
    B001=B(x+0,y+0,z+1),
    B010=B(x+0,y+1,z+0),
    B011=B(x+0,y+1,z+1),
    B100=B(x+1,y+0,z+0),
    B101=B(x+1,y+0,z+1),
    B110=B(x+1,y+1,z+0),
    B111=B(x+1,y+1,z+1);
	Real
		dxBy=(-B000.y()-B001.y()-B010.y()-B011.y()+B100.y()+B101.y()+B110.y()+B111.y())/4,
		dxBz=(-B000.z()-B001.z()-B010.z()-B011.z()+B100.z()+B101.z()+B110.z()+B111.z())/4,
		dyBx=(-B000.x()-B001.x()+B010.x()+B011.x()-B100.x()-B101.x()+B110.x()+B111.x())/4,
		dyBz=(-B000.z()-B001.z()+B010.z()+B011.z()-B100.z()-B101.z()+B110.z()+B111.z())/4,
		dzBx=(-B000.x()+B001.x()-B010.x()+B011.x()-B100.x()+B101.x()-B110.x()+B111.x())/4,
		dzBy=(-B000.y()+B001.y()-B010.y()+B011.y()-B100.y()+B101.y()-B110.y()+B111.y())/4;
	//DEBUG
	Vector3d s=p->sigma;
	psi.x()+=(s.z()*dyBz-s.y()*dzBy)*field_sign-s.y()*s.z()*A.x();
	psi.y()+=(s.x()*dzBx-s.z()*dxBz)*field_sign-s.z()*s.x()*A.y();
	psi.z()+=(s.y()*dxBy-s.x()*dyBx)*field_sign-s.x()*s.y()*A.z();
	A+=psi;
	A.x()/=1+s.y()+s.z();
	A.y()/=1+s.z()+s.x();
	A.z()/=1+s.x()+s.y();
  /*if(DEBUG&&A!=Vector3d::Zero()){//DEBUG
    std::cout<<"DEBUG:psi"<<((field_sign==1)?"E":"H")<<"("<<x<<":"<<y<<":"<<z<<")="<<psi.transpose()<<std::endl;
		std::cout<<"B000="<<B000.transpose()<<std::endl;
		std::cout<<"B001="<<B001.transpose()<<std::endl;
		std::cout<<"B010="<<B010.transpose()<<std::endl;
		std::cout<<"B011="<<B011.transpose()<<std::endl;
		std::cout<<"B100="<<B100.transpose()<<std::endl;
		std::cout<<"B101="<<B101.transpose()<<std::endl;
		std::cout<<"B110="<<B110.transpose()<<std::endl;
		std::cout<<"B111="<<B111.transpose()<<std::endl;
		std::cout<<"dxBy="<<dxBy<<std::endl;
    std::cout<<"dxBz="<<dxBz<<std::endl;
    std::cout<<"dyBx="<<dyBx<<std::endl;
    std::cout<<"dyBz="<<dyBz<<std::endl;
    std::cout<<"dzBx="<<dzBx<<std::endl;
    std::cout<<"dzBy="<<dzBy<<std::endl;
    std::cout<<"sigma="<<s.transpose()<<std::endl;
    std::cout<<"A="<<A.transpose()<<std::endl;
	}//DEBUG*/
  p++;
}

void Solver::applyPML(){
	uint h=pml->thickness;
	PMLpoint* p=pml->points;
	uint mx=size.x();
	uint my=size.y();
	uint mz=size.z();
	for(uint z=0;z<h;z++)
	for(uint y=0;y<my;y++)
	for(uint x=0;x<mx;x++)
		applyPMLat(x,y,z,p);
	for(uint z=mz-h;z<mz;z++)
	for(uint y=0;y<my;y++)
	for(uint x=0;x<mx;x++)
		applyPMLat(x,y,z,p);
	for(uint z=h;z<mz-h;z++)
	for(uint y=0;y<h;y++)
	for(uint x=0;x<mx;x++)
		applyPMLat(x,y,z,p);
	for(uint z=h;z<mz-h;z++)
	for(uint y=my-h;y<my;y++)
	for(uint x=0;x<mx;x++)
		applyPMLat(x,y,z,p);
	for(uint z=h;z<mz-h;z++)
	for(uint y=h;y<my-h;y++)
	for(uint x=0;x<h;x++)
		applyPMLat(x,y,z,p);
	for(uint z=h;z<mz-h;z++)
	for(uint y=h;y<my-h;y++)
	for(uint x=mx-h;x<mx;x++)
		applyPMLat(x,y,z,p);
}

void Solver::handleBoundaryE(){
	if(bcond[5]==CYCLIC)copyEyz();
	else if(bcond[5]==REFLECTIVE)reflectEyz();
	if(bcond[3]==CYCLIC)copyExz();
	else if(bcond[3]==REFLECTIVE)reflectExz();
	if(bcond[1]==CYCLIC)copyExy();
	else if(bcond[1]==REFLECTIVE)reflectExy();
}

void Solver::handleBoundaryH(){
	if(bcond[4]==CYCLIC)copyHyz();
	else if(bcond[4]==REFLECTIVE)reflectEyz();
	if(bcond[2]==CYCLIC)copyHxz();
	else if(bcond[2]==REFLECTIVE)reflectExz();
	if(bcond[0]==CYCLIC)copyHxy();
	else if(bcond[0]==REFLECTIVE)reflectHxy();
}

void Solver::handleBoundary(){
	if(pml)applyPML();
	if(field_sign==1)handleBoundaryE();
	else handleBoundaryH();
}

void Solver::step(){
  mx=size.x(),my=size.y(),mz=size.z();

  //calculate E
  Vector3d* pE=E.data;
  for(uint z=0;z<mz;z++,pE+=E.size.x())
  for(uint y=0;y<my;y++,pE++)
  for(uint x=0;x<mx;x++,pE++){
    (*pE)+=H.rot(Vector3i(x,y,z));
  }

  //apply boundary E
  field_sign=1;
	handleBoundary();

  //calculate H
  Vector3d* pH=H.data+H.size.x()*H.size.y()+H.size.x()+1;
  for(uint z=0;z<mz;z++,pH+=H.size.x())
  for(uint y=0;y<my;y++,pH++)
  for(uint x=0;x<mx;x++,pH++){
    (*pH)-=E.rot(Vector3i(x,y,z));
  }

  //apply boundary H
  field_sign=-1;
	handleBoundary();
	//DEBUG
	Real Emax=0;
    for(uint z=0;z<size.z();z++)
    for(uint y=0;y<size.y();y++)
    for(uint x=0;x<size.x();x++){
      Real L=E(x,y,z).norm();
      if(L>Emax)Emax=L;
    }
    cout<<"Emax="<<Emax<<endl;

}

void Solver::print(){
	cout<<"E:"<<endl;
	for(int z=0;z<size.z()+1;z++)
    for(int y=0;y<size.y()+1;y++){
      for(int x=0;x<size.x()+1;x++){
        if(E(x,y,z)!=Vector3d::Zero())
          cout<<"("<<x<<":"<<y<<":"<<z<<"):"<<E(x,y,z).transpose()<<endl;
      }
		}
	cout<<"H:"<<endl;
	for(int z=0;z<size.z()+1;z++)
    for(int y=0;y<size.y()+1;y++){
      for(int x=0;x<size.x()+1;x++){
        if(H(x,y,z)!=Vector3d::Zero())
          cout<<"("<<x<<":"<<y<<":"<<z<<"):"<<H(x,y,z).transpose()<<endl;
      }
		}
}
