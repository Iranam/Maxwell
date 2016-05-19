Vector3d v=Vector3d::Zero();
if((bc[0]==PML)&&(z<t)){
	Real d=(t-z)/t;
	v.z()=sigma*d*d;
}
if((bc[1]==PML)&&(z>mz-t)){
	Real d=(t+z-mz)/t;
	v.z()=sigma*d*d;
}
if((bc[2]==PML)&&(y<t)){
	Real d=(t-y)/t;
	v.y()=sigma*d*d;
}
if((bc[3]==PML)&&(y>my-t)){
	Real d=(t+y-my)/t;
	v.y()=sigma*d*d;
}
if((bc[4]==PML)&&(x<t)){
	Real d=(t-x)/t;
	v.x()=sigma*d*d;
}
if((bc[5]==PML)&&(x>mx-t)){
	Real d=(t+x-mx)/t;
	v.x()=sigma*d*d;
}
p->sigma=v;
p->E.setZero();
p->H.setZero();
p++;
