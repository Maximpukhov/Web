#include "mc_growth.h" 


void random_displacements(float* ax_,float* ay_,float* az_, unsigned short int config_); // вычислить
                       // случайное смещение (ax_,ay_,az_) атома относительно его равновесного положения
void randn2(float* x1, float* x2); // генерирует 2 нормально распределённых случайных числа



void axyz(int x, int y, int z) { // "обновление" смещения атома (x,y,z), с. 34-38
  float ax_,ay_,az_,ax2,ay2,az2,A_xx,A_yy,A_zz,A_xy,A_xz,A_yz;
  float Bx,By,Bz,Bxx,Bxy,Bxz,Byx,Byy,Byz,Bzx,Bzy,Bzz;
  int dir,factor,x2,y2,z2;
  unsigned short int config_;
  float* p;
  
  if( atoms(x,y,z).type == 0)
    return;
  config_ = atoms(x,y,z).config;
  
  p=AA_[config_]; 
  A_xx=p[0]; A_yy=p[1]; A_zz=p[2]; A_xy=p[3]; A_xz=p[4]; A_yz=p[5];
  
  Bx=atoms(x,y,z).B0.x;
  By=atoms(x,y,z).B0.y;
  Bz=atoms(x,y,z).B0.z;

  neighbors_t nbs;
  factor = atoms.neighbors(x,y,z,nbs);


  // Vectorised version start 
  /*
  float B_[3];
  float a[dir_number*3] __attribute((aligned(64)));

  for (dir=0; dir<dir_number; dir++){
    a[dir*3 + 0] = atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.x;
    a[dir*3 + 1] = atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.y;
    a[dir*3 + 2] = atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.z;
  }

  __m128 Bv[3] = { _mm_set_ps(atoms(x,y,z).B0.x,0,0,0), 
                   _mm_set_ps(atoms(x,y,z).B0.y,0,0,0),
                   _mm_set_ps(atoms(x,y,z).B0.z,0,0,0) };

  int vects = dir_number*3/4;
  for(int i=0; i < 3; i++){
    __m128 *m = (__m128 *)BB_new[config_][i];
    for(int j=0; j < vects; j++){
      __m128 tmp = _mm_mul_ps(m[j],((__m128*)a)[j]);
      Bv[i] = _mm_add_ps(Bv[i],tmp);
    }
    float *p = (float*)(Bv + i);
    B_[i] = 0;
    for(int j = 0; j < 4; j++){
      B_[i] += p[j];
    }
  }
  */
  // Vectorised version end


  for (dir=0; dir<dir_number; dir++) 
  {
    ax2=atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.x;
    ay2=atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.y;
    az2=atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.z;

    p=BB[config_][dir];
    Bxx=p[0]; Bxy=p[1]; Bxz=p[2]; Byx=p[3]; Byy=p[4]; Byz=p[5]; Bzx=p[6]; Bzy=p[7]; Bzz=p[8]; 
    Bx+=Bxx*ax2+Bxy*ay2+Bxz*az2;
    By+=Byx*ax2+Byy*ay2+Byz*az2;
    Bz+=Bzx*ax2+Bzy*ay2+Bzz*az2;
    v_ochered_Edef(nbs.x[dir],nbs.y[dir],nbs.z[dir]);
  }


  v_ochered_Edef(x,y,z);


  ax_=ay_=az_=0;
  random_displacements(&ax_,&ay_,&az_,config_);
  ax_ -= 0.5*( A_xx*Bx+A_xy*By+A_xz*Bz );
  ay_ -= 0.5*( A_xy*Bx+A_yy*By+A_yz*Bz );
  az_ -= 0.5*( A_xz*Bx+A_yz*By+A_zz*Bz );
  
  atoms(x,y,z).a.x = ax_;
  atoms(x,y,z).a.y = ay_;
  atoms(x,y,z).a.z = az_;
}


void random_displacements(float* ax_,float* ay_,float* az_, unsigned short int config_) {
  float a1,a2,a3,a4,coeff;
  float* p;
  coeff=sqrt(0.5*param.T);
  randn2(&a1,&a2);
  randn2(&a3,&a4);
  p=transform_array[config_];
  *ax_ = coeff* (p[0]*a1+p[3]*a2+p[4]*a3);
  *ay_ = coeff* (p[3]*a1+p[1]*a2+p[5]*a3);
  *az_ = coeff* (p[4]*a1+p[5]*a2+p[2]*a3);
}


void randn2(float* x1, float* x2) { // генерирует 2 нормально распределённых случайных числа
  float V1,V2,S,f,t;
  t=1.0/(RAND_MAX+1.);
  while(1) {
    V1= 2*t*(rand()+t*rand())-1;
    V2= 2*t*(rand()+t*rand())-1;
    S= V1*V1+V2*V2;
    if (S<1) break;
  }
  f=sqrt(-2*log(S)/S);
  *x1=V1*f;
  *x2=V2*f;
}


