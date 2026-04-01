#include "mc_growth.h" 


float Epa_term1(int* nb_type, float nb_a[][3], int i);
float Epa_term2(int* nb_type, float nb_a[][3], int i, int j);
float Epa_term3(int* nb_type, float nb_a[][3], int i, int j);


void calc_Edef_pa(int x, int y, int z) {
  int dir,factor;
  float E;
  float ax_,ay_,az_,ax2,ay2,az2,A_xx,A_yy,A_zz,A_xy,A_xz,A_yz;
  float Bx,By,Bz,Bxx,Bxy,Bxz,Byx,Byy,Byz,Bzx,Bzy,Bzz;
  unsigned short int config_;
  float* p;
  int nb_type[dir_number+1];
  float nb_a[dir_number+1][3];
  
  if (atoms(x,y,z).type == 0)  {
    atoms(x,y,z).Edef = 0;
    return;
  }
  
  if( atoms(x,y,z).defect_f.x != 0 ||  atoms(x,y,z).defect_f.y != 0 ||
      atoms(x,y,z).defect_f.z !=0) {
      atoms(x,y,z).Edef = 0;
      return;
  }
  
  config_= atoms(x,y,z).config;
  p = AA_[config_]; 
  A_xx = p[0]; A_yy = p[1]; A_zz = p[2]; A_xy = p[3]; A_xz = p[4]; A_yz = p[5];
  Bx = atoms(x,y,z).B0.x; By = atoms(x,y,z).B0.y; Bz = atoms(x,y,z).B0.z;
  
  neighbors_t nbs;
  factor = atoms.neighbors(x,y,z,nbs);
  

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
    //v_ochered_Edef(x2,y2,z2);
  }
  //v_ochered_Edef(x,y,z);
  ax_=ay_=az_=0;
  //random_displacements(&ax_,&ay_,&az_,config_);
  ax_ -= 0.5*( A_xx*Bx+A_xy*By+A_xz*Bz );
  ay_ -= 0.5*( A_xy*Bx+A_yy*By+A_yz*Bz );
  az_ -= 0.5*( A_xz*Bx+A_yz*By+A_zz*Bz );
  
  if (z%2==0) factor=1; else factor=-1;
  for (dir=0; dir<dir_number; dir++) {
    ax2=atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.x;
    ay2=atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.y;
    az2=atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.z;
    nb_type[dir] = atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).type;
    nb_a[dir][0] = ax2 * factor;
    nb_a[dir][1] = ay2 * factor;
    nb_a[dir][2] = az2 * factor;
  }
  nb_type[dir_number] = atoms(x,y,z).type;
  //nb_a[dir_number][0] = ax_*factor;  // =ax[x][y][z]*factor; 
  //nb_a[dir_number][1] = ay_*factor;  // =ay[x][y][z]*factor; 
  //nb_a[dir_number][2] = az_*factor;  // =az[x][y][z]*factor;
  nb_a[dir_number][0]=atoms(x,y,z).a.x*factor; // actual displacements
  nb_a[dir_number][1]=atoms(x,y,z).a.y*factor; 
  nb_a[dir_number][2]=atoms(x,y,z).a.z*factor;
  
  E=0;
  E+=0.5*Epa_term1(nb_type,nb_a,0);
  E+=0.5*Epa_term1(nb_type,nb_a,1);
  E+=0.5*Epa_term1(nb_type,nb_a,2);
  E+=0.5*Epa_term1(nb_type,nb_a,3);
  E+=1./3*Epa_term2(nb_type,nb_a,0,1);
  E+=1./3*Epa_term2(nb_type,nb_a,0,2);
  E+=1./3*Epa_term2(nb_type,nb_a,0,3);
  E+=1./3*Epa_term2(nb_type,nb_a,1,2);
  E+=1./3*Epa_term2(nb_type,nb_a,1,3);
  E+=1./3*Epa_term2(nb_type,nb_a,2,3);
  E+=1./3*Epa_term3(nb_type,nb_a,0,4);
  E+=1./3*Epa_term3(nb_type,nb_a,0,5);
  E+=1./3*Epa_term3(nb_type,nb_a,0,6);
  E+=1./3*Epa_term3(nb_type,nb_a,1,7);
  E+=1./3*Epa_term3(nb_type,nb_a,1,8);
  E+=1./3*Epa_term3(nb_type,nb_a,1,9);
  E+=1./3*Epa_term3(nb_type,nb_a,2,10);
  E+=1./3*Epa_term3(nb_type,nb_a,2,11);
  E+=1./3*Epa_term3(nb_type,nb_a,2,12);
  E+=1./3*Epa_term3(nb_type,nb_a,3,13);
  E+=1./3*Epa_term3(nb_type,nb_a,3,14);
  E+=1./3*Epa_term3(nb_type,nb_a,3,15);
  if (!(E>-1e100 && E<1e100)) {fprintf(stderr,"Error!!! E(%d,%d,%d)=%f\n",x,y,z,E); exit(0);}
  atoms(x,y,z).Edef_pa = E;
}


float Epa_term1(int* nb_type, float nb_a[][3], int i) {
  float eps, in_brackets;
  if (nb_type[i]==0) return 0;
  eps=0.5*(param.eps[nb_type[i]]+param.eps[nb_type[dir_number]]);
  in_brackets=2*( dx_neighbor[i]*(nb_a[i][0]-nb_a[dir_number][0]) 
                + dy_neighbor[i]*(nb_a[i][1]-nb_a[dir_number][1]) 
                + dz_neighbor[i]*(nb_a[i][2]-nb_a[dir_number][2]) 
                - 3*eps);
  return param.A*in_brackets*in_brackets;
}


float Epa_term2(int* nb_type, float nb_a[][3], int i, int j) {
  float sum_eps, in_brackets;
  int dxi,dyi,dzi,dxj,dyj,dzj;
  if (nb_type[i]==0 || nb_type[j]==0) return 0;
  sum_eps=0.5*(param.eps[nb_type[i]]+param.eps[nb_type[j]]+2*param.eps[nb_type[dir_number]]);
  dxi=dx_neighbor[i];  dxj=dx_neighbor[j];
  dyi=dy_neighbor[i];  dyj=dy_neighbor[j];
  dzi=dz_neighbor[i];  dzj=dz_neighbor[j];
  in_brackets = dxi*nb_a[j][0] + dxj*nb_a[i][0] - (dxi+dxj)*nb_a[dir_number][0] 
              + dyi*nb_a[j][1] + dyj*nb_a[i][1] - (dyi+dyj)*nb_a[dir_number][1] 
              + dzi*nb_a[j][2] + dzj*nb_a[i][2] - (dzi+dzj)*nb_a[dir_number][2] 
              + sum_eps;
  return param.B*in_brackets*in_brackets;
}


float Epa_term3(int* nb_type, float nb_a[][3], int i, int j) {
  float sum_eps, in_brackets;
  int dxi,dyi,dzi,dxj,dyj,dzj;
  if (nb_type[i]==0 || nb_type[j]==0) return 0;
  sum_eps=0.5*(2*param.eps[nb_type[i]]+param.eps[nb_type[j]]+param.eps[nb_type[dir_number]]);
  dxi=dx_neighbor[i];  dxj=dx_neighbor[j];
  dyi=dy_neighbor[i];  dyj=dy_neighbor[j];
  dzi=dz_neighbor[i];  dzj=dz_neighbor[j];
  in_brackets = (2*dxi-dxj)*nb_a[i][0] - dxi*nb_a[j][0] + (dxj-dxi)*nb_a[dir_number][0] 
              + (2*dyi-dyj)*nb_a[i][1] - dyi*nb_a[j][1] + (dyj-dyi)*nb_a[dir_number][1] 
              + (2*dzi-dzj)*nb_a[i][2] - dzi*nb_a[j][2] + (dzj-dzi)*nb_a[dir_number][2] 
              + sum_eps;
  return param.B*in_brackets*in_brackets;
}


