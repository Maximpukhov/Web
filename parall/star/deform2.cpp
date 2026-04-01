#include "mc_growth.h" 



void add_A_term1(float* A, int* nb_type, int i);
void add_A_term2(float* A, int* nb_type, int i, int j);
void add_A_term3(float* A, int* nb_type, int i, int j);
void add_B_term1(float* B, int* nb_type, int i);
void add_B_term2(float* Bi, float* Bj, int* nb_type, int i, int j);
void add_B_term3(float* Bi, float* Bj, int* nb_type, int i, int j);
void add_B0_term1(float* B, int* nb_type, int i);
void add_B0_term2(float* B, int* nb_type, int i, int j);
void add_B0_term3(float* B, int* nb_type, int i, int j);
float E_term1(int* nb_type, float nb_a[][3], int i);
float E_term2(int* nb_type, float nb_a[][3], int i, int j);
float E_term3(int* nb_type, float nb_a[][3], int i, int j);
void jacobi(float a[4][4], float d[4], float v[4][4]); // диагонализация матрицы a методом Якоби




void calc_AA() {
  int i;
  unsigned short int config_;
  int nb_type[dir_number];
  float A_temp[6];
  for (config_=0; ; config_++) {
    fill_nb_type(nb_type,config_);
    for (i=0; i<6; i++) A_temp[i]=0;
    add_A_term1(A_temp,nb_type,0);
    add_A_term1(A_temp,nb_type,1);
    add_A_term1(A_temp,nb_type,2);
    add_A_term1(A_temp,nb_type,3);
    add_A_term2(A_temp,nb_type,0,1);
    add_A_term2(A_temp,nb_type,0,2);
    add_A_term2(A_temp,nb_type,0,3);
    add_A_term2(A_temp,nb_type,1,2);
    add_A_term2(A_temp,nb_type,1,3);
    add_A_term2(A_temp,nb_type,2,3);
    add_A_term3(A_temp,nb_type,0,4);
    add_A_term3(A_temp,nb_type,0,5);
    add_A_term3(A_temp,nb_type,0,6);
    add_A_term3(A_temp,nb_type,1,7);
    add_A_term3(A_temp,nb_type,1,8);
    add_A_term3(A_temp,nb_type,1,9);
    add_A_term3(A_temp,nb_type,2,10);
    add_A_term3(A_temp,nb_type,2,11);
    add_A_term3(A_temp,nb_type,2,12);
    add_A_term3(A_temp,nb_type,3,13);
    add_A_term3(A_temp,nb_type,3,14);
    add_A_term3(A_temp,nb_type,3,15);
    for (i=0; i<6; i++) AA[config_][i]=A_temp[i];
    if (config_==(Nconfig-1)) break;
  }
}


void calc_BB() {
  int i,k;
  unsigned short int config_;
  int nb_type[dir_number];
  float B_temp[dir_number][9];
  for (config_=0; ; config_++) {
    fill_nb_type(nb_type,config_);
    for (i=0; i<dir_number; i++)
      for (k=0; k<9; k++)
        B_temp[i][k]=0;
    add_B_term1(B_temp[0],nb_type,0);
    add_B_term1(B_temp[1],nb_type,1);
    add_B_term1(B_temp[2],nb_type,2);
    add_B_term1(B_temp[3],nb_type,3);
    add_B_term2(B_temp[0],B_temp[1],nb_type,0,1);
    add_B_term2(B_temp[0],B_temp[2],nb_type,0,2);
    add_B_term2(B_temp[0],B_temp[3],nb_type,0,3);
    add_B_term2(B_temp[1],B_temp[2],nb_type,1,2);
    add_B_term2(B_temp[1],B_temp[3],nb_type,1,3);
    add_B_term2(B_temp[2],B_temp[3],nb_type,2,3);
    add_B_term3(B_temp[0],B_temp[ 4],nb_type,0,4);
    add_B_term3(B_temp[0],B_temp[ 5],nb_type,0,5);
    add_B_term3(B_temp[0],B_temp[ 6],nb_type,0,6);
    add_B_term3(B_temp[1],B_temp[ 7],nb_type,1,7);
    add_B_term3(B_temp[1],B_temp[ 8],nb_type,1,8);
    add_B_term3(B_temp[1],B_temp[ 9],nb_type,1,9);
    add_B_term3(B_temp[2],B_temp[10],nb_type,2,10);
    add_B_term3(B_temp[2],B_temp[11],nb_type,2,11);
    add_B_term3(B_temp[2],B_temp[12],nb_type,2,12);
    add_B_term3(B_temp[3],B_temp[13],nb_type,3,13);
    add_B_term3(B_temp[3],B_temp[14],nb_type,3,14);
    add_B_term3(B_temp[3],B_temp[15],nb_type,3,15);
    for (i=0; i<dir_number; i++)
      for (k=0; k<9; k++){
        BB[config_][i][k]=B_temp[i][k];
      }
    if (config_==(Nconfig-1)) break;
  }
}


void calc_Edef(int x, int y, int z) {
  int dir,x2,y2,z2,factor;
  float E;
  float ax_,ay_,az_,ax2,ay2,az2,A_xx,A_yy,A_zz,A_xy,A_xz,A_yz;
  float Bx,By,Bz,Bxx,Bxy,Bxz,Byx,Byy,Byz,Bzx,Bzy,Bzz;
  unsigned short int config_;
  float* p;
  int nb_type[dir_number+1];
  float nb_a[dir_number+1][3];
  
  if (atoms(x,y,z).type == 0){
    atoms(x,y,z).Edef = 0;
    return;
  }
  
  if( atoms(x,y,z).defect_f.x!=0 || atoms(x,y,z).defect_f.y!=0 || 
      atoms(x,y,z).defect_f.z!=0) {
      atoms(x,y,z).Edef = 0;
      return;
  }
  
  config_= atoms(x,y,z).config;
  p=AA_[config_]; 
  A_xx=p[0]; A_yy=p[1]; A_zz=p[2]; A_xy=p[3]; A_xz=p[4]; A_yz=p[5];
  
  Bx = atoms(x,y,z).B0.x; 
  By = atoms(x,y,z).B0.y;
  Bz = atoms(x,y,z).B0.z;
  
  neighbors_t nbs;
  factor = atoms.neighbors(x,y,z,nbs);
  
  for (dir=0; dir<dir_number; dir++) {
    ax2 = atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.x;
    ay2 = atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.y;
    az2 = atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.z;
    p=BB[config_][dir];
    Bxx=p[0]; Bxy=p[1]; Bxz=p[2]; Byx=p[3]; Byy=p[4]; Byz=p[5]; Bzx=p[6]; Bzy=p[7]; Bzz=p[8]; 
    Bx+=Bxx*ax2+Bxy*ay2+Bxz*az2;
    By+=Byx*ax2+Byy*ay2+Byz*az2;
    Bz+=Bzx*ax2+Bzy*ay2+Bzz*az2;
  }
  ax_=ay_=az_=0;
  //random_displacements(&ax_,&ay_,&az_,config_);
  ax_ -= 0.5*( A_xx*Bx+A_xy*By+A_xz*Bz );
  ay_ -= 0.5*( A_xy*Bx+A_yy*By+A_yz*Bz );
  az_ -= 0.5*( A_xz*Bx+A_yz*By+A_zz*Bz );
  
  for (dir=0; dir<dir_number; dir++) {
    nb_type[dir] = atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).type;
    nb_a[dir][0] = atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.x * factor;
    nb_a[dir][1] = atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.y * factor;
    nb_a[dir][2] = atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).a.z * factor;
  }
  nb_type[dir_number]=atoms(x,y,z).type;
  nb_a[dir_number][0]=ax_*factor;  // =ax[x][y][z]*factor; // displacements according to minimal Energy
  nb_a[dir_number][1]=ay_*factor;  // =ay[x][y][z]*factor; 
  nb_a[dir_number][2]=az_*factor;  // =az[x][y][z]*factor;
  //nb_a[dir_number][0]=ax[x][y][z]*factor; 
  //nb_a[dir_number][1]=ay[x][y][z]*factor; 
  //nb_a[dir_number][2]=az[x][y][z]*factor;
  
  E=0;
  E+=E_term1(nb_type,nb_a,0);
  E+=E_term1(nb_type,nb_a,1);
  E+=E_term1(nb_type,nb_a,2);
  E+=E_term1(nb_type,nb_a,3);
  E+=E_term2(nb_type,nb_a,0,1);
  E+=E_term2(nb_type,nb_a,0,2);
  E+=E_term2(nb_type,nb_a,0,3);
  E+=E_term2(nb_type,nb_a,1,2);
  E+=E_term2(nb_type,nb_a,1,3);
  E+=E_term2(nb_type,nb_a,2,3);
  E+=E_term3(nb_type,nb_a,0,4);
  E+=E_term3(nb_type,nb_a,0,5);
  E+=E_term3(nb_type,nb_a,0,6);
  E+=E_term3(nb_type,nb_a,1,7);
  E+=E_term3(nb_type,nb_a,1,8);
  E+=E_term3(nb_type,nb_a,1,9);
  E+=E_term3(nb_type,nb_a,2,10);
  E+=E_term3(nb_type,nb_a,2,11);
  E+=E_term3(nb_type,nb_a,2,12);
  E+=E_term3(nb_type,nb_a,3,13);
  E+=E_term3(nb_type,nb_a,3,14);
  E+=E_term3(nb_type,nb_a,3,15);
  if (!(E>-1e100 && E<1e100)) {fprintf(stderr,"Error!!! E(%d,%d,%d)=%f\n",x,y,z,E); exit(0);}
  atoms(x,y,z).Edef = E;
}


void calc_B0(int x, int y, int z) {
  int dir,x2,y2,z2,factor;
  float B[3];
  int nb_type[dir_number+1];
  if( atoms(x,y,z).type == 0){
    //atoms(x,y,z).B0.x = atoms(x,y,z).B0.y = atoms(x,y,z).B0.z = 0; 
    memset(&atoms(x,y,z).B0, 0, sizeof(atoms(0,0,0).B0) );
    return;
  }

  neighbors_t nbs;
  factor = atoms.neighbors(x,y,z,nbs);
  for (dir=0; dir<dir_number; dir++) {
    nb_type[dir]=atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).type;
  }
  nb_type[dir_number]=atoms(x,y,z).type;
  
  B[0]=B[1]=B[2]=0;
  add_B0_term1(B,nb_type,0);
  add_B0_term1(B,nb_type,1);
  add_B0_term1(B,nb_type,2);
  add_B0_term1(B,nb_type,3);
  add_B0_term2(B,nb_type,0,1);
  add_B0_term2(B,nb_type,0,2);
  add_B0_term2(B,nb_type,0,3);
  add_B0_term2(B,nb_type,1,2);
  add_B0_term2(B,nb_type,1,3);
  add_B0_term2(B,nb_type,2,3);
  add_B0_term3(B,nb_type,0,4);
  add_B0_term3(B,nb_type,0,5);
  add_B0_term3(B,nb_type,0,6);
  add_B0_term3(B,nb_type,1,7);
  add_B0_term3(B,nb_type,1,8);
  add_B0_term3(B,nb_type,1,9);
  add_B0_term3(B,nb_type,2,10);
  add_B0_term3(B,nb_type,2,11);
  add_B0_term3(B,nb_type,2,12);
  add_B0_term3(B,nb_type,3,13);
  add_B0_term3(B,nb_type,3,14);
  add_B0_term3(B,nb_type,3,15);
  //B0x[x][y][z]=B[0]*factor;
  //B0y[x][y][z]=B[1]*factor;
  //B0z[x][y][z]=B[2]*factor;
  atoms(x,y,z).B0.x = B[0]*factor - param.defect_F * atoms(x,y,z).defect_f.x;
  atoms(x,y,z).B0.y = B[1]*factor - param.defect_F * atoms(x,y,z).defect_f.y;
  atoms(x,y,z).B0.z = B[2]*factor - param.defect_F * atoms(x,y,z).defect_f.z;
}


void add_B0_term1(float* B, int* nb_type, int i) {
  float eps,coeff;
  if (nb_type[i]==0) return;
  eps=0.5*(param.eps[nb_type[i]]+param.eps[nb_type[dir_number]]);
  coeff=24*param.A*eps;
  B[0]+=coeff*dx_neighbor[i];
  B[1]+=coeff*dy_neighbor[i];
  B[2]+=coeff*dz_neighbor[i];
}


void add_B0_term2(float* B, int* nb_type, int i, int j) {
  float sum_eps,coeff;
  if (nb_type[i]==0 || nb_type[j]==0) return;
  sum_eps=0.5*(param.eps[nb_type[i]]+param.eps[nb_type[j]]+2*param.eps[nb_type[dir_number]]);
  coeff=-2*param.B*sum_eps;
  B[0]+=coeff*(dx_neighbor[i]+dx_neighbor[j]);
  B[1]+=coeff*(dy_neighbor[i]+dy_neighbor[j]);
  B[2]+=coeff*(dz_neighbor[i]+dz_neighbor[j]);
}


void add_B0_term3(float* B, int* nb_type, int i, int j) {
  float sum_eps,coeff;
  if (nb_type[i]==0 || nb_type[j]==0) return;
  sum_eps=0.5*(2*param.eps[nb_type[i]]+param.eps[nb_type[j]]+param.eps[nb_type[dir_number]]);
  coeff=2*param.B*sum_eps;
  B[0]+=coeff*(dx_neighbor[j]-dx_neighbor[i]);
  B[1]+=coeff*(dy_neighbor[j]-dy_neighbor[i]);
  B[2]+=coeff*(dz_neighbor[j]-dz_neighbor[i]);
}


float E_term1(int* nb_type, float nb_a[][3], int i) {
  float eps, in_brackets;
  if (nb_type[i]==0) return 0;
  eps=0.5*(param.eps[nb_type[i]]+param.eps[nb_type[dir_number]]);
  in_brackets=2*( dx_neighbor[i]*(nb_a[i][0]-nb_a[dir_number][0]) 
                + dy_neighbor[i]*(nb_a[i][1]-nb_a[dir_number][1]) 
                + dz_neighbor[i]*(nb_a[i][2]-nb_a[dir_number][2]) 
                - 3*eps);
  return param.A*in_brackets*in_brackets;
}


float E_term2(int* nb_type, float nb_a[][3], int i, int j) {
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


float E_term3(int* nb_type, float nb_a[][3], int i, int j) {
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


void add_A_term1(float* A, int* nb_type, int i) {
  int dx,dy,dz;
  if (nb_type[i]==0) return;
  dx=2*dx_neighbor[i];
  dy=2*dy_neighbor[i];
  dz=2*dz_neighbor[i];
  A[0]+=param.A*dx*dx;
  A[1]+=param.A*dy*dy;
  A[2]+=param.A*dz*dz;
  A[3]+=param.A*dx*dy;
  A[4]+=param.A*dx*dz;
  A[5]+=param.A*dy*dz;
}


void add_A_term2(float* A, int* nb_type, int i, int j) {
  int dx,dy,dz;
  if (nb_type[i]==0 || nb_type[j]==0) return;
  dx=dx_neighbor[i]+dx_neighbor[j];
  dy=dy_neighbor[i]+dy_neighbor[j];
  dz=dz_neighbor[i]+dz_neighbor[j];
  A[0]+=param.B*dx*dx;
  A[1]+=param.B*dy*dy;
  A[2]+=param.B*dz*dz;
  A[3]+=param.B*dx*dy;
  A[4]+=param.B*dx*dz;
  A[5]+=param.B*dy*dz;
}


void add_A_term3(float* A, int* nb_type, int i, int j) {
  int dx,dy,dz;
  if (nb_type[i]==0 || nb_type[j]==0) return;
  dx=dx_neighbor[j]-dx_neighbor[i];
  dy=dy_neighbor[j]-dy_neighbor[i];
  dz=dz_neighbor[j]-dz_neighbor[i];
  A[0]+=param.B*dx*dx;
  A[1]+=param.B*dy*dy;
  A[2]+=param.B*dz*dz;
  A[3]+=param.B*dx*dy;
  A[4]+=param.B*dx*dz;
  A[5]+=param.B*dy*dz;
}


void add_B_term1(float* B, int* nb_type, int i) {
  int dx,dy,dz;
  if (nb_type[i]==0) return;
  dx=dx_neighbor[i];
  dy=dy_neighbor[i];
  dz=dz_neighbor[i];
  B[0]+= -8*param.A*dx*dx;
  B[1]+= -8*param.A*dx*dy;
  B[2]+= -8*param.A*dx*dz;
  B[3]+= -8*param.A*dy*dx;
  B[4]+= -8*param.A*dy*dy;
  B[5]+= -8*param.A*dy*dz;
  B[6]+= -8*param.A*dz*dx;
  B[7]+= -8*param.A*dz*dy;
  B[8]+= -8*param.A*dz*dz;
}


void add_B_term2(float* Bi, float* Bj, int* nb_type, int i, int j) {
  int dx,dy,dz,dxi,dyi,dzi,dxj,dyj,dzj;
  if (nb_type[i]==0 || nb_type[j]==0) return;
  dx=dx_neighbor[i]+dx_neighbor[j];
  dy=dy_neighbor[i]+dy_neighbor[j];
  dz=dz_neighbor[i]+dz_neighbor[j];
  dxi=dx_neighbor[j];  dxj=dx_neighbor[i];
  dyi=dy_neighbor[j];  dyj=dy_neighbor[i];
  dzi=dz_neighbor[j];  dzj=dz_neighbor[i];
  Bi[0]+= -2*param.B*dx*dxi;   Bj[0]+= -2*param.B*dx*dxj;
  Bi[1]+= -2*param.B*dx*dyi;   Bj[1]+= -2*param.B*dx*dyj;
  Bi[2]+= -2*param.B*dx*dzi;   Bj[2]+= -2*param.B*dx*dzj;
  Bi[3]+= -2*param.B*dy*dxi;   Bj[3]+= -2*param.B*dy*dxj;
  Bi[4]+= -2*param.B*dy*dyi;   Bj[4]+= -2*param.B*dy*dyj;
  Bi[5]+= -2*param.B*dy*dzi;   Bj[5]+= -2*param.B*dy*dzj;
  Bi[6]+= -2*param.B*dz*dxi;   Bj[6]+= -2*param.B*dz*dxj;
  Bi[7]+= -2*param.B*dz*dyi;   Bj[7]+= -2*param.B*dz*dyj;
  Bi[8]+= -2*param.B*dz*dzi;   Bj[8]+= -2*param.B*dz*dzj;
}


void add_B_term3(float* Bi, float* Bj, int* nb_type, int i, int j) {
  int dx,dy,dz,dxi,dyi,dzi,dxj,dyj,dzj;
  if (nb_type[i]==0 || nb_type[j]==0) return;
  dx=dx_neighbor[j]-dx_neighbor[i];
  dy=dy_neighbor[j]-dy_neighbor[i];
  dz=dz_neighbor[j]-dz_neighbor[i];
  dxi=2*dx_neighbor[i]-dx_neighbor[j];  dxj=-dx_neighbor[i];
  dyi=2*dy_neighbor[i]-dy_neighbor[j];  dyj=-dy_neighbor[i];
  dzi=2*dz_neighbor[i]-dz_neighbor[j];  dzj=-dz_neighbor[i];
  Bi[0]+= 2*param.B*dx*dxi;   Bj[0]+= 2*param.B*dx*dxj;
  Bi[1]+= 2*param.B*dx*dyi;   Bj[1]+= 2*param.B*dx*dyj;
  Bi[2]+= 2*param.B*dx*dzi;   Bj[2]+= 2*param.B*dx*dzj;
  Bi[3]+= 2*param.B*dy*dxi;   Bj[3]+= 2*param.B*dy*dxj;
  Bi[4]+= 2*param.B*dy*dyi;   Bj[4]+= 2*param.B*dy*dyj;
  Bi[5]+= 2*param.B*dy*dzi;   Bj[5]+= 2*param.B*dy*dzj;
  Bi[6]+= 2*param.B*dz*dxi;   Bj[6]+= 2*param.B*dz*dxj;
  Bi[7]+= 2*param.B*dz*dyi;   Bj[7]+= 2*param.B*dz*dyj;
  Bi[8]+= 2*param.B*dz*dzi;   Bj[8]+= 2*param.B*dz*dzj;
}


void calc_AA_() {  // заполняем массивы AA_ и transform_array
  int ii,ip,iq,ir;
  unsigned short int config_;
  float a[4][4], d[4], v[4][4]; // метод Якоби: на входе - a, на выходе - d, v, nrot
  float a_[4][4], tr[4][4], d1[4], d2[4], dmax;
  
  for (config_=0; ; config_++) { // обнулим сначала массивы AA_ и transform_array
    for (ii=0; ii<6; ii++) AA_[config_][ii]=transform_array[config_][ii]=0; 
    if (config_==(Nconfig-1)) break;
  }
  
  for (config_=0; ; config_++) { // и заполним их только для "хороших" конфигураций
   if (good_config[config_]) {
    
    a[1][1]=AA[config_][0]; // считываем инф. из массива AA в матрицу a
    a[2][2]=AA[config_][1];
    a[3][3]=AA[config_][2];
    a[1][2]=a[2][1]=AA[config_][3];
    a[1][3]=a[3][1]=AA[config_][4];
    a[2][3]=a[3][2]=AA[config_][5];
    
    // Диагонализуем матрицу a методом Якоби
    jacobi(a,d,v);
    // Теперь в массиве d лежат собственные значения матрицы a, а в v - матрица преобразования.
    // (Матрица a при этом испортилась.)
    
    // Проверим собств. значения на невырожденность, потом положим в массив d1 их обратные величины, 
    // а в массив d2 - корни из обратных величин.
    dmax=0;
    for (ii=1; ii<=3; ii++)  if (d[ii]>dmax)  dmax=d[ii]; // нашли наибольшее собств. значение dmax
    if (dmax==0) {fprintf(stderr,"Error!!! matrix AA[%d] is zero\n",config_); exit(0);}
    for (ii=1; ii<=3; ii++) {
      if (d[ii]<-1e-3*dmax) {fprintf(stderr,"Error!!! negative eigenvalue of AA[%d]\n",config_);exit(0);}
      if (d[ii]<1e-3*dmax)  d1[ii]=d2[ii]=0;
      else {d1[ii]=1.0/d[ii]; d2[ii]=sqrt(1.0/d[ii]);}
    }
    
    // составляем матрицу a_=a^(-1) и tr=a^(-1/2)
    for (ip=1;ip<=3;ip++)
      for (iq=1;iq<=3;iq++)
        a_[ip][iq]=tr[ip][iq]=0;
    for (ip=1;ip<=3;ip++)
      for (iq=1;iq<=3;iq++)
        for (ir=1;ir<=3;ir++) {
          a_[ip][iq]+=v[ip][ir]*d1[ir]*v[iq][ir];  // a_=v*diag(d1)*vT - обратная матрица для a
          tr[ip][iq]+=v[ip][ir]*d2[ir]*v[iq][ir];  // tr=v*diag(d2)*vT - корень обратной матрицы для a
        }
    
    // собираем AA_ из a_,  transform_array из tr
    AA_[config_][0]=a_[1][1];
    AA_[config_][1]=a_[2][2];
    AA_[config_][2]=a_[3][3];
    AA_[config_][3]=a_[1][2];
    AA_[config_][4]=a_[1][3];
    AA_[config_][5]=a_[2][3];
    transform_array[config_][0]=tr[1][1];
    transform_array[config_][1]=tr[2][2];
    transform_array[config_][2]=tr[3][3];
    transform_array[config_][3]=tr[1][2];
    transform_array[config_][4]=tr[1][3];
    transform_array[config_][5]=tr[2][3];
    
   }
   if (config_==(Nconfig-1)) break;
  }
}


#define ROTATE(a,i,j,k,l) g=a[i][j]; h=a[k][l]; a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);
#define n 3
void jacobi(float a[n+1][n+1], float d[n+1], float v[n+1][n+1]) {
  // Диагонализация матрицы a методом Якоби, см. Numerical Recipes, параграф 11.1, стр. 467-468.
  // На входе - матрица a, на выходе - массив собств. значений d, матрица поворота v.
  int nrot,j,iq,ip,i;
  float tresh,theta,tau,t,sm,s,h,g,c;
  float b[n+1],z[n+1];
  int normal_exit; // для проверки нормального выхода из цикла
  
    for (ip=1;ip<=n;ip++) { //Initialize to the identity matrix.
      for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
      v[ip][ip]=1.0;
    }
    for (ip=1;ip<=n;ip++) { //Initialize b and d to the diagonal
      b[ip]=d[ip]=a[ip][ip];  // of a.
      z[ip]=0.0; // This vector will accumulate terms
      // of the form tapq as in equation (11.1.14).
    }
    nrot=0;
    normal_exit=0;
    for (i=1;i<=50;i++) {
      sm=0.0;
      for (ip=1;ip<=n-1;ip++) { // Sum off-diagonal elements.
        for (iq=ip+1;iq<=n;iq++)
          sm += fabs(a[ip][iq]);
      }
      if (sm == 0.0) { // The normal return, which relies
      // on quadratic convergence to machine underflow.
        normal_exit=1; break;
      }
      if (i < 4)
        tresh=0.2*sm/(n*n); // ...on the first three sweeps.
      else
        tresh=0.0; // ...thereafter.
      for (ip=1;ip<=n-1;ip++) {
        for (iq=ip+1;iq<=n;iq++) {
          g=100.0*fabs(a[ip][iq]);
          // After four sweeps, skip the rotation if the off-diagonal element is small.
          if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
            && (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
              a[ip][iq]=0.0;
          else if (fabs(a[ip][iq]) > tresh) {
            h=d[iq]-d[ip];
            if ((float)(fabs(h)+g) == (float)fabs(h))
              t=(a[ip][iq])/h; // t = 1/(2 theta)
            else {
              theta=0.5*h/(a[ip][iq]); // Equation (11.1.10).
              t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
              if (theta < 0.0) t = -t;
            }
            c=1.0/sqrt(1+t*t);
            s=t*c;
            tau=s/(1.0+c);
            h=t*a[ip][iq];
            z[ip] -= h;
            z[iq] += h;
            d[ip] -= h;
            d[iq] += h;
            a[ip][iq]=0.0;
            for (j=1;j<=ip-1;j++) { // Case of rotations 1 <= j < p.
              ROTATE(a,j,ip,j,iq)
            }
            for (j=ip+1;j<=iq-1;j++) { // Case of rotations p < j  < q.
              ROTATE(a,ip,j,j,iq)
            }
            for (j=iq+1;j<=n;j++) { // Case of rotations q < j <= n.
              ROTATE(a,ip,j,iq,j)
            }
            for (j=1;j<=n;j++) {
              ROTATE(v,j,ip,j,iq)
            }
            ++(nrot);
          }
        }
      }
      for (ip=1;ip<=n;ip++) {
        b[ip] += z[ip];
        d[ip]=b[ip]; // Update d with the sum of ta_pq ,
        z[ip]=0.0; // and reinitialize z.
      }
    }
    if (normal_exit==0) {fprintf(stderr,"Error!!! Too many iterations in jacobi\n"); exit(0);}
    // Теперь в массиве d лежат собственные значения матрицы a, а в матрице v - матрица преобразования.
    // Матрица a при этом портится.
}
#undef n
#undef ROTATE



