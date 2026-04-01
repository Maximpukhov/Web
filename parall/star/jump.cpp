#include "mc_growth.h" 


#define N_spisok 1000
int spisok[N_spisok][3];
int n_spisok;


void v_spisok(int x, int y, int z);
void delete_jump(int x, int y, int z, int dir);


void calc_jump_probability(int x, int y, int z) {
  int dir;
  float E_a;
  if (n_jumps(x, y, z) ==0 ) {
    jump_probability(x,y,z)=0;
    return;
  }
  if( unlikely( (atoms(x,y,z).defect_f.x != 0) || 
                (atoms(x,y,z).defect_f.y != 0) || 
                (atoms(x,y,z).defect_f.z != 0) ) ){
	
	char filename[100];
	FILE *f;
	sprintf(filename,"log.txt");
	f=fopen(filename,"at");
  fprintf(f,"	atom (%d,%d,%d) near defect is on the surface: (%hhd, %hhd, %hhd)\n",x,y,z,
        atoms(x,y,z).defect_f.x, atoms(x,y,z).defect_f.y, atoms(x,y,z).defect_f.z);
	fclose(f);
	erase_defect(x,y,z);
  }
/*
	char filename[100];
	FILE *f;
	sprintf(filename,"log1.txt");
	f=fopen(filename,"at");
*/
  E_a=0;
  neighbors_t nbs;
  atoms.neighbors(x,y,z,nbs);
  for (dir=0; dir<dir_number; dir++) {
    if(atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).type!=0){
      if(dir<4){
	  E_a+=0.5*(param.E1[atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).type]+param.E1[atoms(x,y,z).type]);
//	  fprintf(f,"dir=%d  E1[%d]=%7.6f+E1[%d]=%7.6f\n",dir,atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).type,
//	  param.E1[atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).type],atoms(x,y,z).type,param.E1[atoms(x,y,z).type]);
	  }
	  else{
	  E_a+=0.5*(param.E2[atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).type]+param.E2[atoms(x,y,z).type]);
//	  fprintf(f,"dir=%d  E2[%d]=%7.6f+E2[%d]=%7.6f\n",dir,atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).type,
//	  param.E2[atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).type],atoms(x,y,z).type,param.E2[atoms(x,y,z).type]);
	  }
	}
	
  }
  
  E_a-=atoms(x,y,z).Edef;
  
//  fprintf(f,"1 E_a=%7.6f\n",E_a);
/*
  E_a1 = param.E1[1]*n1_config[atoms(x,y,z).config] 
      + param.E2[1]*n2_config[atoms(x,y,z).config] 
      - atoms(x,y,z).Edef;
*/
//  fprintf(f,"2 E_a=%7.6f\n",E_a);	  
//	fclose(f);	  
  if (E_a<0) E_a=0;
  jump_probability(x,y,z) = param.p0*exp(-E_a/param.T)*n_jumps(x, y, z);
}


void set_jump_info(int x, int y, int z) { // устанавливает jumps и n_jumps для данного атома
  int dir,x2,y2,z2,n_jumps_;
  int j_array[dir_number];
  if ( atoms(x,y,z).type == 0 || n1_config[atoms(x,y,z).config]==4 || 
      z<2 || z>=Lz-2) {
    n_jumps(x,y,z) = 0;
    jumps(x,y,z) = 0; 
    return;
  }
  n_jumps_=0;
  for (dir=0; dir<dir_number; dir++) j_array[dir]=0;

  neighbors_t nbs;
  atoms.neighbors(x,y,z,nbs);
  
  for (dir=0; dir<dir_number; dir++) {
    if( atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).type==0 && 
      n1_config[atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).config]>(dir<4)){
      j_array[dir]=1; n_jumps_++;
    }
  }
  n_jumps(x,y,z) = n_jumps_;
  jumps(x,y,z) = massiv_to_config(j_array);
}


int jump(int x, int y, int z, int dir, int* x2_, int* y2_, int* z2_) {
  int xd,yd,zd, x2,y2,z2,dir2, x3,y3,z3,dir3, n,n_spisok2, bad_jump, atom_erase;

  atoms.one_neighbor(x,y,z,dir,xd,yd,zd);
  *x2_=xd; *y2_=yd; *z2_=zd;
	
  if ( unlikely(atoms(x,y,z).type==0) ) {
    fprintf(stderr,"Error!!! jump from empty site (%d,%d,%d): %d ->\n",x,y,z,atoms(x,y,z).type); 
    exit(0);
  }
	
  if ( unlikely(atoms(xd,yd,zd).type!=0) ) {
    fprintf(stderr,"Error!!! jump to occupied site (%d,%d,%d): %d -> (%d,%d,%d): %d\n",x,y,z,atoms(x,y,z).type,xd,yd,zd,atoms(xd,yd,zd).type); 
    exit(0);
  }
 
  if (z<4 || zd<4 || z>=Lz-4 || zd>=Lz-4) {
    delete_jump(x,y,z,dir);
    calc_jump_probability(x,y,z);
    return 0;
  }

  atoms(xd,yd,zd).type = atoms(x,y,z).type;
  atoms(x,y,z).type=0;
  
  n_spisok=0;
  v_spisok(x,y,z);
  v_spisok(xd,yd,zd);

  neighbors_t nbs, nbsd;
  atoms.neighbors(x,y,z,nbs);
  atoms.neighbors(xd,yd,zd,nbsd);

  for (dir2=0; dir2<dir_number; dir2++) {
  /*
    atoms.one_neighbor(x,y,z,dir2,x2,y2,z2);
    v_spisok(x2,y2,z2);
    if( x2 != nbs.x[dir2] || y2 != nbs.y[dir2] || z2 != nbs.z[dir2])
      abort();
      
    atoms.one_neighbor(xd,yd,zd,dir2,x2,y2,z2);
    v_spisok(x2,y2,z2);

    if( x2 != nbsd.x[dir2] || y2 != nbsd.y[dir2] || z2 != nbsd.z[dir2])
      abort();
    */
    v_spisok(nbs.x[dir2],nbs.y[dir2],nbs.z[dir2]);
    v_spisok(nbsd.x[dir2],nbsd.y[dir2],nbsd.z[dir2]);
  }
  bad_jump=0;
  for (n=0; n<n_spisok; n++) {
    x2=spisok[n][0]; y2=spisok[n][1]; z2=spisok[n][2];
    set_config(x2,y2,z2);
    if (atoms(x2,y2,z2).type > 0 && !good_config[atoms(x2,y2,z2).config])
      {bad_jump=1; break;}
  }

  if (bad_jump) {
    atoms(x,y,z).type = atoms(xd,yd,zd).type;
    atoms(xd,yd,zd).type = 0;
    for (n=0; n<n_spisok; n++) {
      x2=spisok[n][0]; y2=spisok[n][1]; z2=spisok[n][2];
      set_config(x2,y2,z2);
      if (atoms(x2,y2,z2).type>0 && !good_config[atoms(x2,y2,z2).config]){
        fprintf(stderr,"Error!!! in jump()\n");
        exit(0);
      }
	  axyz(x2,y2,z2);
      spisok_flag(x2,y2,z2)=0;
    }
    delete_jump(x,y,z,dir);
    calc_jump_probability(x,y,z);
    return 0;
  }

	atom_erase = 0;
  if(zd < 6) {
		atoms(xd,yd,zd).type = 0;
		atom_erase = 1;
	}

  for (n=0; n<n_spisok; n++) {
    x2=spisok[n][0]; y2=spisok[n][1]; z2=spisok[n][2];
	if(atoms(x2,y2,z2).type==3)
	  atoms(x2,y2,z2).type=1;
    calc_B0(x2,y2,z2);
//	axyz(x2,y2,z2);
  }
	
  axyz(xd,yd,zd); // выставим ax,ay,az для прыгнувшего атома
  n_spisok2=n_spisok;
  for (n=0; n<n_spisok2; n++) {
    x2=spisok[n][0]; y2=spisok[n][1]; z2=spisok[n][2];
    atoms.neighbors(x2, y2, z2, nbs);
    for (dir3=0; dir3<dir_number; dir3++) {
      v_spisok(nbs.x[dir3],nbs.y[dir3],nbs.z[dir3]);
    }
  }
  
  for (n=0; n<n_spisok; n++) {
    x3=spisok[n][0]; y3=spisok[n][1]; z3=spisok[n][2];
    set_jump_info(x3,y3,z3);
    if (n_jumps(x3, y3, z3)==0) calc_jump_probability(x3,y3,z3);
    v_ochered_Edef(x3,y3,z3);
    spisok_flag(x3,y3,z3)=0;
  }

		if(atom_erase == 1) move_atom_in_spisok(x,y,z,x,y,z);//remove_atom_from_spisok(x,y,z);
		else move_atom_in_spisok(x,y,z,xd,yd,zd);

  return 1;
}


int deposition(int type_of_new_atom, int* x_, int* y_, int* z_) {
  int x,y,z,i,f,dir, x2,y2,z2,dir2, x3,y3,z3,dir3, n,n_spisok2;
  if (type_of_new_atom<=0) {fprintf(stderr,"Error!!! type_of_new_atom<=0\n"); exit(0);}
  
  z = Lz-1;
  x = (z%2) + 2*random_(Lx/2);
  y = (z%2) + 2*((z/2+x/2)%2) + 4*random_(Ly/4);
  if( unlikely(x<0 || x>=Lx || y<0 || y>=Ly) ){
    fprintf(stderr,"Error!!! in deposition()\n");
    exit(0);
  }
  
  if( unlikely(atoms(x,y,z).type!=0) ){
    fprintf(stderr,"Error!!! in deposition()\n");
    exit(0);
  }
  
  f=0;
  for (i=0; i<1000; i++) {
    if (z%2==0) 
      dir = (random_(2)==0) ? 1 : 2;
    else
      dir = (random_(2)==0) ? 0 : 3;
    atoms.one_neighbor(x,y,z,dir,x2,y2,z2);
    
    if (atoms(x2,y2,z2).type!=0) {
      dir=random_(4);
      atoms.one_neighbor(x,y,z,dir,x2,y2,z2);
    }
		
		if(z < 6) return 1;

    if (atoms(x2,y2,z2).type == 0)
      x=x2, y=y2, z=z2;
    
    if (atoms(x,y,z).type == 0 && good_config[atoms(x,y,z).config]) {
      f=1; break;
    }
  }
  if( unlikely(f==0) ) {
    fprintf(stderr,"Error!!! I can not deposit an atom\n"); 
    exit(0);
  }
  
  if ( unlikely(z<4 || z>=Lz-4) ){
    fprintf(stderr,"Error!!! z<4 || z>=Lz-4 for deposited atom\n"); 
    exit(0);
  }
  
  *x_=x; *y_=y; *z_=z;
  
  atoms(x,y,z).type=type_of_new_atom;
  
  n_spisok=0;
  v_spisok(x,y,z);
  
  neighbors_t nbs;
  atoms.neighbors(x,y,z,nbs);

  for (dir2=0; dir2<dir_number; dir2++) {
    v_spisok(nbs.x[dir2],nbs.y[dir2],nbs.z[dir2]);
    set_config(nbs.x[dir2],nbs.y[dir2],nbs.z[dir2]);
    if ( unlikely(atoms(nbs.x[dir2],nbs.y[dir2],nbs.z[dir2]).type > 0 && 
      !good_config[atoms(nbs.x[dir2],nbs.y[dir2],nbs.z[dir2]).config]) ){
      fprintf(stderr,"Error!!! bad config of atom (%d,%d,%d)\n",
          nbs.x[dir2],nbs.y[dir2],nbs.z[dir2]);
      abort();
    }
  }
  
  for (n=0; n<n_spisok; n++) {
    x2=spisok[n][0]; y2=spisok[n][1]; z2=spisok[n][2];
    calc_B0(x2,y2,z2);
  }
  
  axyz(x,y,z); // выставим ax,ay,az для нового атома
  
  n_spisok2=n_spisok;
  for (n=0; n<n_spisok2; n++) {
    x2=spisok[n][0]; y2=spisok[n][1]; z2=spisok[n][2];
    atoms.neighbors(x2,y2,z2,nbs);
    for (dir3=0; dir3<dir_number; dir3++) {
      v_spisok(nbs.x[dir3],nbs.y[dir3],nbs.z[dir3]);
    }
  }
  
  for (n=0; n<n_spisok; n++) {
    x3=spisok[n][0]; y3=spisok[n][1]; z3=spisok[n][2];
    set_jump_info(x3,y3,z3);
    if (n_jumps(x3, y3, z3) == 0) calc_jump_probability(x3,y3,z3);
    v_ochered_Edef(x3,y3,z3);
    spisok_flag(x3,y3,z3)=0;
  }
  
  add_atom_to_spisok(x,y,z);

  return 1;
}


void v_spisok(int x, int y, int z) {
  if (spisok_flag(x,y,z)) return;
  if ( unlikely(n_spisok>=N_spisok) ) {
    fprintf(stderr,"Error!!! N_spisok is too small\n"); 
    exit(0);
  }
  
  spisok[n_spisok][0]=x;
  spisok[n_spisok][1]=y;
  spisok[n_spisok][2]=z;
  n_spisok++;
  spisok_flag(x,y,z) = 1;
}


void delete_jump(int x, int y, int z, int dir) {
  int d, n_jumps_;
  int j_array[dir_number];
  fill_nb_type(j_array,jumps(x,y,z));
  j_array[dir] = 0;
  jumps(x,y,z) = massiv_to_config(j_array);
  n_jumps_=0;
  for (d=0; d<dir_number; d++)  n_jumps_+=j_array[d];
  n_jumps(x, y, z) = n_jumps_;
}

