#include "mc_growth.h" 

void save_grid(char* filename, float radius) {
  int x,y,z,xg,yg,dx,dy;
  float R2,dz2,z_;
  FILE* f;
  float grid[Lx][Ly];
  R2=radius*radius;
  for (xg=0; xg<Lx; xg++)
    for (yg=0; yg<Ly; yg++)
      grid[xg][yg]=0;
  ZXY
   if (atoms(x,y,z).type>0)
    for (xg=0; xg<Lx; xg++)
     for (yg=0; yg<Ly; yg++) {
      dx=x-xg; dy=y-yg;
      dz2=R2-(dx*dx+dy*dy);
      if (dz2>0) {z_=z+sqrt(dz2); if (z_>grid[xg][yg]) grid[xg][yg]=z_;}
     }
  f=fopen(filename,"wt");
  if (f==NULL) {fprintf(stderr,"!!! cannot open file '%s' for writing\n",filename); return;}
  for (xg=0; xg<Lx; xg++) {
    for (yg=0; yg<Ly; yg++)
      fprintf(f,"%.2f  ",grid[xg][yg]);
    fprintf(f,"\n");
  }
  fclose(f);
}


void read_defects(const char* filename) {
  int x,y,z,i,n;
  FILE *f;
  if (n_defects==0) return;
  if (n_defects>N_defects) {fprintf(stderr,"!!! N_defects is too small\n"); exit(0);}
  f=fopen(filename,"rt");
  if (f==NULL) {fprintf(stderr,"!!! cannot open file '%s' for reading\n",filename); exit(0);}
  i=0;
  while (1) {
    n=fscanf(f,"%d%d%d",&x,&y,&z);
    if (n!=3) break;
    if (i>=N_defects) break;
    defect_x[i]=x; defect_y[i]=y; defect_z[i]=z;
    i++;
  }
  fclose(f);
  if (n_defects>0 && i<n_defects) 
    {fprintf(stderr,"!!! read_defects: there are less than %d entries in file %s\n",n_defects,filename); exit(0);}
  if (n_defects<0) n_defects=i;
}


void write_params_to_log(const char* filename) {
  int i;
  FILE *f;
  f=fopen(filename,"at");
  if (f==NULL) {fprintf(stderr,"!!! cannot open log file '%s'\n",filename); return;}
  fprintf(f,"\n\n---------- parameters: ----------\n");
  
  fprintf(f, "Lx=%d    Ly=%d    Lz=%d\n", Lx,Ly,Lz);
  fprintf(f, "E1[1]=%f K =%f eV\n", param.E1[1], param.E1[1]/11800);
  fprintf(f, "E2[1]=%f K =%f eV\n", param.E2[1], param.E2[1]/11800);
  fprintf(f, "E1[2]=%f K =%f eV\n", param.E1[2], param.E1[2]/11800);
  fprintf(f, "E2[2]=%f K =%f eV\n", param.E2[2], param.E2[2]/11800);
  fprintf(f, "E1[3]=%f K =%f eV\n", param.E1[3], param.E1[3]/11800);
  fprintf(f, "E2[3]=%f K =%f eV\n", param.E2[3], param.E2[3]/11800);
  fprintf(f, "p0=%.2e s^(-1)\n", param.p0);
  fprintf(f, "T=%.2f K\n", param.T);
  fprintf(f, "A=%f K/d^2    B=%f K/d^2\n", param.A,param.B);
  fprintf(f, "eps[1]=%f\n", param.eps[1]);
  fprintf(f, "eps[2]=%f\n", param.eps[2]);
  fprintf(f, "eps[3]=%f\n", param.eps[3]);
  fprintf(f, "dep_type=%d\n", param.dep_type);
  fprintf(f, "p_deposition=%f s^(-1) =%3.2f ML/s\n", param.p_deposition, param.p_deposition*(8.0/Lx/Ly));
  fprintf(f, "moves_percent=%.2f%%\n", param.moves_percent);
  fprintf(f, "time_for_moves=%.2e s\n", param.time_for_moves);
  fprintf(f, "defect_F=%f K/d\n", param.defect_F);
  fprintf(f, "z_surface=%d\n", param.z_surface);
  fprintf(f, "random_seed=%d\n", param.random_seed);
  fprintf(f, "load: %s\n", param.load);
  if(param.show_short!=0) fprintf(f, "show_short every %.2f ML\n", param.show_short);
  else fprintf(f,"no show_short\n");
  if(param.show_long!=0) fprintf(f, "show_long every %.2f ML\n", param.show_long);
  else fprintf(f,"no show_long\n");
  if(param.show_sum!=0) fprintf(f, "show_sum every %.2f ML\n", param.show_sum);
  else fprintf(f,"no show_sum\n");	
  fprintf(f, "experiment_ML=%.2f ML\n", param.experiment_ML);
  
  if (n_defects<=0)  fprintf(f, "no defects\n");
  else {
    fprintf(f, "%d defects\n", n_defects);
  }
  
  fprintf(f,"---------------------------------\n\n");
  fclose(f);
}


void write_to_log(const char* filename) {
  FILE *f;
  f=fopen(filename,"at");
  if (f==NULL) {fprintf(stderr,"!!! cannot open log file '%s'\n",filename); return;}
  fprintf(f, "%.2f ML	t: %.2f s	rate: %.3f ML/s (type %d, T %.0f)	deposited: %6d	jumps: %10d	(bad %.0f%%, back %.0f%%)	moves: %12.0f\n", 
    current.n_deposited/(Lx*Ly/8.),current.t, param.p_deposition*(8.0/Lx/Ly), param.dep_type, param.T,
    current.n_deposited, current.n_jumps, (current.n_bad_jumps*100.0)/current.n_jumps, (current.n_jumps_back*100.0)/current.n_jumps,current.n_moves);
  fclose(f);
}



void read_control(const char* filename) {
  float ML,v,T; int type,n;
  FILE *f;
  f=fopen(filename,"rt");
  if (f==NULL) {fprintf(stderr,"!!! cannot open file '%s' for reading\n",filename); return;}
  while (1) {
    n=fscanf(f,"%f%f%d%f",&ML,&v,&type,&T);
    if (n!=4) break;
    if (ML>current.n_deposited/(Lx*Ly/8.)) break;
    param.p_deposition=v*Lx*Ly/8.0;
    param.dep_type=type;
    param.T=T;
  }
  fclose(f);
}


void read_parameters(const char* filename) {
  int n;
  char buf[110], buf2[110];
  FILE *f;
  f=fopen(filename,"rt");
  if (f==NULL) {fprintf(stderr,"!!! cannot open file '%s' for reading\n",filename); exit(0);}
  while (1) {
    n=fscanf(f,"%99s%99s",buf,buf2); if (n!=2) break;
    
    if (!strcmp(buf,"Lx")) {Lx=atoi(buf2); continue;}
    if (!strcmp(buf,"Ly")) {Ly=atoi(buf2); continue;}
    if (!strcmp(buf,"Lz")) {Lz=atoi(buf2); continue;}
    
    if (!strcmp(buf,"load")) {strcpy(param.load,buf2); continue;}
    if (!strcmp(buf,"show_short")) {param.show_short=atof(buf2); continue;}
    if (!strcmp(buf,"show_long")) {param.show_long=atof(buf2); continue;}
    if (!strcmp(buf,"show_sum")) {param.show_sum=atof(buf2); continue;}
    if (!strcmp(buf,"E1[1]")) {param.E1[1]=atof(buf2); continue;}
    if (!strcmp(buf,"E1[1](eV)")) {param.E1[1]=atof(buf2)*11800; continue;}
    if (!strcmp(buf,"E2[1]")) {param.E2[1]=atof(buf2); continue;}
    if (!strcmp(buf,"E2[1](eV)")) {param.E2[1]=atof(buf2)*11800; continue;}	
	if (!strcmp(buf,"E1[2]")) {param.E1[2]=atof(buf2); continue;}
    if (!strcmp(buf,"E1[2](eV)")) {param.E1[2]=atof(buf2)*11800; continue;}
    if (!strcmp(buf,"E2[2]")) {param.E2[2]=atof(buf2); continue;}
    if (!strcmp(buf,"E2[2](eV)")) {param.E2[2]=atof(buf2)*11800; continue;}
	if (!strcmp(buf,"E1[3]")) {param.E1[3]=atof(buf2); continue;}
    if (!strcmp(buf,"E1[3](eV)")) {param.E1[3]=atof(buf2)*11800; continue;}
    if (!strcmp(buf,"E2[3]")) {param.E2[3]=atof(buf2); continue;}
    if (!strcmp(buf,"E2[3](eV)")) {param.E2[3]=atof(buf2)*11800; continue;}
    if (!strcmp(buf,"p0")) {param.p0=atof(buf2); continue;}
    if (!strcmp(buf,"p_deposition")) {param.p_deposition=atof(buf2); continue;}
    if (!strcmp(buf,"p_deposition(ML/s)")) {param.p_deposition=atof(buf2)*Lx*Ly/8.0; continue;}
    if (!strcmp(buf,"T")) {param.T=atof(buf2); continue;}
    if (!strcmp(buf,"A")) {param.A=atof(buf2); continue;}
    if (!strcmp(buf,"B")) {param.B=atof(buf2); continue;}
    if (!strcmp(buf,"eps[1]")) {param.eps[1]=atof(buf2); continue;}	
	if (!strcmp(buf,"eps[2]")) {param.eps[2]=atof(buf2); continue;}	
	if (!strcmp(buf,"eps[3]")) {param.eps[3]=atof(buf2); continue;}	
    if (!strcmp(buf,"defect_F")) {param.defect_F=atof(buf2); continue;}
    if (!strcmp(buf,"moves_percent")) {param.moves_percent=atof(buf2); continue;}
	if (!strcmp(buf,"time_for_moves")) {param.time_for_moves=atof(buf2)*1e-6; continue;}
    if (!strcmp(buf,"dML_control")) {param.dML_control=atof(buf2); continue;}
    if (!strcmp(buf,"experiment_ML")) {param.experiment_ML=atof(buf2); continue;}
    
    if (!strcmp(buf,"random_seed")) {param.random_seed=atoi(buf2); continue;}
    if (!strcmp(buf,"dep_type")) {param.dep_type=atoi(buf2); continue;}
    if (!strcmp(buf,"z_surface")) {param.z_surface=atoi(buf2); continue;}
    if (!strcmp(buf,"n_defects")) {n_defects=atoi(buf2); continue;}
    if (!strcmp(buf,"z_cap")) {param.z_cap=atoi(buf2); continue;}
    if (!strcmp(buf,"Edef_pa_N")) {param.Edef_pa_N=atoi(buf2); continue;}
    
    fprintf(stderr,"!!! unknown parameter '%s'\n",buf); exit(0);
  }
  fclose(f);
  check_Lxyz();
}



void show_me(char* filename, char key){

    int x,y,z,k,n;
    long n_atoms;

    FILE* f;

    f=fopen(filename,"wt");
    if (f==NULL) {fprintf(stderr,"!!! cannot open file '%s' for writing\n",filename); return;}

    n_atoms=0;
    ZYX_a {
        if (atoms(x,y,z).type!=0/*||defects[i][j][k]==1*/){
          ++n_atoms;
        }
    }
    n_atoms+=n_defects;
	  
    fprintf(f,"%ld\n",n_atoms);
	fprintf(f,"%.2f	%.2f	%d	%d	%d	%d\n",current.n_deposited/(Lx*Ly/8.),current.t,current.n_deposited, current.n_jumps,
												current.n_bad_jumps,current.n_jumps_back);
	if(key == 0) {	//show_short
		ZYX_a {
			if(atoms(x,y,z).type!=0) {
				fprintf(f,"%d	",atoms(x,y,z).type);
				fprintf(f,"%3d %3d %3d\n",x,y,z);
			}
		}
	}
	
	else if (key == 1) {	//show_long
		ZYX_a {
			if(atoms(x,y,z).type!=0) {
				fprintf(f,"%d	",atoms(x,y,z).type);
				fprintf(f,"%3d %.4f	%3d %.4f	%3d %.4f\n",x,atoms(x,y,z).a.x,y,atoms(x,y,z).a.y,z,atoms(x,y,z).a.z);
			}
		}
	}
	
	else if (key == 2) {	//show_sum
		ZYX_a {
			if(atoms(x,y,z).type!=0) {
				fprintf(f,"%d	",atoms(x,y,z).type);
				fprintf(f,"%.4f	%.4f	%.4f\n", x + atoms(x,y,z).a.x, y + atoms(x,y,z).a.y, z + atoms(x,y,z).a.z);
			}
		}
	}
		
    if(n_defects!=0){
        for(n=0;n<n_defects;++n){
        //for(k=0;k<Lz;++k)
        //  for(i=0;i<Lx;++i)	
        //    for(j=0;j<Ly;++j)
        //	  if(i==defect_x[n]&&j==defect_y[n]&&k==defect_z[n]){
            fprintf(f,"4 ");
		    fprintf(f,"%d 0 %d 0 %d 0\n",defect_x[n],defect_y[n],defect_z[n]);
//		    fprintf(f,"%d %d %d\n",defect_x[n],defect_y[n],defect_z[n]);
        //	    }

        }
    }
		
    /*
    for(k=0;k<Lz;++k)
      for(i=0;i<Lx;++i)	
        for(j=0;j<Ly;++j)
          if(defects[i][j][k]==1){
    fprintf(f,"3 ");
    fprintf(f,"%d 0 %d 0 %d 0\n",i,j,k);
           }
    */		
    fclose(f);


}

void show_me_Edef(char* filename){

  int x,y,z,k;
  long n_atoms;

  FILE* f;
  //char filename[100];
  
  //sprintf(filename,"atoms.xyz");
  f=fopen(filename,"wt");
  if ( unlikely(f==NULL) ) {
    fprintf(stderr,"!!! cannot open file '%s' for writing\n",filename); 
    return;
  }

    n_atoms=0;
    ZYX_a {
        if (atoms(x,y,z).type!=0/*||defects[i][j][k]==1*/){
          ++n_atoms;
        }
    }
    n_atoms+=n_defects;
	  
    fprintf(f,"%ld\n",n_atoms);
	fprintf(f,"%.2f	%.2f	%d	%d	%d	%d\n",current.n_deposited/(Lx*Ly/8.),current.t,current.n_deposited, current.n_jumps,
												current.n_bad_jumps,current.n_jumps_back);
	
	    ZYX_a {
      if(atoms(x,y,z).type!=0){
        fprintf(f,"%d	",atoms(x,y,z).type);
		fprintf(f,"%3d %.4f	%3d %.4f	%3d %.4f\n",x,atoms(x,y,z).a.x,y,atoms(x,y,z).a.y,z,atoms(x,y,z).a.z/*,atoms(x,y,z).Edef_pa/param.Edef_pa_N*/);
		}
	}
  
  fclose(f);
}
