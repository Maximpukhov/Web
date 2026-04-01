#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define Lx_max 200
#define Ly_max 200
#define Lz_max 200

int random_(int);

int random_(int n) { // случайное число от 0 до n-1
  return rand() % n;
}

int main(){

int x0[Lx_max/4],y0[Ly_max/4],z0;
int x1[Lx_max/2],y1[Ly_max/2],x_1[Lx_max/2],y_1[Ly_max/2];
int i,j;
int nclasters;

nclasters=0;

for(i=0;i<Lx_max/4;++i){
  x0[i]=0;
  y0[i]=0;
  }

for(i=0;i<Lx_max/2;++i){
  x1[i]=0;
  y1[i]=0;
  x_1[i]=0;
  y_1[i]=0;
  }
  
int random_seed;
//random_seed=19;

printf("\nrandom_seed: ");
scanf("%d",&random_seed);
srand(random_seed);

//for(i=0;i<10;++i)
//printf("\n%d",random_(200));

//z0=34;nclasters=1;


printf("\nz0: ");
scanf("%d",&z0);
if(z0>=Lz_max||z0<4){ printf("\nerror z0=%d",z0); exit(0);}


printf("\nclasters of defects (%d max): ",Lx_max/4);
scanf("%d",&nclasters);
if(nclasters>Lx_max/4){ printf("\nerror nclasters=%d (%d max)",nclasters,Lx_max/4); exit(0);}

printf("\nz0=%d\nclasters of defects: %d",z0,nclasters);
//z0=4*nz+2
//x0=4*nx
//y0=4*ny


// ПЕРЕПИСАТЬ!!!!
for(i=0;i<nclasters;++i){
  x0[i]=4*random_(Lx_max/4);
  if(z0%4==0){ x0[i]+=2; if(x0[i]>Lx_max) x0[i]-=Lx_max;
    }
  y0[i]=4*random_(Ly_max/4);
//  printf("\nx0[%d]=%d, y0[%d]=%d",i,x0[i],i,y0[i]);
  if(z0%2==0){
    x1[2*i]=x0[i]-1;	if(x1[2*i]<0) x1[2*i]+=Lx_max;
	y1[2*i]=y0[i]-1;	if(y1[2*i]<0) y1[2*i]+=Ly_max;
	x1[2*i+1]=x0[i]+1;	if(x1[2*i+1]>Lx_max) x1[2*i+1]-=Lx_max;
	y1[2*i+1]=y0[i]+1;	if(y1[2*i+1]>Ly_max) y1[2*i+1]-=Ly_max;
	x_1[2*i]=x0[i]-1;	if(x_1[2*i]<0) x_1[2*i]+=Lx_max;
	y_1[2*i]=y0[i]+1;	if(y_1[2*i]>Ly_max) y_1[2*i]-=Ly_max;
	x_1[2*i+1]=x0[i]+1;	if(x_1[2*i+1]>Lx_max) x_1[2*i+1]-=Lx_max;
	y_1[2*i+1]=y0[i]-1;	if(y_1[2*i+1]<0) y_1[2*i+1]+=Ly_max;
    }
  else{
    x1[2*i]=x0[i]+1;	if(x1[2*i]>Lx_max) x1[2*i]-=Lx_max;
	y1[2*i]=y0[i]-1;	if(y1[2*i]<0) y1[2*i]+=Ly_max;
	x1[2*i+1]=x0[i]-1;	if(x1[2*i+1]<0) x1[2*i+1]+=Lx_max;
	y1[2*i+1]=y0[i]+1;	if(y1[2*i+1]>Ly_max) y1[2*i+1]-=Ly_max;
	x_1[2*i]=x0[i]-1;	if(x_1[2*i]<0) x_1[2*i]+=Lx_max;
	y_1[2*i]=y0[i]-1;	if(y_1[2*i]<0) y_1[2*i]+=Ly_max;
	x_1[2*i+1]=x0[i]+1;	if(x_1[2*i+1]>Lx_max) x_1[2*i+1]-=Lx_max;
	y_1[2*i+1]=y0[i]+1;	if(y_1[2*i+1]<0) y_1[2*i+1]+=Ly_max;
	}
  }

for(i=0;i<nclasters;++i){
  printf("\n[%d][%d][%d]	[%d][%d][%d]	[%d][%d][%d]	[%d][%d][%d]	[%d][%d][%d]",
  x0[i],y0[i],z0,  x1[2*i],y1[2*i],z0+1,  x1[2*i+1],y1[2*i+1],z0+1,  x_1[2*i],y_1[2*i],z0-1,  x_1[2*i+1],y_1[2*i+1],z0-1);
  }

  FILE *f;
  
  char filename[100];
  sprintf(filename,"defects.txt");
  f=fopen(filename,"wt");
  
for(i=0;i<nclasters;++i){
  fprintf(f,"%d %d %d\n%d %d %d\n%d %d %d\n%d %d %d\n%d %d %d\n",
  x0[i],y0[i],z0,  x1[2*i],y1[2*i],z0+1,  x1[2*i+1],y1[2*i+1],z0+1,  x_1[2*i],y_1[2*i],z0-1,  x_1[2*i+1],y_1[2*i+1],z0-1);
  }
  
  fclose(f);

  
  
return 0;
}