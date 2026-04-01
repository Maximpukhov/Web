#include "mc_growth.h" 
#include "optimization.h"

void set_defect_f(void);



int random_(int n) { // случайное число от 0 до n-1
  return rand() % n;  // it's not very good!!!!!!!!!!!!!!!!!!!!!
}


double rand01(void) {  // случайное число между 0 и 1
  return ((double)rand()/(RAND_MAX+1.)+(double)rand())/(RAND_MAX+1.);
}

void fill_nb_type(int* nb_type, unsigned short int config_) {
  unsigned short int conf;
  int i;
  conf=config_;
  for (i=0; i<dir_number; i++) {
    nb_type[i]=((conf&0x8000)!=0);
    conf<<=1;
  }
}


unsigned short int massiv_to_config(int* nb_type) {
  unsigned short int conf;
  int i;
  conf=0;
  for (i=0; i<dir_number; i++) {
    conf<<=1;
    conf+=(nb_type[i]!=0);
  }
  return conf;
}


void fill_n1_n2_and_good_config() {
  int i,j,n1_,n2_,n2__;
  unsigned short int config_;
  int nb_type[dir_number];
  for (config_=0; ; config_++) {
    fill_nb_type(nb_type,config_);
    n1_=n2_=n2__=0;
    for (i=0; i<dir_number; i++)
      if (pr_neighbor[i]<0 && nb_type[i]>0)
        n1_++;
    for (i=0; i<dir_number; i++)
      if (pr_neighbor[i]>=0) {
        j=pr_neighbor[i];
        if (nb_type[i]>0 && nb_type[j]>0)
          n2_++;
        if (nb_type[i]>0 /*&& nb_type[j]>0*/)
          n2__++;
      }
    n1_config[config_]=n1_;
    n2_config[config_]=n2__;
    good_config[config_]=(n2_>=2);
    if (config_==(Nconfig-1)) break;
  }
}


void v_ochered_Edef(int x,int y,int z) {
  if (n_jumps(x, y, z) ==0 || in_ochered_Edef(x,y,z)) return;
  struct coord c = { x, y, z };
  ochered_Edef.push_back(c);
  in_ochered_Edef(x,y,z) = 1;
}


void update_Edef() {
  int i,x,y,z;
  i=10;
  std::list<coord>::iterator it = ochered_Edef.begin();
  for (; ochered_Edef.size() != 0;) {
    struct coord c = *(ochered_Edef.begin());
    x=c.x;
    y=c.y;
    z=c.z;
    ochered_Edef.pop_front();
    calc_Edef(x,y,z);
    calc_jump_probability(x,y,z);
    in_ochered_Edef(x,y,z)=0;
  }
}


void create_initial_structure(int *initial_layers) 
{
  int x,y,z;

  ZXY{
    memset(&atoms(x,y,z), 0, sizeof(atom_t));
    jump_probability(x,y,z) = n_jumps(x, y, z) = 0;
    spisok_flag(x,y,z) = in_ochered_Edef(x,y,z) = jumps(x,y,z) = 0;
    number_in_spisok_atomov(x,y,z) = -1;
  }

  ochered_Edef.clear();
  spisok_atomov.clear();
  
  for (z=0; z<Lz; z++)
    ALL_X ALL_Y  atoms(x,y,z).type = initial_layers[z]; // расставим атомы в решётке

// some geometric structures
//***********************************************************
/*
// ball in the center of structure
  int br;
  br=40;	// radius
  for (z=0; z<Lz; z++)	
	ALL_X ALL_Y
  if((x-Lx/2.)*(x-Lx/2.)+(y-Ly/2.)*(y-Ly/2.)+(z-param.z_surface/2.)*(z-param.z_surface/2.)<br*br)
    atoms(x,y,z).type=2;
*/	
//***********************************************************
/*
// cube
  int cr;
  cr=20;	// side
  for (z=0; z<Lz; z++)	
  ALL_X ALL_Y
  if(x>(Lx/2-cr/2)&&x<(Lx/2+cr/2)&&y>(Ly/2-cr/2)&&y<(Ly/2+cr/2)&&z>(param.z_surface/2-cr/2)&&z<(param.z_surface/2+cr/2))
    atoms(x,y,z).type=2;
*/
//***********************************************************
/*
// inverted prmd
// 110
int a;
float alpha;

a = 50;
alpha = 45;

for (z=0; z<Lz; z++)	
	ALL_X ALL_Y
		if(z<=param.z_surface&&
//		   +(x-Lx/2)+(y-Ly/2)<a*sin(M_PI/4)&&
//		   -(x-Lx/2)+(y-Ly/2)<a*sin(M_PI/4)&&
//		   -(x-Lx/2)-(y-Ly/2)<a*sin(M_PI/4)&&
//		   +(x-Lx/2)-(y-Ly/2)<a*sin(M_PI/4))
		   
		   -(x-Lx/2)*a*a/2*sin(M_PI*alpha/180)*tan(M_PI*alpha/180)-(y-Ly/2)*a*a/2*sin(M_PI*alpha/180)*tan(M_PI*alpha/180)+z*a*a*pow(sin(M_PI*alpha/180),2)+a*a*a/2*pow(sin(M_PI*alpha/180),2)*tan(M_PI*alpha/180)-param.z_surface*a*a*pow(sin(M_PI*alpha/180),2)>0)
		   atoms(x,y,z).type = 2;
*/
//***********************************************************
/*
//tablet
float r, R, h, alpha;
int i;
r = 60;	// tablet top radus
R = 90;	// tablet bottom radius
h = 9;		// tablet height
alpha = 14;	// angle of side
// 200x200: tr=80	0,-2,-4,-6,-8,-10,-12,-13	
// 300x300: tr=120	
// 400x400: tr=160	0,-2,-6,-8,-12,-14,-17,-20

  for (z=0; z<Lz; z++)	
	ALL_X ALL_Y {
		for(i = 0; i < h; i++) {
			// r(alpha) | R, h - constant
//			if(z>param.z_surface&&z<=param.z_surface+1+i&&(x-Lx/2)*(x-Lx/2)+(y-Ly/2)*(y-Ly/2)<pow(R-i/tan(M_PI*alpha/180),2))
			// R(alpha) | r, h - constant
			if(z>param.z_surface&&z<=param.z_surface+1+i&&(x-Lx/2)*(x-Lx/2)+(y-Ly/2)*(y-Ly/2)<pow(r+(h-i)/tan(M_PI*alpha/180),2))
				atoms(x,y,z).type=2;	
// WL			
//			if(z==param.z_surface||z==param.z_surface-1) atoms(x,y,z).type=1;	
		}
	}

*/	

//***********************************************************
/*
// V-trench 30dgr
  for (z=0; z<Lz; z++)	
  ALL_X ALL_Y
  if(z<=param.z_surface&&((x-Lx/4)/50./sqrt(2)+(y-Ly/4)/50./sqrt(2)-(z-param.z_surface-5)/50./sqrt(2)/0.41)<1
					   &&((Lx/4-x)/50./sqrt(2)+(Ly/4-y)/50./sqrt(2)-(z-param.z_surface-5)/50./sqrt(2)/0.41)<1)
    atoms(x,y,z).type=0;

  for (z=0; z<Lz; z++)		
  ALL_X ALL_Y
  if(z<=param.z_surface&&((x-3*Lx/4)/50./sqrt(2)+(y-3*Ly/4)/50./sqrt(2)-(z-param.z_surface-5)/50./sqrt(2)/0.41)<1
					   &&((3*Lx/4-x)/50./sqrt(2)+(3*Ly/4-y)/50./sqrt(2)-(z-param.z_surface-5)/50./sqrt(2)/0.41)<1)
    atoms(x,y,z).type=0;					   
//atom_type[0][0][120]=0;
*/
//***********************************************************
/*
// U-trench 30dgr
  for (z=0; z<Lz; z++)	
  ALL_X ALL_Y
  if(z<=param.z_surface&&((x+y-100)*(x+y-100)/160.-(z-100)<0||(x+y-300)*(x+y-300)/160.-(z-100)<0))
    atoms(x,y,z).type=0;
*/
//***********************************************************
/*
// 
  for (z=0; z<Lz; z++)	
  ALL_X ALL_Y
  if(z<=param.z_surface&&((x+y<358&&x+y>42&&(x+y<158||x+y>242))&&z>40))
	if(atoms(x,y,z).type!=0)
	  if(random_(4)==0) atoms(x,y,z).type=2;
*/
//***********************************************************
/*
// sin
  int a;
//	15	30dgr
//	25	45dgr
  a = 25;
  for (z=0; z<Lz; z++)	
  ALL_X ALL_Y{
  if(((x+y)*sqrt(2)>(Lx+Ly)/8*sqrt(2)&&(x+y)*sqrt(2)<(Lx+Ly)/8*3*sqrt(2))||((x+y)*sqrt(2)>(Lx+Ly)/8*5*sqrt(2)&&(x+y)*sqrt(2)<(Lx+Ly)/8*7*sqrt(2)))
  {
    if(z<=param.z_surface&&a*sin(3.14*sqrt(2)/(Lx/2*sqrt(2))*(x+y)+(Lx+Ly)/8*sqrt(2))+param.z_surface-2*a<z)
	  atoms(x,y,z).type=0;
//defects	  
	else if(z<=param.z_surface&&40*sin(3.14*sqrt(2)/(Lx/2*sqrt(2))*(x+y)+(Lx+Ly)/8*sqrt(2))+param.z_surface-2*31<z)
	  if(random_(4)==0) atoms(x,y,z).type=3;
  }
  else
    if(z<=param.z_surface&&a*sin(3.14*sqrt(2)/(Lx/2*sqrt(2))*(x+y)+(Lx+Ly)/8*sqrt(2))+param.z_surface-2*a<z)
	  atoms(x,y,z).type=0;
//defects	  
	else if(z<=param.z_surface&&40*sin(3.14*sqrt(2)/(Lx/2*sqrt(2))*(x+y)+(Lx+Ly)/8*sqrt(2))+param.z_surface-2*31<z)
	  if(random_(4)==0) atoms(x,y,z).type=3;	  
  }

//erase defects near surface
ZXY 
  if(atoms(x,y,z).type!=0){
    set_config(x,y,z);
    if(n2_config[atoms(x,y,z).config]<12)
      atoms(x,y,z).type=1;
  }
*/
//***********************************************************
// paraboloid
/*
  for (z=0; z<Lz; z++)	
  ALL_X ALL_Y{
  if(z<=param.z_surface&&(z-50)/10.>(x-Lx/2)*(x-Lx/2)/20/20.+(y-Ly/2)*(y-Ly/2)/20/20.)
    atoms(x,y,z).type=2;
  }
*/ 
//***********************************************************

// pits

float r = 50;
float alpha = 30;			
																// cone
alpha = 90 - alpha;																		// cone
float a_c = sin(M_PI*alpha/180)/cos(M_PI*alpha/180);	// cone
float h = r/a_c; 																			// cone
// param.z_surface should be at least r/a_c						// cone
//float p = 3*r;

for (z=0; z<Lz; z++)
	ALL_X ALL_Y{

/*
		if(z<=param.z_surface&&z-37>1.41*sqrt((x-Lx/2+0.5)*(x-Lx/2+0.5)+(y-Ly/2+0.5)*(y-Ly/2+0.5))+random_(4)
			||z-37>1.41*sqrt((x+0.5)*(x+0.5)+(y+0.5)*(y+0.5))+random_(4)
			||z-37>1.41*sqrt((x-Lx+0.5)*(x-Lx+0.5)+(y+0.5)*(y+0.5))+random_(4)
			||z-37>1.41*sqrt((x+0.5)*(x+0.5)+(y-Ly+0.5)*(y-Ly+0.5))+random_(4)
			||z-37>1.41*sqrt((x-Lx+0.5)*(x-Lx+0.5)+(y-Ly+0.5)*(y-Ly+0.5))+random_(4)
		)
		atoms(x,y,z).type=0;
*/		

// hexagon / square finite pits
//	z^2/c^2 > x^2/a^2 + y^2/a^2
//	a^2/c^2	= sin(alpha)/cos(alpha), alpha - angle of pit 
/*
		if( z > param.z_surface - h&&(
//				(z+h-param.z_surface)*(z+h-param.z_surface)*a_c*a_c>(x-Lx/2)*(x-Lx/2)+(y-Ly/2)*(y-Ly/2)||	// hexa
				(z+h-param.z_surface)*(z+h-param.z_surface)*a_c*a_c>(x)*(x)+(y)*(y)
			||(z+h-param.z_surface)*(z+h-param.z_surface)*a_c*a_c>(x-Lx)*(x-Lx)+(y)*(y)
			||(z+h-param.z_surface)*(z+h-param.z_surface)*a_c*a_c>(x)*(x)+(y-Ly)*(y-Ly)
			||(z+h-param.z_surface)*(z+h-param.z_surface)*a_c*a_c>(x-Lx)*(x-Lx)+(y-Ly)*(y-Ly)
				)
		)
			atoms(x,y,z).type=0;
*/

/*
// Hexagon along <100> / and square along <100> if remove Lx/2_Ly/2 1st line
		if(z>-1&&
//				(x-Lx/2)*(x-Lx/2)+(y-Ly/2)*(y-Ly/2)<pow(r,2)||	// central pit for hexa
				(x)*(x)+(y)*(y)<pow(r,2)
			||(x-Lx)*(x-Lx)+(y)*(y)<pow(r,2)
			||(x)*(x)+(y-Ly)*(y-Ly)<pow(r,2)
			||(x-Lx)*(x-Lx)+(y-Ly)*(y-Ly)<pow(r,2)				
		)
		atoms(x,y,z).type=0;
*/		
		
/*	
// Hexagon along <110>
// Lx = Ly = p*5/2^0.5	~ 1592
// 											~ 1856	(p = 3.5*r, r = 150)
		if(z>-1
			&&(x - Lx *  0/15) * (x - Lx *  0/15) + (y - Ly *  0/15) * (y - Ly *  0/15) < pow(r,2)	//	1
			||(x - Lx *  0/15) * (x - Lx *  0/15) + (y - Ly * 15/15) * (y - Ly * 15/15) < pow(r,2)	//	2
			||(x - Lx *  1/15) * (x - Lx *  1/15) + (y - Ly * 11/15) * (y - Ly * 11/15) < pow(r,2)	//	3
			||(x - Lx *  2/15) * (x - Lx *  2/15) + (y - Ly *  7/15) * (y - Ly *  7/15) < pow(r,2)	//	4
			||(x - Lx *  3/15) * (x - Lx *  3/15) + (y - Ly *  3/15) * (y - Ly *  3/15) < pow(r,2)	//	5
			||(x - Lx *  4/15) * (x - Lx *  4/15) + (y - Ly * 14/15) * (y - Ly * 14/15) < pow(r,2)	//	6
			||(x - Lx *  5/15) * (x - Lx *  5/15) + (y - Ly * 10/15) * (y - Ly * 10/15) < pow(r,2)	//	7
			||(x - Lx *  6/15) * (x - Lx *  6/15) + (y - Ly *  6/15) * (y - Ly *  6/15) < pow(r,2)	//	8
			||(x - Lx *  7/15) * (x - Lx *  7/15) + (y - Ly *  2/15) * (y - Ly *  2/15) < pow(r,2)	//	9
			||(x - Lx *  8/15) * (x - Lx *  8/15) + (y - Ly * 13/15) * (y - Ly * 13/15) < pow(r,2)	//	10
			||(x - Lx *  9/15) * (x - Lx *  9/15) + (y - Ly *  9/15) * (y - Ly *  9/15) < pow(r,2)	//	11
			||(x - Lx * 10/15) * (x - Lx * 10/15) + (y - Ly *  5/15) * (y - Ly *  5/15) < pow(r,2)	//	12
			||(x - Lx * 11/15) * (x - Lx * 11/15) + (y - Ly *  1/15) * (y - Ly *  1/15) < pow(r,2)	//	13
			||(x - Lx * 12/15) * (x - Lx * 12/15) + (y - Ly * 12/15) * (y - Ly * 12/15) < pow(r,2)	//	14
			||(x - Lx * 13/15) * (x - Lx * 13/15) + (y - Ly *  8/15) * (y - Ly *  8/15) < pow(r,2)	//	15
			||(x - Lx * 14/15) * (x - Lx * 14/15) + (y - Ly *  4/15) * (y - Ly *  4/15) < pow(r,2)	//	16
			||(x - Lx * 15/15) * (x - Lx * 15/15) + (y - Ly *  0/15) * (y - Ly *  0/15) < pow(r,2)	//	17
			||(x - Lx * 15/15) * (x - Lx * 15/15) + (y - Ly * 15/15) * (y - Ly * 15/15) < pow(r,2)	//	18
			
			||(x + Lx *  1/15) * (x + Lx *  1/15) + (y - Ly *  4/15) * (y - Ly *  4/15) < pow(r,2)	//	a
			||(x - Lx * 11/15) * (x - Lx * 11/15) + (y - Ly * 16/15) * (y - Ly * 16/15) < pow(r,2)	//	b
			||(x - Lx * 16/15) * (x - Lx * 16/15) + (y - Ly * 11/15) * (y - Ly * 11/15) < pow(r,2)	//	c
			||(x - Lx *  4/15) * (x - Lx *  4/15) + (y + Ly *  1/15) * (y + Ly *  1/15) < pow(r,2)	//	d
		)
		atoms(x,y,z).type=0;
*/
/*
// Square along <110>
// Lx = Ly = p*2*2^0.5	~ 1272	(p = 3*r, r = 150)
// 											~ 1484	(p = 3.5*r, r = 150)
		if(z>-1
			&&(x - Lx * 0/4) * (x - Lx * 0/4) + (y - Ly * 0/4) * (y - Ly * 0/4) < pow(r,2)	//	1
			||(x - Lx * 2/4) * (x - Lx * 2/4) + (y - Ly * 0/4) * (y - Ly * 0/4) < pow(r,2)	//	2
			||(x - Lx * 4/4) * (x - Lx * 4/4) + (y - Ly * 0/4) * (y - Ly * 0/4) < pow(r,2)	//	3
			||(x - Lx * 1/4) * (x - Lx * 1/4) + (y - Ly * 1/4) * (y - Ly * 1/4) < pow(r,2)	//	4
			||(x - Lx * 3/4) * (x - Lx * 3/4) + (y - Ly * 1/4) * (y - Ly * 1/4) < pow(r,2)	//	5
			||(x - Lx * 0/4) * (x - Lx * 0/4) + (y - Ly * 2/4) * (y - Ly * 2/4) < pow(r,2)	//	6
			||(x - Lx * 2/4) * (x - Lx * 2/4) + (y - Ly * 2/4) * (y - Ly * 2/4) < pow(r,2)	//	7
			||(x - Lx * 4/4) * (x - Lx * 4/4) + (y - Ly * 2/4) * (y - Ly * 2/4) < pow(r,2)	//	8
			||(x - Lx * 1/4) * (x - Lx * 1/4) + (y - Ly * 3/4) * (y - Ly * 3/4) < pow(r,2)	//	9
			||(x - Lx * 3/4) * (x - Lx * 3/4) + (y - Ly * 3/4) * (y - Ly * 3/4) < pow(r,2)	//	10
			||(x - Lx * 0/4) * (x - Lx * 0/4) + (y - Ly * 4/4) * (y - Ly * 4/4) < pow(r,2)	//	11
			||(x - Lx * 2/4) * (x - Lx * 2/4) + (y - Ly * 4/4) * (y - Ly * 4/4) < pow(r,2)	//	12
			||(x - Lx * 4/4) * (x - Lx * 4/4) + (y - Ly * 4/4) * (y - Ly * 4/4) < pow(r,2)	//	13
		)
		atoms(x,y,z).type=0;
*/			

  }

/*
for (z = param.z_surface-2; z <= param.z_surface; z++)	
	ALL_X ALL_Y{
		if(atoms(x,y,z).type == 1) atoms(x,y,z).type = 2;
	}
*/

/*
//cap
for (z=0; z<param.z_surface+1+3; z++){
    ALL_X ALL_Y  
		if(atoms(x,y,z).type == 0)
			atoms(x,y,z).type = 2;
	}


for (z=0; z<Lz; z++)
  ALL_X ALL_Y{
	  if(z<=param.z_surface+3&&z-37-4>1.41*sqrt((x-Lx/2+0.5)*(x-Lx/2+0.5)+(y-Ly/2+0.5)*(y-Ly/2+0.5))+random_(4))
		  atoms(x,y,z).type=0;
  }


//+4  	
for (z=0; z<Lz; z++)	
  ALL_X ALL_Y{
	  if(z<=param.z_surface+3&&z-37-4>1.41*sqrt((x+0.5)*(x+0.5)+(y+0.5)*(y+0.5))+random_(4))
		  atoms(x,y,z).type=0;
	  if(z<=param.z_surface+3&&z-37-4>1.41*sqrt((x-Lx+0.5)*(x-Lx+0.5)+(y+0.5)*(y+0.5))+random_(4))
		  atoms(x,y,z).type=0;
	  if(z<=param.z_surface+3&&z-37-4>1.41*sqrt((x+0.5)*(x+0.5)+(y-Ly+0.5)*(y-Ly+0.5))+random_(4))
		  atoms(x,y,z).type=0;
	  if(z<=param.z_surface+3&&z-37-4>1.41*sqrt((x-Lx+0.5)*(x-Lx+0.5)+(y-Ly+0.5)*(y-Ly+0.5))+random_(4))
		  atoms(x,y,z).type=0;
  }  
*/



/*
for (z=0; z<50; z++){
    ALL_X ALL_Y  
		if(atoms(x,y,z).type == 0)
			atoms(x,y,z).type = 2;
	}
*/
/*
for (z=0; z<Lz; z++)	
  ALL_X ALL_Y{
//  if(z>100-((x-Lx/2)*(x-Lx/2)+(y-Ly/2)*(y-Ly/2))-20*20+2*20*sqrt((x-Lx/2)*(x-Lx/2)+(y-Ly/2)*(y-Ly/2)));
	if(z>=param.z_surface+1+6.5+a+a*cos(sqrt((x-Lx/2)*(x-Lx/2)+(y-Ly/2)*(y-Ly/2))/b)&&(x-Lx/2)*(x-Lx/2)+(y-Ly/2)*(y-Ly/2)<3.14*b*3.14*b)
    atoms(x,y,z).type=0;
  }
*/  
  
  
 
	
  if(param.z_cap>-100)
    for (z=0; z<param.z_surface+1+param.z_cap; z++){
      ALL_X ALL_Y  
        if(atoms(x,y,z).type == 0)
          atoms(x,y,z).type = 2; // закрываем структуру param.z_cap слоями кремния
    }
	
// //cap 2
//   if(param.z_cap>-100)
//     for (z=param.z_surface+1+param.z_cap; z<param.z_surface+1+param.z_cap+3; z++){
//       ALL_X ALL_Y  
//         if(atoms(x,y,z).type == 0)
//           atoms(x,y,z).type = 2; // закрываем структуру param.z_cap слоями Ge
//     }
	
																								// remove atoms with bad config
																								bool bad_config_key = false;
																								int count_bad_cycles = 0;
																								int count_bad_atoms = 0;
																								label_check_config:	{
  ZXY set_config(x,y,z);
  ZXY if ( unlikely(atoms(x,y,z).type>0 && !good_config[atoms(x,y,z).config]) ){
  	fprintf(stderr,"Error!!! atom (%d,%d,%d) has wrong initial configuration: t=%d\n",
  	    x,y,z,atoms(x,y,z).type); 
																								atoms(x,y,z).type = 0;
																								count_bad_atoms ++;
																								bad_config_key = true;
//  	abort();
  }
																									if (bad_config_key) {
																										bad_config_key = false;
																										count_bad_cycles ++;
																										goto label_check_config;
																									}
																								}
																								if(count_bad_cycles !=0)
																									fprintf(stderr,"count_bad_cycles: %d	count_bad_atoms: %d\n",
																													count_bad_cycles, count_bad_atoms);
  
  ZXY set_jump_info(x,y,z);
  set_defect_f();
  ZXY calc_B0(x,y,z);
  ZXY v_ochered_Edef(x,y,z);
  update_Edef();
  ZXY if (atoms(x,y,z).type>0) add_atom_to_spisok(x,y,z);
}





void load_initial_structure(char* filename) {

//char is_filename[100];
//sprintf(is_filename,"load.xyz");


int i,j,k,i1,j1,k1,At,/*z_min,*/z_max,ii,n;
float ai,aj,ak;
long n_atoms;
//float step_time;
int x,y,z;
int defect_x1[N_defects], defect_y1[N_defects], defect_z1[N_defects];

  for(i==0;i<N_defects;i++){
	defect_x1[i]=0;
	defect_y1[i]=0;
	defect_z1[i]=0;
  }  
  
  ZXY{
    memset(&atoms(x,y,z), 0, sizeof(atom_t));
    jump_probability(x,y,z) = n_jumps(x, y, z) = jumps(x,y,z) = 0;
    spisok_flag(x,y,z) = in_ochered_Edef(x,y,z) = 0;
    number_in_spisok_atomov(x,y,z) = -1;
  }

  ochered_Edef.clear();
  spisok_atomov.clear();
  
  ii=0;
  //z_min=10000;
  z_max=0;

  FILE* f;
  f=fopen(/*is_*/filename,"r");

  if (f==NULL) {
    fprintf(stderr,"!!! load_initial: cannot open file '%s' for reading\n",filename);
    exit(0);
  }

  if( fscanf(f,"%ld",&n_atoms) != 1 ) {
    fprintf(stderr,"!!! load_initial: bad file format: '%s'\n",filename);
    exit(0);
  }
  
  if( fscanf(f,"%lf%lf%d%d%d%d",&current.t,&current.t,&current.n_deposited,&current.n_jumps,&current.n_bad_jumps,&current.n_jumps_back) != 6 ) {
    fprintf(stderr,"!!! load_initial: bad file format: '%s'\n",filename);
    exit(0);
  }		
												
												
//  if( fscanf(f,"%*[^\n]") < 0 ){
//    fprintf(stderr,"!!! load_initial: bad file format: '%s'\n",filename);
//    exit(0);
//  }

    while(1){
      n=fscanf(f,"%d%d%f%d%f%d%f",&At,&i1,&ai,&j1,&aj,&k1,&ak);
      if (n!=7) break;
//      n=fscanf(f,"%d%d%d%d",&At,&i1,&j1,&k1);
//      if (n!=4) break;
	  i=i1;j=j1;k=k1;
      if(z_max<=k) z_max=k;
	  
	  if(At==4){
		defect_x1[ii]=i; defect_y1[ii]=j; defect_z1[ii]=k;
		ii++;
	  }
	  
	  else{
		atoms(i,j,k).type = At;
		atoms(i,j,k).a.x=ai;
		atoms(i,j,k).a.y=aj;
		atoms(i,j,k).a.z=ak;
	  }


  	  if (feof(f)) 
  	    break;
      }

      fclose(f);

  
  read_defects("defects.txt");
  
  int n_defects1=0, n_defects2=0;
  int exists=0;      
  
  for(i=0;i<ii;i++){			//По всем прочитанным из *.xyz файла междоузлиям
	exists=0;
	for(j=0;j<n_defects;j++){	//По всем прочитанным из defects.txt файла междоузлиям
	  if(defect_x1[i]==defect_x[j]&&defect_y1[i]==defect_y[j]&&defect_z1[i]==defect_z[j]){
		exists=1;		// Если в такое междоузлие уже есть
		n_defects1++;	// Количество одинаковых междоузлий
	  }
	}
	if(exists==0){		// Если есть не совпадающие междоузлия
	  n_defects2++;		// Количество добавляемых междоузлий
	  defect_x[n_defects-1+n_defects2]=defect_x1[i];	// Добавляем их в список междоузлий
	  defect_y[n_defects-1+n_defects2]=defect_y1[i];
	  defect_z[n_defects-1+n_defects2]=defect_z1[i];
	}
  }
  
  n_atoms-=ii;				// Исключаем междоузлия, прочитанные из *.xyz файла
  n_defects+=n_defects2;	// Обновляем количество междоузлий

  
  if(param.z_cap>-100)
    for (z=0; z<z_max+1+param.z_cap; z++){
      ALL_X ALL_Y  
        if(atoms(x,y,z).type == 0)
          atoms(x,y,z).type = 1; // закрываем структуру param.z_cap слоями кремния
    }
	  
  ZXY set_config(x,y,z);
  ZXY if (atoms(x,y,z).type>0 && !good_config[atoms(x,y,z).config]){
  	fprintf(stderr,"Error!!! atom (%d,%d,%d) has wrong initial configuration\n",x,y,z);
  	abort();
  }


  ZXY set_jump_info(x,y,z);
  set_defect_f();
  ZXY calc_B0(x,y,z);
  ZXY v_ochered_Edef(x,y,z);
  update_Edef();
  ZXY if (atoms(x,y,z).type>0) add_atom_to_spisok(x,y,z);
 
//current.t+=0.000001;
 
}






void set_config(int x, int y, int z) {
  int dir,x2,y2,z2,factor;
  int nb_type[dir_number];
  
  neighbors_t nbs;
  factor = atoms.neighbors(x,y,z,nbs);
  
  for (dir=0; dir<dir_number; dir++) {
    nb_type[dir]=atoms(nbs.x[dir],nbs.y[dir],nbs.z[dir]).type;
  }
  atoms(x,y,z).config = massiv_to_config(nb_type);
}


void set_defect_f(void) {
  int x,y,z,n,dx,dy,dz,x2,y2,z2,sum_dx,sum_dy,sum_dz,sum;
  for (x=0; x<Lx; x++)  // сначала очистим всё
   for (y=0; y<Ly; y++)
    for (z=0; z<Lz; z++){
      atoms(x,y,z).defect_f.x = atoms(x,y,z).defect_f.y = atoms(x,y,z).defect_f.z = 0;
    }

  for (n=0; n<n_defects; n++) {  // затем пробежимся по всем дефектам и добавим вклады от каждого из них
    if ( unlikely(n>=N_defects) ) {fprintf(stderr,"Error!!! in set_defect_f(): n>=N_defects\n"); exit(0);}
    x=defect_x[n]; y=defect_y[n]; z=defect_z[n];
    if ( unlikely(x%2!=y%2 || y%2!=z%2 || (x/2+y/2+z/2)%2!=1) || unlikely(x<0||x>=Lx||x<0||y>=Ly||z<0||z>=Lz)) 
      {fprintf(stderr,"Error!!! bad defect position (%d,%d,%d)\n",x,y,z); exit(0);}
    sum_dx=sum_dy=sum_dz=sum=0;
    for (dx=-1; dx<=1; dx+=2)
     for (dy=-1; dy<=1; dy+=2)
      for (dz=-1; dz<=1; dz+=2) {
        x2=x+dx; y2=y+dy; z2=z+dz;
        boundary(x2,y2,z2);
        if ( unlikely(x2<0 || x2>=Lx || y2<0 || y2>=Ly || z2<0 || z2>=Lz) ) 
          {fprintf(stderr,"Error!!! bad defect position (%d,%d,%d).\n",x,y,z); exit(0);}
        if (x2%2==y2%2 && y2%2==z2%2 && (x2/2+y2/2+z2/2)%2==0) { // если x2,y2,z2 принадл. решётке алмаза
//          if ( unlikely(atoms(x2,y2,z2).type==0) )
//            {fprintf(stderr,"Error!!! bad defect position (%d,%d,%d)..\n",x,y,z); exit(0);}
		  if ( unlikely(atoms(x2,y2,z2).type!=0) ){
			atoms(x2,y2,z2).defect_f.x += dx; // добавляем вклад от дефекта (x,y,z) в силу, действ-ю на (x2,y2,z2)
			atoms(x2,y2,z2).defect_f.y += dy;
			atoms(x2,y2,z2).defect_f.z += dz;
		  }
          sum_dx+=dx; sum_dy+=dy; sum_dz+=dz; sum++;
        }
      }
    if ( unlikely(sum_dx!=0 || sum_dy!=0 || sum_dz!=0 || sum!=4) )
      {fprintf(stderr,"Error!!! bad defect position (%d,%d,%d)...\n",x,y,z); exit(0);}
  }
}


void erase_defect(int x3, int y3, int z3) {	// поиск и удаление дефектов рядом с (x3,y3,z3) и его влияния на окружающие атомы решетки
  int i,j,dx,dy,dz,dx2,dy2,dz2,x,y,z,x2,y2,z2;
  char filename[100];
  
  for (dx2=-1; dx2<=1; dx2+=2)			// поиск дефектов рядом с атомом (x3,y3,z3)
	for (dy2=-1; dy2<=1; dy2+=2)
	  for(dz2=-1; dz2<=1; dz2+=2) {
	  x=x3+dx2; y=y3+dy2; z=z3+dz2;
	  boundary(x,y,z);
	  for(i=0;i<n_defects;i++){		// поиск дефекта (x,y,z) в списке дефектов
		if(x==defect_x[i]&&y==defect_y[i]&&z==defect_z[i]){	// если нашли дефект (x,y,z)

		  for (dx=-1; dx<=1; dx+=2)
			for (dy=-1; dy<=1; dy+=2)
			  for (dz=-1; dz<=1; dz+=2) {
				x2=x+dx; y2=y+dy; z2=z+dz;
				boundary(x2,y2,z2);
				if (x2%2==y2%2 && y2%2==z2%2 && (x2/2+y2/2+z2/2)%2==0 && unlikely(atoms(x2,y2,z2).type!=0)) {
				  atoms(x2,y2,z2).defect_f.x -= dx; // удаляем вклад от дефекта (x,y,z) в силу, действ-ю на окружающие атомы решетки
				  atoms(x2,y2,z2).defect_f.y -= dy;
				  atoms(x2,y2,z2).defect_f.z -= dz;
				}
			  }
			  
		  for (dx=-1; dx<=1; dx+=2)
			for (dy=-1; dy<=1; dy+=2)
			  for (dz=-1; dz<=1; dz+=2) {
				x2=x+dx; y2=y+dy; z2=z+dz;
				boundary(x2,y2,z2);
				if (x2%2==y2%2 && y2%2==z2%2 && (x2/2+y2/2+z2/2)%2==0 && unlikely(atoms(x2,y2,z2).type!=0)) {
				  calc_B0(x2,y2,z2);
				}
			  }
		  	
		  FILE *f;
		  sprintf(filename,"log.txt");
		  f=fopen(filename,"at");
		  fprintf(f,"	-- defect (%d,%d,%d) removed\n",x,y,z);
		  fclose(f);
		  
		  for(j=i;j<n_defects-1;j++){
			defect_x[j]=defect_x[j+1];	// сдвигаем список дефектов на место удаляемого (x2,y2,z2)
			defect_y[j]=defect_y[j+1];
			defect_z[j]=defect_z[j+1];
		  }
		n_defects--;	// сокращаем список каждый раз, как нашли дефект рядом с атомом (x,y,z)
		}
	  }
   }

//   ZXY calc_B0(x,y,z);	// обновляем влияние всех дефектов на все атомы кристалла
}


void add_atom_to_spisok(int x,int y,int z) 
{
  struct coord c = {x, y, z};
  number_in_spisok_atomov(x,y,z) = spisok_atomov.size();
  spisok_atomov.push_back(c);
}

void remove_atom_from_spisok(int x,int y,int z)
{
	int n;
  n = number_in_spisok_atomov(x,y,z);
	spisok_atomov.erase(spisok_atomov.begin() + n);
}


void move_atom_in_spisok(int x_old, int y_old, int z_old, int x_new, int y_new, int z_new) {
  int n;
  n=number_in_spisok_atomov(x_old,y_old,z_old);
  struct coord c = { x_new, y_new, z_new };
  spisok_atomov[n] = c;
  number_in_spisok_atomov(x_new,y_new,z_new) = n;
  number_in_spisok_atomov(x_old,y_old,z_old) = -1;
}


void check_Lxyz(void) {
  if ( unlikely(Lx<=0 || Lx%4!=0) ) {fprintf(stderr,"Error!!! wrong Lx=%d\n",Lx); exit(0);}
  if ( unlikely(Ly<=0 || Ly%4!=0) ) {fprintf(stderr,"Error!!! wrong Ly=%d\n",Ly); exit(0);}
  if ( unlikely(Lz<=0 || Lz%4!=0) ) {fprintf(stderr,"Error!!! wrong Lz=%d\n",Lz); exit(0);}
}
