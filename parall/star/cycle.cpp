#include "mc_growth.h" 



void mc_step(void);   // шаг Монте-Карло
double calc_P_jump_sum(void); // сосчитать суммарную вероятность всех прыжков
void choose_jump(double P_jump_sum, int* x, int* y, int* z, int* dir); // выбрать прыжок
void do_many_axyz(void/*int x, int y, int z*/);  // сделать "шевеления" атомов после прыжка или осаждения

void main_loop(void) { // бесконечный цикл, состоящий из шагов Монте-Карло
  int i, n = 0;
  int x, y, z;
  float ML_next_control;
  char filename[100];
  float ML_deposited;
  double t;
	float T;

  ML_deposited = current.n_deposited / (Lx * Ly / 8.);

  ML_next_control = param.dML_control * ((int) (ML_deposited / param.dML_control));// the last control check

  // если исходная структура без смещений и планируется моделирование роста, то 1000 раз пересчитываем смещения всех атомов
  if(((ML_deposited==0)&&param.experiment_ML!=0)){
    for(n=0; n<1000; n++){
			std::random_shuffle(spisok_atomov.begin(),spisok_atomov.end());
      for(i=0;i<spisok_atomov.size();i++){
        struct coord c;
        c = spisok_atomov[i];
        axyz(c.x,c.y,c.z);
      }
    }
		ZXY v_ochered_Edef(x,y,z);
		update_Edef();
	}

  // Thermal annealing
  if (param.Edef_pa_N>0) {
  
	read_control("control.txt");
	float Edef_sum = 0;
	   
	FILE *f;

    for(n=0; n<param.Edef_pa_N; n++) {
	  
		std::random_shuffle(spisok_atomov.begin(),spisok_atomov.end());
		for (i=0; i<spisok_atomov.size(); i++) {
			struct coord c = spisok_atomov[i];
			axyz(c.x,c.y,c.z);
			calc_Edef_pa(c.x,c.y,c.z);
			Edef_sum += atoms(c.x,c.y,c.z).Edef_pa;
			atoms(c.x,c.y,c.z).Edef_pa = 0;
		}
	  
		sprintf(filename,"Esum_log.txt");
		f=fopen(filename,"at");
		if (f==NULL) {fprintf(stderr,"!!! cannot open log file '%s'\n",filename); return;}
		fprintf(f,"%d	%.2f	%.5f\n",n,param.T,Edef_sum); 
		fclose(f);
		Edef_sum = 0;	  
    }
	
	T = param.T;
	
	do {
		T -= 1.;
		if(T < 0) T = 0;
		
		for(n = 0; n < (T==0 ? 100 : 10); n++){
			
			std::random_shuffle(spisok_atomov.begin(),spisok_atomov.end());
			for (i=0; i<spisok_atomov.size(); i++) {
				struct coord c = spisok_atomov[i];
				axyz(c.x,c.y,c.z);
			}
		}

		for (i=0; i<spisok_atomov.size(); i++) {
			struct coord c = spisok_atomov[i];
			calc_Edef_pa(c.x,c.y,c.z);
			Edef_sum += atoms(c.x,c.y,c.z).Edef_pa;
			atoms(c.x,c.y,c.z).Edef_pa = 0;
		}

		sprintf(filename,"Esumcool_log.txt");
		f=fopen(filename,"at");
		if (f==NULL) {fprintf(stderr,"!!! cannot open log file '%s'\n",filename); return;}
		fprintf(f,"%.2f	%.5f\n",T,Edef_sum); 
		fclose(f);
		Edef_sum = 0;
		
	} while(T > 0);
	
	sprintf(filename,"aaT0N%d.xyz",param.Edef_pa_N);
	show_me(filename,true);
    comp.notify(filename);
  }

	float show_short;
  float show_long;
  float show_sum;
	
	if(param.show_short != 0)	show_short = param.show_short * ((int) (ML_deposited / param.show_short));
	else show_short = 0;
	
	if(param.show_long != 0) show_long = param.show_long * ((int) (ML_deposited / param.show_long));
	else show_long = 0;
	
  if(param.show_sum != 0) show_sum = param.show_sum * ((int) (ML_deposited / param.show_sum));
	else show_sum = 0;
	
  t=current.t;
  
// main loop
  while (1) { // и запустим цикл
	
  if (ML_deposited >= ML_next_control) {
      read_control("control.txt");
      write_to_log("log.txt");
      ML_next_control += param.dML_control;
    }

	if(param.show_short != 0)
		if(ML_deposited >= show_short) {
			sprintf(filename,"A%.2f.xyz",show_short);
			show_me(filename,0);
			comp.notify(filename);
			show_short += param.show_short;		
		}
	
	if(param.show_long != 0)
		if(ML_deposited >= show_long) {
			sprintf(filename,"_A%.2f.xyz",show_long);
			show_me(filename,1);
			comp.notify(filename);
			show_long += param.show_long;		
		}
		
	if(param.show_sum != 0)
		if(ML_deposited >= show_sum) {
			sprintf(filename,"__A%.2f.xyz",show_sum);
			show_me(filename,2);
			comp.notify(filename);
			show_sum += param.show_sum;
		}

  if (ML_deposited>=param.experiment_ML) break;
		
    mc_step();
	
//	was 1 per 1e-6 s (mjus) at 670K
	if(current.t-t>=param.time_for_moves){	// calculate displacements per param.time_for_moves seconds
  	  do_many_axyz();
	  t=current.t;
	}

/*
    i++;
    
    if (i%10000==0) {
      save_grid("M.txt",5);
      printf("step %d  t=%e  deposited=%d\n",i,current.t,current.n_deposited);
    }
    if (current.t>=t_next) {
      sprintf(filename,"M%.0f.txt",t_next);
      save_grid(filename,5);
      t_next+=dt;
    }
*/

  ML_deposited=current.n_deposited/(Lx*Ly/8.);
  }
/*
  write_to_log("log.txt");
  if(param.experiment_ML!=0){
	sprintf(filename,"_A%.2f.xyz",ML_deposited);
	show_me(filename,1);
	comp.notify(filename);
  }
*/
}

void mc_step(void) { // шаг Монте-Карло
  double P_jump_sum,P_total,Px,dt;
  int x,y,z,dir,select,res,x2,y2,z2;
  
    update_Edef(); // приводим энергии деформ. и вероятности прыжков в соответствие с текущим положением
    
    P_jump_sum=calc_P_jump_sum(); // сосчитаем суммарную вероятность всех прыжков
    P_total=P_jump_sum+param.p_deposition;
    Px=P_total*rand01();  // выберем, что делать - прыжок или осаждение
    if (Px<param.p_deposition) { // делаем осаждение
      select=1;
      res=deposition(param.dep_type,&x,&y,&z);
    }
    else {                       // или делаем прыжок
      select=2;
      choose_jump(P_jump_sum,&x,&y,&z,&dir); // выберем, какой прыжок
      res=jump(x,y,z,dir,&x2,&y2,&z2);       // и выполним его
    }
    dt=-log(rand01())/P_total;   // обновим время
    current.t+=dt;

    // сделаем "шевеления" атомов после прыжка или осаждения
//    if (select==1 && res==1) do_many_axyz(x,y,z); 
//    if (select==2 && res==1) do_many_axyz((x+x2)/2,(y+y2)/2,(z+z2)/2);
    
    // и напоследок обновим счётчики событий
    if (select==1 && res==1) current.n_deposited++;
    if (select==2 && res==0) current.n_bad_jumps++;
    if (select==2 && res==1) {
      current.n_jumps++;
      if (x2==current.prev_x && y2==current.prev_y && z2==current.prev_z) current.n_jumps_back++;
      current.prev_x=x; current.prev_y=y; current.prev_z=z;
    }
    
//    if (select==1 && res==1) printf("                            %d   %d   %d.\n",x,y,z);
    //if (select==2 && res==0) printf("(%d %d %d) -> (%d %d %d) - otmenen\n",x,y,z,x2,y2,z2);
    //if (select==2 && res==1) printf("(%d %d %d) -> (%d %d %d)\n",x,y,z,x2,y2,z2);
    
}

double calc_P_jump_sum(void) { // сосчитать суммарную вероятность всех прыжков
  double real_sum = 0;
  int x, y, z;
    auto start = std::chrono::steady_clock::now();
  
  #pragma omp parallel reduction (+: real_sum) num_threads(4)
  {
  # pragma omp for
    for (z=2; z<Lz-2; z++) 
    {
      for (y=((int)z%2); y<Ly; y+=2)
      {

      
        for (x=((int)z%2)+2*(((int)z/2+(int)y/2)%2); x<Lx; x+=4) 
        {
          real_sum += jump_probability(x,y,z);
        }
      }
    }
  }
  auto end = std::chrono::steady_clock::now();
  sumtime+= std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
  return real_sum;
}

void choose_jump(double P_jump_sum, int* x_, int* y_, int* z_, int* dir_) { // выбрать прыжок
  int x,y,z,f,dir,n,nx;
  int j_array[dir_number];
  double P,Px;
  
  if ( unlikely(P_jump_sum<=0) ) {fprintf(stderr,"Error!!! in choose_jump(): P_jump_sum<=0\n"); exit(0);}
  
  Px=P_jump_sum*rand01();
  P=0; f=0;
  ZYX {
    P+=jump_probability(x,y,z);
    if (P>Px) {f=1; goto choose_jump_exit;} // нашли нужный узел (x,y,z) - выходим из цикла
  }
  choose_jump_exit: // выходим из цикла сюда
  
  if ( unlikely(f==0) )  {fprintf(stderr,"Error!!! in choose_jump() P=%f Px=%f P_jump_sum=%f\n",P,Px,P_jump_sum); exit(0);}
  if ( unlikely(x<0||x>=Lx||y<0||y>=Ly||z<2||z>=Lz-2) )  {fprintf(stderr,"Error!!! in choose_jump().\n"); exit(0);}
  if ( unlikely(n_jumps(x, y, z) <=0) ) {fprintf(stderr,"Error!!! in choose_jump()..\n"); exit(0);}
  
  fill_nb_type(j_array, jumps(x,y,z) );
  nx=1+random_( n_jumps(x,y,z) ); // которое по счёту направление прыжка мы выбираем
  n=0; f=0;
  for (dir=0; dir<dir_number; dir++) {
    n+=j_array[dir];
    if (n==nx) {f=1; break;}
  }
  
  if ( unlikely(f==0) )  {fprintf(stderr,"Error!!! in choose_jump()...\n"); exit(0);}
  if ( unlikely(dir<0 || dir>=dir_number) ) {fprintf(stderr,"Error!!! in choose_jump()....\n"); exit(0);}
  
  *x_=x; *y_=y; *z_=z; *dir_=dir;
}

/*
void do_many_axyz(int x0, int y0, int z0) {  // сделать "шевеления" атомов после прыжка или осаждения
  int x,y,z, i,i_, I2,Imax;
  I2=500; //1000;
  Imax=1000000;
  i=0;
  for (i_=0; i_<Imax; i_++) {
    z = 2 + random_(Lz-4);
    x = (z%2) + 2*random_(Lx/2);
    y = (z%2) + 2*((z/2+x/2)%2) + 4*random_(Ly/4);
    if (x<0 || x>=Lx || y<0 || y>=Ly || z<2 || z>=Lz-2)
      {fprintf(stderr,"Error!!! in do_many_axyz()\n"); exit(0);}
    if (atom_type[x][y][z]!=0) {
      axyz(x,y,z);
      i++;
    }
    if (i>=I2) break;
  }
}
*/

void do_many_axyz(void/*int x0, int y0, int z0*/) {  // сделать "шевеления" атомов после прыжка или осаждения
  int i,I,n;
  I=(int) (param.moves_percent/100.*spisok_atomov.size());
  for (i=0; i<I; i++) {
    n=random_(spisok_atomov.size());
    struct coord c = spisok_atomov[n];
    axyz(c.x,c.y,c.z);
	current.n_moves++;
  }
}

