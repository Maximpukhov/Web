#include <ctype.h>
#include "mc_growth.h"
#include "compress.h"

//------- глобальные переменные и массивы -------------

int Lx=40, Ly=40, Lz=100;

struct param_struct param;
struct current_struct current;

//char_array atom_type;
//float_array ax,ay,az;
//float_array Edef;
//u_s_int_array config; // конфигурация соседей
//float_array B0x,B0y,B0z; // величины B при нулевых смещениях атомов окружения
//char_array defect_fx, defect_fy, defect_fz; // силы, действующие на атомы со стороны дефектов
std::list<struct coord> ochered_Edef;
std::vector<struct coord> spisok_atomov;

dlattice_t <atom_t> atoms;

//float_array jump_probability;
dlattice_t <float> jump_probability;

//char_array n_jumps; //  кол-во разрешённых прыжков из данного узла
dlattice_t <char> n_jumps;

//u_s_int_array jumps; // разрешённые прыжки из данного узла (каждый бит соответствует направлению)
dlattice_t <unsigned short int> jumps;

//char_array in_ochered_Edef; //  =1 для атомов, поставленных в очередь на обновление Edef
dlattice_t <char> in_ochered_Edef; //  =1 для атомов, поставленных в очередь на обновление Edef

//int_array number_in_spisok_atomov; // используется вместе с массивом spisok_atomov
dlattice_t <int> number_in_spisok_atomov; // используется вместе с массивом spisok_atomov


//char_array spisok_flag; //  используется в jump()
dlattice_t <char> spisok_flag; //  используется в jump()
unsigned long sumtime = 0; 
//#define Nconfig 65536 // 2^16
//#define dir_number 16 // число 1-х и 2-х соседей
float AA[Nconfig][6]; // Axx,Ayy,Azz,Axy,Axz,Ayz
float AA_[Nconfig][6]; // (A^-1)xx ... (A^-1)yz
float BB[Nconfig][dir_number][9]; // Bxx,Bxy,Bxz,Byx,Byy,Byz,Bzx,Bzy,Bzz
//float BB_new[Nconfig][3][dir_number*3] __attribute((aligned(64)));
float transform_array[Nconfig][6]; // для random_displacements
char good_config[Nconfig]; //  =1, если конфигурация "правильная"; =0, если нет
char n1_config[Nconfig]; //  число 1-х соседей для данной конфигурации
char n2_config[Nconfig]; //  число 2-х соседей для данной конфигурации

int n_defects; // число дефектов
int defect_x[N_defects], defect_y[N_defects], defect_z[N_defects]; // координаты дефектов



//------- main -------------

ckpt_compressor comp;

int main(int argc, char** argv) {
  int z;


  if( argc > 1 && !strcmp(argv[1], "-tar") ){
      int freq = 10;
      if( argc > 2 && isdigit(argv[2][0])){
        freq = atoi(argv[2]);
      }
      comp.init("mc_sim_ckpt", freq);
  }
  
 // сначала обнулим все счётчики
  current.n_deposited=current.n_jumps=current.n_bad_jumps=current.n_jumps_back=current.n_moves=0;
  current.prev_x=current.prev_y=current.prev_z=-1;  
  current.t=0;

  //  назначим значения параметров по умолчанию
 
  sprintf(param.load,"0");
  param.show_short = 0.1; 		// сохранять короткий *.xyz файл (без смещений) каждые param.show_short монослоев
  param.show_long = 1.0; 		// сохранять подробный *.xyz файл (со смещениями) каждые param.show_long монослоев
  param.show_sum = 0.0; 		// сохранять *.xyz файл с актуальными координатами атомов (x+ax) каждые param.show_sum монослоев
  param.E1[0]=0;
  param.E2[0]=0;
  param.E1[1]=(0.35+0.35*0.25)*11800;
  param.E2[1]=(0.15+0.15*0.25)*11800;
  param.E1[2]=0.35*11800;
  param.E2[2]=0.15*11800;
  param.E1[3]=0.35*11800*10;
  param.E2[3]=0.15*11800*10;
  param.p0=1e13;
  param.p_deposition=0.1*Lx*Ly/8.0;
  param.T=670;
  param.A=48.5 * 1.36e-10*1.36e-10/1.6e-19 * 11800 / 8.0;
  param.B=13.8   * 1.36e-10*1.36e-10/1.6e-19 * 11800 / 8.0;
  param.eps[0]=0;
  param.eps[1]=0;
  param.eps[2]=0.04;
  param.eps[3]=0.1;
  param.defect_F= (5 * 1.5 * 11800) / sqrt(3.0);
  param.moves_percent=1;	// percent of ALL atoms to axyz()
  param.time_for_moves=1e-6;	// calculate atom displacements each param.time_for_moves seconds (1e-6 at 670 K in previous versions)
  param.dML_control=0.01;
  param.experiment_ML=6;
  param.random_seed=19;
  param.dep_type=2;
  param.z_surface=20; // Lz/2;
  param.z_cap = -100; // количество закрывающих монослоев (-100 = закрывающего слоя нет)
  n_defects=-1; // отрицательное n_defects означает, что берём их столько, сколько есть в файле defects.txt
  // А если n_defects положительно либо 0, то берём из файла defects.txt только указанное число дефектов.

  param.Edef_pa_N = 0;  // количество смещений по которым усреднять деформацию на атом, 0 - не считать  
  
  read_parameters("parameters.txt");  // а теперь прочитаем значения параметров из parameters.txt

  // Tune lattice once we know parameters 
  atoms.tune(Lx,Ly,Lz);
  jump_probability.tune(Lx,Ly,Lz);
  n_jumps.tune(Lx,Ly,Lz);
  jumps.tune(Lx,Ly,Lz);
  in_ochered_Edef.tune(Lx,Ly,Lz);
  number_in_spisok_atomov.tune(Lx,Ly,Lz);
  spisok_flag.tune(Lx,Ly,Lz);
  
  
  // ----  здесь заканчивается назначение величин параметров и начинается собственно расчёт
  
  
  srand(param.random_seed);
//  read_defects("defects.txt");
  
  fill_n1_n2_and_good_config();
  calc_AA();
  calc_AA_();
  calc_BB();


  if(strcmp(param.load,"0")==0){
    read_defects("defects.txt");
    int initial_layers[Lz];

    for (z=0; z<Lz; z++) {
      if (z<=param.z_surface) initial_layers[z]=1;
      else                    initial_layers[z]=0;
    }
    create_initial_structure(initial_layers);
    printf("create_initial_structure - OK.\n");
    char filename[100];	
    sprintf(filename,"_initial.xyz");
    show_me(filename,true);
    comp.notify(filename);
    printf("initial show_me - OK\n\n");
  }else{  
    load_initial_structure(param.load);
    printf("load_initial_structure - OK.\n");
  }

	if(current.n_deposited==0)  write_params_to_log("log.txt");
  auto start = std::chrono::steady_clock::now();
  main_loop();
  auto end = std::chrono::steady_clock::now();
  printf("End progrramm = %lu ns \n", std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
  printf("End summ time = %lu ns \n",sumtime);
  comp.finish();
  printf("Finished successfully.\n");

  return 1;
}

