#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <list>
#include <vector>
#include "compress.h"
#include "optimization.h"
#include <algorithm>    // std::shuffle
#include <chrono>
#include "dlattice.h"
#include <stdint.h>


static inline uint32_t rng_u32(uint32_t &s) {
  s ^= s << 13;
  s ^= s >> 17;
  s ^= s << 5;
  return s;
}

static inline double rng01(uint32_t &s) {
  return (double)(rng_u32(s) >> 8) * (1.0 / 16777216.0);
}

static inline int rng_int(uint32_t &s, int n) {
  // [0; n)
  return (n <= 0) ? 0 : (int)(rng_u32(s) % (uint32_t)n);
}


//------- глобальные переменные и массивы -------------

extern int Lx, Ly, Lz;

struct atom_t{
  float Edef;
  float Edef_pa;
  struct {
    float x;
    float y;
    float z;
  } a;
  struct {
    float x;
    float y;
    float z;
  } B0;
  unsigned short config;
  char type;
  struct {
    char x;
    char y;
    char z;
  } defect_f;
};

struct param_struct {
  float E1[4];
  float E2[4];
  float p0;
  float p_deposition;
  float T;
  float A;
  float B;
  float eps[4];
  float defect_F;
  float moves_percent;
  double time_for_moves;
  float dML_control;
  float experiment_ML;
  int random_seed;
  int dep_type;
  int z_surface;
  int z_cap;
  char load[100];
  float show_short;
  float show_long;
  float show_sum;
  int Edef_pa_N;
  int run_mode;
};
extern struct param_struct param;

enum run_mode_enum {
  RUN_MODE_PERF = 0,
  RUN_MODE_VALIDATION = 1
};

struct current_struct {
  double t;
  int n_deposited;
  int n_jumps;
  int n_bad_jumps;
  int n_jumps_back;
  double n_moves;
  int prev_x, prev_y, prev_z;
};
extern struct current_struct current;

/* Change independent definitions to grouping one to acheive better cache hit rate
extern char_array atom_type;
extern float_array ax,ay,az;
extern float_array Edef;
extern char_array defect_fx, defect_fy, defect_fz; // силы, действующие на атомы со стороны дефектов
extern float_array B0x,B0y,B0z; // величины B при нулевых смещениях атомов окружения
extern u_s_int_array config; // конфигурация соседей
*/
//extern atom_array atoms;

extern dlattice_t <atom_t> atoms;

//extern float_array jump_probability;
extern dlattice_t <float> jump_probability;

//extern char_array n_jumps; //  кол-во разрешённых прыжков из данного узла
extern dlattice_t <char> n_jumps; //  кол-во разрешённых прыжков из данного узла
//extern u_s_int_array jumps; // разрешённые прыжки из данного узла (каждый бит соответствует направлению)
extern dlattice_t <unsigned short int> jumps; // разрешённые прыжки из данного узла (каждый бит соответствует направлению)
//extern char_array in_ochered_Edef; //  =1 для атомов, поставленных в очередь на обновление Edef
extern dlattice_t <char> in_ochered_Edef; //  =1 для атомов, поставленных в очередь на обновление Edef

//extern int_array number_in_spisok_atomov; // используется вместе с массивом spisok_atomov
extern dlattice_t <int> number_in_spisok_atomov; // используется вместе с массивом spisok_atomov

//extern char_array spisok_flag; //  используется в jump()
extern dlattice_t <char> spisok_flag; //  используется в jump()

extern std::list <struct coord> ochered_Edef;
extern std::vector <struct coord> spisok_atomov;

extern ckpt_compressor comp;
extern unsigned long sumtime;

#define Nconfig 65536 // 2^16

extern float AA[Nconfig][6]; // Axx,Ayy,Azz,Axy,Axz,Ayz
extern float AA_[Nconfig][6]; // (A^-1)xx ... (A^-1)yz
extern float BB[Nconfig][dir_number][9]; // Bxx,Bxy,Bxz,Byx,Byy,Byz,Bzx,Bzy,Bzz
//extern float BB_new[Nconfig][3][dir_number*3] __attribute((aligned(64)));
extern float transform_array[Nconfig][6]; // для random_displacements
extern char good_config[Nconfig]; //  =1, если конфигурация "правильная"; =0, если нет
extern char n1_config[Nconfig]; //  число 1-х соседей для данной конфигурации
extern char n2_config[Nconfig]; //  число 2-х соседей для данной конфигурации


extern int n_defects; // число дефектов
#define N_defects 50000
extern int defect_x[N_defects], defect_y[N_defects], defect_z[N_defects]; // координаты дефектов


//------- макросы -------------

#define ZXY       ALL_Z ALL_X ALL_Y   // перебор атомов решётки алмаза, кроме 2 верхних и 2 нижних слоёв
#define ALL_Z     for (z=2; z<Lz-2; z++)       // z может быть любым
#define ALL_X     for (x=((int)z%2); x<Lx; x+=2)  // x должно иметь ту же чётность, что и z
#define ALL_Y     for (y=((int)z%2)+2*(((int)z/2+(int)x/2)%2); y<Ly; y+=4)  // y идёт через 4

#define ZYX       ALL_Z ALL_Y1 ALL_X1   // перебор атомов решётки алмаза, кроме 2 верхних и 2 нижних слоёв
#define ALL_Y1     for (y=((int)z%2); y<Ly; y+=2)  // x должно иметь ту же чётность, что и z
#define ALL_X1     for (x=((int)z%2)+2*(((int)z/2+(int)y/2)%2); x<Lx; x+=4)  // y идёт через 4

#define ZYX_a       ALL_Z_a ALL_Y1 ALL_X1   // перебор всех атомов решётки алмаза
#define ALL_Z_a     for (z=0; z<Lz; z++)       // z может быть любым


//------- объявления функций -------------

// из deform.c
void axyz(int x, int y, int z); // "обновление" смещения атома (x,y,z), с. 34-38

// из deform2.c
void calc_Edef(int x, int y, int z); // обновить Edef для атома (x,y,z)
void calc_B0(int x, int y, int z);
void calc_AA();
void calc_AA_();
void calc_BB();

// из deform3.c
void calc_Edef_pa(int x,int y,int z); // посчитать деформацию на атом (x,y,z)

// из geometry.c
int random_(int n); // случайное число от 0 до n-1
double rand01(void);  // случайное число между 0 и 1
void calc_x2y2z2(int x, int y, int z, int dir, int* x2, int* y2, int* z2);
void erase_defect(int x, int y, int z);	// удалить дефект, действующий на атом (x,y,z) и его влияние на окружающие атомы

//void boundary(int* x, int* y, int* z); // изменить (x,y,z), если нужно, так чтобы попасть внутрь ящика
static inline void boundary(int &x, int &y, int &z) { // изменить (x,y,z), если нужно, так чтобы попасть внутрь ящика
  if ( unlikely(x<0) )   x+=Lx;
  if ( unlikely(x>=Lx) ) x-=Lx;
  if ( unlikely(y<0) )   y+=Ly;
  if ( unlikely(y>=Ly) ) y-=Ly;
  if ( unlikely((z<0) || (z>=Lz)) ) {
    fprintf(stderr,"Error!!! z out of range\n");
    exit(0);
  }
}


void v_ochered_Edef(int x, int y, int z); // поставить атом (x,y,z) в очередь на обновление Edef
void update_Edef(); // обновить Edef для атомов, поставленных в очередь
void fill_nb_type(int* nb_type, unsigned short int config_);
unsigned short int massiv_to_config(int* nb_type);
void fill_n1_n2_and_good_config();
void create_initial_structure(int *);
void load_initial_structure(char* filename);
void set_config(int x, int y, int z);
void add_atom_to_spisok(int x,int y,int z);
void remove_atom_from_spisok(int x,int y,int z);
void move_atom_in_spisok(int x_old, int y_old, int z_old, int x_new, int y_new, int z_new);
void check_Lxyz(void);

// из jump.c
int jump(int x, int y, int z, int dir, int* x2, int* y2, int* z2);
int deposition(int type_of_new_atom, int* x, int* y, int* z);
// Планирование/применение осаждения без использования глобального rand().
// Нужно для параллельной генерации событий.
int deposition_plan(int type_of_new_atom, uint32_t &seed, int* x, int* y, int* z);
int deposition_at(int type_of_new_atom, int x, int y, int z);
void set_jump_info(int x, int y, int z); // устанавливает jumps и n_jumps для данного атома
void calc_jump_probability(int x, int y, int z);

// из cycle.c
void main_loop(void); // бесконечный цикл, состоящий из шагов Монте-Карло

// из output.c
void save_grid(char* filename, float radius);
void read_defects(const char* filename);
void write_params_to_log(const char* filename);
void read_control(const char* filename);
void write_to_log(const char* filename);
void read_parameters(const char* filename);
void show_me(char* filename, char key);
