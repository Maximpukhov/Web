#include <ctype.h>
#include "mc_growth.h"
#include <mpi.h>
#include "compress.h"
#include <chrono>
#include <cstring>
#include <vector>
#include <list>

//------- глобальные переменные -------------
int Lx=40, Ly=40, Lz=100;

// Переменные для MPI, которые будут использоваться в других файлах
int mpi_rank = 0;
int mpi_size = 1;

struct param_struct param;
struct current_struct current;

std::list<struct coord> ochered_Edef;
std::vector<struct coord> spisok_atomov;

dlattice_t <atom_t> atoms;
dlattice_t <float> jump_probability;
dlattice_t <char> n_jumps;
dlattice_t <unsigned short int> jumps;
dlattice_t <char> in_ochered_Edef;
dlattice_t <int> number_in_spisok_atomov;
dlattice_t <char> spisok_flag;
unsigned long sumtime = 0;

float AA[Nconfig][6];
float AA_[Nconfig][6];
float BB[Nconfig][dir_number][9];
float transform_array[Nconfig][6];
char good_config[Nconfig];
char n1_config[Nconfig];
char n2_config[Nconfig];

int n_defects;
int defect_x[N_defects], defect_y[N_defects], defect_z[N_defects];

ckpt_compressor comp;

//------- main -------------

int main(int argc, char** argv) {
  int z;

  // --- Инициализация MPI ---
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  if (mpi_rank == 0) {
      printf("Starting MPI simulation (Replicated State) on %d processes.\n", mpi_size);
  }

  // Обнуляем счетчики
  current.n_deposited=current.n_jumps=current.n_bad_jumps=current.n_jumps_back=current.n_moves=0;
  current.prev_x=current.prev_y=current.prev_z=-1;
  current.t=0;

  // Дефолтные параметры
  sprintf(param.load,"0");
  param.show_short = 0.1; param.show_long = 1.0; param.show_sum = 0.0;
  param.E1[0]=0; param.E2[0]=0;
  param.E1[1]=(0.35+0.35*0.25)*11800; param.E2[1]=(0.15+0.15*0.25)*11800;
  param.E1[2]=0.35*11800; param.E2[2]=0.15*11800;
  param.E1[3]=0.35*11800*10; param.E2[3]=0.15*11800*10;
  param.p0=1e13; param.p_deposition=0.1*Lx*Ly/8.0;
  param.T=670;
  param.A=48.5 * 1.36e-10*1.36e-10/1.6e-19 * 11800 / 8.0;
  param.B=13.8   * 1.36e-10*1.36e-10/1.6e-19 * 11800 / 8.0;
  param.eps[0]=0; param.eps[1]=0; param.eps[2]=0.04; param.eps[3]=0.1;
  param.defect_F= (5 * 1.5 * 11800) / sqrt(3.0);
  param.moves_percent=1; param.time_for_moves=1e-6;
  param.dML_control=0.01; param.experiment_ML=6;
  param.random_seed=19;
  param.run_mode=RUN_MODE_PERF;
  param.dep_type=2;
  param.z_surface=20; param.z_cap = -100;
  n_defects=-1; param.Edef_pa_N = 0;

  read_parameters("parameters.txt");


  srand(param.random_seed);

  if (mpi_rank == 0) {
      printf("Global Random Seed: %d (Sync Mode)\n", param.random_seed);
      printf("Lattice: %dx%dx%d\n", Lx, Ly, Lz);
      printf("Run Mode: %s\n", (param.run_mode == RUN_MODE_VALIDATION) ? "validation" : "perf");
  }

  atoms.tune(Lx,Ly,Lz);
  jump_probability.tune(Lx,Ly,Lz);
  n_jumps.tune(Lx,Ly,Lz);
  jumps.tune(Lx,Ly,Lz);
  in_ochered_Edef.tune(Lx,Ly,Lz);
  number_in_spisok_atomov.tune(Lx,Ly,Lz);
  spisok_flag.tune(Lx,Ly,Lz);

  fill_n1_n2_and_good_config();
  calc_AA(); calc_AA_(); calc_BB();

  if(strcmp(param.load,"0")==0){
    read_defects("defects.txt");
    int initial_layers[Lz];
    for (z=0; z<Lz; z++) {
      if (z<=param.z_surface) initial_layers[z]=1;
      else                    initial_layers[z]=0;
    }
    create_initial_structure(initial_layers);

    // Только Rank 0 пишет файлы
    if (mpi_rank == 0) {
        printf("create_initial_structure - OK.\n");
        char filename[100];
        sprintf(filename,"_initial.xyz");
        show_me(filename,true);
        comp.init("mc_sim_ckpt", 10);
        comp.notify(filename);
    }
  } else {
    load_initial_structure(param.load);
    if(mpi_rank==0) printf("load_initial_structure - OK.\n");
  }

  if (mpi_rank == 0) {
    FILE* logf = fopen("log.txt", "wt");
    if (logf == NULL) {
      fprintf(stderr, "!!! cannot open log file 'log.txt' for truncating\n");
    } else {
      fclose(logf);
    }
  }

  if(current.n_deposited==0 && mpi_rank == 0) write_params_to_log("log.txt");

  // Синхронизация перед стартом замера времени
  MPI_Barrier(MPI_COMM_WORLD);
  double start_time = MPI_Wtime();

  main_loop(); // ЗАПУСК

  MPI_Barrier(MPI_COMM_WORLD);
  double end_time = MPI_Wtime();

  if (mpi_rank == 0) {
    printf("Simulation finished.\n");
    printf("Total Execution Time: %f seconds\n", end_time - start_time);
    // ИСПРАВЛЕНИЕ: (long) перед переменной
    printf("Total Jumps: %ld\n", (long)current.n_jumps);
    comp.finish();
  }

  MPI_Finalize();
  return 0;
}
