#include "mc_growth.h"
#include <mpi.h>
#include <vector>
#include <deque>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <utility>

extern int mpi_rank;
extern int mpi_size;

struct task_t {
    int id;
    int select;
    int x, y, z;
    uint32_t seed;
    double dt;
};

struct result_t {
    int id;
    int valid;
    int dir;
    int src_x, src_y, src_z;
    int dst_x, dst_y, dst_z;
};

static task_t make_empty_task(int id) {
    task_t task = {};
    task.id = id;
    task.select = 0;
    task.x = task.y = task.z = -1;
    task.seed = 0;
    task.dt = 0.0;
    return task;
}

static bool same_site(int ax, int ay, int az, int bx, int by, int bz) {
    return ax == bx && ay == by && az == bz;
}

static bool source_already_selected(const std::vector<task_t>& tasks, int upto, int x, int y, int z) {
    for (int i = 0; i < upto; i++) {
        if (tasks[i].select == 2 && tasks[i].x >= 0 && same_site(tasks[i].x, tasks[i].y, tasks[i].z, x, y, z)) {
            return true;
        }
    }
    return false;
}

void mc_step(void);
double calc_P_jump_sum(void);
void do_many_axyz(void);

static int choose_dir_readonly(int x, int y, int z, uint32_t seed) {
    int j_array[dir_number];
    fill_nb_type(j_array, jumps(x, y, z));
    int nj = n_jumps(x, y, z);
    if (nj <= 0) return -1;

    int nx = 1 + rng_int(seed, nj);
    int n = 0;
    for (int dir = 0; dir < dir_number; dir++) {
        n += j_array[dir];
        if (n == nx) return dir;
    }
    return -1;
}

static result_t evaluate_task_readonly(const task_t& task, int id) {
    result_t res = {};
    res.id = id;
    res.valid = 0;
    res.src_x = task.x;
    res.src_y = task.y;
    res.src_z = task.z;

    if (task.select == 1 && task.x >= 0) {
        res.valid = 1;
        res.dst_x = task.x;
        res.dst_y = task.y;
        res.dst_z = task.z;
        return res;
    }

    if (task.select == 2 && task.x >= 0) {
        uint32_t seed = task.seed;
        int dir = choose_dir_readonly(task.x, task.y, task.z, seed);
        if (dir >= 0) {
            int x2, y2, z2;
            atoms.one_neighbor(task.x, task.y, task.z, dir, x2, y2, z2);
            res.valid = 1;
            res.dir = dir;
            res.dst_x = x2;
            res.dst_y = y2;
            res.dst_z = z2;
        }
    }

    return res;
}

static bool check_conflict(const result_t& a, const result_t& b) {
    if (a.valid == 0 || b.valid == 0) return false;

    if (a.dst_x == b.dst_x && a.dst_y == b.dst_y && a.dst_z == b.dst_z) return true;

    int d_src = abs(a.src_x - b.src_x) + abs(a.src_y - b.src_y) + abs(a.src_z - b.src_z);
    if (d_src <= 4) return true;

    int d_dst = abs(a.dst_x - b.dst_x) + abs(a.dst_y - b.dst_y) + abs(a.dst_z - b.dst_z);
    if (d_dst <= 4) return true;

    if (a.dst_x == b.src_x && a.dst_y == b.src_y && a.dst_z == b.src_z) return true;
    if (b.dst_x == a.src_x && b.dst_y == a.src_y && b.dst_z == a.src_z) return true;

    return false;
}

void main_loop(void) {
    float ML_next_control;
    char filename[100];
    float ML_deposited;
    double t_local;

    ML_deposited = current.n_deposited / (Lx * Ly / 8.0f);
    ML_next_control = param.dML_control * ((int)(ML_deposited / param.dML_control));

    if ((ML_deposited == 0) && param.experiment_ML != 0) {
        int x, y, z;
        for (int n = 0; n < 1000; n++) {
            std::random_shuffle(spisok_atomov.begin(), spisok_atomov.end());
            for (int i = 0; i < (int)spisok_atomov.size(); i++) {
                coord c = spisok_atomov[i];
                axyz(c.x, c.y, c.z);
            }
        }
        ZXY v_ochered_Edef(x, y, z);
        update_Edef();
    }

    t_local = current.t;

    while (1) {
        if (ML_deposited >= ML_next_control) {
            if (mpi_rank == 0) {
                read_control("control.txt");
                write_to_log("log.txt");
                sprintf(filename, "conf_ML%.2f.xyz", ML_deposited);
                show_me(filename, 2);
            }
            ML_next_control += param.dML_control;
        }

        if (ML_deposited >= param.experiment_ML) break;

        mc_step();

        if (current.t - t_local >= param.time_for_moves) {
            do_many_axyz();
            t_local = current.t;
        }

        ML_deposited = current.n_deposited / (Lx * Ly / 8.0f);
    }
}

void mc_step(void) {
    static unsigned long step_id = 0;
    static std::deque<task_t> deferred_tasks;
    step_id++;

    update_Edef();

    double P_jump_sum = calc_P_jump_sum();
    double P_total = P_jump_sum + param.p_deposition;
    if (unlikely(P_total <= 0)) return;

    int n_events = (param.run_mode == RUN_MODE_VALIDATION) ? 1 : mpi_size;
    if (n_events <= 0) n_events = 1;

    std::vector<task_t> tasks(n_events);
    std::vector<result_t> results(n_events);

    if (mpi_rank == 0) {
        for (int i = 0; i < n_events; i++) {
            tasks[i] = make_empty_task(i);

            if (mpi_rank == 0 && param.run_mode == RUN_MODE_PERF && !deferred_tasks.empty()) {
                tasks[i] = deferred_tasks.front();
                deferred_tasks.pop_front();
                tasks[i].id = i;
                continue;
            }

            tasks[i].seed = param.random_seed + step_id * 1000 + i;

            double u = rng01(tasks[i].seed);
            if (u <= 0) u = 1e-12;
            tasks[i].dt = (-log(u) / P_total);

            double Px = P_total * rng01(tasks[i].seed);
            if (Px < param.p_deposition) {
                tasks[i].select = 1;
                if (!deposition_plan(param.dep_type, tasks[i].seed, &tasks[i].x, &tasks[i].y, &tasks[i].z)) {
                    tasks[i] = make_empty_task(i);
                }
            } else {
                tasks[i].select = 2;
            }
        }

        for (int i = 0; i < n_events; i++) {
            if (tasks[i].select == 2 && tasks[i].x < 0) {
                double available_sum = 0.0;
                int x, y, z;

                ZYX {
                    double pnode = jump_probability(x, y, z);
                    if (pnode <= 0) continue;
                    if (source_already_selected(tasks, i, x, y, z)) continue;
                    available_sum += pnode;
                }

                if (available_sum <= 0.0) {
                    tasks[i] = make_empty_task(i);
                    continue;
                }

                double Px = available_sum * rng01(tasks[i].seed);
                double P = 0.0;
                bool found = false;

                ZYX {
                    double pnode = jump_probability(x, y, z);
                    if (pnode <= 0) continue;
                    if (source_already_selected(tasks, i, x, y, z)) continue;
                    P += pnode;
                    if (P > Px) {
                        tasks[i].x = x;
                        tasks[i].y = y;
                        tasks[i].z = z;
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    tasks[i] = make_empty_task(i);
                }
            }
        }
    }

    MPI_Bcast(tasks.data(), n_events * (int)sizeof(task_t), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (param.run_mode == RUN_MODE_VALIDATION) {
        result_t res = {};
        if (mpi_rank == 0) {
            res = evaluate_task_readonly(tasks[0], 0);
        }
        MPI_Bcast(&res, (int)sizeof(result_t), MPI_BYTE, 0, MPI_COMM_WORLD);
        results[0] = res;
    } else {
        int my_idx = mpi_rank;
        result_t my_res = {};
        if (my_idx >= 0 && my_idx < n_events) {
            my_res = evaluate_task_readonly(tasks[my_idx], my_idx);
        }
        MPI_Allgather(&my_res, (int)sizeof(result_t), MPI_BYTE,
                      results.data(), (int)sizeof(result_t), MPI_BYTE, MPI_COMM_WORLD);
    }

    std::vector<bool> executed(n_events, false);
    double dt_sum = 0.0;
    int executed_count = 0;

    for (int i = 0; i < n_events; i++) {
        if (results[i].valid == 0) {
            if (tasks[i].select == 2) current.n_bad_jumps++;
            continue;
        }

        bool conflict = false;
        for (int j = 0; j < i; j++) {
            if (executed[j] && check_conflict(results[i], results[j])) {
                conflict = true;
                break;
            }
        }
        if (conflict) {
            if (mpi_rank == 0 && param.run_mode == RUN_MODE_PERF) {
                deferred_tasks.push_back(tasks[i]);
            }
            continue;
        }

        if (tasks[i].select == 1) {
            if (atoms(results[i].dst_x, results[i].dst_y, results[i].dst_z).type != 0) {
                if (mpi_rank == 0 && param.run_mode == RUN_MODE_PERF) {
                    deferred_tasks.push_back(tasks[i]);
                }
                continue;
            }
            if (deposition_at(param.dep_type, results[i].dst_x, results[i].dst_y, results[i].dst_z)) {
                current.n_deposited++;
                executed[i] = true;
                dt_sum += tasks[i].dt;
                executed_count++;
            }
        } else if (tasks[i].select == 2) {
            if (atoms(results[i].dst_x, results[i].dst_y, results[i].dst_z).type != 0) {
                if (mpi_rank == 0 && param.run_mode == RUN_MODE_PERF) {
                    deferred_tasks.push_back(tasks[i]);
                } else {
                    current.n_bad_jumps++;
                }
                continue;
            }
            int x2, y2, z2;
            int res = jump(results[i].src_x, results[i].src_y, results[i].src_z, results[i].dir, &x2, &y2, &z2);
            if (res == 1) {
                current.n_jumps++;
                executed[i] = true;
                dt_sum += tasks[i].dt;
                executed_count++;
            } else {
                current.n_bad_jumps++;
            }
        }
    }

    if (executed_count > 0) {
        current.t += dt_sum / executed_count;
    }
}

double calc_P_jump_sum(void) {
    double real_sum = 0;
    int x, y, z;
    for (z = 2; z < Lz - 2; z++)
        for (y = (z % 2); y < Ly; y += 2)
            for (x = (z % 2) + 2 * (((z / 2 + y / 2) % 2)); x < Lx; x += 4)
                real_sum += jump_probability(x, y, z);
    return real_sum;
}

void do_many_axyz(void) {
    int I = (int)(param.moves_percent / 100.0 * spisok_atomov.size());
    for (int i = 0; i < I; i++) {
        int n = random_(spisok_atomov.size());
        coord c = spisok_atomov[n];
        axyz(c.x, c.y, c.z);
        current.n_moves++;
    }
}
