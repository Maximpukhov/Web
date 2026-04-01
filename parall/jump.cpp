#include "mc_growth.h"

#define N_spisok 1000

int spisok[N_spisok][3];
int n_spisok;

// ------------------------------------------------------------
// Прототипы
// ------------------------------------------------------------

void v_spisok(int x, int y, int z);
void delete_jump(int x, int y, int z, int dir);

// ------------------------------------------------------------
// Расчёт вероятности прыжка для атома (x,y,z)
// ------------------------------------------------------------

void calc_jump_probability(int x, int y, int z)
{
    int dir;
    float E_a;

    if (n_jumps(x, y, z) == 0)
    {
        jump_probability(x, y, z) = 0;
        return;
    }

    // проверка дефектов
    if (unlikely(
        atoms(x, y, z).defect_f.x != 0 ||
        atoms(x, y, z).defect_f.y != 0 ||
        atoms(x, y, z).defect_f.z != 0))
    {
        FILE* f = fopen("log.txt", "at");
        fprintf(f,
            "\tatom (%d,%d,%d) near defect is on the surface: (%hhd,%hhd,%hhd)\n",
            x, y, z,
            atoms(x, y, z).defect_f.x,
            atoms(x, y, z).defect_f.y,
            atoms(x, y, z).defect_f.z);
        fclose(f);

        erase_defect(x, y, z);
    }

    E_a = 0;

    neighbors_t nbs;
    atoms.neighbors(x, y, z, nbs);

    for (dir = 0; dir < dir_number; dir++)
    {
        if (atoms(nbs.x[dir], nbs.y[dir], nbs.z[dir]).type != 0)
        {
            if (dir < 4)
            {
                E_a += 0.5f * (
                    param.E1[atoms(nbs.x[dir], nbs.y[dir], nbs.z[dir]).type] +
                    param.E1[atoms(x, y, z).type]);
            }
            else
            {
                E_a += 0.5f * (
                    param.E2[atoms(nbs.x[dir], nbs.y[dir], nbs.z[dir]).type] +
                    param.E2[atoms(x, y, z).type]);
            }
        }
    }

    // учитываем деформацию
    E_a -= atoms(x, y, z).Edef;

    if (E_a < 0) E_a = 0;

    jump_probability(x, y, z) =
        param.p0 * exp(-E_a / param.T) * n_jumps(x, y, z);
}

// ------------------------------------------------------------
// Установка допустимых направлений прыжка
// ------------------------------------------------------------

void set_jump_info(int x, int y, int z)
{
    int dir;
    int j_array[dir_number];

    if (atoms(x, y, z).type == 0 ||
        n1_config[atoms(x, y, z).config] == 4 ||
        z < 2 || z >= Lz - 2)
    {
        n_jumps(x, y, z) = 0;
        jumps(x, y, z) = 0;
        return;
    }

    if (z <= param.z_surface) {
        n_jumps(x, y, z) = 0;
        jumps(x, y, z) = 0;
        return;
    }

    int n_jumps_ = 0;
    for (dir = 0; dir < dir_number; dir++)
        j_array[dir] = 0;

    neighbors_t nbs;
    atoms.neighbors(x, y, z, nbs);

    for (dir = 0; dir < dir_number; dir++)
    {
        if (atoms(nbs.x[dir], nbs.y[dir], nbs.z[dir]).type == 0 &&
            n1_config[atoms(nbs.x[dir], nbs.y[dir], nbs.z[dir]).config] > (dir < 4))
        {
            j_array[dir] = 1;
            n_jumps_++;
        }
    }

    n_jumps(x, y, z) = n_jumps_;
    jumps(x, y, z) = massiv_to_config(j_array);
}

// ------------------------------------------------------------
// Прыжок атома
// ------------------------------------------------------------

int jump(int x, int y, int z, int dir, int* x2_, int* y2_, int* z2_)
{
    int xd, yd, zd;
    int x2, y2, z2;
    int dir2, dir3;
    int n, n_spisok2;
    int bad_jump;
    int atom_erase;

    atoms.one_neighbor(x, y, z, dir, xd, yd, zd);
    *x2_ = xd;
    *y2_ = yd;
    *z2_ = zd;

    if (unlikely(atoms(x, y, z).type == 0))
    {
        fprintf(stderr, "Error!!! jump from empty site (%d,%d,%d)\n", x, y, z);
        exit(0);
    }

    if (unlikely(atoms(xd, yd, zd).type != 0))
    {
        fprintf(stderr, "Error!!! jump to occupied site\n");
        exit(0);
    }

    if (z < 4 || zd < 4 || z >= Lz - 4 || zd >= Lz - 4)
    {
        delete_jump(x, y, z, dir);
        calc_jump_probability(x, y, z);
        return 0;
    }

    // выполняем прыжок
    atoms(xd, yd, zd).type = atoms(x, y, z).type;
    atoms(x, y, z).type = 0;

    n_spisok = 0;
    v_spisok(x, y, z);
    v_spisok(xd, yd, zd);

    neighbors_t nbs, nbsd;
    atoms.neighbors(x, y, z, nbs);
    atoms.neighbors(xd, yd, zd, nbsd);

    for (dir2 = 0; dir2 < dir_number; dir2++)
    {
        v_spisok(nbs.x[dir2], nbs.y[dir2], nbs.z[dir2]);
        v_spisok(nbsd.x[dir2], nbsd.y[dir2], nbsd.z[dir2]);
    }

    // проверка корректности конфигураций
    bad_jump = 0;
    for (n = 0; n < n_spisok; n++)
    {
        x2 = spisok[n][0];
        y2 = spisok[n][1];
        z2 = spisok[n][2];

        set_config(x2, y2, z2);
        if (atoms(x2, y2, z2).type > 0 &&
            !good_config[atoms(x2, y2, z2).config])
        {
            bad_jump = 1;
            break;
        }
    }

    if (bad_jump)
    {
        // откат прыжка
        atoms(x, y, z).type = atoms(xd, yd, zd).type;
        atoms(xd, yd, zd).type = 0;

        for (n = 0; n < n_spisok; n++)
        {
            x2 = spisok[n][0];
            y2 = spisok[n][1];
            z2 = spisok[n][2];

            set_config(x2, y2, z2);
            axyz(x2, y2, z2);
            spisok_flag(x2, y2, z2) = 0;
        }

        delete_jump(x, y, z, dir);
        calc_jump_probability(x, y, z);
        return 0;
    }

    atom_erase = 0;
    if (zd < 6)
    {
        atoms(xd, yd, zd).type = 0;
        atom_erase = 1;
    }

    // обновление окружения
    for (n = 0; n < n_spisok; n++)
    {
        x2 = spisok[n][0];
        y2 = spisok[n][1];
        z2 = spisok[n][2];

        if (atoms(x2, y2, z2).type == 3)
            atoms(x2, y2, z2).type = 1;

        calc_B0(x2, y2, z2);
    }

    axyz(xd, yd, zd);

    n_spisok2 = n_spisok;
    for (n = 0; n < n_spisok2; n++)
    {
        x2 = spisok[n][0];
        y2 = spisok[n][1];
        z2 = spisok[n][2];

        atoms.neighbors(x2, y2, z2, nbs);
        for (dir3 = 0; dir3 < dir_number; dir3++)
            v_spisok(nbs.x[dir3], nbs.y[dir3], nbs.z[dir3]);
    }

    for (n = 0; n < n_spisok; n++)
    {
        int x3 = spisok[n][0];
        int y3 = spisok[n][1];
        int z3 = spisok[n][2];

        set_jump_info(x3, y3, z3);
        if (n_jumps(x3, y3, z3) == 0)
            calc_jump_probability(x3, y3, z3);

        v_ochered_Edef(x3, y3, z3);
        spisok_flag(x3, y3, z3) = 0;
    }

    if (atom_erase)
        move_atom_in_spisok(x, y, z, x, y, z);
    else
        move_atom_in_spisok(x, y, z, xd, yd, zd);

    return 1;
}

// ------------------------------------------------------------
// Осаждение (deposition)
//
// В текущем состоянии проекта прототипы deposition/deposition_plan/deposition_at
// были объявлены в mc_growth.h, но реализация отсутствовала (ошибка линковки).
// Реализация ниже сделана минимально инвазивной и совместимой с текущими
// структурами данных:
//  - deposition_plan(): детерминированно (по seed) выбирает кандидат (x,y,z)
//    вблизи поверхности.
//  - deposition_at(): пытается добавить атом в (x,y,z), проверяет корректность
//    конфигураций и обновляет локальное окружение (config/B0/axyz/jumps/Edef).
//
// Важно: это не «идеальная» физическая модель адсорбции, но она корректно
// поддерживает рост структуры и валидную выдачу в формате .xyz для OVITO.
// ------------------------------------------------------------

static inline void random_lattice_xy_for_z(int z, uint32_t &seed, int *x, int *y)
{
    // Сэмплируем только допустимые узлы алмазной решётки для фиксированного z.
    // x: z%2, z%2+2, ...
    int xp = (z & 1);
    int xi = rng_int(seed, Lx / 2);
    int xx = xp + 2 * xi;

    // y: z%2 + 2*((z/2 + x/2)%2) + 4*k
    int y0 = xp + 2 * (((z / 2 + xx / 2) & 1));
    int yi = rng_int(seed, Ly / 4);
    int yy = y0 + 4 * yi;

    *x = xx;
    *y = yy;
}

static inline bool is_reasonable_deposition_site(int x, int y, int z)
{
    if (z < 2 || z >= Lz - 2) return false;
    if (atoms(x, y, z).type != 0) return false;

    // Хотим осаждать на «поверхность»: у пустого узла должен быть хотя бы один
    // ближайший сосед (1-й) в кристалле.
    set_config(x, y, z);
    if (n1_config[atoms(x, y, z).config] <= 0) return false;

    return true;
}

int deposition_plan(int type_of_new_atom, uint32_t &seed, int *x, int *y, int *z)
{
    (void)type_of_new_atom;

    // Диапазон поиска возле поверхности. Берём небольшой «коридор», чтобы
    // поддерживать шероховатость, но не улетать в потолок.
    int z0 = param.z_surface + 1;
    if (z0 < 2) z0 = 2;

    int z1 = param.z_surface + 6;
    if (z1 > Lz - 3) z1 = Lz - 3;
    if (z1 < z0) z1 = z0;

    // Сначала пробуем случайными попытками (быстро).
    for (int it = 0; it < 300; it++)
    {
        int zz = z0 + rng_int(seed, (z1 - z0 + 1));
        int xx, yy;
        random_lattice_xy_for_z(zz, seed, &xx, &yy);

        if (is_reasonable_deposition_site(xx, yy, zz))
        {
            *x = xx;
            *y = yy;
            *z = zz;
            return 1;
        }
    }

    // Если случайные попытки не сработали — линейный поиск в самом нижнем
    // слое над подложкой.
    int zz = z0;
    int xx, yy;
    for (yy = (zz & 1); yy < Ly; yy += 2)
        for (xx = (zz & 1) + 2 * (((zz / 2 + yy / 2) & 1)); xx < Lx; xx += 4)
            if (is_reasonable_deposition_site(xx, yy, zz))
            {
                *x = xx;
                *y = yy;
                *z = zz;
                return 1;
            }

    return 0;
}

static inline void push_unique(std::vector<coord> &v, int x, int y, int z)
{
    for (const auto &c : v)
        if (c.x == x && c.y == y && c.z == z) return;
    v.push_back(coord{ x, y, z });
}

int deposition_at(int type_of_new_atom, int x, int y, int z)
{
    if (z < 2 || z >= Lz - 2) return 0;
    if (atoms(x, y, z).type != 0) return 0;

    // Временно ставим атом и проверяем корректность конфигураций локально.
    atoms(x, y, z).type = (char)type_of_new_atom;

    std::vector<coord> upd;
    upd.reserve(1 + dir_number);
    push_unique(upd, x, y, z);

    neighbors_t nbs;
    atoms.neighbors(x, y, z, nbs);
    for (int d = 0; d < dir_number; d++)
        push_unique(upd, nbs.x[d], nbs.y[d], nbs.z[d]);

    // Проверяем, что после осаждения ни один из затронутых атомов
    // не попадает в «плохую» конфигурацию.
    bool bad = false;
    for (const auto &c : upd)
    {
        set_config(c.x, c.y, c.z);
        if (atoms(c.x, c.y, c.z).type > 0 && !good_config[atoms(c.x, c.y, c.z).config])
        {
            bad = true;
            break;
        }
    }

    if (bad)
    {
        atoms(x, y, z).type = 0;
        // Восстановим config (на всякий случай) для окружения.
        for (const auto &c : upd)
            set_config(c.x, c.y, c.z);
        return 0;
    }

    // Добавляем атом в список атомов структуры.
    add_atom_to_spisok(x, y, z);

    // Обновляем окружение: B0/axyz/jumps/вероятности/Edef.
    for (const auto &c : upd)
    {
        if (atoms(c.x, c.y, c.z).type == 3)
            atoms(c.x, c.y, c.z).type = 1;

        calc_B0(c.x, c.y, c.z);

        if (atoms(c.x, c.y, c.z).type != 0)
            axyz(c.x, c.y, c.z);

        set_jump_info(c.x, c.y, c.z);
        if (n_jumps(c.x, c.y, c.z) == 0)
            jump_probability(c.x, c.y, c.z) = 0;
        else
            calc_jump_probability(c.x, c.y, c.z);

        v_ochered_Edef(c.x, c.y, c.z);
    }

    return 1;
}

int deposition(int type_of_new_atom, int *x, int *y, int *z)
{
    // Старый интерфейс: используется глобальный rand(). Оставляем для совместимости.
    // Пытаемся несколько раз найти позицию.
    uint32_t seed = (uint32_t)rand();
    for (int it = 0; it < 300; it++)
    {
        int xx, yy, zz;
        if (!deposition_plan(type_of_new_atom, seed, &xx, &yy, &zz))
            return 0;
        if (deposition_at(type_of_new_atom, xx, yy, zz))
        {
            *x = xx;
            *y = yy;
            *z = zz;
            return 1;
        }
    }
    return 0;
}

// ------------------------------------------------------------
// Добавление координаты в список обновлений
// ------------------------------------------------------------

void v_spisok(int x, int y, int z)
{
    if (spisok_flag(x, y, z)) return;

    if (unlikely(n_spisok >= N_spisok))
    {
        fprintf(stderr, "Error!!! N_spisok is too small\n");
        exit(0);
    }

    spisok[n_spisok][0] = x;
    spisok[n_spisok][1] = y;
    spisok[n_spisok][2] = z;
    n_spisok++;

    spisok_flag(x, y, z) = 1;
}

// ------------------------------------------------------------
// Удаление направления прыжка
// ------------------------------------------------------------

void delete_jump(int x, int y, int z, int dir)
{
    int d;
    int j_array[dir_number];

    fill_nb_type(j_array, jumps(x, y, z));
    j_array[dir] = 0;

    jumps(x, y, z) = massiv_to_config(j_array);

    int n_jumps_ = 0;
    for (d = 0; d < dir_number; d++)
        n_jumps_ += j_array[d];

    n_jumps(x, y, z) = n_jumps_;
}
