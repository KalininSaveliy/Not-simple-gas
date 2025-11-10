#define _USE_MATH_DEFINES  // to get PI value
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
// #include <filesystem>
#include "cnpy.h"



// config value
bool isDebug;
std::string save_folder;

size_t n_time;      // число итераций по времени
size_t save_time;   // период сохранений
size_t n_x, n_y;    // размер секти по координате
size_t n_v;         // размер сетки по скоростям
double Knudsen;  // число Кнудсена
double Mach;     // число Маха (v0 / sqrt(kT/m))
double T1;       // температура пластины сверху относительная
double T2;       // температура пластины снизу относительная
double v_cut;  // максимальная рассматривая скорость

/*
XY(ij) space:
                      ^ Y
        j : n_y-1     |
          : 5         |
          : 4         |
plate_j   : 3 ________0_____plate_____,_______________> X in lambda where lambda - длина свободного пробега
plate_j-1 : 2         |
          : 1         |
          : 0   1   2 | 3   4   5   6   7   8   n_x-1
          """"""""""""^"""""""""""""""^""""""""""""""" i
          |   alpha   |     beta      |    gamma     |  where (alpha, beta, gamma in (0.0, 1.0)) and (alpha + beta + gammma = 1)
                          plate_len
        plate_beg_i = 3
        plate_end_i = 7
*/
double alpha, beta;
double real_plate_len;  // длина пластины в метрах
// config value end

// extra value
size_t plate_beg_i, plate_end_i;  // исключая конец
size_t plate_j;  // координата пластины по Y (пластина находится между plate_j и plate_j - 1)
double h;  // шаг по координатной сетке
double dv;  // шаг по скоростной сетке
double tau;  // шаг по временной сетке
double denom_up, denom_down;  // denominator(T)
// extra value end

// return coord x in metre
double real_x(size_t i) {
    return (0.5 + i - plate_beg_i) * h * Knudsen * real_plate_len;
}
double real_y(size_t j) {
    return (0.5 + j - plate_j) * h * Knudsen * real_plate_len;
}
double speed(size_t i) {
    return (i - n_v / 2.0) * dv;
}

double Maxwell2(double vx, double vy, double v0_x, double v0_y, double n=1.0, double temp=1.0) {
    double c = n / (2 * M_PI * temp);
    return c * std::exp(-((vx - v0_x) * (vx - v0_x) + (vy - v0_y) * (vy - v0_y)) / (2 * temp));
}

double denominator(double T) {
    double total = 0;
    double v;
    for (size_t ii = n_v / 2 + 1; ii < n_v - 1; ++ii) {
        v = speed(ii);
        total += v * std::exp(-v * v / (2 * T));
    }
    return total;
}

// retrun 1dim index in distribution array
int idx(size_t x, size_t y, size_t vx, size_t vy) {
    return x + y * n_x + vx * (n_x * n_y) + vy * (n_x * n_y * n_v);
}

void initDistribution(double* f) {
    double n_bound;
    for (size_t ii = 0; ii < n_v; ++ii) {
        double vx = speed(ii);
        for (size_t jj = 0; jj < n_v; ++jj) {
            double vy = speed(jj);
            // n_0 = Maxwell2(vx, vy, 0.0, 0.0);
            n_bound = Maxwell2(vx, vy, Mach, 0.0);
            for (size_t j = 0; j < n_y; ++j) {
                f[idx(0, j, ii, jj)] = n_bound;
                for (size_t i = 1; i < n_x; ++i) {
                    f[idx(i, j, ii, jj)] = n_bound;
                }
            }
        }
    }
}

void make_iteration_x(double* f_old, double* f_new) {
    // double g;
    #pragma omp parallel for
    for (size_t ii = 0; ii < n_v; ++ii) {
        double g = speed(ii) * tau / h;
        for (size_t jj = 0; jj < n_v; ++jj) {
            for (size_t j = 0; j < n_y; ++j) {
                f_new[idx(0, j, ii, jj)] = f_old[idx(0, j, ii, jj)];  // const value
                if (g > 0) {
                    for (size_t i = 1; i < n_x; ++i) {
                        f_new[idx(i, j, ii, jj)] = f_old[idx(i, j, ii, jj)] - g * (f_old[idx(i,   j, ii, jj)] -
                                                                                   f_old[idx(i-1, j, ii, jj)]);
                    }
                } else {
                    f_new[idx(n_x - 1, j, ii, jj)] = f_old[idx(n_x - 1, j, ii, jj)];
                    for (size_t i = 1; i < n_x - 1; ++i) {
                        f_new[idx(i, j, ii, jj)] = f_old[idx(i, j, ii, jj)] - g * (f_old[idx(i+1, j, ii, jj)] -
                                                                                   f_old[idx(i,   j, ii, jj)]);
                    }
                }
            }
        }
    }
}

void make_iteration_y(double* f_old, double* f_new) {
    // double g;
    #pragma omp parallel for
    for (size_t jj = 0; jj < n_v; ++jj) {
        double g = speed(jj) * tau / h;
        for (size_t ii = 0; ii < n_v; ++ii) {
            for (size_t j = 0; j < n_y - 1; ++j) {
                f_new[idx(0, j, ii, jj)] = f_old[idx(0, j, ii, jj)];  // const value
            }
            for (size_t i = 1; i < n_x; ++i) {
                if (g > 0) {
                    f_new[idx(i, 0, ii, jj)] = f_old[idx(i, 0, ii, jj)];
                    for (size_t j = 1; j < n_y; ++j) {
                        f_new[idx(i, j, ii, jj)] = f_old[idx(i, j, ii, jj)] - g * (f_old[idx(i, j,   ii, jj)] -
                                                                                   f_old[idx(i, j-1, ii, jj)]);
                    }
                } else {
                    f_new[idx(i, n_y - 1, ii, jj)] = f_old[idx(i, n_y - 1, ii, jj)];
                    for (size_t j = 0; j < n_y - 1; ++j) {
                        f_new[idx(i, j, ii, jj)] = f_old[idx(i, j, ii, jj)] - g * (f_old[idx(i, j+1, ii, jj)] -
                                                                                   f_old[idx(i, j,   ii, jj)]);
                    }
                }
            }
        }
    }
    // Diffusion reflection
    double nom_up, nom_down;  // up - above plate, down - below plate
    double v_y;
    double liftForce = 0.0;
    for (size_t i = plate_beg_i; i < plate_end_i; ++i) {
        for (size_t ii = 0; ii < n_v; ++ii) {
            nom_down = 0;
            nom_up = 0;
            for (size_t jj = 0; jj < n_v; ++jj) {
                v_y = speed(jj);
                if (v_y > 0)
                    nom_down += v_y * f_old[idx(i, plate_j - 1, ii, jj)];
                else
                    nom_up   -= v_y * f_old[idx(i, plate_j, ii, jj)];  // (v_y <= 0) => ( -= )
            }
            liftForce += nom_down - nom_up;  // поток, для получения силы надо *2 / tau
            for (size_t jj = 0; jj < n_v; ++jj) {
                v_y = speed(jj);
                if (v_y > 0)
                    f_new[idx(i, plate_j,     ii, jj)] = nom_up   / denom_up   * std::exp(-v_y * v_y / (2 * T1));
                else if (v_y < 0)
                    f_new[idx(i, plate_j - 1, ii, jj)] = nom_down / denom_down * std::exp(-v_y * v_y / (2 * T2));
            }
        }
    }
    std::ofstream fout("data/main/Lift.txt", std::ios_base::app);
    if (fout.is_open())
        fout << liftForce << '\n';
    fout.close();
}

void make_iteration(double* f_old, double* f_new) {
    make_iteration_x(f_old, f_new);
    std::swap(f_old, f_new);
    make_iteration_y(f_old, f_new);
}


bool load_config(const std::string& config_name) {
    bool isLoaded = false;
    std::ifstream fin;
    std::cout << "I'm here\n";
    fin.open(config_name);
    if (fin.is_open()) {
        fin >> n_time
            >> save_time
            >> n_x
            >> n_y
            >> n_v
            >> v_cut
            >> alpha
            >> beta
            >> Knudsen
            >> Mach
            >> T1
            >> T2
            >> real_plate_len
            >> isDebug
            >> save_folder;

        double plate_len = 1 / Knudsen;  // длина пластины в длинах свободного пробега
        // TODO: add check for alpha, beta, n_x value
        int n_alpha = std::floor(n_x * alpha);  // целое число ячеек сетки для области до пластины
        int n_beta  = std::floor(n_x * beta);    // целое число ячеек сетки для пластины
        // int n_gamma = n_x - n_alpha - n_beta;   // целое число ячеек сетки для области после пластины

        plate_beg_i = n_alpha;
        plate_end_i = n_alpha + n_beta;
        plate_j     = n_y / 2;
        h           = plate_len / n_beta;
        dv          = v_cut / (n_v / 2);
        tau         = h / v_cut;
        denom_up    = denominator(T1);
        denom_down  = denominator(T2);

        if (isDebug) {
            std::cout
                << n_x << ' ' << n_y << '\n' << n_v << ' ' << v_cut << '\n'
                << alpha << ' ' << beta << '\n' << Knudsen << ' ' << Mach << '\n'
                << T1 << ' ' << T2 << ' ' << real_plate_len << '\n'
                << plate_beg_i << ' ' << plate_end_i << ' ' << plate_j << '\n'
                << h << ' ' << dv << ' ' << n_time << '\n';
        }
        isLoaded = true;
        fin.close();
    } else {
        std::cout << "Error: Can not open config file with name: " << std::quoted(config_name) << '\n';
    }
    return isLoaded;
}

bool save_grid(const std::string& folder) {
    bool isSaved = false;
    std::ofstream fout(folder + "grid.csv");
    // std::ofstream fout(folder + "grid.csv", std::ios::out | std::ios::trunc);
    if (fout.is_open()) {
        fout << "x,y\n";
        size_t n = std::max(n_x, n_y);
        for (size_t i = 0; i < n; ++i) {
            if (i < n_x) {fout << real_x(i);}
            // (i < n_x) ? (fout << real_x(i)) : (fout << " ");
            fout << ',';
            if (i < n_y) {fout << real_y(i);}
            fout << '\n';
        }
        fout.close();
        isSaved = true;
    } else {
        std::cout << "Error: can't open file to save grid << \n";
    }
    return isSaved;
}

/** @brief Save distribuiton f(x, y, vx, vy) to file
 * 
 * @param f is a distribution to save
 * @param filename name of file to save without extension
 * @param saveAll if true save all distribution. This work slower.
 *                if false save only concentration in each point (x, y). This work faster.
 * @return Void
 */
void save(double* f, std::string& filename, bool saveAll=false) {
    std::cout << "pre_save\n";
    cnpy::npy_save(filename + ".npy", f, {n_v, n_v, n_y, n_x}, "w");
    std::cout << "post_save\n\n";
    // filename += saveAll ? ".txt" : ".csv";
    // std::ofstream fout;
    // fout.open(filename);
    // if (fout.is_open()) {
    //     if (saveAll) {
    //         fout << "# sum(f(vx,vy)) for each (x, y) point in relative units\n";
    //         fout << "# next line is the concentration n = sum of f(vx, vy) in point (x,y)\n";
    //         fout << n_v << ' ' << dv << ' ' << v_cut << "  # n_v, dv and max_velocity\n\n";
    //     }
    //     double n, sum;
    //     for (int j = 0; j < n_y; ++j) {
    //         for (int i = 0; i < n_x; ++i) {
    //             sum = 0.0;
    //             if (saveAll) {fout << '\n';}
    //             for (int ii = 0; ii < n_v; ++ii) {
    //                 for (int jj = 0; jj < n_v; ++jj) {
    //                     n = f[idx(i, j, ii, jj)];
    //                     sum += n;
    //                     if (saveAll) {fout << n << ' ';}
    //                 }
    //                 if (saveAll) {fout << '\n';}
    //             }
    //             fout << sum * dv * dv;  // \int f(vx, vy) dvx dvy
    //             if (i != n_x - 1) {fout << ',';}
    //         }
    //         if (!saveAll) {fout << '\n';}
    //     }
    // } else {
    //     std::cout << "Error: Could not save distribution to file " << std::quoted(filename) << "\n";
    // }
    // fout.close();
}

bool make_simulation(const std::string& config_name) {
    if (!load_config(config_name))
        return false;

    std::string folder = "./data/" + save_folder + '/';
    // std::filesystem::create_directories(folder);
    if (!save_grid(folder))
        return false;
    
    double* distribute_buffer = new double[n_x * n_y * n_v * n_v];
    double* distribution      = new double[n_x * n_y * n_v * n_v];
    initDistribution(distribution);

    std::string file = folder + "0";  // TODO: maybe should add name of file to config
    save(distribution, file, isDebug);  // TODO: add time to save function

    for (size_t t_i = 1; t_i <= n_time; ++t_i) {
        make_iteration(distribution, distribute_buffer);
        if (t_i % save_time == 0) {
            file = folder + std::to_string(t_i);
            save(distribution, file, isDebug);
            std::cout << "Iter " << t_i << " / " << n_time << " was done.\n";
        }
        // std::swap(distribution, distribute_buffer);
    }
    delete[] distribution;
    delete[] distribute_buffer;
    return true;
}


int main(int argc, char** argv) {  // TODO: create dir to save file
    if (argc < 2) {
        std::cout << "Sorry, you forgot add config file\n";
        return 0;
    }
    std::string config_name = argv[1];

    if (make_simulation(config_name)) {
        std::cout << "Simulation " << std::quoted(config_name)<< " was done\n";
    } else {
        std::cout << "Error: Simultaion " << std::quoted(config_name) << " end with false\n";
    }

    return 0;
}
