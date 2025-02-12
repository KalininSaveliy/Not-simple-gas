#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

// #define _USE_MATH_DEFINES  // to get PI value


// config value
bool isDebug;

int n_time;      // число итераций по времени
int n_x, n_y;    // размер секти по координате
int n_v;         // размер сетки по скоростям
double Knudsen;  // число Кнудсена
double Mach;     // число Маха
double T1;       // температура пластины сверху относительная
double T2;       // температура пластины снизу относительная
double v_cut;  // максимальная рассматривая скорость

/*
XY space:
|                            |
|          ________          | plate
|                            |
   alpha  |  beta  |  gamma    where (alpha, beta, gamma in (0.0, 1.0)) and (alpha + beta + gammma = 1)
 ---------0--------^-------> X in lambda where lambda - длина свободного пробега
                plate_len
*/
double alpha, beta;
double real_plate_len;  // длина пластины в метрах
// config value end

// extra value
int plate_beg_i, plate_end_i;  // исключая конец
int plate_y;  // координата пластины по Y
double h;  // шаг по координатной сетке
double dv;  // шаг по скоростной сетке
double tau;  // шаг по временной сетке
// extra value end


bool load_config(const std::string& config_name) {
    bool isLoaded = true;
    std::string full_cfg = "config/" + config_name + ".txt";
    std::ifstream fin;
    fin.open(full_cfg);
    if (fin.is_open()) {
        std::cout << "config: " << full_cfg << '\n';
        std::string line;
        char delim = '=';

        std::getline(fin, line, delim);
        fin >> n_time;
        std::getline(fin, line, delim);
        fin >> n_x >> n_y;
        std::getline(fin, line, delim);
        fin >> n_v >> v_cut;
        if (n_v % 2 == 0)
            n_v++;
        std::getline(fin, line, delim);
        fin >> alpha >> beta;
        std::getline(fin, line, delim);
        fin >> Knudsen >> Mach;
        std::getline(fin, line, delim);
        fin >> T1 >> T2;
        std::getline(fin, line, delim);
        fin >> real_plate_len;

        std::getline(fin, line, delim);
        fin >> isDebug;

        double plate_len = 1 / Knudsen;  // длина пластины в длинах свободного пробега
        // TODO: add check for alpha, beta, n_x value
        int n_alpha = std::floor(n_x * alpha);  // целое число ячеек сетки для области до пластины
        int n_beta = std::floor(n_x * beta);    // целое число ячеек сетки для пластины
        // int n_gamma = n_x - n_alpha - n_beta;   // целое число ячеек сетки для области после пластины

        plate_beg_i = n_alpha;
        plate_end_i = n_alpha + n_beta;
        plate_y = n_y / 2;
        h = plate_len / n_beta;
        dv = v_cut / (n_v / 2);
        tau = h / v_cut;

        if (isDebug) {
            std::cout << n_x << ' ' << n_y << '\n' << n_v << ' ' << v_cut << '\n'
                    << alpha << ' ' << beta << '\n' << Knudsen << ' ' << Mach << '\n'
                    << T1 << ' ' << T2 << ' ' << real_plate_len << '\n'
                    << h << ' ' << dv << ' ' << n_time << '\n';
        }
    } else {
        std::cout << "Error: Can not open config file with name: " << full_cfg << '\n';
        isLoaded = false;
    }
    fin.close();
    return isLoaded;
}

double coord(int i) {
    return (i - plate_beg_i) * h;
}
double speed(int i) {
    return (i - n_v / 2) * dv;
}

double Maxwell2(double vx, double vy, double v0_x, double v0_y, double n=1.0, double temp=1.0) {
    double c = n / (2 * M_PI * temp);
    return c * std::exp(-((vx - v0_x) * (vx - v0_x) + (vy - v0_y) * (vy - v0_y)) / (2 * temp));
}

// retrun 1dim index in distribution array
int idx(int x, int y, int vx, int vy) {
    return x + y * n_x + vx * (n_x * n_y) + vy * (n_x * n_y * n_v);
}
void initDistribution(double* f) {
    double n_0, n_bound;
    for (int ii = 0; ii < n_v; ++ii) {
        double vx = speed(ii);
        for (int jj = 0; jj < n_v; ++jj) {
            double vy = speed(jj);
            n_0 = Maxwell2(vx, vy, 0.0, 0.0);
            n_bound = Maxwell2(vx, vy, Mach, 0.0);
            for (int j = 0; j < n_y; ++j) {
                f[idx(0, j, ii, jj)] = n_bound;
                for (int i = 1; i < n_x; ++i) {
                    f[idx(i, j, ii, jj)] = n_0;
                }
            }
        }
    }
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
    filename += saveAll ? ".txt" : ".csv";
    std::ofstream fout;
    fout.open(filename);
    if (fout.is_open()) {
        if (saveAll) {
            fout << "# sum(f(vx,vy)) for each (x, y) point in relative units\n";
            fout << "# next line is the concentration n = sum of f(vx, vy) in point (x,y)\n";
            fout << n_v << ' ' << dv << ' ' << v_cut << "  # n_v, dv and max_velocity\n\n";
        } else {
            fout << "# x,y,sum(f(vx | vy)) in relative units\n";
        }
        double x, y, n, sum;
        for (int i = 0; i < n_x; ++i) {
            x = coord(i);
            for (int j = 0; j < n_y; ++j) {
                y = coord(j);
                sum = 0.0;
                fout << x << ',' << y;
                saveAll ? (fout << '\n') : (fout << ',');
                for (int ii = 0; ii < n_v; ++ii) {
                    for (int jj = 0; jj < n_v; ++jj) {
                        n = f[idx(i, j, ii, jj)];
                        sum += n;
                        if (saveAll) {fout << n << ' ';}
                    }
                    if (saveAll) {fout << '\n';}
                }
                fout << sum * dv * dv << '\n';  // \int f(vx, vy) dvx dvy
            }
        }
    } else {
        std::cout << "Error: Could not save distribution to file " << filename << '\n';
    }
    fout.close();
}

void make_iteration_x(double* f_old, double* f_new) {
    double g;
    for (int ii = 0; ii < n_v; ++ii) {
        g = speed(ii) * tau / h;
        for (int jj = 0; jj < n_v; ++jj) {
            for (int j = 0; j < n_y; ++j) {
                if (g > 0) {
                    f_new[idx(0, j, ii, jj)] = f_old[idx(0, j, ii, jj)];  // const value
                    for (int i = 1; i < n_x; ++i) {
                        f_new[idx(i, j, ii, jj)] = f_old[idx(i, j, ii, jj)] - g * (f_old[idx(i,   j, ii, jj)] -
                                                                                   f_old[idx(i-1, j, ii, jj)]);
                    }
                } else {
                    for (int i = 0; i < n_x - 1; ++i) {
                        f_new[idx(i, j, ii, jj)] = f_old[idx(i, j, ii, jj)] - g * (f_old[idx(i+1, j, ii, jj)] -
                                                                                   f_old[idx(i,   j, ii, jj)]);
                    }
                }
            }
        }
    }
}

void make_iteration_y(double* f_old, double* f_new) {
    double g;
    for (int jj = 0; jj < n_v; ++jj) {
        g = speed(jj) * tau / h;
        for (int ii = 0; ii < n_v; ++ii) {
            for (int i = 0; i < n_x; ++i) {
                if (g > 0) {
                    for (int j = 1; j < n_y; ++j) {
                        f_new[idx(i, j, ii, jj)] = f_old[idx(i, j, ii, jj)] - g * (f_old[idx(i, j,   ii, jj)] -
                                                                                   f_old[idx(i, j-1, ii, jj)]);
                    }
                } else {
                    for (int j = 0; j < n_y; ++j) {
                        f_new[idx(i, j, ii, jj)] = f_old[idx(i, j, ii, jj)] - g * (f_old[idx(i, j-1, ii, jj)] -
                                                                                   f_old[idx(i, j, ii, jj)]);
                    }
                }
            }
        }
    }
    // Diffusion reflection
}

void make_iteration(double* f_old, double* f_new) {
    make_iteration_x(f_old, f_new);
    std::swap(f_old, f_new);
    make_iteration_y(f_old, f_new);
}

int main(int argc, char** argv) {  // TODO: create dir to save file;
    std::string config_name;
    if (argc < 2) {
        std::cout << "Sorry, you forgot add config file\n";
        return 0;
    } else {
        config_name = argv[1];
    }
    if (!load_config(config_name))
        return 0;

    double* prev_distribution = new double[n_x * n_y * n_v * n_v];
    double* distribution      = new double[n_x * n_y * n_v * n_v];
    initDistribution(distribution);

    std::string folder = "data/" + config_name + '/';
    std::string file = folder + "0";  // TODO: maybe should add name of file to config
    save(distribution, file, isDebug);  // TODO: add time to save function

    for (int t = 1; t < n_time; ++t) {
        make_iteration(prev_distribution, distribution);
        if (t % 10 == 0) {  // TODO: maybe should add 10 to config
            file = folder + std::to_string(t);
            save(distribution, file, isDebug);
            std::cout << "Iter " << t << " / " << n_time << " was done.\n";
        }
        std::swap(prev_distribution, distribution);
    }

    return 0;
}
