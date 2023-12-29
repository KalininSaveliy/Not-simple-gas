// #include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <tuple>
#include <unordered_map>  // TODO: удалить после дебага

#include "source/KorobovGrid.h"
#include "source/Vec3d.h"

using v_t = vec3d_t<double>;
using ind_t = vec3d_t<int>;


std::vector<int> neg_f_table; // TODO: 

const int space = 8;
const int n_vx(12), n_vy(12), n_vz(12);

const double v_cut = 5.0;
const double dv_x = 2 * v_cut / n_vx;
const double dv_y = 2 * v_cut / n_vy;
const double dv_z = 2 * v_cut / n_vz;
const double dv_xyz[3] = {dv_x, dv_y, dv_z};
const double time = 20;
const double tau = 0.02;
const double square_ef = 1.0;  // эффективаная площадь сечения (для "тведых шариков" - d^2)
const double f_max = 1 / (2 * M_PI) / std::sqrt(2 * M_PI);


double calc_v(int i, double dv) {
    return -v_cut + (i + 0.5) * dv;
}

v_t calc_v(const ind_t& ind) {
        v_t v;
        v.x = calc_v(ind.x, dv_x);
        v.y = calc_v(ind.y, dv_y);
        v.z = calc_v(ind.z, dv_z);
        return v;
}

int find_n0() {  // число узлов внутри сферы радиуса v_cut
    int n0 = 0;
    for (int i = 0; i < n_vx; ++i) {
        for (int j = 0; j < n_vy; ++j) {
            for (int k = 0; k < n_vz; ++k) {
                if (calc_v({i, j, k}).norm() < v_cut) {
                    n0++;
                }
            }
        }
    }
    return n0;
}

int next_prime(int n) {  // if n is prime then func return n
    int size = 2 * n;
    bool is_prime[size];
    for (int i = 0; i < size; ++i)
        is_prime[i] = true;

    for (int i = 2; i * i < size; ++i) {
        if (is_prime[i]) {
            for (int j = i * i; j < size; j += i) {
                is_prime[j] = false;
            }
        }
    }
    int next_n = n;
    while (!is_prime[next_n]) {
        next_n++;
    }
    return next_n;
}

const int n0 = find_n0();
int n_col;
// const int n_col = next_prime(std::ceil(4 * tau * n0 * v_cut * square_ef * n0 * dv_x * dv_y * dv_z * f_max / 2 / std::sqrt(2)));


struct Distribution {
    static constexpr double coef = 0.063493635934240969786;  // (2*pi)^(-3/2) = 0.063493635934240969786
    std::array<std::array<std::array<double, n_vz>, n_vy>, n_vz> f;

    // double (*)(v_t) fixed_temperature(double (*func_with_T)(v_t, double), double temp) {
    //     return *func_with_T(v_t, temp);
    // }

    static double Maxwell (const v_t& v, double temp = 1.0) {
        return  coef * std::exp(v.squared() / (-2.0 * temp)) / sqrt(std::pow(temp, 3));
    }

    static double TwoGauss (const v_t& v, double temp = 1.0) {  // TODO: сделать v_t u параметром
        v_t u(2, 0, 0);
        return coef * (0.5 * (std::exp((v - u).squared() / (-2.0 * temp)) + std::exp((v + u).squared() / (-2.0 * temp)))) / sqrt(std::pow(temp, 3));
    }

    Distribution(double (*VelocityDistFunc)(const v_t&, double), double temp) {  // TODO: think about T as a second variable
        for (int i = 0; i < n_vx; ++i) {
            for (int j = 0; j < n_vy; ++j) {
                for (int k = 0; k < n_vz; ++k) {
                    v_t v = calc_v({i, j, k});
                    if (v.norm() < v_cut) {
                        f[i][j][k] = VelocityDistFunc(v, temp);
                    } else {
                        f[i][j][k] = 0;
                    }
                }
            }
        }
    }

    Distribution() : Distribution(Maxwell, 1.0) {}

    Distribution(double temp) : Distribution(Maxwell, temp) {}

    double& operator() (const ind_t& ind) {
        return f[ind.x][ind.y][ind.z];
    }

    double moment_2(int alpha, int beta) const {  // ind = 1, 2 or 3
        double dv_alpha = dv_xyz[alpha - 1];
        double dv_beta = dv_xyz[beta - 1];
        double m = 0.0;

        int ijk[3];  // TODO: hard constant 3

        for (int i = 0; i < n_vx; ++i) {
            ijk[0] = i;
            for (int j = 0; j < n_vy; ++j) {
                ijk[1] = j;
                for (int k = 0; k < n_vz; ++k) {
                    ijk[2] = k;
                    m += calc_v(ijk[alpha - 1], dv_alpha) * calc_v(ijk[beta - 1], dv_beta) * f[i][j][k];
                }
            }
        }
        return m * dv_x * dv_y * dv_z;
    }

    double get_concentration() const {
        double n = 0.0;
        for (auto& col : f) {
            for (auto& row : col) {
                for (auto& el : row) {
                    n += el;
                }
            }
        }
        return n * dv_x * dv_y * dv_z;
    }

    v_t get_momentum() const {
        v_t momentum(0, 0, 0);
        for (int i = 0; i < n_vx; ++i) {
            for (int j = 0; j < n_vy; ++j) {
                for (int k = 0; k < n_vz; ++k) {
                    momentum += calc_v({i, j, k}) * f[i][j][k];
                }
            }
        }
        return momentum * dv_x * dv_y * dv_z;
    }

    double get_energy() const {
        double energy = 0.0;
        for (int i = 0; i < n_vx; ++i) {
            for (int j = 0; j < n_vy; ++j) {
                for (int k = 0; k < n_vz; ++k) {
                    energy += calc_v({i, j, k}).squared() * f[i][j][k];
                }
            }
        }
        return energy * dv_x * dv_y * dv_z;
    }

    double get_temp() const {
        return this->get_energy() / this->get_concentration() / 3;
    }

    double max_difference_from_Maxwell() const {
        double temp = this->get_temp();
        double max_eps = Maxwell(calc_v({0, 0, n_vz - 1})) - f[0][0][n_vz - 1];
        double eps;
        double maxwell_val;
        for (int i = 0; i < n_vx; ++i) {
            for (int j = 0; j < n_vy; ++j) {
                for (int k = 0; k < n_vz - 1; ++k) {
                    if (f[i][j][k] > 0) {
                        maxwell_val = Maxwell(calc_v({i, j, k}), temp);
                        eps = std::abs(maxwell_val - f[i][j][k]) / maxwell_val;
                        if (eps > max_eps)
                            max_eps = eps;
                    }
                }
            }
        }
        return max_eps;
    }

    void save_to_file(const std::string& filename) const {
        double conc = this->get_concentration();
        double temp = this->get_temp();
        double ener = this->get_energy();
        v_t moment = this->get_momentum();
        std::cout << "n = " << conc << " , T = " << temp << " , E = " << ener << " , vec_p = " << moment << '\n';

        std::ofstream out(filename);
        out << "# f(v_x, v_y, v_z), v_z = const\n";
        out << "# nv_x = " << n_vx << ", nv_y = " << n_vy << ", nv_z = " << n_vz << '\n';
        out << "# n = " << conc << " , T = " << temp << " , E = " << ener << " , vec_p = " << moment << '\n';

        for (int i = 0; i < n_vx; ++i) {
            // for (int j = 0; j < n_vy; ++j) {
                // for (int k = 0; k < n_vz; ++k) {
                    int j = n_vy / 2;
                    int k = n_vz / 2;  
                    v_t v = calc_v({i, j, k});  
                    out << v.x << " " << f[i][j][k] << std::endl;
                    // out << v.x << " " << v.y << " " << f[i][j][k] << std::endl;
                    // out << v.x << " " << v.y << " "  << v.z << " " << f[i][j][k] << std::endl;
                // }
                // out << std::endl;
            // }
            // out << std::endl;
        }
        out.close();
    }
};

struct collision_t {  // TODO: полное говно, прям говнище
    ind_t ind_1;
    ind_t ind_2;
    double s;
    double e;
    bool isGood;

    private:
        int find_near_ind(double v, double dv) {
            int i = std::lround((v + v_cut) / dv - 0.5);
            // if (i < 0) {i = 0;}
            return i;
        }

        ind_t find_near_ind(const v_t& vec) {
            ind_t ind;
            ind.x = find_near_ind(vec.x, dv_x);
            ind.y = find_near_ind(vec.y, dv_y);
            ind.z = find_near_ind(vec.z, dv_z);
            return ind;
        }

        ind_t calc_allowed_ind(double x, double y, double z) {
            v_t v = v_t(x, y, z) * 2 * v_cut - v_t(v_cut);
            ind_t ind = find_near_ind(v);
            return ind;
        }

        double get_relative_energy(v_t vec, v_t v_centre) {
            return (vec - v_centre).squared();
        }

        int find_min_ind(double v, double dv) {
            int i = std::floor((v + v_cut) / dv - 0.5);  // TODO: устранить повтор в find_near_i
            // if (i < 0) {i = 0;}
            return i;
        }

        ind_t find_min_ind(const v_t vec) {
            ind_t min_ind;
            min_ind.x = find_min_ind(vec.x, dv_x);
            min_ind.y = find_min_ind(vec.y, dv_y);
            min_ind.z = find_min_ind(vec.z, dv_z);
            return min_ind;
        }

        ind_t find_add_ind(v_t u1, v_t v_cm, ind_t sum_ind, double energy_0, bool isEnergy0_bigger) {
            ind_t centre = find_min_ind(u1);  // cenetre index
            ind_t best_ind;
            double big_number = 12.0 * v_cut * v_cut;
            double min = big_number;

            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int k = 0; k < 2; k++) {
                        ind_t cur_ind = centre + ind_t(i, j, k);
                        ind_t paired_ind = sum_ind - cur_ind;
                        v_t cur_vec = calc_v(cur_ind);
                        v_t paired_vec = calc_v(paired_ind);
                        if ((cur_vec.squared() < v_cut * v_cut) && (paired_vec.squared() < v_cut * v_cut)) {
                            double cur_energy = get_relative_energy(cur_vec, v_cm);
                            // want: (E_0 > E and isEnergy0_bigger = true) or (E_0 < E and isEnergy0_bigger = false)
                            if (((energy_0 > cur_energy) == isEnergy0_bigger) && (energy_0 != cur_energy)) {
                                double distance = (u1 - cur_vec).squared();
                                if (distance < min) { // хотим найти min
                                    min = distance;
                                    best_ind = cur_ind;
                                }
                            }
                        }
                    }
                }
            }
            if (min == big_number) {
                this->isGood = false;
            }
            return best_ind;
        }

        double find_r(double energy_0, ind_t lam_mu, ind_t lam_mu_add, v_t v_cm) {
            double energy_1 = get_relative_energy(calc_v(lam_mu), v_cm);
            double energy_2 = get_relative_energy(calc_v(lam_mu_add), v_cm);
            return (energy_0 - energy_1) / (energy_2 - energy_1);
        }

    public:
        collision_t() {}
        collision_t(std::vector<double>& vec) {  // TODO: а что делать, если в векторе меньше 8 элементов ???
            ind_1 = calc_allowed_ind(vec[0], vec[1], vec[2]);
            ind_2 = calc_allowed_ind(vec[3], vec[4], vec[5]);
            if (calc_v(ind_1).norm() < v_cut && calc_v(ind_2).norm() < v_cut) {
                isGood = true;
            } else {
                isGood = false;
            }
            s = vec[6] * square_ef;
            e = vec[7] * 2 * M_PI;
        }

        std::tuple<ind_t, ind_t, ind_t, ind_t, double> get_new_v() {
            ind_t lam, lam_add, mu, mu_add;
            double r;
            if (isGood) {
                // Находим скорости полсе столкновения
                double theta = 2.0 * std::acos(std::sqrt(s));
                v_t g, new_g, v1, v2;
                v1 = calc_v(ind_1);
                v2 = calc_v(ind_2);
                g = v2 - v1;
                double g_norm = g.norm();
                double g_xy = g.xy();
                if (g_xy == 0) {
                    new_g.x = g_norm * std::sin(e) * std::sin(theta);
                    new_g.y = g_norm * std::cos(e) * std::sin(theta);
                    new_g.z = g_norm * std::cos(theta);
                } else {
                    new_g.x = g.x * std::cos(theta) - (g.x * g.z / g_xy) * std::cos(e) * std::sin(theta) + g_norm * g.y / g_xy * std::sin(e) * std::sin(theta);
                    new_g.y = g.y * std::cos(theta) - (g.y * g.z / g_xy) * std::cos(e) * std::sin(theta) - g_norm * g.x / g_xy * std::sin(e) * std::sin(theta);
                    new_g.z = g.z * std::cos(theta) + g_xy * std::cos(e) * std::sin(theta);
                }
                v_t u1 = (v1 + v2) / 2 - new_g / 2;
                v_t u2 = (v1 + v2) / 2 + new_g / 2;

                // Находим аппроксимипующие скорости (индексы) и r
                ind_t near1 = find_near_ind(u1);
                ind_t near2 = find_near_ind(u2);
                v_t u_near1 = calc_v(near1);
                v_t u_near2 = calc_v(near2);
                if ((u_near1.norm() >= v_cut) || (u_near2.norm() >= v_cut)) {
                    isGood = false;
                } else {
                    v_t v_cm = (v1 + v2) / 2;  // скорость центра масс
                    double energy_0 = get_relative_energy(v1, v_cm);
                    double energy_near = get_relative_energy(u_near1, v_cm);
                    if (energy_0 == energy_near) {
                        lam = near1;
                        lam_add = near1;
                        mu = near2;
                        mu_add = near2;
                        r = 1;
                    } else {
                        if (energy_0 < energy_near) {
                            lam_add = near1;
                            mu_add = near2;
                            lam = this->find_add_ind(u1, v_cm, near1 + near2, energy_0, true);
                            mu = lam_add + mu_add - lam;
                        } else {
                            lam = near1;
                            mu = near2;
                            lam_add = this->find_add_ind(u1, v_cm, near1 + near2, energy_0, false);
                            mu_add = lam + mu - lam_add;
                        }
                        if (isGood) {
                            r = find_r(energy_0, lam, lam_add, v_cm);
                        }
                    }
                }
                // TODO: cтереть после проверки
                // if (isGood) {
                //     double energy_0 = v1.squared() + v2.squared();
                //     double energy_1 = calc_v(lam).squared() + calc_v(mu).squared();
                //     double energy_2 = calc_v(lam_add).squared() + calc_v(mu_add).squared();
                //     double r_add = (energy_0 - energy_1) / (energy_2 - energy_1);
                //     v_t v_lam = calc_v(lam);
                //     v_t v_mu = calc_v(mu);
                //     v_t v_lam_add = calc_v(lam_add);
                //     v_t v_mu_add = calc_v(mu_add);
                //     std::cout << energy_0 << " " << energy_1 << " " << energy_2 << '\n';
                //     std::cout << v_lam + v_mu - v1 - v2 << v_lam_add + v_mu_add - v1 - v2<< '\n';
                //     std::cout << r << " " << r_add << "\n\n";
                // }  
            }
            return std::make_tuple(lam, mu, lam_add, mu_add, r);
        }

        bool operator== (const collision_t& col) const {  // just velocity
            // if ((this->ind_1 == col.ind_1) && (this->ind_2 == col.ind_2))
            if (this->ind_2 == col.ind_2)
                return true;
            return false;
        }

        friend std::ostream& operator<< (std::ostream& out, const collision_t& col) {
            out << calc_v(col.ind_1) << "\t" << calc_v(col.ind_2) << "\t" << col.s << "\t" << col.e;
            return out;
        }
};

int magic_hash(int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

template<>
struct std::hash<collision_t> {
    std::size_t operator() (const collision_t& col) const {
        int hash = 0;
        ind_t ind;
        // for (int i = 0; i < 1; ++i) {
                // ind = col.ind_1;
                ind = col.ind_2;
            hash += magic_hash(ind.x) / 6 + magic_hash(ind.y) / 3 + magic_hash(ind.z) / 2;
        // }
        return hash;
    }
};

// std::vector<collision_t> make_collision_array() {
//     KorobovGrid kor_grid(n_col, space);
//     std::vector<collision_t> arr_col;
//     int i = 0;
//     int j = 0;
//     for (auto& row : kor_grid) {
//         arr_col.push_back(collision_t(row));
//         collision_t col(row);
//         if (calc_v(col.ind_1).norm() < v_cut && calc_v(col.ind_2).norm() < v_cut) {
//             // std::cout << i  << " " << j << '\t' << col << '\n';
//             // std::cout << calc_v(col.ind_1).norm() << " " << calc_v(col.ind_2).norm() << "\n\n";
//             ++j;
//         }
//         ++i;
//     }
//     std::cout << "Korobov grid: p = " << i << " , v1 and v2 < v_cut: " << j << "\n\n";

//     // TODO: стереть после дебага
//     // std::unordered_map<collision_t, int> table;
//     // table.reserve(arr_col.size());
//     // for (auto& col : arr_col) {
//     //     // if (col.ind_2.x == arr_col[5000].ind_2.x && col.ind_2.y == arr_col[5000].ind_2.y && col.ind_2.z == arr_col[5000].ind_2.z)
//     //     // if (col.isGood)
//     //         table[col]++;
//     // }
//     // int a = 0;
//     // for (auto& el : table) {
//     //     std::cout << el.first << "\t" << el.second << "\n";
//     //     a += el.second;
//     // }
//     // std::cout << table.size() << " " << a << std::endl;
//     //
//     std::random_shuffle(arr_col.begin(), arr_col.end());
//     return arr_col;
// }


bool change_positive_value(double** value, double* delta, int n) {
    // if (n == 0) {
    //     return true;
    // } else {
    //     if (*(value[n-1]) + *(delta[n-1]) < 0) {
    //         return false;
    //     } else {
    //         if (change_positive_value(value, delta, n - 1)) {
    //             *(value[n-1]) += *(delta[n-1]);
    //             return true;
    //         } else {
    //             return false;
    //         }
    //     }
    // }
    bool isPositive = true;
    double max_change = 0.0;
    int i = 0;
    
    while (i < n) {
        if (*(value[i]) + delta[i] < 0) {
            isPositive = false;
            break;
        } else {
            ++i;
        }
    }
    if (isPositive) {
        for (int j = 0; j < n; ++j) {
            *(value[j]) += delta[j];
            // if (std::abs(delta[j]) > max_change)
            //     max_change = std::abs(delta[j]);
        }
    }
    return isPositive;
}

double make_iteration(Distribution& f) {
    static KorobovGrid kor_grid(n_col, space, true);
    // static auto arr_col = make_collision_array();
    const double c = 1 / (4 * std::sqrt(2)) * square_ef * dv_x * dv_y * dv_z * tau * n0 * n0 / n_col;
    double max_delta = 0.0;

    int count(0), count_neg(0);
    for (auto& row : kor_grid) {
        collision_t col(row);
        if (col.isGood) {
            double r;
            ind_t lam, mu, lam_add, mu_add;
            std::tie(lam, mu, lam_add, mu_add, r) = col.get_new_v();

            // TODO: cтереть после проверки
            // if (col.isGood && count == 56) {
            //     v_t v1 = calc_v(col.ind_1);
            //     v_t v2 = calc_v(col.ind_2);
            //     double energy_0 = v1.squared() + v2.squared();
            //     double energy_1 = calc_v(lam).squared() + calc_v(mu).squared();
            //     double energy_2 = calc_v(lam_add).squared() + calc_v(mu_add).squared();
            //     double r_add = (energy_0 - energy_1) / (energy_2 - energy_1);
            //     v_t v_lam = calc_v(lam);
            //     v_t v_mu = calc_v(mu);
            //     v_t v_lam_add = calc_v(lam_add);
            //     v_t v_mu_add = calc_v(mu_add);
            //     std::cout << energy_0 << " " << energy_1 << " " << energy_2 << '\n';
            //     std::cout << (v_lam + v_mu).x  << " " << (v_lam_add + v_mu_add).x << " " << (v1 + v2).x <<'\n';
            //     std::cout << r << " " << r_add << "\n\n";
            // }  

            if (col.isGood) {
                double omega = (f(lam) * f(mu) * std::pow(f(lam_add) * f(mu_add) / f(lam) / f(mu), r) - f(col.ind_1) * f(col.ind_2))
                               * (calc_v(col.ind_1) - calc_v(col.ind_2)).norm();
                double delta = c * omega;
                
                if (omega != omega) {
                    std::cout << count  << ' ' << omega << '\t' << col << '\n';
                    std::cout << col.ind_1 << ' ' << col.ind_2 << ' ' << lam << ' ' << mu << lam_add << ' ' << mu_add << '\n';
                    std::cout << f(col.ind_1) << ' ' << f(col.ind_2) << ' ' << f(lam) << ' ' << f(mu) << f(lam_add) << ' ' << f(mu_add) << '\n';
                }

                double* val[]{&f(col.ind_1), &f(col.ind_2), &f(lam), &f(mu), &f(lam_add), &f(mu_add)};
                double del[]{delta, delta, (r - 1) * delta, (r - 1) * delta, -r * delta, -r * delta};
                bool wasChanged = change_positive_value(val, del, 6);
                // if (change > max_delta)
                //     max_delta = change;
                count++;
                if (!wasChanged) {
                    count_neg++;
                }
                // std::cout << '\n';
                // std::cout << count << " " << f.get_concentration() << " " << f.get_momentum() << " " << f.get_energy() << '\n';
                // count++;
            }
        }
    }
    kor_grid.random_shift();

    // TODO: если не нужна таблица, то эти строчки не нужны
    neg_f_table.push_back(count_neg);
    neg_f_table.push_back(count);

    return max_delta;
    // TODO: удалить после дебага
    // std::unordered_map<collision_t, int> table;
    // table.reserve(arr_col.size());
    // for (auto& col : arr_col) {
    //     // if (col.ind_2.x == arr_col[5000].ind_2.x && col.ind_2.y == arr_col[5000].ind_2.y && col.ind_2.z == arr_col[5000].ind_2.z)
    //     if (col.isGood)
    //         table[col]++;
    // }
    // int a = 0;
    // for (auto& el : table) {
    //     std::cout << el.first << "\t" << el.second << "\n";
    //     a += el.second;
    // }
    // std::cout << table.size() << " " << a << std::endl;
    
    // std::cout << count << "\n";
}


int main() {
    int n_iter = time / tau;
    n_col = next_prime(250000);
    // std::vector<int> p_arr{500009, 1000003, 2000003, 4000037, 8000009};

    Distribution f(Distribution::TwoGauss, 1.0);
    f.save_to_file("data/Distribution/p_" + std::to_string(n_col) + "_init.txt");
    // Distribution g(Distribution::Maxwell, f.get_energy() / f.get_concentration() / 3);
    

    double n = f.get_concentration();
    v_t p = f.get_momentum();
    double eng = f.get_energy();
    std::cout << "Check " << f.moment_2(2, 2) - eng << '\n';
    std::cout << n << " " << p << " " << eng << "\n\n";

    std::ofstream out("data/relax_time_" + std::to_string(n_col) + ".txt");
    out << "# Max change (t), [t] - время свободного пробега\n";

    int i = 0;
    double eps = 10e-18;
    double delta = 2 * eps;
    while (i < n_iter) {
        delta = make_iteration(f);
        double momemnt_xx = f.moment_2(1, 1);
        out << i * tau << '\t' << momemnt_xx << '\t' << f.moment_2(2,2) << '\t' << f.moment_2(3,3) << '\n';
        if (i % 1 == 0) {
            // delta = f.max_difference_from_Maxwell();
            std::string filename = "data/Distribution/p_" + std::to_string(n_col) + "_f_" + std::to_string(i) + ".txt";
            f.save_to_file(filename);
            std::cout << "n_col " << n_col << ", step " << i + 1 << ", delta " << momemnt_xx << "\n\n";
        } else {
            std::cout << i + 1 << '\n';
        }
        ++i;
    }
    out.close();
    std::cout << "Relaxtion time: " << i * tau << '\n';

    // int n = 4;
    // double a(10), b(8), c(6), d(4);
    // double w(-10.000), x(-6), y(-4), z(-5);
    // double* val[n]{&a, &b, &c, &d};
    // double* del[n]{w, x, y, z};
    // change_positive_value(val, del, n);
    // std::cout << a << " " << b << " " << c << " " << d <<std::endl;

    std::vector<int> p_arr = {n_col};
    std::ofstream out2("data/neg_f.txt");
    out2 << "N of iter\t ";
    for (const auto prime : p_arr) {
        out2 << "p = " << prime << " N of neg_f and N of total_f\t";
    }
    out2 << '\n';
    // Переписываем одномерные данные в читаемую таблицу
    for (int i = 0; i < n_iter; ++i) {
        out2 << i + 1;
        for (int j = 0; j < p_arr.size(); ++j) {
            int ind = j * 2 * n_iter + 2 * i;
            out2 << '\t' << neg_f_table[ind] << '\t' << neg_f_table[ind + 1];
            if (neg_f_table[ind] > 0) {
                out2 << "\t here";
            }
        }
        out2 << '\n';
    }
    out2.close();

    return 0;
}
