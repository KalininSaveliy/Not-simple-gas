#include <array>
// #include <chrono>  // for srand(time)
#include <ctime>  // for srand(time)
#include <cmath>  // for pi, round, acos
#include <cstdlib>
#include <iostream>
#include <tuple>


const int dim = 3;
const int n_col = 20 * 1000;  // количество столкновений
const int v_cut = 5;
const int n_vx = 20, n_vy = 20, n_vz = 20;

const double d_vx = 2.0 * v_cut / n_vx;
const double d_vy = 2.0 * v_cut / n_vy;
const double d_vz = 2.0 * v_cut / n_vz;

// union collision_t {
//     double val[8];
//     struct {
//         vec3d v1;
//         vec3d v2;
//         double b;
//         double e;
//     };
// };

// int get_static() {
//     static int seed = 4;
//     seed += 3;
//     return seed;
// }

double rand_0_1() {  // return number from [0, 1)
    // std::srand(std::time(nullptr) * get_static());  // use current time as seed
    // std::srand(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count() % RAND_MAX);
    return 1.0 * std::rand() / RAND_MAX;
}

double rand_sym(double max) {
    return rand_0_1() * 2.0 * max - max;
}

int find_near_i(double v, double dv) {
    int i = std::lround((v + v_cut) / dv - 0.5);
    if (i < 0) {
        i = 0;
    }
    return i;
}

int find_min_i(double v, double dv) {
    int i = std::floor((v + v_cut) / dv - 0.5);  // TODO: устранить повтор в find_near_i
    if (i < 0) {i = 0;}
    return i;
}

double calc_v(int i, double dv) {
    return -v_cut + (i + 0.5) * dv;
}

double calc_allowed_v(double v, double dv) {
    int i = find_near_i(v, dv);
    return calc_v(i, dv);
}

template <typename Type>
struct vec3d {
    Type x;
    Type y;
    Type z;

    vec3d(): x(0), y(0), z(0) {}
    vec3d(Type x0, Type y0, Type z0): x(x0), y(y0), z(z0) {}

    double squared() {
        return this->x * this->x + this->y * this->y + this->z * this->z; 
    }

    double norm() {
        return std::sqrt(this->squared());
    }

    double xy() {
        return std::sqrt(this->x * this->x + this->y * this->y);
    }

    bool operator== (const vec3d& vec) {
        return (this->x == vec.x) && (this->y == vec.y) && (this->z == vec.z);
    }

    vec3d operator- (const vec3d& vec) {
        vec3d res;
        res.x = this->x - vec.x;
        res.y = this->y - vec.y;
        res.z = this->z - vec.z;
        return res;
    }

    vec3d operator+ (const vec3d& vec) {
        vec3d res;
        res.x = this->x + vec.x;
        res.y = this->y + vec.y;
        res.z = this->z + vec.z;
        return res;
    }

    vec3d operator+= (const vec3d& vec) {
        this->x += vec.x;
        this->y += vec.y;
        this->z += vec.z;
        return *this;
    }

    vec3d operator* (Type a) {
        vec3d res;
        res.x = this->x * a;
        res.y = this->y * a;
        res.z = this->z * a;
        return res;
    }

    vec3d operator/ (Type a) {
        vec3d res;
        res.x = this->x / a;
        res.y = this->y / a;
        res.z = this->z / a;
        return res;
    }

    friend std::ostream& operator<< (std::ostream& out, const vec3d& vec) {
        out << "{" << vec.x << " " << vec.y << " " << vec.z << "}";
        return out;
    }
};

struct collision_t {
    using v_t = vec3d<double> ;
    using indexes_t = vec3d<int>;
    v_t v1;
    v_t v2;
    double b;
    double e;
    bool isGood;

    private:
        double init_rand_v(double dv) {
            return calc_allowed_v(rand_sym(v_cut), dv);
        }

        void init_right_v(v_t& vec) {
            vec.x = init_rand_v(d_vx);
            vec.y = init_rand_v(d_vy);
            vec.z = init_rand_v(d_vz);
        }

        indexes_t find_near_ijk(const v_t vec) {
            indexes_t indexes;
            indexes.x = find_near_i(vec.x, d_vx);
            indexes.y = find_near_i(vec.y, d_vy);
            indexes.z = find_near_i(vec.z, d_vz);
            return indexes;
        }

        indexes_t find_min_ijk(const v_t vec) {
            indexes_t min_ind;
            min_ind.x = find_min_i(vec.x, d_vx);
            min_ind.y = find_min_i(vec.y, d_vy);
            min_ind.z = find_min_i(vec.z, d_vz);
            return min_ind;
        }

        v_t calc_v3d(const indexes_t ind) {
            v_t v;
            v.x = calc_v(ind.x, d_vx);
            v.y = calc_v(ind.y, d_vy);
            v.z = calc_v(ind.z, d_vz);
            return v;
        }

        double get_relative_energy(v_t vec, v_t v_centre) {
            return (vec - v_centre).squared();
        }

        indexes_t find_add_ijk(v_t u1, v_t v_cm, indexes_t sum_ind, double energy_0, bool isEnergy0_bigger) {
            indexes_t centre = find_min_ijk(u1);  // cenetre index
            indexes_t best_ind;
            double big_number = 12.0 * v_cut * v_cut;
            double min = big_number;

            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int k = 0; k < 2; k++) {
                        indexes_t cur_ind = centre + indexes_t(i, j, k);
                        indexes_t paired_ind = sum_ind - cur_ind;
                        v_t cur_vec = calc_v3d(cur_ind);
                        v_t paired_vec = calc_v3d(paired_ind);
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
        
        double find_r(double energy_0, indexes_t lam_mu, indexes_t lam_mu_add, v_t v_cm) {
            double energy_1 = get_relative_energy(calc_v3d(lam_mu), v_cm);
            double energy_2 = get_relative_energy(calc_v3d(lam_mu_add), v_cm);
            return (energy_0 - energy_1) / (energy_2 - energy_1);
        }

    public:
        collision_t() {
            isGood = true;
            b = rand_0_1();
            e = 2.0 * M_PI * rand_0_1();
            init_right_v(v1);
            init_right_v(v2);
            if ((v1 == v2) || (v1.norm() > v_cut) || (v2.norm() > v_cut)) {
                isGood = false;
            }
        }

        std::tuple<indexes_t, indexes_t, indexes_t, indexes_t, double> get_new_v() {
            indexes_t lam, lam_add, mu, mu_add;
            double r;
            if (isGood) {
                // Находим скорости полсе столкновения
                v_t new_g;
                double theta = 2.0 * std::acos(b);
                v_t g = v2 - v1;
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
                indexes_t near1 = find_near_ijk(u1);
                indexes_t near2 = find_near_ijk(u2);
                v_t u_near1 = calc_v3d(near1);
                v_t u_near2 = calc_v3d(near2);
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
                            lam = this->find_add_ijk(u1, v_cm, near1 + near2, energy_0, true);
                            mu = lam_add + mu_add - lam;
                        } else {
                            lam = near1;
                            mu = near2;
                            lam_add = this->find_add_ijk(u1, v_cm, near1 + near2, energy_0, false);
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
                //     double energy_1 = calc_v3d(lam).squared() + calc_v3d(mu).squared();
                //     double energy_2 = calc_v3d(lam_add).squared() + calc_v3d(mu_add).squared();
                //     double r_add = (energy_0 - energy_1) / (energy_2 - energy_1);
                //     v_t v_lam = calc_v3d(lam);
                //     v_t v_mu = calc_v3d(mu);
                //     v_t v_lam_add = calc_v3d(lam_add);
                //     v_t v_mu_add = calc_v3d(mu_add);
                //     std::cout << energy_0 - energy_1 << " " << energy_2 - energy_0 << '\n';
                //     std::cout << v_lam + v_mu - v1 - v2 << v_lam_add + v_mu_add - v1 - v2<< '\n';
                //     std::cout << r << " " << r_add << "\n\n";
                // }
                
            }
            return std::make_tuple(lam, mu, lam_add, mu_add, r);
        }

        friend std::ostream& operator<< (std::ostream& out, const collision_t& col) {
            out << col.v1 << " " << col.v2 << " " << col.b << " " << col.e;
            return out;
        }
};

using v_grid_t = std::array<std::array<std::array<std::array<double, dim>, n_vz>, n_vy>, n_vx>;
using arr_col_t = std::array<collision_t, n_col>;


v_grid_t create_v_grid() {
    v_grid_t v_grid;
    int i = 0;
    for (auto& slice : v_grid) {
        int j = 0;
        for (auto& row : slice) {
            int k = 0;
            for (auto& el : row) {
                double vx = calc_v(i, d_vx);
                double vy = calc_v(j, d_vy);
                double vz = calc_v(k, d_vz);
                el = {vx, vy, vz};
                ++k;
            }
            ++j;
        }
        ++i;
    }
    return v_grid;
}



int main() {
    std::srand(std::time(nullptr));
    // std::srand(1);
    // auto v_grid = create_v_grid();
    arr_col_t arr_col;

    int count = 0;
    for (auto& col : arr_col) {
        if (col.isGood) {
            double r;
            collision_t::indexes_t a, b, c, d;
            std::tie(a, b, c, d, r) = col.get_new_v();
            if (col.isGood) {
                std::cout << a + b << c + d  << " " << r << '\n';
                count++;
            }
        }
    }
    std::cout << count << '\n';
    
    return 0;
}
