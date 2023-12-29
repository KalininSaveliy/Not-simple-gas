#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "KorobovGrid.h"
#include "Iterator3dim.h"
#include "vec3d.h"
using v_t = vec3d_t<double>;
using ind_t = vec3d_t<int>;

namespace gas_relaxation {

    const double NORM_COEF = 1 / (2 * M_PI) / sqrt(2 * M_PI);  // нормировочный коэффициент для распределения Максвелла по скоростям

    struct DistributionFunc {
        double temp = 1.0;
        DistributionFunc(double temperature) : temp(temperature) {}
        double operator() (const v_t&);
    };

    struct Maxwell : DistributionFunc {
        double operator() (const v_t& v) {
            return NORM_COEF * std::exp(v.squared() / (-2.0 * temp)) / sqrt(std::pow(temp, 3));
        }  
    };
    struct TwoGauss : DistributionFunc {
        v_t u;
        TwoGauss(double temperature = 1.0, v_t shift = {2.0, 0.0, 0.0}) : DistributionFunc(temperature), u(shift) {}
        double operator() (const v_t& v) {
            return NORM_COEF * 0.5 * (std::exp((v - u).squared() / (-2.0 * temp)) + std::exp((v + u).squared() / (-2.0 * temp))) / sqrt(std::pow(temp, 3));
        }
    };


    template<std::size_t N_vx, std::size_t N_vy, std::size_t N_vz>
    struct IntegralSpace {
        public:
            double b_max_squared = 1.0;  // в площадях d^2, где d - диаметр молекулы
        private:
            double v_cut = 5.0;  // в скоростях \ksi = sqrt(T_0 / m)
            v_t dv = calc_dv();

        public:
            double calc_dv(int n) const {
                return 2 * v_cut / n;
            }
            v_t calc_dv() const {
                return v_t(calc_dv(N_vx), calc_dv(N_vy), calc_dv(N_vz));
            }
            double calc_v(int i, double dv) const {
                return -v_cut + (i + 0.5) * dv;
            }
            v_t calc_v(const ind_t& ind) const {
                v_t v;
                v.x = calc_v(ind.x, dv.x);
                v.y = calc_v(ind.y, dv.y);
                v.z = calc_v(ind.z, dv.z);
                return v;
            }

            double get_v_cut() const {
                return this->v_cut;
            }
            v_t get_dv() const {
                return this->dv;
            }
            void set_v_cut(double v_max) {
                this->v_cut = v_max;
                this->dv = calc_dv();
            }

            IntegralSpace() {}
            IntegralSpace(double v_max) : v_cut(v_max) {}
            IntegralSpace(double v_max, double _b_max_squared) : v_cut(v_max), b_max_squared(_b_max_squared) {}
    };

    template<std::size_t N_vx = 10, std::size_t N_vy = 10, std::size_t N_vz = 10>
    class Collision : IntegralSpace<N_vx, N_vy, N_vz> {
        public:
            ind_t ind_1;
            ind_t ind_2;
            double s;
            double e;
            bool isGood;

        private:
            int find_near_ind(double v, double dv) {
                int i = std::lround((v + this->get_v_cut()) / dv - 0.5);
                // if (i < 0) {i = 0;}
                return i;
            }

            ind_t find_near_ind(const v_t& vec) {
                ind_t ind;
                ind.x = find_near_ind(vec.x, this->get_dv().x);
                ind.y = find_near_ind(vec.y, this->get_dv().y);
                ind.z = find_near_ind(vec.z, this->get_dv().z);
                return ind;
            }

            ind_t calc_allowed_ind(double x, double y, double z) {
                v_t v = v_t(x, y, z) * 2 * this->get_v_cut() - v_t(this->get_v_cut());
                ind_t ind = find_near_ind(v);
                return ind;
            }

            double get_relative_energy(v_t vec, v_t v_centre) {
                return (vec - v_centre).squared();
            }

            int find_min_ind(double v, double dv) {
                int i = std::floor((v + this->get_v_cut()) / dv - 0.5);  // TODO: устранить повтор в find_near_i
                // if (i < 0) {i = 0;}
                return i;
            }

            ind_t find_min_ind(const v_t vec) {
                ind_t min_ind;
                min_ind.x = find_min_ind(vec.x, this->get_dv().x);
                min_ind.y = find_min_ind(vec.y, this->get_dv().y);
                min_ind.z = find_min_ind(vec.z, this->get_dv().z);
                return min_ind;
            }

            ind_t find_add_ind(v_t u1, v_t v_cm, ind_t sum_ind, double energy_0, bool isEnergy0_bigger) {
                ind_t centre = find_min_ind(u1);  // cenetre index
                ind_t best_ind;
                double big_number = 12.0 * this->get_v_cut() * this->get_v_cut();
                double min = big_number;

                // Iterator3dim<2, 2, 2> it({0,0,0});
                // for (Iterator3dim<2,2,2> it({0,0,0}); *it != vec3d_t<int>(2,0,0); ++it){}
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        for (int k = 0; k < 2; k++) {
                            ind_t cur_ind = centre + ind_t(i, j, k);
                            ind_t paired_ind = sum_ind - cur_ind;
                            v_t cur_vec = this->calc_v(cur_ind);
                            v_t paired_vec = this->calc_v(paired_ind);
                            if ((cur_vec.squared() < this->get_v_cut() * this->get_v_cut()) && (paired_vec.squared() < this->get_v_cut() * this->get_v_cut())) {
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
                double energy_1 = get_relative_energy(this->calc_v(lam_mu), v_cm);
                double energy_2 = get_relative_energy(this->calc_v(lam_mu_add), v_cm);
                return (energy_0 - energy_1) / (energy_2 - energy_1);
            }

        public:
            Collision() {}
            Collision(std::vector<double>& vec, double v_max, double b_max_sq = 1.0) : IntegralSpace<N_vx, N_vy, N_vz>(v_max, b_max_sq) {  // TODO: а что делать, если в векторе меньше 8 элементов ???
                ind_1 = calc_allowed_ind(vec[0], vec[1], vec[2]);
                ind_2 = calc_allowed_ind(vec[3], vec[4], vec[5]);
                if (this->calc_v(ind_1).norm() < this->get_v_cut() && this->calc_v(ind_2).norm() < this->get_v_cut()) {
                    isGood = true;
                } else {
                    isGood = false;
                }
                s = vec[6] * this->b_max_squared;
                e = vec[7] * 2 * M_PI;
            }

            std::tuple<ind_t, ind_t, ind_t, ind_t, double> get_new_v() {
                ind_t lam, lam_add, mu, mu_add;
                double r;
                if (isGood) {
                    // Находим скорости полсе столкновения
                    double theta = 2.0 * std::acos(std::sqrt(s));
                    v_t g, new_g, v1, v2;
                    v1 = this->calc_v(ind_1);
                    v2 = this->calc_v(ind_2);
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
                    v_t u_near1 = this->calc_v(near1);
                    v_t u_near2 = this->calc_v(near2);
                    if ((u_near1.norm() >= this->get_v_cut()) || (u_near2.norm() >= this->get_v_cut())) {
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

            bool operator== (const Collision& col) const {  // just velocity
                // if ((this->ind_1 == col.ind_1) && (this->ind_2 == col.ind_2))
                if (this->ind_2 == col.ind_2)
                    return true;
                return false;
            }

            friend std::ostream& operator<< (std::ostream& out, const Collision& col) {
                out << col->calc_v(col.ind_1) << "\t" << col->calc_v(col.ind_2) << "\t" << col.s << "\t" << col.e;
                return out;
            }
    };

    template<std::size_t N_vx = 10, std::size_t N_vy = 10, std::size_t N_vz = 10>
    class Distribution : public IntegralSpace<N_vx, N_vy, N_vz>, public std::array<std::array<std::array<double, N_vz>, N_vy>, N_vx> {
        private:
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
                int next_prime = n;
                while (!is_prime[next_prime]) {  // от [n до 2n-1] точно есть простое, так что цикл не вечный
                    next_prime++;
                }
                return next_prime;
            }
            bool change_positive_value(double** value, double* delta, int n) {
                // if (n == 0) {  // рекурсивное решение
                //     return true;
                // } else {
                //     if (*(value[n-1]) + *(delta[n-1]) < 0)
                //         return false;
                //     else
                //         return change_positive_value(value, delta, n - 1);
                // }
                bool isPositive = true;
                double max_change = 0.0;
                int i = 0;
                
                while (i < n) {
                    if (*(value[i]) + delta[i] < 0) {
                        isPositive = false;
                        break;
                    } else
                        ++i;
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
            int find_n0() {  // число узлов внутри сферы радиуса v_cut
                int n0 = 0;
                for (auto& ind : *this) {
                    if (this->calc_v(ind).norm() < this->get_v_cut())
                        n0++;
                }
                return n0;
            }
            int n0 = find_n0();

        public:
            void set_v_cut(double v_cut) {
                this->IntegralSpace<N_vx, N_vy, N_vz>::set_v_cut(v_cut);
                n0 = find_n0();
            }

            typedef Iterator3dim<N_vx, N_vy, N_vz, ind_t> iterator;
            typedef Iterator3dim<N_vx, N_vy, N_vz, ind_t> const_iterator;  // TODO: подумать над разницей между iterator и const_iterator
            iterator begin() {
                return iterator({0, 0, 0});
            }
            iterator end() {
                return iterator({N_vx, 0, 0});
            }
            const_iterator begin() const {
                return const_iterator({0, 0, 0});
            }
            const_iterator end() const {
                return const_iterator({N_vx, 0, 0});
            }

            double& operator() (const ind_t& ind) {
                return (*this)[ind.x][ind.y][ind.z];
            }
            double operator() (const ind_t& ind) const {  // TODO: is it ok? (two overload for operator() )
                return (*this)[ind.x][ind.y][ind.z];
            }

            double d3v() const {
                v_t d_v = this->get_dv();
                return  d_v.x * d_v.y * d_v.z;
            }
            double concentration() const {
                double n = 0.0;
                for (auto& ind : *this)
                    n += (*this)(ind);
                return n * d3v();
            }
            v_t momentum() const {
                v_t momentum(0.0, 0.0, 0.0);
                for (auto& ind : *this) 
                    momentum += this->calc_v(ind) * (*this)(ind);
                return momentum * d3v();
            }
            double energy() const {
                double energy = 0.0;
                for (auto& ind : *this)
                    energy += this->calc_v(ind).squared() * (*this)(ind);
                return energy * d3v();
            }
            double temperature() const {
                return energy() / concentration() / 3;
            }
            // M_{i,j} = integral(v_i * v_j * f * d^3v)
            double M_2(int alpha, int beta) const {   // ind = 1, 2 or 3
                double m = 0.0;
                for (auto& ind : *this)
                    m += this->calc_v(ind[alpha], this->get_dv()[alpha]) * this->calc_v(ind[beta], this->get_dv()[beta]) * (*this)(ind);
                return m * d3v();
            }

            void rationing() {
                double n = concentration();
                for (auto& ind : *this)
                    (*this)(ind) /= n;
            }

            template<typename Functional>
            void init(Functional obj_or_func, bool haveRationing = true) {  // obj_or_func должен быть вызываемым от v_t, т.е. obj_or_func(v) должно выдавать double
                for (ind_t& ind : *this) {
                    v_t v = this->calc_v(ind);
                    if (v.norm() < this->get_v_cut())
                        (*this)(ind) = obj_or_func(v);
                    else
                        (*this)(ind) = 0.0;
                }
                if (haveRationing)
                    rationing();
            }
            void init(double (*VelosityDistFuncT)(const v_t&, double), double temp = 1.0, bool haveRationing = true) {  // TODO: delete or make constructor
                for (ind_t& ind : *this) {
                    v_t v = this->calc_v(ind);
                    if (v.norm() < this->get_v_cut())
                        (*this)(ind) = VelosityDistFuncT(v, temp);
                    else
                        (*this)(ind) = 0.0;
                }
                if (haveRationing)
                    rationing();
            }

            Distribution() {}
            Distribution(double v_max) : IntegralSpace<N_vx, N_vy, N_vz>(v_max) {}
            Distribution(double v_max, double b_max_sq) : IntegralSpace<N_vx, N_vy, N_vz>(v_max, b_max_sq) {}
            template<typename Functional>
            Distribution(Functional VelocityDistFunc, double v_max = 5.0, double b_max_sq = 1.0, bool haveRationing = true) : IntegralSpace<N_vx, N_vy, N_vz>(v_max, b_max_sq) {
                init(VelocityDistFunc, haveRationing);
            }

            void add_header(std::ofstream& out) const {
                out << "# n_vx = " << N_vx << ", n_vy = " << N_vy << ", n_vz = " << N_vz << '\n';
            }
            void add_param(std::ofstream& out) const {
                out << "n = " << concentration() << "\nT = " << temperature() << "\nE = " << energy() << "\nvec_p = " << momentum() << '\n';
            }
            void save_to_file_3dim(std::string& filename, bool haveParam = false) const {
                std::ofstream outfile(filename);
                outfile << "# f(v_x, v_y, v_z)\n";
                add_header(outfile);
                if (haveParam)
                    add_param(outfile);
                for (auto& ind : *this) {
                    v_t v = calc_v(ind);
                    outfile << v.x << ' ' << v.y << ' ' << v.z << ' ' << (*this)(ind) << '\n';
                }
            }
            void save_to_file_1dim(std::string filename, ind_t print = ind_t(-1, N_vy / 2, N_vz / 2), bool haveParam = false) const {
                std::ofstream outfile(filename);
                add_header(outfile);
                outfile << "# f(v_x, v_y, v_z)\n";
                int out_i = 0;
                for (int i : {1, 2, 3}) {
                    if (print[i] >= 0) {
                        out_i += i;
                        outfile << ", where v_";
                        switch (i)
                        {
                        case 1:
                            outfile << "x";
                            break;
                        case 2:
                            outfile << "y";
                            break;
                        case 3:
                            outfile << "z";
                            break;
                        
                        default:
                            break;
                        }
                        outfile << " = " << this->calc_v(print[i], this->dv[i]) << ' ';
                    }
                }
                outfile << "\n\n";
                if (haveParam)
                    add_param(outfile);

                out_i = 1 + 2 + 3 - out_i;
                int not_out1 = out_i % 3 + 1;  // 1, 2 or 3
                int not_out2 = (out_i + 1) % 3 + 1;
                for (auto& ind : (*this)) {
                    if (ind[not_out1] == print[not_out1] && ind[not_out2] == print[not_out2])
                        outfile << this->calc_v(ind[out_i], this->dv[out_i]) << '\t' << (*this)(ind) << '\n';
                }
            }
            
            void make_iteration(double tau, double n_col, bool haveShuffle = false, bool haveShifft = false) {
                n_col = next_prime(n_col);
                KorobovGrid kor_grid(n_col, 8, haveShuffle, haveShifft);
                // static auto arr_col = make_collision_array();
                const double c = 1 / (4 * std::sqrt(2)) * this->b_max_squared * d3v() * tau * n0 * n0 / n_col;
                double max_delta = 0.0;

                int count(0), count_neg(0);
                for (auto& row : kor_grid) {
                    Collision<N_vx, N_vy, N_vz> col(row, this->get_v_cut(), this->b_max_squared);
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
                            double omega = ((*this)(lam) * (*this)(mu) * std::pow((*this)(lam_add) * (*this)(mu_add) / (*this)(lam) / (*this)(mu), r) - (*this)(col.ind_1) * (*this)(col.ind_2))
                                        * (this->calc_v(col.ind_1) - this->calc_v(col.ind_2)).norm();
                            double delta = c * omega;
                            
                            // if (omega != omega) {
                            //     std::cout << count  << ' ' << omega << '\t' << col << '\n';
                            //     std::cout << col.ind_1 << ' ' << col.ind_2 << ' ' << lam << ' ' << mu << lam_add << ' ' << mu_add << '\n';
                            //     std::cout << f(col.ind_1) << ' ' << f(col.ind_2) << ' ' << (*this)(lam) << ' ' << (*this)(mu) << (*this)(lam_add) << ' ' << (*this)(mu_add) << '\n';
                            // }

                            double* val[]{&(*this)(col.ind_1), &(*this)(col.ind_2), &(*this)(lam), &(*this)(mu), &(*this)(lam_add), &(*this)(mu_add)};
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
                if (haveShifft)
                    kor_grid.random_shift();

                // TODO: если не нужна таблица, то эти строчки не нужны
                // neg_f_table.push_back(count_neg);
                // neg_f_table.push_back(count);

                // return max_delta;
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
            void make_iteration(double tau, bool haveShuffle = false, bool haveShift = false) {  // TODO: неоптимально, каждый раз надо считать n_col
                int n_col = next_prime(std::ceil(2 * tau * n0 * this->v_cut * this->b_max_squared * n0 * d3v() * NORM_COEF / std::sqrt(2)));
                make_iteration(tau, n_col, haveShuffle, haveShift);
            }

            void simulate(int n_iter, double tau, double n_col) {
                for (int i = 0; i < n_iter; ++i) {
                    make_iteration(tau, n_col, true, true);
                }
            }
            void simulate(double time, double tau, double n_col) {
                int n_iter = time / tau;
                simulate(n_iter, tau, n_col);
            }

            double max_difference_from_Maxwell() const;  // TODO: дописать
    };
}
