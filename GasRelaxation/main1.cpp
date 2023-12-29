#include <iostream>
#include <fstream>
#include <string>

#include "source\GasRelaxation.h"


const size_t n_v(30);

void save_to_file(std::string& filename, gas_relaxation::Distribution<n_v, n_v, n_v>& f) {
    std::ofstream out(filename);
    out << "# Nv = " << n_v << ", v_cut = " << f.get_v_cut() << '\n';
    out << '#' << f.concentration() << '\t' << f.energy() << '\n';
    out << "# v_x   f(v_x) = int f(v) dv_y dv_z\n";
    for (int i = 0; i < n_v; ++i) {
        double n = 0.0;
        for (int j = 0; j < n_v; ++j) {
            for (int k = 0; k < n_v; ++k)
                n += f[i][j][k];
        }
        out << f.calc_v(i, f.get_dv().x) << '\t' << n * f.get_dv().y * f.get_dv().z << '\n';
    }
    out.close();
}

int main() {
    int n_col(250000);
    double tau(0.02), time(20);
    double v_cut(5.0);

    using namespace gas_relaxation;
    Distribution<n_v, n_v, n_v> f(v_cut);
    f.init(TwoGauss(1.0, {2.0, 0.0, 0.0}), false);

    std::string file_f = "./data/sample_f/";
    std::string file_m = "./data/sample_m.txt";
    std::ofstream out(file_m);
    out << "# N_v = " << n_v << ", N_col = " << n_col << '\n';
    out << "# time\t M_xx\t M_yy\t M_zz\n";
    out << 0.0 << '\t' << f.M_2(1,1) << '\t' << f.M_2(2,2) << '\t' << f.M_2(3,3) << '\n';

    for (int i = 0; i < time / tau; ++i) {
        f.make_iteration(tau, n_col,  true, true);
        out << (i+1) * tau << '\t' << f.M_2(1,1) << '\t' << f.M_2(2,2) << '\t' << f.M_2(3,3) << '\n';
        std::string filename = file_f + 'f' + std::to_string(i) + ".txt";
        save_to_file(filename, f);
        std::cout << "Step " << i+1 << " / " << time/tau << '\n';
    }
    out.close();

    return 0;
}
