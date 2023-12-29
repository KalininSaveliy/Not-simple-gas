#include <iostream>
#include <string>
#include <vector>

#include ".\source\GasRelaxation.h"

int main() {
    using namespace gas_relaxation;
    std::string data_dir = ".\\data\\";
    
    double v_cut(5.0);
    const std::size_t N1(4), N2(6), N3(8), N4(10), N5(30);
    double temperature(1.0);
    v_t u(2.0, 0.0, 0.0);
    Distribution<N1, N1, N1> f1(TwoGauss(temperature, u), v_cut, 1.0, true);
    Distribution<N2, N2, N2> f2(TwoGauss(temperature, u), v_cut, 1.0, true);
    Distribution<N3, N3, N3> f3(TwoGauss(temperature, u), v_cut, 1.0, true);
    Distribution<N4, N4, N4> f4(TwoGauss(temperature, u), v_cut, 1.0, true);
    Distribution<N5, N5, N5> f5(TwoGauss(temperature, u), v_cut, 1.0, true);

    int n_col1(1000), n_col2(5000), n_col3(10000), n_col4(125000), n_col5(250000);
    double time(20.0), tau(0.02);

    std::cout << "1) n = " << f1.concentration() << "\tT = " << f1.temperature() << "\tE = " << f1.energy() << "\tvec_p = " << f1.momentum() << '\n';
    std::cout << "2) n = " << f2.concentration() << "\tT = " << f2.temperature() << "\tE = " << f2.energy() << "\tvec_p = " << f2.momentum() << '\n';
    std::cout << "3) n = " << f3.concentration() << "\tT = " << f3.temperature() << "\tE = " << f3.energy() << "\tvec_p = " << f3.momentum() << '\n';
    std::cout << "4) n = " << f4.concentration() << "\tT = " << f4.temperature() << "\tE = " << f4.energy() << "\tvec_p = " << f4.momentum() << '\n';
    std::cout << "5) n = " << f5.concentration() << "\tT = " << f5.temperature() << "\tE = " << f5.energy() << "\tvec_p = " << f5.momentum() << '\n';


    std::string filename = data_dir + "several_f.txt";
    std::ofstream out(filename);
    out << "# v_cut = " << v_cut << ", tau = " << tau << '\n';
    out << "#1) N = " << N1 << ", n_col = " << n_col4 << '\n';
    out << "#2) N = " << N2 << ", n_col = " << n_col4 << '\n';
    out << "#3) N = " << N3 << ", n_col = " << n_col4 << '\n';
    out << "#4) N = " << N4 << ", n_col = " << n_col4 << '\n';
    out << "#5) N = " << N5 << ", n_col = " << n_col5 << '\n';

    out << "0.0\t" << f1.M_2(1,1) << '\t'
                   << f2.M_2(1,1) << '\t'
                   << f3.M_2(1,1) << '\t'
                   << f4.M_2(1,1) << '\t'
                   << f5.M_2(1,1) << '\n';

    std::cout << "Start\n";
    for (int i = 0; i < time / tau; ++i) {
        f1.make_iteration(tau, n_col4, true, true);
        f2.make_iteration(tau, n_col4, true, true);
        f3.make_iteration(tau, n_col4, true, true);
        f4.make_iteration(tau, n_col4, true, true);
        f5.make_iteration(tau, n_col5, true, true);
        
        out << tau * (i + 1) << '\t' << f1.M_2(1,1) << '\t'
                                     << f2.M_2(1,1) << '\t'
                                     << f3.M_2(1,1) << '\t'
                                     << f4.M_2(1,1) << '\t'
                                     << f5.M_2(1,1) << '\n';

        std::cout << "Step " << i + 1 << " was done\n";
    }
    out.close();

    return 0;
}
