#include "KorobovGrid.h"

#include <iostream>
#include <algorithm>

double KorobovGrid::func_H (int p, int s, int z) {  // p is prime
    double sum = 0;
    for (int i = 1; i <= (p - 1) / 2; ++i) {
        double product = 1;
        int r = (i * z) % p;
        for (int j = 1; j <= s; ++j) {
            double trash;
            double m = std::modf(static_cast<double>(r) / p, &trash);
            product *= (1 - 2 * m) * (1 - 2 * m);
            r = (r * z) % p;
        }
        sum += product;
    }
    return sum;
}

KorobovGrid::coef_t KorobovGrid::korobov_coef(int prime, int n) {
    if (n == 8) {  // TODO: make mith switch
        if (prime == 1009) {
            coef_t a{1, 253, 442, 836, 627, 218, 668, 501};
            return a;
        }
        if (prime == 5003) {
            coef_t a{1, 2362, 699, 48, 3310, 3534, 2304, 3787};
            return a;
        }
        if (prime == 10007) {
            coef_t a{1, 436, 9970, 3882, 1369, 6471, 9389, 741};
            return a;
        }
        if (prime == 11621) {  // nv = 10^3
            coef_t a{1, 4541, 5027, 3963, 6675, 3607, 5398, 3629};
            return a;
        }
        if (prime == 25013) {
            coef_t a{1, 2647, 2969, 4861, 10385, 24821, 17049, 5251};
            return a;
        }
        if (prime == 50021) {
            coef_t a{1, 11281, 7537, 39218, 32534, 11977, 5816, 32765};
            return a;
        }
        if (prime == 85009) {
            coef_t a{1, 18332, 21647, 10792, 23001, 9692, 4934, 512};
            return a;
        }
        if (prime == 85049) {  // nv = 20^3
            coef_t a{1, 25190, 70560, 52398, 30189, 37801, 83635, 16871};
            return a;
        }
        if (prime == 100003) {
            coef_t a{1, 20285, 68883, 49739, 25348, 68757, 93907, 46351};
            return a;
        }
        if (prime == 125003) {
            coef_t a{1, 14571, 58947, 21124, 40418, 41545, 87669, 19342};
            return a;
        }
        if (prime == 200003) {
            coef_t a{1, 47369, 188507, 54145, 156036, 158419, 37051, 42494};
            return a;
        }
        if (prime == 250007) {
            coef_t a{1, 5532, 102170, 188620, 166629, 15819, 8258, 181982};
            return a;
        }
        if (prime == 289951) {  // nv = 30^3
            coef_t a{1, 5861, 137103, 106462, 289181, 126246, 262805, 80393};
            return a;
        }
        if (prime == 300017) {
            coef_t a{1, 81575, 103565, 136172, 101475, 54078, 262899, 170731};
            return a;
        }
        if (prime == 500009) {
            coef_t a{1, 42535, 193663, 307439, 182488, 487373, 37415, 418387};
            return a;
        }
        if (prime == 1000003) {
            coef_t a{1, 417564, 171019, 163483, 410620, 615303, 611111, 188073};
            return a;
        }
        if (prime == 1500007) {
            coef_t a{1, 413996, 388189, 943278, 1496508, 434758, 733031, 985685};
            return a;
        }
        if (prime == 2000003) {
            coef_t a{1, 832685, 1269182, 1228431, 532894, 174792, 458201, 527381};
            return a;
        }
        if (prime == 3000017) {
            coef_t a{1, 368334, 166765, 2871452, 407635, 979274, 1865572, 2703215};
            return a;
        }
        if (prime == 4000037) {
            coef_t a{1, 72362, 210611, 92212, 583028, 681897, 2974319, 1680656};
            return a;
        }
        if (prime == 6000011) {
            coef_t a{1, 1323844, 1723313, 1392620, 251332, 5750225, 908859, 5328166};
            return a;
        }
        if (prime == 8000009) {
            coef_t a{1, 93973, 6914802, 3957321, 907968, 4380879, 3879127, 4791477};
            return a;
        }
        if (prime == 10000019) {
            coef_t a{1, 1833663, 3609180, 7252140, 5522715, 8914182, 4652083, 6262402};
            return a;
        }
    }

    std::vector<double> h((prime - 1) / 2 + 1);
    double h_min = prime;
    int ind_min = 0;
    for (int i = 1; i < h.size(); ++i) {
        h[i] = func_H(prime, n, i);
        if (h_min > h[i]) {
            h_min = h[i];
            ind_min = i;
        }
    }
    coef_t coef(n);
    coef[0] = 1;
    coef[1] = ind_min;
    int r = ind_min;
    for (int i = 2; i < coef.size(); ++i) {
        r = (r * ind_min) % prime;
        coef[i] = r;
    }
    std::cout << "I recommend saving these coefficients so as not to waste time calculating them again (it's a long time)\n";
    std::cout << "You can use another consturct with this coeficients\n";
    std::cout << "Pime number is " << prime << ", dimension of space is " << coef.size() << '\n' << '{';
    for (int i = 0; i < coef.size() - 1; ++i) {
        std::cout << coef[i] << ", ";
    }
    std::cout << coef[coef.size() - 1] << "}\n\n";
    return coef;
}

void KorobovGrid::random_shift() {
    if (this->size() > 0) {
        int n = (*this)[0].size();  // длина строки в сетке 
        double rand_shift_ar[n];
        for (int i = 0; i < n; ++i) {
            rand_shift_ar[i] = random<double>(0.0, 1.0);
            // std::cout << "Random shift: " << rand_shift_ar[i] << '\n';  // TODO: delete after debug
        }

        double trash;
        for (int k = 0; k < this->size(); ++k) {
            for (int i = 0; i < n; ++i) {
                (*this)[k][i] = std::modf(rand_shift_ar[i] + (*this)[k][i], &trash);
            }
        }
    }
}

void KorobovGrid::shuffle() {
    std::random_shuffle(this->begin(), this->end());
}

KorobovGrid::KorobovGrid (const KorobovGrid& grid) : grid_t(grid) {}
KorobovGrid::KorobovGrid (KorobovGrid&& grid) : grid_t(std::move(grid)) {}

KorobovGrid::KorobovGrid (const int prime, const coef_t& coef, bool haveShuffle, bool haveShift)
            : grid_t(prime, std::vector<double>(coef.size())) {

    double trash;
    for (int k = 0; k < prime; ++k) {
        for (int i = 0; i < coef.size(); ++i) {
            (*this)[k][i] = std::modf(static_cast<double>(coef[i]) * (k + 1) / prime, &trash);
        }
    }
    if (haveShuffle)
        this->shuffle();
    if (haveShift)
        this->random_shift();
}

KorobovGrid::KorobovGrid (const int prime, const int n, bool haveShaffle, bool haveShift)
            : KorobovGrid(prime, korobov_coef(prime, n), haveShaffle, haveShift) {}

const KorobovGrid& KorobovGrid::operator= (const KorobovGrid& grid) {
    if (this != &grid)
        *this = grid;
    return *this;
}
const KorobovGrid& KorobovGrid::operator= (KorobovGrid&& grid) {
    if (this != &grid) {
        *this = std::move(grid);
    }
    return *this;
}
