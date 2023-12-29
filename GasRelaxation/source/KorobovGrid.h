#pragma once

#include <vector>
#include <random>

// Korobov's grid
class KorobovGrid : public std::vector<std::vector<double>> {
    private:
        using grid_t = std::vector<std::vector<double>>;
        using coef_t = std::vector<int>;
        // using row_t = std::vector<double>;

        template<typename Numeric, typename Generator = std::mt19937>
        Numeric random(Numeric from, Numeric to) {
            thread_local static Generator gen(std::random_device{}());

            using dist_type = typename std::conditional
            <
                std::is_integral<Numeric>::value
                , std::uniform_int_distribution<Numeric>
                , std::uniform_real_distribution<Numeric>
            >::type;

            thread_local static dist_type dist;

            return dist(gen, typename dist_type::param_type{from, to});
        }

    public:
        static double func_H (int prime, int s, int z);

        static coef_t korobov_coef(int prime, int n);

        void random_shift();

        void shuffle();

        KorobovGrid () = default;
        KorobovGrid (const KorobovGrid& grid);
        KorobovGrid (KorobovGrid&& grid);

        KorobovGrid (const int prime, const coef_t& coef, bool haveShuffle = false, bool haveShift = false);
        KorobovGrid (const int prime, const int n, bool haveShaffle = false, bool haveShift = false);

        const KorobovGrid& operator= (const KorobovGrid& grid);

        const KorobovGrid& operator= (KorobovGrid&& grid);

        ~KorobovGrid() = default;  // TODO: is it OK?
};
