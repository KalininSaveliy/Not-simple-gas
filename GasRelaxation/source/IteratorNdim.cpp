// # pragma once

#include <cstddef>

// debug include
#include <iostream>

template<typename ValueType, std::size_t... dims>
class IteratorNdim {
    public:
    int dd[sizeof...(dims)] = {dims...};
    IteratorNdim() {
        const int size = sizeof...(dims);
        for (int i = 0; i < size; ++i) {
            std::cout << dd[i] << ' ';
        }
        std::cout << "size: " << size << '\n';
        // std::cout << (dims) << '\n';
    }
};


int main() {
    std::cout << __cplusplus << '\n';
    int x = 10;
    IteratorNdim<char, 2, 2, 3, 4, 5, 0> it;
    x++;

    return 0;
}
