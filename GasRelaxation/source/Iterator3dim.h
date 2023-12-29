#pragma once

#include <cstddef>
#include "Vec3d.h"

// Итератор для трехмерного массива. *it = {i_x, i_y, i_z}. Чтобы получить само значение используйте Сontainer(*it) или Container[*it] 
template<std::size_t N_x, std::size_t N_y, std::size_t N_z, typename ValueType>
// template<typename ValueType>
class Iterator3dim : public std::iterator<std::input_iterator_tag, ValueType> {
    private:
        ValueType ind;

    public:
        Iterator3dim(const ValueType& other_ind) : ind(other_ind) {} // TODO: лучше это отпрваить в private,наваерное
        Iterator3dim(const Iterator3dim& it) : ind(it.ind) {}
        bool operator== (const Iterator3dim& other_it) const {
            return this->ind == other_it.ind;
        }
        bool operator!= (const Iterator3dim& other_it) const {
            return this->ind != other_it.ind;
        }
        ValueType& operator* () {
            return ind;
        }
        Iterator3dim& operator++ () {  // end = {Nx, 0, 0}
            if (ind == ValueType(N_x, 0, 0))
                return *this;

            ind.z++;
            if (ind.z == N_z) {
                ind.z = 0;
                ind.y++;
            }
            if (ind.y == N_y) {
                ind.y = 0;
                ind.x++;
            }
            return *this;
        }
};