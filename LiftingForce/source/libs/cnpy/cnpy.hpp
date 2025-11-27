#pragma once

#include <vector>
#include <string>

namespace npy {
    
    template<typename T>
    void save_npy(
        const std::string& filename,
        const T* data,
        const std::vector<size_t>& shape,
        bool fortran_order = false
    );
    
}

#include "cnpy.tpp"
