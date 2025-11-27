#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <type_traits>

namespace npy {

    template<typename T>
    std::string descr();

    template<> inline std::string descr<double>() { return "<f8"; }
    template<> inline std::string descr<float>()  { return "<f4"; }
    template<> inline std::string descr<int>()    { return "<i4"; }
    template<> inline std::string descr<uint8_t>(){ return "|u1"; }


    template<typename T>
    std::string GetPreambule(
        const std::vector<size_t>& shape,
        bool fortran_order
    ){
        // magic string
        const std::string magic_str("\x93NUMPY");

        // version 1.0
        uint8_t major(1), minor(0);
        std::string version;
        version.push_back(major);
        version.push_back(minor);

        // header (python dict)
        std::string header;
        header =  "{'descr': '" + descr<T>() + "', 'fortran_order': " + (fortran_order ? "True" : "False");
        header += ", 'shape': (";
        if (shape.size() == 1) {
            header += std::to_string(shape[0]) + ",)";
        } else {
            for (size_t i = 0; i < shape.size(); ++i) {
                header += std::to_string(shape[i]);
                if (i + 1 != shape.size())
                    header += ',';
            }
            header += "),}";
        }

        // space padding
        const size_t divider = 8;
        size_t preambule_len = magic_str.size() + version.size() + 2 + header.size() + 1; // 2 bytes = header_len after padding, 1 byte = '\n'
        size_t remains = preambule_len % divider;
        if (remains != 0) {
            size_t padding_len = divider - remains;
            header.append(padding_len, ' '); // padding with spaces before '\n'
        }
        header.push_back('\n');

        // length
        std::string len_str; // in little endian
        len_str.push_back(header.size() & 0xFF);        // low byte
        len_str.push_back((header.size() >> 8) & 0xFF); // hi  byte

        return magic_str + version + len_str + header;
    }


    template<typename T>
    void save_npy(
        const std::string& filename,
        const T* data,
        const std::vector<size_t>& shape,
        bool fortran_order
    ){
        static_assert(std::is_trivially_copyable<T>::value, "T must be trivially copyable");

        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open()) throw std::runtime_error("Cannot open file: " + filename);

        std::string preambule = GetPreambule<T>(shape, fortran_order);
        out.write(preambule.data(), preambule.size());

        size_t total = 1;
        for (size_t d : shape)
            total *= d;
        out.write(reinterpret_cast<const char*>(data), total * sizeof(T));
    }
}

