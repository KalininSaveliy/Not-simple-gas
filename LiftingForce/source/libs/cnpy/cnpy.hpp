#pragma once
#include <fstream>
#include <vector>
#include <string>
#include <type_traits>
#include <cstdint>
#include <sstream>

namespace npy {

template<typename T>
std::string type_descr();

template<> inline std::string type_descr<double>() { return "<f8"; }
template<> inline std::string type_descr<float>()  { return "<f4"; }
template<> inline std::string type_descr<int>()    { return "<i4"; }
template<> inline std::string type_descr<uint8_t>(){ return "|u1"; }

// ----------------------------------------------------------

template<typename T>
void save_npy(
    const std::string& filename,
    const T* data,
    const std::vector<size_t>& shape,
    bool fortran_order = false
){
    static_assert(std::is_trivially_copyable<T>::value, "T must be trivially copyable");

    // ------------------- build shape string -------------------
    std::ostringstream oss_shape;
    oss_shape << "(";
    if (shape.size() == 1) {
        oss_shape << shape[0] << ",)";
    } else {
        for (size_t i = 0; i < shape.size(); ++i) {
            oss_shape << shape[i];
            if (i + 1 != shape.size()) oss_shape << ",";
        }
        oss_shape << ",)"; // trailing comma TODO: double comma
    }
    std::string shape_str = oss_shape.str();

    // ------------------- python dict header -------------------
    std::ostringstream oss;
    oss << "{'descr': '" << type_descr<T>()
        << "', 'fortran_order': "
        << (fortran_order ? "True" : "False")
        << ", 'shape': " << shape_str << ", }";

    std::string header = oss.str();

    // Must end with '\n'
    header.push_back('\n');

    // ------------------- padding to 16 bytes -------------------
    const size_t magic_and_version = 10; // 6 magic + 2 ver + 2 header_len
    size_t header_len = header.size();
    size_t rem = (magic_and_version + header_len) % 16;
    if (rem != 0) {
        size_t pad = 16 - rem;
        // pad with spaces before '\n'
        header.insert(header.end() - 1, pad, ' ');
        header_len = header.size();
    }

    // ------------------- open file -------------------
    std::ofstream out(filename, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot open file: " + filename);

    // ------------------- write magic -------------------
    const char magic[] = "\x93NUMPY";
    out.write(magic, 6);

    // version 1.0
    uint8_t major = 1, minor = 0;
    out.write(reinterpret_cast<char*>(&major), 1);
    out.write(reinterpret_cast<char*>(&minor), 1);

    // header length (uint16, little-endian)
    uint16_t hlen = static_cast<uint16_t>(header_len);
    out.write(reinterpret_cast<char*>(&hlen), 2);

    // header text
    out.write(header.data(), header_len);

    // ------------------- write raw data -------------------
    // contiguous write
    size_t total = 1;
    for (size_t d : shape) total *= d;
    out.write(reinterpret_cast<const char*>(data),
              total * sizeof(T));
}

} // namespace npy
