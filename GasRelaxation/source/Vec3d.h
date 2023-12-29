#pragma once

#include <cmath>  // TODO: include guard?
#include <iostream>

template <typename Type>
struct vec3d_t {
    Type x, y, z;

    vec3d_t(): x(), y(), z() {}
    vec3d_t(const Type& a): x(a), y(a), z(a) {}
    vec3d_t(const Type& x0, const Type& y0, const Type& z0): x(x0), y(y0), z(z0) {}

    vec3d_t(const vec3d_t& vec) : x(vec.x), y(vec.y), z(vec.z) {}
    vec3d_t(vec3d_t&& vec) : x(std::move(vec.x)), y(std::move(vec.y)), z(std::move(vec.z)) {}

    const vec3d_t& operator= (const vec3d_t& vec) {
        if (this != &vec) {
            x = vec.x;
            y = vec.y;
            z = vec.z;
        }
        return *this;
    }

    const vec3d_t& operator= (vec3d_t&& vec) {
        if (this != &vec) {
            x = std::move(vec.x);
            y = std::move(vec.y);
            z = std::move(vec.z);
        } 
        return *this;
    }

    bool operator== (const vec3d_t& vec) const {
        return (x == vec.x) && (y == vec.y) && (z == vec.z);
    }

    bool operator != (const vec3d_t& vec) const {
        return !(*this == vec);
    }

    const vec3d_t operator+ (const vec3d_t& vec) const {
        vec3d_t res(x + vec.x, y + vec.y, z + vec.z);
        return res;
    }

    vec3d_t& operator+= (const vec3d_t& vec) {
        x += vec.x;
        y += vec.y;
        z += vec.z;
        return *this;
    }

    const vec3d_t operator- (const vec3d_t& vec) const {
        vec3d_t res(x - vec.x, y - vec.y, z - vec.z);
        return res;
    }

    vec3d_t& operator-= (const vec3d_t& vec) {
        x -= vec.x;
        y -= vec.y;
        z -= vec.z;
        return *this;
    }

    const vec3d_t operator* (const Type& a) const {
        vec3d_t res(x * a, y * a, z * a);   
        return res;
    }

    const vec3d_t operator/ (const Type& a) const {
        vec3d_t res(x / a, y / a, z / a);
        return res;
    }

    friend std::ostream& operator<< (std::ostream& out, const vec3d_t& vec) {
        out << "{" << vec.x << " " << vec.y << " " << vec.z << "}";
        return out;
    }

    Type& operator[] (int i) {  // TODO: может и нужна эта функия???
        switch (i)
        {
            case 1:
                return this->x;
                // break;
            case 2:
                return this->y;
            case 3:
                return this->z;
            default:
                std::cout << "Wrong index " << i << " in Type& vec3d_t operator[]. Use 1, 2 or 3\n";
                return this->x;
        }
    }
    Type operator[] (int i) const {
        switch (i)
        {
        case 1:
            return x;
        case 2:
            return y;
        case 3:
            return z;
        default:
            std::cout << "Wrong index " << i << " in Type vec3d_t operator[] const. Use 1, 2 or 3\n";
            return x;
        }
    }

    Type squared() const {
        return x * x + y * y + z * z; 
    }

    Type norm() const {
        return std::sqrt(this->squared());
    }

    Type xy() const {
        return std::sqrt(x * x + y * y);
    }

    ~vec3d_t() {}
};

template<typename Type> 
const vec3d_t<Type> operator* (const Type& a, const vec3d_t<Type>& vec) {
    return vec * a;
}

template<typename Type>
const vec3d_t<Type> operator- (const vec3d_t<Type>& vec) {
    vec3d_t<Type> res(-vec.x, -vec.y, -vec.z);
    return res;
}
