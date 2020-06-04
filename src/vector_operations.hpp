#ifndef VECTOR_EXTENSIONS_H
#define VECTOR_EXTENSIONS_H

#include "types.hpp"

#include <array>
#include <algorithm>

using std::array;

vector3 operator+(const vector3& a, const vector3& b);
Jacobian operator+(const Jacobian& a, const Jacobian& b);
vector3 operator-(const vector3& a, const vector3& b);
vector3 operator*(const vector3& a, const double& b);
vector3 operator/(const vector3& a, const double& b);
vector3 crossproduct(const vector3& a, const vector3& b);

double dot(const vector3& a, const vector3& b);
double norm(const vector3& a);
double min(const vector3& a);

#endif