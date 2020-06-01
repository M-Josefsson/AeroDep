/*!
* @file 
* @brief This file contains linear algebra vector operations. 
*/

#include "vector_operations.hpp"

#include <assert.h>
#include <algorithm> //transform
#include <cmath>   //pow
#include <array>

#include <iostream>

using std::array;

vector3 operator+(const vector3& a, const vector3& b){

    vector3 result;
    std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::plus<double>());
    return result;
}

Jacobian operator+(const Jacobian& a, const Jacobian& b){

    array<double, 9> result;
    std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::plus<double>());
    return result;
}


vector3 operator-(const vector3& a, const vector3& b)
{
    vector3 result;
    std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::minus<double>());
    return result;
}

vector3 operator*(const vector3& a, const double& b){  
    
    vector3 result;

    for(size_t i=0; i<3; ++i){
        result[i] = a[i]*b;
    }
    return result;
}

vector3 operator/(const vector3& a, const double& b)
{  
    vector3 result;
    result = a * (1.0/b);
    return result;
}

vector3 crossproduct(const vector3& a, const vector3& b){
    vector3 r;

    r[0] = a[1]*b[2] - a[2]*b[1];
    r[1] = a[2]*b[0] - a[0]*b[2];
    r[2] = a[0]*b[1] - a[1]*b[0];    
    return r;
}

double dot(const vector3& a, const vector3& b){

    double r = 0.0;
    for(size_t i=0; i<3; ++i){
        r+= a[i]*b[i];
    }
    return r;
}

double norm(const vector3& a){       
    double res = 0.0;
    for(size_t i = 0; i < 3; ++i){
        res += a[i]*a[i];
    }
    res = sqrt(res);
    return res;
}