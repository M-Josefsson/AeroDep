/*!
* @file
* @brief This file contains custom type names. 
*/

#ifndef TYPES_H
#define TOPES_H

#include<array>

typedef std::array<double, 3> vector3;
//!< @brief Array with three components representing a 3D vector.

typedef std::array<double, 9> Jacobian;
//!< @brief Array with nine components used to represent the Jacobian of a vector field.

#endif