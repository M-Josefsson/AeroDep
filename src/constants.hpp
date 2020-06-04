/*!
* @file
*
* @brief This file contains all physical and numerical constants.
*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

const double ELEM = 1.6021766208e-19;
//!< @brief Elementary charge.

const double PI = 3.141592653589793;
//!< @brief pi

const double H = 6.62607015e-34;
//!< @brief Planck's constant.

const double KB = 1.380649e-23;
//!< @brief Boltzmann's constant.

const double MU0 = 4.0*PI*1e-7;
//!< @brief Magnetic permeability of vacuum.

const double NU = 3e15;
//!< @brief Main absorption frequency of the gas in UV range.

const double MIN_DIST = 0.5e-9;
//!< @brief Particles closer than this value will be considered to have collided.

const int MAX_ITER = 1000000;
//!< @brief Maximum number of iterations (of the time evolution) for a single particle.

#endif