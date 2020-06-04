/*****************************************************************************************************//**
* @file 
*
* @brief This file contains the class InputData.
*********************************************************************************************************/

/*****************************************************************************************************//**
*
* @class InputData 
*
* @brief This class stores the values of all parameters used in the simulation.
*
* Upon instantiation default values for each parameter is used, which can then be overwritten.
* The default values are usually based upon Iron particles, in a Nitrogen carrier gas that will be 
* deposited on a Silicon substrate.
*
*********************************************************************************************************/

#ifndef INPUTDATA
#define INPUTDATA

#include "constants.hpp"
#include "types.hpp"

#include <array>
#include <math.h>

struct InputData {
    
    int nbr_of_particles {100};
    //!< Total number of particles to be generated.

    double q {-1.0}; 
    //!< Charge of the particles.

    double z_start {150e-9}; 
    //!< Height at which new particles are generated.

    double control_l {1.5e-6}; 
    //!< Defines the simulation area in the x-y plane (control_^2).

    double T {300.0}; 
    //!< Temeprature.

    double diameter {30e-9}; 
    //!< Particle (mean) diamater.

    double dt {1e-9}; 
    //!< Time step.

    double E0 {300e3}; 
    //!< Strength of external electric field in the z-direction.

    double n_substrate {1.4585}; 
    //!< Refractive index of the substrate. Default: Silicon.
    
    double n_gas {1.00}; 
    //!< Refractive index of the gas. Default: Nitrogen.

    double eta_g {18.13e-6}; 
    //!< Dynamic viscocity(?) of the gas.

    double mfp {66.5e-9};
    //!< particle mean free path.

    double Xi {-2.2e-5}; 
    //!< The particle's magnetic susceptibility.
    
    double interaction_length {500e-9}; 
    //!< Particle-particle forces are only calculated for distances below this value.

    double density {7310.0}; 
    //!< Particle density. Default: Iron.

    double m_saturation {1.713e6}; 
    //!< Particle's (saturation) magnetisation. Default: Iron (bulk).

    double diameter_std {0.0}; 
    //!< Standard deviation used for a log-norm distribution for diameters.

    double diameter_std2 {0.0}; 
    //!< Standard deviation used for a log-norm distribution for diameters of doubly charged particles.
    
    double double_charge_fraction {0.0}; 
    //!< Fraction of particles with charge 2*q.

    double alignment_field_strength {1e-5}; 
    //!< H-field strength above which the particles magnetization aligns with the local field.

    std::array<double, 4> eps {8.82e-12, 3.9, 3.9, 1.0}; 
    //!< Dielectric permittivities. 0 - vacuum, 1 - substrate, 2 - substrate, 3 - gas.

    vector3 B {0.0, 0.0, 0.0}; 
    //!< External magnetic field.

    vector3 v_g {0.0, 0.0, 0.0}; 
    //!< Gas velocity.


    bool print_trajectory {false}; 
    //!< Print trajectories for every particle?
    
    bool remove_surface_charge {true};
    //!< Remove a particle's charge upon collision?

    bool calcMagnetic {true}; 
    //!< Include magnetic forces?

    bool magnetic_ferro {true}; 
    //!< Use ferromagnetic force calculation? If false, paramegnetic force calculation is used.

    double AH131;
    //!< Hamaker constant of the substrate material.

    double AH132;
    //!< Hamaker constant between the substrate (1) and a particle (2) in a gas 3.

    double AH232 {4e-19}; 
    //!< Hamaker constant for two metal particles in a gas.

    /**************************************************************************************************//**
    *
    * @brief Calculates and sets the values of the Hamaker constants (AH131 and AH132). 
    *
    * The values are based on refractive indicies, relative permittivities, and termperature. This function 
    * must be called for the Hamaker constants to be defined. 
    *
    ******************************************************************************************************/
    void SetHamakerConstants(){  
        double n1 = n_substrate;
        double n3 = n_gas;  
        AH131 = 3.0/4.0*KB*T*pow((eps[1] - eps[3])/(eps[1] + eps[3]), 2) + 
                        3.0*H*NU/(16.0*sqrt(2.0))* pow(n1*n1 - n3*n3, 2)/pow(n1*n1 + n3*n3, 3.0/2.0);    
        AH132 = sqrt(AH232*AH131); 
    }
};

#endif