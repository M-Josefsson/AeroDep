#ifndef FIELDS_AND_FORCES
#define FIELDS_AND_FORCES

#include "types.hpp"
#include "particle.hpp"
#include "data_struct.cpp"

#include <vector>

using std::vector;

vector3 Get_total_force( Particle& particle, 
                                  const vector<Particle>& frozen_particles, 
                                  const InputData& data);

vector3 E_field_deposited_particle( const Particle& frozen_particle, 
                                    const vector3& p_to_p, 
                                    const double& r, 
                                    const double& E0, 
                                    const array<double, 4>& eps);

Jacobian Gradient_E_field_deposited_particle( const Particle& frozen_particle, 
                                              const vector3& p_to_p , 
                                              const double& r, 
                                              const double& E0, 
                                              const array<double, 4>& eps);

vector3 F_image_particle_substrate( const Particle& p, const array<double, 4>& eps);

vector3 F_image_particle_particle( const Particle& p, 
                                   const Particle& frozen_particle, 
                                   const vector3& p_to_p, 
                                   const double& S, 
                                   const array<double, 4>& eps);

vector3 F_dipole( const Particle& par, 
                  const vector3& E_field, 
                  const Jacobian& dE, 
                  const array<double, 4>& eps);

vector3 F_vdW_particle_substrate( const Particle& p, const double& AH);

vector3 F_vdW_particle_particle( const Particle& p, 
                                 const Particle& frozen_particle, 
                                 const vector3& p_to_p, 
                                 const double& r, 
                                 const double& AH);

vector3 F_ferromagnetic_dipole_dipole( const Particle& p,  
                                       const Particle& frozen_particle, 
                                       const vector3& p_to_p, 
                                       const double& dist);

vector3 F_paramagnetic( const Particle& p, 
                        const vector3& r_diff, 
                        const double& dist, 
                        const vector3& B_ext);  
                                        
vector3 H_field_dipole( const Particle& frozen_particle, 
                        const vector3& p_to_p, 
                        const double& dist);  
#endif