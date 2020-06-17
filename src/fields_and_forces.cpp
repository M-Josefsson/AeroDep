/*****************************************************************************************************//**
* @file
* 
* @brief This file contains functions for calculating the two-particle forces and for calculating the 
*  electric and magnetic fields around already deposited particles.
*
*********************************************************************************************************/

#include "fields_and_forces.hpp"
#include "vector_operations.hpp"

#include <cmath>
#include <algorithm> //std::min
#include <iostream>

//Uncomment for parallelized code
//#include <omp.h>
//#pragma omp declare reduction(vector3plus: vector3: omp_out = omp_out + omp_in)
//#pragma omp declare reduction(Jacobianplus: Jacobian: omp_out = omp_out + omp_in)

using std::vector;

/*****************************************************************************************************//**
* @brief Returns the total external force acting on the particle.
*
* The total force includes all two particle forces and all forces between the incoming particle and 
* its surroundings, such as the force due to interactions with the substrate or the external electric 
* field. \n                 
*
* @param particle The incoming particle
* @param frozen_particles Vector containing the already deposited (frozen) particles.
* @param data Data structure containing all the input and environment variables 
* 
* 
* @return  The total force acting on the incoming particle
*
*********************************************************************************************************/
vector3 Get_total_force( Particle& particle, const vector<Particle>& frozen_particles, 
                         const InputData& data){

    Jacobian gradient_E_field;
    gradient_E_field.fill(0.0);
    vector3 E_field{0.0, 0.0, data.E0}, H_field_dipoles{0.0, 0.0, 0.0}, F_ext{0.0, 0.0, 0.0};    

    /*Uncomment for parallelized code  
    //Determine how many threads to use. Based on empirical testing.
    int N = std::min(static_cast<int> ((frozen_particles.size()+2500-1) / 2500), omp_get_max_threads());  
    #pragma omp parallel for if (frozen_particles.size()>3000) shared(frozen_particles, data) \ 
    firstprivate(particle) num_threads(N) \
    reduction(vector3plus: F_ext, E_field, H_field_dipoles) reduction(Jacobianplus: gradient_E_field)    
    */

    for (size_t i = 0; i< frozen_particles.size() ; ++i ){

        const Particle &frozen_particle = frozen_particles[i];
                    
        vector3 frozen_pos{0.0, 0.0, 0.0};
        vector3 p_to_p{0.0, 0.0, 0.0};
        frozen_pos[2] = frozen_particle.pos[2];

        //Include particles in neighboring boxes
        for(double x_offset : {-1.0, 0.0, 1.0}){ 
            for(double y_offset : {-1.0, 0.0, 1.0}){ 
                            
                frozen_pos[0] = frozen_particle.pos[0] + x_offset * data.control_l;
                frozen_pos[1] = frozen_particle.pos[1] + y_offset * data.control_l;                

                p_to_p = particle.pos - frozen_pos;
                double dist = norm(p_to_p);                

                if(dist < data.interaction_length){  

                    E_field = E_field + E_field_deposited_particle( frozen_particle, p_to_p, dist, data.E0, data.eps);                    
                    gradient_E_field = gradient_E_field + Gradient_E_field_deposited_particle( frozen_particle, p_to_p, dist, data.E0, data.eps);

                    F_ext = F_ext + F_image_particle_particle(particle, frozen_particle, p_to_p, dist, data.eps);
                    F_ext = F_ext + F_vdW_particle_particle(particle, frozen_particle, p_to_p, dist, data.AH232); 

                    if(data.calcMagnetic){

                        if (data.magnetic_ferro){

                            H_field_dipoles = H_field_dipoles + H_field_dipole(frozen_particle, p_to_p, dist);
                            F_ext = F_ext + F_ferromagnetic_dipole_dipole(particle, frozen_particle, p_to_p, dist);  

                        }else{
                            F_ext = F_ext + F_paramagnetic(particle, p_to_p, dist, data.B);                                    
                          
                        }
                    }                                    
                }            
            }            
        }
    } 
   
    F_ext = F_ext + F_dipole( particle, E_field, gradient_E_field, data.eps);

    //Electric field contribution to Lorentz force
    F_ext = F_ext + E_field*(particle.q*ELEM);
    
    //Particle-substrate forces
    F_ext = F_ext + F_image_particle_substrate(particle, data.eps);
    F_ext = F_ext + F_vdW_particle_substrate(particle, data.AH132); 
    
    if (data.calcMagnetic){

        //Magnetic field contribution to Lorentz force
        F_ext = F_ext + crossproduct(particle.vel, data.B)*(particle.q*ELEM); 

        if (data.magnetic_ferro){
            
            vector3 H_tot{0.0, 0.0, 0.0};

            for (size_t i=0; i<3; ++i){                
                H_tot[i] = data.B[i]/MU0 + H_field_dipoles[i];                            
            }

            double H_norm = norm(H_tot);       

            if (H_norm > data.alignment_field_strength){                         
                particle.magnetization = H_tot / H_norm * data.m_saturation; 
            }
        }              
    }
    return F_ext;
}


/*****************************************************************************************************//**
* @brief Calculates the electric field around a deposited particle.
*
* The deposited particle is assumed to be a conductive sphere in an external field E0.
* The field is then given by
* @image html E_field.png width=330px
*
* The notation can be found in Krinke et al.
*
* @param frozen_particle An already deposited particle.
* @param p_to_p Vector pointing from the center of an incoming particle to the center of the frozen 
* particle.
* @param r Distance between particle centers (length of p_to_p).  
* @param E0 The strength of the external electric field.
* @param eps Vector containing dielectric permittivities.
*
* @return The electric field at the point p_to_p relative to the origin.
*
*********************************************************************************************************/
vector3 E_field_deposited_particle( const Particle& frozen_particle,
                                    const vector3& p_to_p , const double& r, 
                                    const double& E0, const array<double, 4>& eps){
    double r3 = pow(r, 3);
    double r5 = pow(r, 5);
    double rp3 = pow(frozen_particle.diameter/2.0, 3);

    vector3 E{0.0, 0.0, 0.0};
    E[0] += 3*E0*rp3 * p_to_p[0]*p_to_p[2]/r5 + frozen_particle.q*ELEM/(4*PI*eps[0]*eps[3])*p_to_p[0]/r3;
    E[1] += 3*E0*rp3 * p_to_p[1]*p_to_p[2]/r5 + frozen_particle.q*ELEM/(4*PI*eps[0]*eps[3])*p_to_p[1]/r3;
    E[2] += - E0*rp3 * (pow(p_to_p[0],2) + pow(p_to_p[1],2)-2*pow(p_to_p[2],2))/r5 + 
                                    frozen_particle.q*ELEM/(4*PI*eps[0]*eps[3])*p_to_p[2]/r3;

    return E; 
}


/*****************************************************************************************************//**
* @brief Calculates Jacobian of the electric field around a frozen particle.
*
* The jacobian is used in the first order term in the Taylor expansion of the electric field.
*
* @param frozen_particle An already deposited particle.
* @param p_to_p Vector pointing from the center of an incoming particle to the center of the frozen 
* particle.
* @param r Distance between particle centers (length of p_to_p).  
* @param E0 The strength of the external electric field.
* @param eps Vector containing dielectric permittivities. 
*
* @return Jacobian of the electric field around the frozen particle at point p_to_p relative to the 
* origin.
*
*********************************************************************************************************/
Jacobian Gradient_E_field_deposited_particle( const Particle& frozen_particle, const vector3& p_to_p , 
                                              const double& r, const double& E0, 
                                              const array<double, 4>& eps){
    
    double r2 = pow(r,2);    
    double const1 = 3.0*E0*pow(frozen_particle.diameter/2, 3) / pow(r,5);
    double const2 = frozen_particle.q*ELEM/(4.0*PI*eps[0]*eps[3]) / pow(r,3);    

    double x = p_to_p[0];
    double y = p_to_p[1];
    double z = p_to_p[2];

    Jacobian dE;
    dE.fill(0.0);

    dE[0] = const1*z*(1.0 - 5.0*x*x/r2) + const2*(1.0 - 3.0*x*x/r2);  //dE_x/dx    
    dE[1] = -5*const1*x*y*z/r2 - 3.0*const2*x*y/r2;  //dE_x/dy    
    dE[2] = const1*x*(1.0 - 5.0*z*z/r2) - 3.0*const2*x*z/r2;  //dE_x/dz
    
    dE[3] = dE[1];  //dE_y/dx
    dE[4] = const1*z*(1.0 - 5.0*y*y/r2) + const2*(1.0 - 3.0*y*y/r2);  //dE_y/dy
    dE[5] = const1*y*(1.0 - 5.0*z*z/r2) - 3.0*const2*y*z/r2; //dE_y/dz

    dE[6] = dE[2] + 2.0*const1*x*(1.0 - 5.0*z*z/r2);  //dE_z/dx
    dE[7] = dE[5] + 2.0*const1*y*(1.0 - 5.0*z*z/r2);  //dE_z/dy
    dE[8] = const1/3.0*(4.0*z + 5.0*z*(1.0 - 3.0*z*z/r2)) + const2*(1.0 - 3.0*z*z/r2);  //dE_z/dz

    return dE; 
}


/*****************************************************************************************************//**
*
* @brief Calculates the Image force between the incoming particle and the substrate.
*
* The force is calculated assuming a conducting substrate and an incoming charge, which gives 
* @image html F_image_PS.png width=230px
*
* The notation can be found in Krinke et al.
*
* @param p The incoming particle
* @param eps Vector containing dielectric permittivities.
*
* @return The particle-substrate Image force.
*
*********************************************************************************************************/
vector3 F_image_particle_substrate(const Particle& p, const array<double, 4>& eps){

    vector3 F{0.0, 0.0, 0.0};
    F[2] = -pow(p.q*ELEM, 2)/(4.0*PI*eps[0]*eps[3]*pow(2.0*p.pos[2], 2)) * (eps[2]-eps[3])/(eps[2]+eps[3]);
   
    return F;    
}


/*****************************************************************************************************//**
*
* @brief Calculates the Image force between the incoming particle and a frozen particle.
*
* The particles are assumed to be conductive spheres, potentially with a charge in its center,
* which gives the force
* @image html F_image_PP.png width=650px
*
* The notation can be found in Krinke et al.
*
* @param p The incoming particle.
* @param frozen_particle The frozen (already deposited) particle.
* @param p_to_p Vector pointing from the center of an incoming particle to the center of the frozen 
* particle.
* @param S Distance between particle centers (length of p_to_p).
* @param eps Vector containing dielectric permittivities.
*
* @return The particle-particle Image force.
*
*********************************************************************************************************/
vector3 F_image_particle_particle( const Particle& p, const Particle& frozen_particle, 
                                   const vector3& p_to_p, const double& S, 
                                   const array<double, 4>& eps){

    double S2 = pow(S, 2);
    double dp2 = pow(frozen_particle.diameter, 2);
    vector3 m = p_to_p*(-1.0)/S;    
    
    vector3 F = m * ( pow(p.q*ELEM, 2)*pow(frozen_particle.diameter/S, 3) * 
                            ((8.0*S2 - dp2) / ( pow(4.0*S2 - dp2,2)*2.0)) );

    if (frozen_particle.q != 0){

        dp2 = pow(p.diameter, 2);
        F = F + m * ( pow(frozen_particle.q*ELEM,2)*pow(p.diameter/S, 3) * 
                      ((8.0*S2 - dp2) / ( pow(4.0*S2 - dp2, 2)*2.0)) );
    }

    F = F * (1.0 / (4.0*PI*eps[0]*eps[3]));    
    return F;
}


/*****************************************************************************************************//**
*
* @brief Calculates the dipole force acting on the incoming particle.
* 
* The calculation is based on a first order Taylor expansion of the local electric field. 
* This provides the force
* @image html F_dipole.png width=170px
*
* The notation can be found in Krinke et al.
*
* @param particle The incoming particle 
* @param E_field Local electric field at the center of the incoming particle.
* @param dE Jacobian of the local E_field at the center of the incoming particle.
* @param eps Vector containing dielectric permittivities.
*
* @return The dipole force.
*
*********************************************************************************************************/
vector3 F_dipole( const Particle& particle, const vector3& E_field, 
                  const Jacobian& dE, const array<double, 4>& eps){  

    vector3 p = E_field * 4*PI*eps[0]*eps[3]*pow(particle.diameter/2.0, 3);

    vector3 F{0.0, 0.0, 0.0};
    F[0] = p[0]*dE[0] + p[1]*dE[3] + p[2]*dE[6];
    F[1] = p[0]*dE[1] + p[1]*dE[4] + p[2]*dE[7];
    F[2] = p[0]*dE[2] + p[1]*dE[5] + p[2]*dE[8];

    return F;
}


/*****************************************************************************************************//**
*
* @brief Calculates the van der Waals force between a particle and the substrate.
*
* The force is given by 
* @image html F_vdW_PS.png width=220px
*
* @param particle The incoming particle
* @param AH Hamaker constant AH132 for the substrate (1), the particle (2) and the gas (3).
*
* @return The particle-substrate van der Waals force. 
*
*********************************************************************************************************/
vector3 F_vdW_particle_substrate(const Particle& particle, const double& AH){

    vector3 F{0.0, 0.0, 0.0};
    double S = particle.pos[2] - particle.diameter/2;
    F[2]  = -2.0*AH/3.0*pow(particle.diameter/2.0, 3)/(S*S*pow(S + 2.0*particle.diameter/2.0, 2));
    
    return F;
}


/*****************************************************************************************************//**
*
* @brief Calculates the van der Waals force between two spherical particles.                
* 
* The force is given by
* @image html F_vdW_PP.png width=470px
*
* @param particle The incoming particle.
* @param frozen_particle The frozen particle.
* @param p_to_p Vector pointing from the center of an incoming particle to the center of the frozen 
* particle.
* @param r Distance between particle centers (length of p_to_p).
* @param AH Hamaker constant AH232 for two particles (2) in a gas (3).
* 
* @return The particle-particle van der Waals force.
*
*********************************************************************************************************/
vector3 F_vdW_particle_particle( const Particle& particle, const Particle& frozen_particle, 
                                 const vector3& p_to_p, const double& r, const double& AH){    

    double R1 = particle.diameter/2;
    double R2 = frozen_particle.diameter/2;            
    double S = r - particle.diameter/2 - frozen_particle.diameter/2;

    vector3 F = p_to_p * (-1.0)/r * (32.0*AH/3.0 * pow( R1*R2/(2.0*R1+ S)/(2.0*R2 + S), 2)
                            * (R1*R2 / pow(S, 2)) * ( (R1 + R2 + S) / pow(2.0*R1 + 2.0*R2 + S, 2) ) );  

    F = F + p_to_p * (-1.0)/r * 4.0/3.0*AH*R1*R2*( R1 + R2 + S ) / 
                                        ( S*(2*R1 + S)*(2*R2 + S)*(2*R1 + 2*R2 + S));                                       
    
    return F;
}


/*****************************************************************************************************//**
*
* @brief Calculates the force between two magnetic dipoles.
*
* The force is given by 
* @image html F_ferro.png width=650px
* where r refers to a particle's position.
*
* @param particle The incoming particle
* @param frozen_particle The frozen particle.
* @param p_to_p Vector pointing from the center of an incoming particle to the center of the frozen 
* particle.
* @param dist Distance between particle centers (length of p_to_p)
*
* @return The magnetic dipole-dipole force.
*
*********************************************************************************************************/
vector3 F_ferromagnetic_dipole_dipole( const Particle& particle,  const Particle& frozen_particle, 
                                       const vector3& p_to_p, const double& dist){

    vector3 m1 = particle.magnetization*particle.V;
    vector3 m2 = frozen_particle.magnetization * frozen_particle.V;

    double m1r = dot(m1, p_to_p);
    double m2r = dot(m2, p_to_p);

    vector3 F = m2*m1r + m1*m2r + p_to_p*dot(m1,m2) - p_to_p *(5.0*m1r*m2r/pow(dist, 2));
    F = F*(3.0*MU0/(4.0*PI*pow(dist, 5)));

    return F;
}


/*****************************************************************************************************//**
*
* @brief Experimental feature! Calculates the magnetic force on a weakly magnetic particle.
*
* Warning: this is an experimental feature that is not properly tested.
*
* Calculates the magnetic force on an incoming particle with susceptibility Xi.
* The calculation is based on evaluating the gradient in the magnetic field that arises from
* corrections due to the magnetic moment of the already deposited particle. The particles are assumed to 
* be single domained with magnetization susceptibility Xi. The particles are also assumed to be small 
* enough such that the magnetic field gradient is homogeneous through the particle. The force calculation
* is based on Mikkelsen et al. JoMaMM 293, 578 (2005).
*
* The force is given by
* @image html F_para.png width=650px
*
* @param particle The incoming particle.
* @param p_to_p Vector pointing from the center of an incoming particle to the center of the frozen 
* particle.
* @param dist Distance between particle centers (length of p_to_p).
* @param B_ext The external B-field.
*
* @return The magnetic force.  
*
*********************************************************************************************************/
vector3 F_paramagnetic( const Particle& particle, const vector3& p_to_p, 
                        const double& dist, const vector3& B_ext){

    vector3 H = B_ext/MU0;   
    double HdotR = dot(H, p_to_p);
    double HdotH = dot(H, H);    

    vector3 F{0.0, 0.0, 0.0}, f1{0.0, 0.0, 0.0}, f2{0.0, 0.0, 0.0};

    //M1DdH2
    f1 = f1 - p_to_p*(15.0*pow(HdotR, 2)/pow(dist, 3));        
    f1 = f1 + p_to_p*(3.0*HdotH/dist);    
    f1 = f1 + H*(3.0*HdotR/dist);
    f1 = f1 * 3.0*pow(particle.Xi/(particle.Xi+3.0), 2) * pow(particle.diameter/2.0, 3) / pow(dist, 4);    
    F = F + f1;    

    //dM1DdH2
    f2 = f2 - p_to_p*(12.0*pow(HdotR, 2)/pow(dist, 3));
    f2 = f2 + H*(HdotR*3.0/dist);
    f2 = f2 - p_to_p*(HdotH*3.0/dist);    
    f2 = f2 * 3.0 * pow(particle.Xi/(particle.Xi+3.0), 3) * pow(particle.diameter/(2.0*dist), 6) / dist;
    F = F + f2;   

    F = F * MU0 * particle.V;
    return F;    
}


/*****************************************************************************************************//**
*
* @brief Calculates the magnetic field around a magnetic dipole.
*
* The field at position r around a magnetic dipole with magnetic moment m is 
* @image html H_dipole.png width=270px
*
* @param frozen_particle The frozen particle.
* @param p_to_p Vector pointing from the center of an incoming particle to the center of the frozen 
* particle.
* @param dist Distance between particle centers (length of p_to_p).
*
* @return The local H-field around a magnetic dipole. 
*
*********************************************************************************************************/
vector3 H_field_dipole(const Particle& frozen_particle, const vector3& p_to_p, const double& dist){

    vector3 H_field = p_to_p * 3.0 * dot(frozen_particle.magnetization, p_to_p)/pow(dist, 2) 
                                                                    - frozen_particle.magnetization;    
    H_field = H_field * frozen_particle.V / (4.0*PI*pow(dist, 3) );

    return H_field;
}