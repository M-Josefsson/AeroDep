/*!
* @file 
* @brief This file contains the Particle class.
*
*/

/*!
*
* @class Particle 
* @brief This class represents a single particle. 
*
* It is used for storing the particle's 
* parameters and for integrating the equations of motion. 
*/

#include "particle.hpp"
#include "vector_extensions.hpp"
#include "fields_and_forces.hpp"
#include "constants.hpp"

#include <fstream>
#include <iostream>
#include <math.h>

using std::vector;
using std::endl;

/*! 
* @brief Creates a Particle at height z_start_. 
* 
* The particle's position is randomized in the xy-plane within the control area (data.control_l^2). 
* 
* @param z_start z-component of the initial position
* @param d_p the particle's diamater
* @param data Data sutrcture containing the user input data and environment data.
* @param r Array contining three uniformely sampled random number used for the y- and z-component of the initial position.
* @param rm Array representing a 3D point uniformely sampled inside a sphere. User for the initial random magnetization.
*/

Particle::Particle (const double& z_start_, const double& d_p, const InputData& data, 
                    const vector3& r, const vector3& rm){
        z_start = z_start_;                                                    
        diameter = d_p;
        q = data.q;
        Xi = data.Xi;        
        V = 4.0/3.0*PI*pow(diameter/2.0, 3);
        mass = data.density*V;    
        Cc = 1.0 + data.mfp/diameter*(2.514 + 0.8*exp(-0.55*diameter*0.5/diameter)); 

        magnetization = rm/norm(rm) * data.m_saturation;
        Set_initial_conditions(data.T, data.v_g, r, data.control_l);                    
}


/*!
*
* @brief Sets the initial position and velocity.
* Sets the starting position of a Particle. z=z_start and the x,y coordinates are set to random position
* using the random numbers provided in r.
*
* @param T Temperature
* @param v_g Carrier gas velocity.
* @param r Vector contining three uniformely sampled random number used for the y- and z-component of the initial position.
* @param control_l Defines the side-length of the simulation volume.
*/
void Particle::Set_initial_conditions(const double& T, const vector3& v_g, const vector3& r, 
                                    const double control_l){
        
    for(int i = 0; i<2; ++i){                
        pos[i] = (r[i]-0.5)*control_l;
    }

    pos[2] = z_start; 

    for(int i = 0; i<3; ++i){        
        vel[i] = v_g[i] + r[2] * sqrt(2*KB*T/mass);     
    }

    iteration = 0;
}


/*!
* @brief Make a forward timestep
*
* Make a forward by timestep increasing time with dt (Euler's method) using the Chandrasekar procedure for particle motion including
* Brownian motion. This routine is based on the one described in the appendix of Zarutskaya and Shapiro JoAS 31, 907 (2000). After the time step dt it checks
* wether the particle has collided with the substrate or another particle.
*
* @param data Data sutrcture containing the user input data and environment data.
* @param frozen_particles Vector containing the already deposited (frozen) particles.
* @param r Array containing 9 random numbers sampled from a normal distribution. Used for Brownian motion.
* @param os Ostream where output (only info) should be written. This is normally std::cout.
*
* @return False if the simulation should be stopped (due to collision or exceeding max number of iterations), otherwise true.
*                        
*/

bool Particle::Step_time(const InputData& data, const vector<Particle>& frozen_particles, const std::array<double,9>& r, ostream& os){    

    double Beta = (3.0*PI*data.eta_g*diameter)/(Cc*mass);
    double eb = exp(-Beta*data.dt);    
    double Bm = Beta*mass;

    double sigma_v = sqrt(KB*data.T/mass*(1.0-eb*eb));
    double sigma_vr = KB*data.T/(mass*Beta)*pow(1.0-eb, 2);
    double sigma_r = 1.0/Beta*sqrt(KB*data.T/mass*(2.0*Beta*data.dt - 3.0 + 4.0*eb-pow(eb,2)));

    vector3 v_new{0.0, 0.0, 0.0};
    vector3 F_ext = Get_total_force(*this, frozen_particles, data); 

    v_new[0] = vel[0]*eb + (data.v_g[0] + F_ext[0])/(Bm)*(1.0-eb) + sigma_v*r[0];
    v_new[1] = vel[1]*eb + (data.v_g[1] + F_ext[1])/(Bm)*(1.0-eb) + sigma_v*r[1];
    v_new[2] = vel[2]*eb + (data.v_g[2] + F_ext[2])/(Bm)*(1.0-eb) + sigma_v*r[2];

    double sigma_sqrt = sqrt( (pow(sigma_r,2) - pow(sigma_vr/sigma_v,2))); //sqrt(abs(sigma_sqrt));
    double sigma_frac = sigma_vr/sigma_v;

    pos[0] = pos[0] + vel[0]/Beta*(1.0-eb)+(data.v_g[0] + F_ext[0]/(Bm)) * (data.dt-1.0/Beta*(1.0-eb)) + r[0]*sigma_frac + r[3]*sigma_sqrt;
    pos[1] = pos[1] + vel[1]/Beta*(1.0-eb)+(data.v_g[1] + F_ext[1]/(Bm)) * (data.dt-1.0/Beta*(1.0-eb)) + r[1]*sigma_frac + r[4]*sigma_sqrt;
    pos[2] = pos[2] + vel[2]/Beta*(1.0-eb)+(data.v_g[2] + F_ext[2]/(Bm)) * (data.dt-1.0/Beta*(1.0-eb)) + r[2]*sigma_frac + r[5]*sigma_sqrt;

    vel = v_new;            
    Check_boundary_conditions(data.control_l, os);

    iteration++;
    if (iteration > MAX_ITER){
        return true;
    }else{ 
        return Check_collision(frozen_particles, data.remove_surface_charge, data.control_l, os);
    }
}

/*!
* Applies periodic boundary conditions in x- and y-direction and 
* resets the z-xomponent if the height exceeds 2*z_start.
*/
void Particle::Check_boundary_conditions(const double& control_l, ostream& os){

    for(int i=0; i<2; ++i){

        if (pos[i] > control_l/2.0){
            pos[i] = pos[i] - control_l;           
        }
        if(pos[i] < -control_l/2.0){
            pos[i] = pos[i] + control_l;          
        }        
    }
    if(pos[2] > 2*z_start){
        pos[2] = z_start;
        os << "z reset" << endl;
    }
}


/*!
* Checks if the Particle is within 'MIN_DIST' from the substrate or from another Particle.
* If it is, its position is frozen and the charge is removed if 'remove_surface_charge' is true.
*
* @return True is the particle has collided with something, otherwise false.
*/      
bool Particle::Check_collision(const vector<Particle>& frozen_particles, const bool& remove_surface_charge, const double& control_l, ostream& os){
    double dist;
    vector3 p_to_p{0.0, 0.0, 0.0}, frozen_pos{0.0, 0.0, 0.0};

    
    for (size_t i = 0; i < frozen_particles.size(); ++i){
        const Particle* frozen_particle = &frozen_particles[i];

        if(fabs(frozen_particle->pos[2]-pos[2]) < (frozen_particle->diameter/2.0 + diameter/2.0)*1.2 ){
            frozen_pos[2] = frozen_particle->pos[2];
            //Include particles in neighbouring boxes

            for(double x_offset : {-1,0,1}){ 
                for(double y_offset : {-1,0,1}){ 
                
                    frozen_pos[0] = frozen_particle->pos[0] + x_offset*control_l;
                    frozen_pos[1] = frozen_particle->pos[1] + y_offset*control_l;                

                    p_to_p = pos - frozen_pos;
                    dist = norm(p_to_p);                
                    
                    if(dist < diameter*0.5 + frozen_particle->diameter*0.5 + MIN_DIST){
                        frozen = true;
                        collision_object = 'p';
                        os << "Collided with a Particle - Distance between particles: " << dist << endl;
                        os << "Steps: " << iteration << endl;

                        if(remove_surface_charge){
                            q = 0;                          
                        }
                        return frozen;                        
                    }            
                }            
            }
        }
    }

    if( (pos[2] - diameter/2.0) < MIN_DIST){
        pos[2] = diameter/2.0;
        frozen = true;
        collision_object = 's';
        os << "Collided with the substrate - Distance to substrate : " << pos[2] << endl;
        os << "Steps: " << iteration << endl;

        if(remove_surface_charge){
            q = 0;            
        }        
    }
    return frozen;
}