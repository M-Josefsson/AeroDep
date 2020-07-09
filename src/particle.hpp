#ifndef PARTICLE
#define PARTICLE

#include "types.hpp"
#include "data_struct.cpp"

#include <vector>

using std::vector;
using std::ostream;
using std::array;

class Particle{
public:
    Particle ( const double& z_start, const double& d_p, const InputData& data, 
               const vector3& r, const vector3& rm);        

    bool Step_time( const InputData& data, const vector<Particle>& frozen_particles, 
                    const array<double, 6>& r, ostream& os);           

    char Get_collision_object(); 

    int Get_iterations() const {return iteration;} 
    //!< @brief Returns the current number of iterations in the time evolution.

    double diameter; 
    //!< @brief The particle's diameter.
    
    double V; 
    //!< @brief The particle's volume.
    
    double mass; 
    //!< @brief The particle's mass.
    
    double q; 
    //!<@brief The particle's charge.
    
    double Xi; 
    //!<@brief The particle's magnetic susceptibility.
    
    vector3 pos; 
    //!<@brief The particle's position.
    
    vector3 vel; 
    //!<@brief The particle's velocity.
    
    vector3 magnetization; 
    //!<@brief The particle's magnetization.

private:
    void Set_initial_conditions( const double& T, const vector3& v_g, 
                                 const vector3& r, const double control_l);

    bool Check_collision( const vector<Particle>& frozen_particles, const bool& remove_surface_charge, 
                          const double& control_l, const bool verbose, ostream& os);

    void Check_boundary_conditions(const double& control_l, const bool verbose, ostream& os);
    
    char collision_object;   
    //!< @brief 'p': particle-particle collision. 's': particle-substrate collision.

    double z_start;
    //!< @brief The particle's start height.

    double Cc;
    //!< @brief The particle's Cunningham factor.

    bool frozen{false};             
    //!< @brief Has the particle collided with something?

    int iteration{0};
    //!< @brief Current number of iteration in the time evolution for this particle.
};

#endif