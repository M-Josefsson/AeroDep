/*****************************************************************************************************//**
* @file 
*
* @brief This file contains the Deposition class.
*
*********************************************************************************************************/

/*****************************************************************************************************//**
*
* @class Deposition 
*
* @brief This class controls the generation of new particles and is the class that the 
* main program (or the user) typically talks to.
*
* After creating a new particle the Deposition class repeatedly calls the time evolution 
* step in the Particle class until collision. 
* A Deposition object keeps track of all the previously generated (now frozen) particles.
* This class also writes the outputs (particle positions and trajectories) to the respective files.
*
*********************************************************************************************************/

#include "deposition.hpp"
#include "constants.hpp"

#include <math.h>  //pow, sqrt, log
#include <ctime>   //time(0)
#include <fstream>
#include <iostream>
#include <stdexcept>

using std::vector;
using std::ostream;
using std::endl;
using std::string;

/*****************************************************************************************************//**
*
* @brief Creates the deposition object responsible for adding and storing particles as 
* well as prints the results.
*
* The generators for the different distributions used are instantiated here. They are all 
* seeded using time(0). 
*
* @param data Data structure containing the user input data and environment data. 
*
*********************************************************************************************************/
Deposition::Deposition(const InputData& data){

    z_start = data.z_start;
    z_start_original = data.z_start;
    diameter = data.diameter;
    frozen_particles.reserve(data.nbr_of_particles);
    
    //Check if log-norm distribution and/or doubly charged particles should be included.
    rand_size = data.diameter_std > 0.0;
    rand_size2 = data.diameter_std2 > 0.0;
    double_charge = data.double_charge_fraction > 1e-4;        

    if ( rand_size ) {
        distributions.Setup_lognorm((data.diameter), data.diameter_std, Random::CHARGE::SINGLE);
    }

    //Setup for sampling for doubly charged particles
    if ( double_charge ){

        d_p2 = Double_charge_diameter(data.diameter,  data.mfp);
        d_p2 *= 1e-9;

        if (rand_size2){
            distributions.Setup_lognorm((d_p2), data.diameter_std2, Random::CHARGE::DOUBLE);
        }
    }
}


/*****************************************************************************************************//**
*
* @brief Adds the already deposited particles supplied as input to the vector of frozen particles.
* 
* @param input_particles Vector of arrays of length 7 where elements 0-2 contains the position, 3-5 the 
* magnetization and 6 the diameter of a particle.  
* @param data Data structure containing the user input data and environment data. 
* @param os Ostream where output (only info) will be written. This is normally std::cout.
*
*********************************************************************************************************/
void Deposition::Add_input_particles( const vector<array<double, 7>>& input_particles, 
                                      const InputData& data, ostream& os){
   
    vector3 rm{0.0, 0.0, 0.0}, r{0.0, 0.0, 0.0};
    double q;

    if(data.remove_surface_charge){
        q = 0.0;
    }else{
        q = data.q;
    }    

    for(size_t i=0; i < input_particles.size(); ++i){
        
        Particle P( data.z_start, input_particles[i][6], data, r, rm);

        P.pos = {input_particles[i][0], input_particles[i][1], input_particles[i][2]};
        P.magnetization = {input_particles[i][3], input_particles[i][4], input_particles[i][5]};
        P.q = q;

        frozen_particles.push_back(P);
        Update_z_start(os);
    }
}


/*****************************************************************************************************//**
*
* @brief Spawns a new particle and performs the time evolution until collision.
*
* Creates a new particle with attributes based on the input in data. Evolves the time (by 
* calling Particle.step_time()) until the particle collided with another particle or the substrate. 
* If data.print_trajectories is True the trajectory of the particle is printed to a new file
* in ./trajectories/.
*
* @param os Ostream where output (only info) will be written. This is normally std::cout.
* @param data Data structure containing the user input data and environment data. 
*
* @return True if the time evolution was successful, otherwise returns False.
*
*********************************************************************************************************/
bool Deposition::Add_particle(ostream& os, const InputData& data){   

    vector3 r1, r3;
    array<double, 6> r2;
    vector<vector3> Pos, magnetization;
    bool collided = false;

    //Generate random numbers for the inial position (x-y-plane) and magnetization.
    distributions.Fill_uniform(r1);
    r3 = distributions.Generate_in_point_sphere();

    //Calculate diamater based on user settings
    double current_q  = data.q;
    double diameter_new = Get_diameter(current_q, data.double_charge_fraction);

    //Create particle
    Particle P( z_start, diameter_new, data, r1, r3);
    os << "Particle created" << endl;

    //Perform the time evolution until collision.
    while(!collided){
        
        distributions.Fill_normal(r2);

        collided = P.Step_time( data, frozen_particles, r2, os);

        if (data.print_trajectory){   

            Pos.push_back(P.pos);
            magnetization.push_back(P.magnetization);
        }
    }    

    if (P.Get_iterations() > MAX_ITER){
        os << "Warning: iterations exceeded one million." << endl;
        return false;
    }
    
    if (data.print_trajectory){
        print_trajectory(Pos, magnetization);   
    }

    frozen_particles.push_back(P);
    Update_z_start(os);
    
    return true;
}


/*****************************************************************************************************//**
*
* @brief Returns the diameter of a new particle. 
*
* The diameter is calculated based on the user input (fraction of double charged particles
* and standard deviations for particle diameters). 
*
* @param current_q The user defined charge of the particles.
* @param double_charge_fraction Fraction of particles having charge 2*current_2.
*
* @return The diameter for the new particle.
*
*********************************************************************************************************/
double Deposition::Get_diameter(double& current_q, const double& double_charge_fraction){

    if (double_charge){

        double r = distributions.uniform();

        if ( double_charge_fraction > r ){

            current_q *= 2.0;

            if (rand_size2){
                return distributions.log_norm(Random::CHARGE::DOUBLE);
            }else{
                return d_p2;
            }
        }
    }

    if (rand_size){
        return distributions.log_norm(Random::CHARGE::SINGLE);
    }else{
        return diameter;
    }
}


/*****************************************************************************************************//**
*
* @brief Calculates the diameter of a doubly charged particle.
*
* Iteratively solves for the diameter of a doubly charged particle. The solution is based in 
* all particles having the same mobility.
*
* @param d The user defined (mean) diameter for a particle with single charge.
* @param mfp the particles' mean free path. 
*
* @return The diameter.
*
*********************************************************************************************************/
double Deposition::Double_charge_diameter(double d, double mfp){
    
    double f, df, Cc1, Z, x, dx;

    d = d/1e-9;
    x = 2.0*d;
    dx = 1e-3;
    mfp = mfp/1e-9;

    Cc1 = 1.0 + mfp/d*(2.514 + 0.8*exp(-0.55*d*0.5/mfp));
    Z = Cc1/d;
    
    f = Mobility(x, mfp, Z); 

    while(fabs(f) > 1e-4){

        df = (Mobility(x+dx, mfp, Z) - f)/dx;        
        x = x - f/df;        
        f = Mobility(x, mfp, Z);
    }    

    return x;

}


/*****************************************************************************************************//**
*
* @brief Calculates the mobility of a particle.
*
* @param d The particle diameter.
* @param mfp the particles' mean free path. 
* @param Z Cc/diameter
*
* @return The mobility.
*
*********************************************************************************************************/
double Deposition::Mobility(const double& d, const double& mfp, const double& Z){

    double Cc2 = 1.0 + mfp/d*(2.514 + 0.8*exp(-0.55*d*0.5/mfp));    
    return 2.0*Cc2/d - Z;
}


/*****************************************************************************************************//**
*
* @brief Update the starting height for new particles.
*
* Updates the starting height of the particles such that it always is a+z_start_ where a is the center 
* of the Particle with largest z-component. 
*
* @param os Ostream where output (only info) will be written. This is normally std::cout.
*
*********************************************************************************************************/
void Deposition::Update_z_start(ostream& os){

    Particle *Particle = &frozen_particles.back();

    if( (Particle->pos[2] + Particle->diameter/2.0) + z_start_original > z_start ){

        z_start = Particle->pos[2] + Particle->diameter/2.0 + z_start_original;
        os << "new z_start = " << z_start << endl;
    }
}


/*****************************************************************************************************//**
*
* @brief Prints the final particle positions (of the frozen particles) to a file. This function 
* is typically called when all particles have been generated.
*
* @param filename The name of the file where the final positions will be printed.
* @param os Ostream where output (only info) will be written. This is normally std::cout.
*
*********************************************************************************************************/
void Deposition::Print_final_positions(const string& filename, ostream& os){

    std::ofstream file;
    file.open(filename);

    if (!file.is_open()){
        os << "Error: Could not open " << filename << endl;
        return;
    }

    for(size_t i = 0; i < frozen_particles.size(); ++i){

        file << frozen_particles[i].pos[0] << "\t"; 
        file << frozen_particles[i].pos[1] << "\t";
        file << frozen_particles[i].pos[2] << "\t";
        file << frozen_particles[i].magnetization[0] << "\t"; 
        file << frozen_particles[i].magnetization[1] << "\t";
        file << frozen_particles[i].magnetization[2] << "\t";
        file << frozen_particles[i].diameter << endl;
    }

    file.close();
}


/*****************************************************************************************************//**
*
* @brief Prints the contents of Pos and magnetization to the file "./trajectories/particle_N"
* where N is the current particle number.
*
* @param Pos vector of vector3 containing the positions at all time steps for the current particle.
* @param magnetization vector of vector3 containing the magnetization at all time steps for the current particle.
*
*********************************************************************************************************/
void Deposition::print_trajectory(const vector<vector3>& Pos, const vector<vector3>& magnetization){
    
    std::ofstream file;
    string path = "./trajectories/particle_" + std::to_string(frozen_particles.size());
    file.open(path);

    if (!file.is_open()){
        throw std::invalid_argument( "Error: Could not open " + path + ". Does the dir exist? " );
    }

    for(size_t i = 0; i < Pos.size(); ++i){    

        file << Pos[i][0] << "\t";
        file << Pos[i][1] << "\t";
        file << Pos[i][2] << "\t";
        file << magnetization[i][0] << "\t";
        file << magnetization[i][1] << "\t";
        file << magnetization[i][2];          
        file << endl;
    }

    file.close();
}

/*****************************************************************************************************//**
*
* @brief Counts how many particles collided with the substrate and how many collided with another
* particle. Prints the result.
*
* @param os Ostream where output (only info) will be written. This is normally std::cout.
*
*********************************************************************************************************/
void Deposition::Finalize(ostream& os){

    int p_count = 0;
    int s_count = 0;

    for(auto p : frozen_particles){

        if(p.Get_collision_object()=='p'){
            p_count++;
        }else if(p.Get_collision_object()=='s'){
            s_count++;
        }
    }

    os << "Total number of particles: " << frozen_particles.size() << endl;
    os << "Particles collided with substrate: " << s_count << endl;
    os << "Particles collided with other particle: " << p_count << endl;
}