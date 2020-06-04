#include <vector>
#include <array>
#include <iostream>

#include "deposition.hpp"
#include "inputreader.hpp" 


/*****************************************************************************************************//**
*
* @file 
*
* @brief This file contains the driver program. 
*
* The functionality of the main function is, in short: \n
* 1) Use an instance of InputReader to read input files. \n
* 2) Instantiate a Deposition object to set up the simulation. \n
* 3) Continue to spawn particles and evolve time until the user defined stopping condition 
* (nbr_of_particles) in reached. \n
* \n
* All particle position are written to the output files for every 50th generated particle. 
*
*********************************************************************************************************/

/*****************************************************************************************************//**
*
* @mainpage 
* 
* This program performs molecular dynamics-type calculations in order to simulate the final steps of the
* depositon process of (magnetic) nanoparticles in an aerosol phase onto an substrate. The particle 
* concentration in the gas is assumed to be very low such that only one particle is in the aerosol phase 
* in the simulation volume at any given time. When the particle collides with the substrate, or another 
* particle, its properties (such as position and magnetization) are frozen, and a new particle is spawned.
*
* The calculations are based on Euler's method for solving Newton's force equation by taking many small 
* successive time steps. The forces included are of electrostratic, magnetic, and van der Waals nature. 
* Interactions between the frozen particles and the incoming particle is taken into account using 
* pair-wise interactions. In addition, stochastic motion governed by Brownian motion is included as it has
* a significat effect on the particles' trajectories for small nano particles.
* 
* AeroDep is mainly meant to be used as a standalone program, but its different source code parts can
* also be used on their own in other projects. This documentations aims describing the idividual parts of 
* the code, thus lowering the bar for understanding the machinery behind the program.
*
*********************************************************************************************************/

int main(int argc, char *argv[]){

    //Start of setup
    InputReader reader;
    const std::string outfile = "./particle_positions";
    ostream& os = std::cout;

    //Read input file and default values
    if(argc>1){

        reader = InputReader(argv[1], os);

        if (!reader.Success()) {
            os << "Failed to read input file: " << argv[1] << std::endl;
            return 1;
        }
        
    }else{
        reader = InputReader();    
    }
    
    InputData data = reader.Get_data();
    Deposition deposition(data);

    //Read optional file with already deposited particles    
    if (argc>2){

        std::vector<std::array<double, 7>> input_particles = reader.Read_particles(argv[2], os);

        if (!reader.Success()) {
            os << "Failed to read Particle infile." << std::endl;
            return 1;
        }
        deposition.Add_input_particles(input_particles, data, os);
    }      
    
    //End of setup and start of the actual simulation
    bool cont;
    for (int i=0; i<data.nbr_of_particles; ++i){

        os << std::endl << i << std::endl;
        cont = deposition.Add_particle(os, data);    

        if (!cont) continue;

        if (i%50 == 0 ) {
            deposition.Print_final_positions(outfile, os);
        }            
    }

    deposition.Print_final_positions(outfile, os);
    deposition.Finalize(os);    
}