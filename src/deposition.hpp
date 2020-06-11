#ifndef DEPOSITION_H
#define DEPOSITION_H

#include "types.hpp"
#include "particle.hpp"
#include "inputreader.hpp"
#include "random.hpp"

class Deposition{
public:
    Deposition(const InputData& data);
    bool Add_particle(ostream& os, const InputData& data);
    void Add_input_particles(const std::vector<std::array<double, 7>>& input_particles, const InputData& data, ostream& os);
    void Print_final_positions(const std::string& filename, std::ostream& os);
    void Finalize(std::ostream& os);

    friend int test_deposition(Deposition deposition);

private:
    void Update_z_start(ostream& os);        
    double Mobility(const double& d1, const double& mfp, const double& Z);
    double Double_charge_diameter(double diameter, double mfp);
    double Get_diameter(double& current_q, const double& double_charge_fraction);
    void print_trajectory(const std::vector<vector3>& Pos, const std::vector<vector3>& magnetization);

    Random distributions;

    std::vector<Particle> frozen_particles;
    //!< @brief Vector containing the frozen (already deposited) particles.

    double z_start;
    //!< @brief Current start height for new particles.

    double z_start_original;
    //!< @brief Original start height for new particles.    

    double diameter; 
    //!< @brief Diameter (or mean diameter if rand_size=True) of single charged particles.

    double d_p2;
    //!< @brief Diameter (or mean diameter if rand_size2=True) of double charged particles.

    bool rand_size;
    //!< @brief Generate particles according to a Log-norm distribution?

    bool rand_size2;
    //!< @brief Generate double charged particles according to a Log-norm distribution?

    bool double_charge;
    //!< @brief Include double charged particles?
    
};
#endif