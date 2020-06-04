#ifndef DEPOSITION_H
#define DEPOSITION_H

#include "types.hpp"
#include "particle.hpp"
#include "inputreader.hpp"

#include <random>

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
    void Generate_random_normal(std::array<double, 9>& r);
    void Generate_random_uniform(vector3& r);
    void Generate_random_point_sphere(vector3& r);
    double Mobility(const double& d1, const double& mfp, const double& Z);
    double Double_charge_diameter(double diameter, double mfp);
    double Get_diameter(double& current_q, const double& double_charge_fraction);
    void print_trajectory(const std::vector<vector3>& Pos, const std::vector<vector3>& magnetization);

    std::vector<Particle> frozen_particles;
    //!< @brief Vector containing the frozen (already deposited) particles.

    double z_start;
    //!< @brief Current start height for new particles.

    double z_start_original;
    //!< @brief Original start height for new particles.    

    double diameter; 
    //!< @brief Diameter (or mean diameter if rand_size=True) of single charegd particles.

    double d_p2;
    //!< @brief Diameter (or mean diameter if rand_size2=True) of double charegd particles.

    bool rand_size;
    //!< @brief Generate particles according to a Log-norm distribution?

    bool rand_size2;
    //!< @brief Generate double charged particles according to a Log-norm distribution?

    bool double_charge;
    //!< @brief Include double charged particles?
    

    std::subtract_with_carry_engine<std::uint_fast64_t, 48, 5, 12> generator;
    std::uniform_real_distribution<double> uniform;
    std::normal_distribution<double> normal;     
    std::lognormal_distribution<double> normal_diameter, normal_diameter_double;
};
#endif