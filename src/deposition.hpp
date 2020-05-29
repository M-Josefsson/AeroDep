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
    void print_trajectory(std::vector<vector3> Pos, std::vector<vector3> magnetization);

    std::vector<Particle> frozen_particles;
    bool rand_size, rand_size2, double_charge;
    double z_start, z_start_original, AH132, AH131, AH232, diameter, d_p2;

    //std::default_random_engine generator;
    std::subtract_with_carry_engine<std::uint_fast64_t, 48, 5, 12> generator;
    std::uniform_real_distribution<double> uniform;
    std::normal_distribution<double> normal;     
    std::lognormal_distribution<double> normal_diameter, normal_diameter_double;
};
#endif