#include "../src/vector_operations.hpp"
#include "../src/particle.hpp"
#include "../src/constants.hpp"
#include "../src/fields_and_forces.hpp"
#include "../src/data_struct.cpp"
#include "../src/inputreader.hpp"
#include "../src/deposition.hpp"
#include "../src/random.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <chrono>
#include <assert.h> 

using std::vector;
using std::string;
using std::cout;
using std::endl;

using vec = vector<double>;
using arr = array<double, 3>;

static double REL_TOL = 1e-12; 

template <typename T>
void print_arr(const T& A){
    for (auto a: A){
        cout << a << " ";
    }
    cout << endl;
}

bool test_close_double(double a, double b){
    if (a==0 && b==0){
        return true;
    }
    return fabs((a-b)/a) < REL_TOL;
}

template <typename T>
bool test_close_arr(const T& a, const T& b){    
    assert(a.size()==b.size());
    bool ret = true;
    for(size_t i = 0; i<a.size(); ++i){
        if(!test_close_double(a[i], b[i])){ret = false;}
    }
    if (ret){
        cout << "passed." << endl;
    }else{
        cout << "failed." << endl << " Excpected:  ";
        print_arr(a);
        cout << endl << " Got:    ";
        print_arr(b);
        cout << endl;
    }
    return ret;
}

int test_AH(const InputData& data){
    bool a = test_close_double( 1.5702527177791352e-19, data.AH132 );
    cout << "AH132 ";
    if (a){
        cout << "passed" << endl;
        return 0;
    }else{
        cout << "failed." << endl;
        return 1;
    }
}

int test_F_vdW_particle_substrate(const Particle& particle, const double& AH){
    cout << "F_vdW_PS ";
    arr t = F_vdW_particle_substrate(particle, AH);
    return !test_close_arr({0,0,-3.6975881790444115e-15}, t);
}

int test_F_image_particle_substrate(const Particle& particle, const array<double, 4>& eps){
    cout << "F_Image_PS ";
    arr t = F_image_particle_substrate(particle, eps);
    return !test_close_arr({0,0,-3.42676448878678e-15}, t);
}

int test_E_field(const Particle& p_inc, const Particle& p2, const InputData& data){
    int c = 0;
    cout << "E_field (1/2) ";
    arr p_to_p = p_inc.pos - p2.pos;
    double dist = norm(p_to_p);
    arr t = E_field_deposited_particle( p2, p_to_p, dist, data.E0, data.eps);
    c += !test_close_arr({0.0,0.0, -135179.60822369}, t);

    cout << "E_field (2/2) ";    
    Particle p3(p2);
    p3.q = 0.0;
    p_to_p = p_inc.pos - p3.pos;
    dist = norm(p_to_p);
    t = E_field_deposited_particle( p3, p_to_p, dist, data.E0, data.eps);
    c += !test_close_arr({0.0,0.0, 9375.000000000002}, t);
    return c;
}

int test_F_image_particle_particle(const Particle& p_inc, const Particle& p2, const InputData& data){
    int c = 0;
    cout << "F_Image_PP (1/2) ";
    arr p_to_p = p_inc.pos - p2.pos;
    double dist = norm(p_to_p);
    arr t = F_image_particle_particle(p_inc, p2, p_to_p, dist, data.eps);
    c += !test_close_arr({0.0,0.0, -9.595106728538234e-16}, t);

    cout << "F_Image_PP (2/2) ";    
    Particle p3(p2);
    p3.q = 0.0;
    p_to_p = p_inc.pos - p3.pos;
    dist = norm(p_to_p);
    t = F_image_particle_particle(p_inc, p3, p_to_p, dist, data.eps);
    c += !test_close_arr({0.0,0.0, -7.977402694968764e-16}, t);
    return c;
}

int test_F_vdW_particle_particle(const Particle& p_inc, const Particle& p2, const InputData& data){
    int c = 0;
    cout << "F_vdW_PP (1/2) ";
    cout << "AH132 " << data.AH132 << endl;
    cout << "AH232 " << data.AH232 << endl;
    arr p_to_p = p_inc.pos - p2.pos;
    double dist = norm(p_to_p);
    arr t = F_vdW_particle_particle(p_inc, p2, p_to_p, dist, data.AH232);
    c += !test_close_arr({0.0, 0.0, -2.43753761069778364e-13}, t);

    cout << "F_vdW_PP (2/2) ";
    arr null{0.0,0.0,0.0}; 
    arr m {-1.0, 1.0, 0.0};
    Particle p3(0.0, 30e-9, data, null, m);
    p3.pos[0] = 20e-9;
    p3.pos[1] = 80e-9;
    p3.pos[2] = 0.0;    
    p_to_p = p_inc.pos - p3.pos;
    dist = norm(p_to_p);
    t = F_vdW_particle_particle(p_inc, p3, p_to_p, dist, data.AH232);
    c += !test_close_arr({8.99834836034805674e-15, 3.5993393441392227e-14, -4.49917418017402664e-14}, t);
    return c;
}

int test_H_field_dipole(const Particle& p_inc, const Particle& p2, const InputData& data){
    int c = 0;
    cout << "H_field_dipole (1/2) ";
    arr p_to_p = p_inc.pos - p2.pos;
    double dist = norm(p_to_p);
    arr t = H_field_dipole(p2, p_to_p, dist);
    c += !test_close_arr({0.0, 0.0, 17843.75}, t);

    cout << "H_field_dipole (2/2) ";
    array<double, 3> m{-1.0, 1.0, 0.0};
    arr null{0.0,0.0,0.0}; 
    Particle p3(0.0, 30e-9, data, null, m);
    p3.pos[0] = 20e-9;
    p3.pos[1] = 80e-9;
    p_to_p = p_inc.pos - p3.pos;
    dist = norm(p_to_p);
    t = H_field_dipole(p3, p_to_p, dist);
    c += !test_close_arr({759.8919620809551, -89.3990543624652, -670.4929077184898}, t);
    return c;
}

int test_F_ferro(const Particle& p_inc, const Particle& p2, const InputData& data){
    int c = 0;
    cout << "F_ferro (1/4) ";
    arr p_to_p = p_inc.pos - p2.pos;
    double dist = norm(p_to_p);
    arr t = F_ferromagnetic_dipole_dipole(p_inc, p2, p_to_p, dist);
    c += !test_close_arr({0.0, 0.0, -1.6290596923211355e-11}, t);   


    arr null{0.0,0.0,0.0}; 
    arr m{0.0,0.0,1.0};
    Particle p3(100.0e-9, 30e-9, data, null, m);
    p3.pos[1] = 80e-9;
    p3.pos[0] = 80e-9;
    cout << "F_ferro (2/4) ";
    p_to_p = p3.pos - p2.pos;
    dist = norm(p_to_p);
    t = F_ferromagnetic_dipole_dipole(p3, p2, p_to_p, dist);
    c += !test_close_arr({-9.903619105253994e-13, -9.903619105253994e-13, 8.37438380223683e-13}, t);   

    cout << "F_ferro (3/4) ";
    arr new_mag = {-1.0, 1.0, 0.0};
    p3.magnetization = new_mag * data.m_saturation;
    t = F_ferromagnetic_dipole_dipole(p3, p2, p_to_p, dist);
    c += !test_close_arr({-1.037695384190216e-12, 1.037695384190216e-12, 0.0}, t);    

    cout << "F_ferro (4/4) ";
    p3.pos[0] = 0.0;
    p3.pos[1] = 20.0e-9;
    double s3 = sqrt(1.0/3.0);
    new_mag = {s3,s3,s3};
    p3.magnetization = new_mag * data.m_saturation;
    p_to_p = p3.pos - p2.pos;
    dist = norm(p_to_p);
    t = F_ferromagnetic_dipole_dipole(p3, p2, p_to_p, dist);
    c += !test_close_arr({4.263468188022102e-12, 1.9677545483178844e-13, -1.0953833652302943e-11}, t);   
    return c;
}

int test_E_field_gradient(const Particle& p_inc, const Particle& p2, const InputData& data){
    int c = 0;
    cout << "E_field_gradient (1/3) ";
    arr p_to_p = p_inc.pos - p2.pos;
    double dist = norm(p_to_p);
    array<double, 9> t = Gradient_E_field_deposited_particle( p2, p_to_p , dist, data.E0, data.eps);
    array<double, 9> ref = {-1304921082236.9004, 0.0, 0.0, 0.0, -1304921082236.9004, 0.0,
                0.0, 0.0, 2609842164473.801};
    c += !test_close_arr(ref, t);


    cout << "E_field_gradient (2/3) ";
    Particle p3(p2);
    p3.q = 0.0;
    p_to_p = p_inc.pos - p3.pos;
    dist = norm(p_to_p);
    t = Gradient_E_field_deposited_particle( p3, p_to_p , dist, data.E0, data.eps);
    ref = {140624999999.99997, 0.0, 0.0, 0.0,140624999999.99997, 0.0,
                0.0, 0.0, -281249999999.99994};
    c += !test_close_arr(ref, t);

    cout << "E_field_gradient (3/3) ";
    p3.pos = {20e-9, 80e-9, 0.0};
    p3.q = -1.0;
    p3.diameter = 30e-9;
    p_to_p = p_inc.pos - p3.pos;
    dist = norm(p_to_p);
    t = Gradient_E_field_deposited_particle( p3, p_to_p , dist, data.E0, data.eps);
    ref = {-609114704473.9462, 185716698033.7088, -233806503301.45377, 185716698033.7088,
        87322913152.46185, -935226013205.8151,-227243057919.38818, -908972231677.5527,
         521791791321.4843};
    c += !test_close_arr(ref, t);
    return c;
}

int test_F_dipole(const Particle& p_inc, const Particle& p2, const InputData& data){
    int c = 0;
    cout << "F_dipole (1/2) ";
    arr p_to_p = p_inc.pos - p2.pos;
    double dist = norm(p_to_p);
    array<double,9> dE = Gradient_E_field_deposited_particle(p2, p_to_p, dist, data.E0, data.eps);
    arr E = E_field_deposited_particle( p2, p_to_p, dist, data.E0, data.eps);
    E[2] += data.E0;
    arr t = F_dipole( p_inc, E, dE, data.eps);
    c += !test_close_arr({0.0, 0.0, 1.60907916600849404e-16}, t);

    cout << "F_dipole (2/2) ";
    arr null{0.0, 0.0, 0.0}; 
    arr m{-1.0, 1.0, 0.0}; 
    Particle p3(0.0, 30e-9, data, null, m);
    p3.pos[0] = 20e-9;
    p3.pos[1] = 80e-9;
    p_to_p = p_inc.pos - p3.pos;
    dist = norm(p_to_p);
    dE = Gradient_E_field_deposited_particle(p3, p_to_p, dist, data.E0, data.eps);
    E = E_field_deposited_particle( p3, p_to_p, dist, data.E0, data.eps);
    E[2] += data.E0;
    t = F_dipole( p_inc, E, dE, data.eps);
    c += !test_close_arr({-1.92334788243923037e-17, -7.69339152975692147e-17, 2.61763435749594285e-17}, t);          
    return c;
}


int test_total_force(Particle& p_inc, const Particle& p2, InputData& data){
    int c = 0;
    cout << "Total_force (1/4): ";
    vector<Particle> fp;
    arr t = Get_total_force(p_inc, fp, data);
    c += !test_close_arr({0.0, 0.0, -5.51896512918311993e-14}, t);

    cout << "Total_force (2/4): ";
    fp.push_back(p2);
    data.calcMagnetic = false;
    t = Get_total_force(p_inc, fp, data);
    c += !test_close_arr({0.0, 0.0, -2.78083854327372518e-13}, t);

    Particle p3(0.0, 30e-9, data, {0, 0, 0}, {-1.0, 1.0, 0.0});
    p3.pos[0] = 20e-9;
    p3.pos[1] = 80e-9;
    fp.push_back(p3);
    
    cout << "Total_force (3/4): ";
    data.calcMagnetic = false;
    t = Get_total_force(p_inc, fp, data);
    c += !test_close_arr({6.89706952827778644e-15, 2.75882781131111458e-14, -3.1263007761534491e-13}, t);

    cout << "Total_force (4/4): ";
    data.calcMagnetic = true;
    t = Get_total_force(p_inc, fp, data);
    c += !test_close_arr({-4.54630997167884765e-13, -1.18157427159361326e-13, -1.61999972162395245e-11}, t);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (size_t i =0; i< 10000; ++i){    
        t = t +  Get_total_force(p_inc, fp, data);
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " [µs]" << std::endl;
    return c;
}   

bool test_input(double value, double expected_value, string key){
    if (!test_close_double(value, expected_value)){
        cout << key << " failed" << endl;
        cout << "Expected: " << expected_value << "    Got: " << value << endl;
        return true;
    }else{
        return false;
    }
}

int test_inputreader_infile(){
    InputReader reader("../tests/infile_test", std::cout);
    InputData data = reader.Get_data();

    int c = 0;

    c += test_input(data.diameter, 45e-9, "diameter");
    c += test_input(data.nbr_of_particles, 37, "nbr_of_particles");
    c += test_input(data.control_l, 1.4e-6, "box_size");
    c += test_input(data.dt, 5.5e-10, "dt");
    c += test_input(data.B[0], 0.1, "Bx");
    c += test_input(data.B[1], 0.2, "By");
    c += test_input(data.B[2], 0.3, "Bz");
    c += test_input(data.z_start, 5.4e-7, "start_height");
    c += test_input(data.interaction_length, 0.8e-6, "interaction_length");
    c += test_input(data.alignment_field_strength, 1e-5, "alignment_field_strength");
    c += test_input(data.q, 0.0, "charge");
    c += test_input(data.T, -230, "temperature");
    c += test_input(data.E0, -200e3, "E0");
    c += test_input(data.n_substrate, 1.1, "n_substrate");
    c += test_input(data.n_gas, 1.05, "n_gas");
    c += test_input(data.mfp, 60e-9, "mean_free_path");
    c += test_input(data.Xi, -5e-5, "particle_susceptibility");
    c += test_input(data.density, 1000, "density");
    c += test_input(data.eps[1], 4, "dielectric_substrate");
    c += test_input(data.eps[3], 0.75, "dielectric_substrate");
    c += test_input(data.m_saturation, 1500, "m_saturation");
    c += test_input(data.diameter_std, 5e-9, "diameter_std");
    c += test_input(data.diameter_std2, 10e-9, "diameter_std2");
    c += test_input(data.double_charge_fraction, 0.9, "double_charge_fraction");

    if(data.print_trajectory){ cout << "print_trajectory is true, expected false" << endl; c+=1; }
    if(data.remove_surface_charge){ cout << "remove_surface_charge is true, expected false"<< endl; }
    if(!data.calcMagnetic){ cout << "magnetic is false, expected true"<< endl; }
    if(!data.magnetic_ferro){ cout << "magnetic_ferro is false, expected true"<< endl; }

    cout << c << " errors found. " << endl;    

    cout << endl << "reading extrernal particles from file." << endl;
    vector<array<double, 7>> in_particle;
    in_particle = reader.Read_particles("../tests/particle_positions_test", std::cout);
    try{
        in_particle = reader.Read_particles("../tests/particle_positions_test_neg_z", std::cout);
        cout << "Error: Read_particles did not catch negative z-value!" << endl;
    }catch(const std::invalid_argument& e ){
        cout << "Gaurd against negative z-values passed" << endl;
    }
    try{
        in_particle = reader.Read_particles("../tests/particle_positions_test_neg_d", std::cout);
        cout << "Error: Read_particles did not catch negative diamater!" << endl;
    }catch(const std::invalid_argument& e ){
        cout << "Gaurd against negative diamater passed" << endl;
    }

    return c;

}

int test_inputreader_default(){
    InputReader reader;
    InputData data = reader.Get_data();

    int c = 0;

    c += test_input(data.diameter, 30e-9, "diameter");
    c += test_input(data.nbr_of_particles, 100, "nbr_of_particles");
    c += test_input(data.control_l, 1.5e-6, "box_size");
    c += test_input(data.dt, 10e-10, "dt");
    c += test_input(data.B[0], 0.0, "Bx");
    c += test_input(data.B[1], 0.0, "By");
    c += test_input(data.B[2], 0.0, "Bz");
    c += test_input(data.z_start, 150e-9, "start_height");
    c += test_input(data.interaction_length, 0.5e-6, "interaction_length");
    c += test_input(data.alignment_field_strength, 1e-5, "alignment_field_strength");
    c += test_input(data.q, -1.0, "charge");
    c += test_input(data.T, 300, "temperature");
    c += test_input(data.E0, 300e3, "E0");
    c += test_input(data.n_substrate, 1.4585, "n_substrate");
    c += test_input(data.n_gas, 1.0, "n_gas");
    c += test_input(data.mfp, 66.5e-9, "mean_free_path");
    c += test_input(data.Xi, -2.2e-5, "particle_susceptibility");
    c += test_input(data.density, 7874, "density");
    c += test_input(data.eps[1], 3.9, "dielectric_substrate");
    c += test_input(data.eps[3], 1.0, "dielectric_gas");
    c += test_input(data.m_saturation, 1.707e6, "m_saturation");
    c += test_input(data.diameter_std, 0.0, "diameter_std");
    c += test_input(data.diameter_std2, 0.0, "diameter_std2");
    c += test_input(data.double_charge_fraction, 0.0, "double_charge_fraction");

    if(data.print_trajectory){ cout << "print_trajectory is true, expected false" << endl; c += 1;}
    if(!data.remove_surface_charge){ cout << "remove_surface_charge is false, expected true"<< endl; c += 1;}
    if(!data.calcMagnetic){ cout << "magnetic is false, expected true"<< endl; c += 1;}
    if(!data.magnetic_ferro){ cout << "magnetic_ferro is false, expected true"<< endl;c += 1;}

    cout << c << " errors found. " << endl;    
    return c;
}

void test_random(){
    Random R;

    //Calculate mean values of the distributions (normal, uniform, and 3D uniform in sphere)  
    double avg_n{0.0}, avg_u{0.0}, x{0.0}, y{0.0}, z{0.0};
    double count_n{0.0}, count_u{0.0}, count_s{0.0};
    array<double, 6> rn;
    arr ru, rs;

    for (size_t i = 0; i<10000; i++){
        R.Fill_normal(rn);
        R.Fill_uniform(ru);
        rs = R.Generate_in_point_sphere();
        for (auto j:rn){
            avg_n += j;
            count_n += 1.0;
        }    
        for(auto j:ru){
            avg_u += j;
            count_u += 1.0;
        }
        x+= rs[0]; 
        y+= rs[1];
        z+= rs[2];
        count_s += 1.0;
    }
    cout << "Testing mean values of the distributions used for random numbers."<< endl;
    cout << "Mean of normal distribution (should be close to 0): " << avg_n/count_n << endl;
    cout << "Mean of normal uniform (should be close to 0.5): " << avg_u/count_u << endl;
    cout << "Means of positions for points generated in sphere (should be zeros): " << endl;
    cout << x/count_s << endl << y/count_s << endl << z/count_s << endl;

}

int test_deposition(Deposition dep){

    int c = 0;
    cout << "\nTesting mobility (1/1) " ;
    if( !test_close_double(dep.Mobility(30e-9, 66.5e-9, 2.0), 542608844.2321931)){
        cout << "Failed.\n Expected " << 542608844.2321931 << " got " << dep.Mobility(30e-9, 66.5e-9, 2.0) << endl;
        c += 1;
    }else{
        cout <<"passed." << endl;
    }

    double d_double_charge = dep.Double_charge_diameter(30e-9, 66.5e-9);
    cout << "Testing Double_charge_diameter (1/1) ";
    if (fabs(d_double_charge - 43.3524) < 1e-3){
        cout << "passed." << endl;
    }else{
        cout << "Failed.\n Expected " << 43.3524 << " got " << d_double_charge << endl;
        c += 1;
    }

    cout << "\nAdding particles from file (new z_start should be 2.4694e-7" << endl;
    InputData data;
    InputReader reader;
    vector<array<double, 7>> in_particle;
    in_particle = reader.Read_particles("../tests/particle_positions_test", std::cout);
    dep.Add_input_particles(in_particle, data, cout);
    cout << "Checking particle position (1/2): ";
    c += !test_close_arr({2.99411e-07, -4.86038e-07, 1.5e-08}, dep.frozen_particles[0].pos);
    cout << "Checking particle position (2/2): ";
    c += !test_close_arr({-3.78074e-07, -7.75986e-09, 6.35748e-08}, dep.frozen_particles[4].pos);
    return c;
}


int main(){
    //setup
    cout.precision(18);
    InputData data;
    data.SetHamakerConstants();
    data.alignment_field_strength = 1e10;
    data.density = 7310.0;
    data.m_saturation = 1.713e6;

    arr r{0.0,0.0,0.0};
    arr m{0.0, 0.0, 1.0};
    Particle p_inc(0.0, 30e-9, data, r, m);
    p_inc.pos[0] = 0.0;
    p_inc.pos[1] = 0.0;
    p_inc.pos[2] = 100e-9;

    Particle p2(0.0, 50e-9, data, r, m);
    p2.pos[0] = 0.0;
    p2.pos[1] = 0.0;
    p2.pos[2] = 0.0;

    //Test forces
    cout << "\n-------------------------------------------------------\n";
    cout << "Testing force calculations." << endl;
    cout << "-------------------------------------------------------\n";

    int nbr_errors = 0;

    nbr_errors += test_AH(data);
    nbr_errors += test_F_vdW_particle_substrate(p_inc, data.AH132);
    nbr_errors += test_F_image_particle_substrate(p_inc, data.eps);
    nbr_errors += test_E_field(p_inc, p2, data);
    nbr_errors += test_F_image_particle_particle(p_inc, p2, data);
    nbr_errors += test_F_vdW_particle_particle(p_inc, p2, data);
    nbr_errors += test_H_field_dipole(p_inc, p2, data);
    nbr_errors += test_F_ferro(p_inc, p2, data);
    nbr_errors += test_E_field_gradient(p_inc, p2, data);
    nbr_errors += test_F_dipole(p_inc, p2, data);
    nbr_errors += test_total_force( p_inc, p2, data);
    //Test paramagnetic

    //Test time step

    //Test inputs
    cout << "\n-------------------------------------------------------\n";
    cout << "Testing inputreader." << endl;
    cout << "-------------------------------------------------------\n";
    cout << "Reading default values." << endl;
    nbr_errors += test_inputreader_default();

    cout <<  "\nReading values from infile (expecting one warning about negative temperature)." << endl;
    nbr_errors += test_inputreader_infile();
    
    //Test deposition
    cout << "\n-------------------------------------------------------\n";
    cout << "Testing functionality of the Deposition class." << endl;
    cout << "-------------------------------------------------------\n";
    Deposition dep(data);
    nbr_errors += test_deposition(dep);

    //Test random number generation
    cout << "\n-------------------------------------------------------\n";
    cout << "Testing Random number generation." << endl;
    cout << "-------------------------------------------------------\n";
    test_random();

    cout << "-------------------------------------------------------\n";
    cout << "\nTesting complete. In total " << nbr_errors <<" errors were found.\n\n" << endl;
    assert(nbr_errors==0);

}