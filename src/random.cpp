/*****************************************************************************************************//**
* @file 
*
* @brief This file contains the Random class.
*
*********************************************************************************************************/

/*****************************************************************************************************//**
*
* @class Random
*
* @brief This class is responsible for generating random numbers from the different distributions used in 
* the simulation.
*
* It contains functionality to generate numbers from a normal, uniform as well as log-norm distributions.
*
*********************************************************************************************************/

#include "random.hpp"

#include <stdexcept>
#include <algorithm>

/*****************************************************************************************************//**
*
* @brief Create a Random object.
*
* Sets up the uniform and normal distributions. 
*********************************************************************************************************/
Random::Random(){

    uniform_engine = std::uniform_real_distribution<double> (0.0, 1.0);
    normal_engine = std::normal_distribution<double> (0.0, 1.0); 

    generator = std::subtract_with_carry_engine<std::uint_fast64_t, 48, 5, 12> (time(0));
    srand(time(0));
}

/*****************************************************************************************************//**
*
* @brief Sets the mean value and standard deviation of one of the two log-norm distributions used in the 
* simulation.
* 
* @param mean Logarithm of the mean value.
* @param std Standard deviation.
* @param option Determines which distribution to change. Alternatives are Random::CHARGE::SINGLE and
*  Random::CHARGE::DOUBLE. 
* 
* @param 
*********************************************************************************************************/
void Random::Setup_lognorm(const double mean, const double std, const int option){

    double mu = log( pow(mean,2) / sqrt(pow(mean,2) + pow(std,2)) );
    double s = log ( 1 + pow(std, 2) / pow(mean, 2) );

    if(option == CHARGE::SINGLE){

        log_norm_1 = std::lognormal_distribution<double> (mu, s);

    }else if(option == CHARGE::DOUBLE){

        log_norm_2 = std::lognormal_distribution<double> (mu, s);

    } else {

        throw std::invalid_argument( "Inavlid option in Setup_lognorm");
    }
}


/*****************************************************************************************************//**
*
* @brief Fills the input vector3 with random numbers from a normal distribution with mean 0 and std 1.
*
*********************************************************************************************************/
void Random::Fill_normal(vector3& r){
    for( auto& i : r) i = normal_engine(generator);
}


/*****************************************************************************************************//**
*
* @brief Fills the input array with random numbers from a normal distribution with mean 0 and std 1.
*
*********************************************************************************************************/
void Random::Fill_normal(std::array<double, 6>& r){
    for( auto& i : r) i = normal_engine(generator);
}


/*****************************************************************************************************//**
*
* @brief Fills the input array with random numbers from a uniform distribution in the range [0,1].
*
*********************************************************************************************************/
void Random::Fill_uniform(vector3& r){
    for( auto& i : r) i = uniform_engine(generator);
}


/*****************************************************************************************************//**
*
* @brief Generates a random number from a log-norm distribution.
*
* The log-norm distributions are used for generating the diameters of particles. Both single and doubly
* charged ones. The means and standard deviations are set up using Setup_lognorm().
* 
* @param option Determines which log-norm distribution to use. Alternatives are Random::CHARGE::SINGLE 
* and Random::CHARGE::DOUBLE.
*
*********************************************************************************************************/
double Random::log_norm(int option){
    
    if (option == Random::CHARGE::SINGLE){

        return log_norm_1(generator);

    }else if(option == Random::CHARGE::DOUBLE){

        return log_norm_2(generator);

    }else{

        throw std::invalid_argument( "Inavlid option in log_norm.");   
    }
}


/*****************************************************************************************************//**
*
* @brief Returns a random 3D point sampled uniformly from inside a sphere with radius 1. 
*
* All coordinates are in the range [-1,1].
*
* @return vector3 containing a random point.
*
*********************************************************************************************************/
vector3 Random::Generate_in_point_sphere(){    

    vector3 r;

    while(true){ 

        Fill_uniform(r);
        
        double n = 0.0;
        
        for(int i=0; i<3; ++i){

            r[i] = (r[i] - 0.5)*2.0;
            n += pow(r[i], 2);
        }
        if (pow(n,0.5)<1.0) return r;    
    }
}