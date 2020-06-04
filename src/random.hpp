#ifndef RANDOM
#define RANDOM

#include <ctime>   //time(0)
#include <random>
#include <array>

#include "types.hpp"


class Random{
public:
    
    Random();
    void Setup_lognorm( const double mean, const double std, const int option);

    //! Single and double charge option for log-norm distributions.
    enum CHARGE { SINGLE = 1, DOUBLE = 2};
    

    double uniform(){ return uniform_engine(generator);}
    //!< Returns a single random number uniformely sampled in the range [0,1].

    double log_norm(int option);

    void Fill_normal( vector3& r);
    void Fill_uniform( vector3& r );
    void Fill_normal( std::array<double, 9>& r );

    vector3 Generate_in_point_sphere();

private:

    std::subtract_with_carry_engine<std::uint_fast64_t, 48, 5, 12> generator;
    //!< Generator for random numbers.

    std::uniform_real_distribution<double> uniform_engine;
    //!< Engine for the uniform distribution.

    std::normal_distribution<double> normal_engine;     
    //!< Engine for the normal distribution.

    std::lognormal_distribution<double> log_norm_1;
    //!< Engine for the log-norm distribution for singly charged particles.

    std::lognormal_distribution<double> log_norm_2;
    //!< Engine for the log-norm distribution for doubly charged particles.
};


#endif