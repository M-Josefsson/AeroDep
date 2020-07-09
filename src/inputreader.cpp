/*****************************************************************************************************//**
* @file 
*
* @brief This file contains the class InputReader.
*
*********************************************************************************************************/

/*****************************************************************************************************//**
*
* @class InputReader 
*
* @brief This class is responsible for reading infiles, both containing input/environment parameters
* as well as frozen (already deposited) particles from a previous run.
*
*********************************************************************************************************/

#include "inputreader.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

using std::endl;
using std::ifstream;
using std::string;
using std::vector;
using std::array;



/*****************************************************************************************************//**
*
* @brief Creates the InputReader object and reads input parameters from the file 'infile'. 
*
* First, the default values for all variables are loaded. Then those values are overwritten, with the
* corresponding values in the input file. Ignores lines beginning with '#'.
*
* @param infile Filename of the infile containing deposition parameters. Syntax must be 'key = value' 
* (spaces and tabs are ignored).
* @param os Ostream where output (only info) will be written. This is normally std::cout.
*
*********************************************************************************************************/
InputReader::InputReader(const std::string infile, std::ostream& os){

    errors = "";    
    ifstream fin( infile, ifstream::in );
    fin.precision( 15 );
    
    if ( ! fin.is_open() ) {

      errors += "Unable to open file \"" + infile + "\"\n";
      success = false;
      os << errors << endl;
      return;

    }

    //Read keys from input file
    string line;
    while (getline(fin, line)){

        line.erase(std::remove_if (line.begin(), line.end(), ::isblank), line.end());        
        std::istringstream is_line(line);        
        string key;

        if (getline(is_line, key, '=')){

            string value;

            if (key[0] == '#') continue;

            if (getline(is_line, value)){

                Read_key(key, value);
                keys_read++;
            }
        }
    }

    fin.close();     
    os << "Infile read. In total " << keys_read << " keys were found." << endl;
    os << errors << endl;

    data.SetHamakerConstants();
}


/*****************************************************************************************************//**
*
* @brief Reads the value for 'key' and converts 'value' to the appropriate data type for key.
*
*********************************************************************************************************/
void InputReader::Read_key(const string& key, const string& value){

    if (key=="diameter"){
        try{ data.diameter = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }

    }else if(key=="particle_number"){
        try{ data.nbr_of_particles = check_pos(stoi(value), key); }catch(...){ Add_error(key, value); }

    }else if(key=="temperature"){
        try{ data.T = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }

    }else if(key=="start_height"){
        try{ data.z_start = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }

    }else if(key=="box_size"){
        try{ data.control_l = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }

    }else if(key=="charge"){
        try{ data.q = stod(value); }catch(...){ Add_error(key, value); }

    }else if(key=="density"){
        try{ data.density = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }

    }else if(key=="dt"){
        try{ data.dt = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }

    }else if(key=="E"){
        try{ data.E0 = stod(value); }catch(...){ Add_error(key, value); }

    }else if(key=="n_substrate"){
        try{ data.n_substrate = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }

    }else if(key=="n_gas"){
        try{ data.n_gas = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }

    }else if(key=="dynamic_viscocity"){
        try{ data.eta_g = stod(value); }catch(...){ Add_error(key, value); }

    }else if(key=="mean_free_path"){
        try{ data.mfp = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }

    }else if(key=="particle_susceptibility"){
        try{ data.Xi = stod(value); }catch(...){ Add_error(key, value); }

    }else if(key=="interaction_length"){
        try{ data.interaction_length = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }

    }else if(key=="v_x"){
        try{ data.v_g[0] = stod(value);}catch(...){ Add_error(key, value); }

    }else if(key=="v_y"){
        try{ data.v_g[1] = stod(value);}catch(...){ Add_error(key, value); }

    }else if(key=="v_z"){
        try{ data.v_g[2] = stod(value);}catch(...){ Add_error(key, value); }

    }else if(key=="dielectric_substrate"){
        try{ data.eps[1] = check_pos(stod(value), key); data.eps[2] = stod(value); }catch(...){ Add_error(key, value); }

    }else if(key=="dielectric_gas"){
        try{ data.eps[3] = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }

    }else if(key=="Bx"){
        try{ data.B[0] = stod(value); }catch(...){ Add_error(key, value); }

    }else if(key=="By"){
        try{ data.B[1] = stod(value); }catch(...){ Add_error(key, value); }

    }else if(key=="Bz"){
        try{ data.B[2] = stod(value); }catch(...){ Add_error(key, value); }

    }else if(key=="m_saturation"){
        try{ data.m_saturation = stod(value); }catch(...){ Add_error(key, value); }

    }else if(key=="alignment_field_strength"){
        try{ data.alignment_field_strength = stod(value); }catch(...){ Add_error(key, value); }

    }else if(key=="diameter_std"){
        try{ data.diameter_std = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }  

    }else if(key=="diameter_std2"){
        try{ data.diameter_std2 = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }   

    }else if(key=="double_charge_fraction"){
        try{ data.double_charge_fraction = check_pos(stod(value), key); }catch(...){ Add_error(key, value); }                

    }else if(key=="print_trajectory"){
        try{ data.print_trajectory = To_bool(value); }catch(...){ Add_error(key, value); }

    }else if(key=="remove_surface_charge"){
        try{ data.remove_surface_charge = To_bool(value); }catch(...){ Add_error(key, value); }

    }else if(key=="verbose"){
        try{ data.verbose = To_bool(value); }catch(...){ Add_error(key, value); }

    }else if(key=="magnetic"){
        try{ data.calcMagnetic = To_bool(value); }catch(...){ Add_error(key, value); }

    }else if(key=="magnetic_type"){        
        if (value == "ferro"){            
            data.magnetic_ferro = true;
        }else if (value == "para"){
            data.magnetic_ferro = false;
        }else{
            Add_error(key, value);
        }   

    }else{
      errors += "Unknown key or value " + key + "," + value + "\n";
      keys_read--;
    }    
}


/*****************************************************************************************************//**
* @brief Adds errors to the error string.
*
* @param key The key in the input file that caused the error.
* @param value The value in the input file that caused the error.
*
*********************************************************************************************************/
void InputReader::Add_error(const string& key, const string& value){

    errors += "Invalid value " + value + " for key " + key + "\n"; 
    success = false;
    keys_read--;
}


/*****************************************************************************************************//**
* @brief Checks whether the input value is positive. If not an error message is added to errors.
*
* @param value The value to be checked.
* @param key The key that is associated with the value. Only used in the error message.
*
*********************************************************************************************************/
double InputReader::check_pos(double value, string key){

    if(value < 0.0){
        errors += "Negative value (" + std::to_string(value) + ")  not allowed for key " + key;
        success = false;
    }
    return value;
}


/*****************************************************************************************************//**
*
* @brief Reads already deposited particles from a file.
*
* The input file must have seven columns where each row corresponds to one particle. Elements 0-2 represent 
* a particle's position, 3-5 its magnetization, and element 6 its diameter.
*
* @param filename The name of the file to be read.
* @param os Ostream where output (only info) will be written. This is normally std::cout.
*
* @return A vector of arrays (length 7) where elements 0-2 represent a particle's position, 
* 3-5 its magnetization, and element 6 its diameter. 
*
*********************************************************************************************************/
vector<array<double, 7>> InputReader::Read_particles(string filename, std::ostream& os){

    vector<array<double, 7>> input_particles;        
    double x, y, z, mx, my, mz, d;
    errors = "";    

    ifstream fin( filename, ifstream::in );
    fin.precision( 15 );
    
    if ( ! fin.is_open() ) {

      errors += "Unable to open file \"" + filename + "\"\n";
      success = false;
      os << errors << endl;
      return input_particles;

    }

    int count = 0;

    while ( (fin >> x) && (fin >> y) && (fin >> z) &&
            (fin >> mx) && (fin >> my) && (fin >> mz) &&
            (fin >> d) ){

        if (z < -d){
            throw std::invalid_argument( "Negative z-position not allowed for input particles." );
        }else if(d < 0){
            throw std::invalid_argument( "Negative diameter not allowed for input particles." );
        }

        array<double, 7> temp;
        temp[0] = x;
        temp[1] = y;        
        temp[2] = z;
        temp[3] = mx;
        temp[4] = my;
        temp[5] = mz;
        temp[6] = d;

        input_particles.push_back(temp);
        count++;
    }

    fin.close(); 

    os << count << " particles read from file " << filename << endl;
    return input_particles;
}


/*****************************************************************************************************//**
*
* @brief Converts the input string to a Boolean value.
*
*********************************************************************************************************/
bool InputReader::To_bool(const string& str) {

    string str_lower(str);
    std::transform (str_lower.begin(), str_lower.end(), str_lower.begin(),
                                [](unsigned char c) -> unsigned char { return std::tolower(c); });

    std::istringstream is(str_lower);
    bool b;

    if (str_lower != "false" && str_lower != "true"){        
        throw std::invalid_argument( str + " cannot be converted to a boolean value." );
    }

    is >> std::boolalpha >> b;
    return b;
}