#ifndef INPUTREADER_H
#define INPUTREADER_H

#include "data_struct.cpp"

#include <vector>
#include <array>
#include <string>

class InputReader{
public:
    InputReader() {data.SetHamakerConstants();}
    InputReader(std::string infile, std::ostream& os);

    bool Success() const {return success;} /*!< @brief Returns True if no error were encountered when reading the files.*/
    InputData Get_data() const {return data;} /*!< @brief Returns an instance of InputData with values as defined in the input file, or default values if not defined in the input.*/
    std::vector<std::array<double, 7>> Read_particles(std::string filename, std::ostream& os);

private:    
    void Read_key(const std::string& key, const std::string& value);
    void Add_error(const std::string& key, const std::string& value);
    bool To_bool(const std::string& str);
    double check_pos(double value, std::string key);

    bool success{true};
    int keys_read{0};
    std::string errors;
    InputData data;
};
#endif