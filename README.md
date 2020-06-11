# AeroDep
Molecular dynamics-type simulations of the deposition process of aerosol nanoparticles.

### Readme table of contents:

1. [Overview](#overview)

2. [Installation and requirements](#installation-and-requirements)

3. [Usage](#usage)

4. [Infiles](#infiles)

5. [Output](#output)

6. [Documentation](#documentation)


### Overview
This program performs molecular dynamics-type calculations to simulate the final steps of the deposition process of (magnetic) aerosol nanoparticles onto an substrate. The particle concentration in the gas is assumed to be very low (<$10^6$ particles per $cm^{-3}$)such that particles in the gas does not interact, and only one particle is in the aerosol phase in the simulation volume at any given time. When the particle collides with the substrate, or another particle, its properties (such as position and magnetization)
are frozen, and a new particle is spawned.

The calculations are based on Euler's method for solving Newton's force equation by taking many small successive time steps. The forces included are of electrostatic, magnetic, and van der Waals nature. Interactions between the frozen particles and the incoming particle are taken into account using pair-wise interactions. In addition, stochastic motion governed by Brownian motion is included as it has a significant effect on the particles' 
trajectories for small nanoparticles.

The core of the program builds upon the work in:
Krinke et al. "[Microscopic aspects of the deposition of nanoparticles from the gas phase](https://doi.org/10.1016/S0021-8502(02)00074-5)" Journal of Aerosol Science 33.10 (2002), but with added functionality (such as magnetism).

The integration scheme is described in Zarutskaya and Shapiro ["Capture of nanoparticles by magnetic filters"](https://doi.org/10.1016/S0021-8502(99)00567-4)  Journal of Aerosol Science 31.8 (2000).


### Installation and requirements 
The program is entirely run and installed using the command line. The only requirement for compiling the program is to use the gcc compiler. This requirement can be relaxed if the flag `-ffast-math` is removed from `src/makefile`.

For compilation `cd` to AeroDep/src and type 
```bash
make
```
which will put the compiled executable in `AeroDep/bin/.`

### Usage
AeroDep is run from `AeroDep/bin` by typing 
```bash
./AeroDep [infile] [particle_file]
``` 
where the commands in `[]` are optional. If no optional arguments are provided the program will execute using its default parameters (see [Infiles](#infiles)). The first optional argument is the file name of the input file containing the simulation parameters (see [Infiles](#infiles)). The second parameter is the file name for a file containing already deposited particles (i.e. an output file from a previous run). If the second option is provided, the particles specified in the file are added to the simulation volume as frozen particles before the simulation is executed.

### Infiles
There are two kinds of input files to AeroDep: parameter input files and particle input files. The parameter input file contains the simulation specific parameters. The syntax in that file have to be `key = value`. Lines starting `#` are ignored, and all  spaces and tabs are always ignored. An example of the contents of a minimal infile is
```
particle_number = 10
```
which sets the total number of particles to 10. All other parameters take their default values. For more examples see the provided example infile `AeroDep/bin/example_infile`.

All available keys are (default values within parenthesis):

- ***particle_number*** - total number of particles in the deposition (100)
- ***diameter*** - particle diameter in m (30e-9)
- ***start_height*** - distance between a newly generated particle and the highest frozen particle or substrate in m (150e-9)
- ***temperature*** - chamber/gas temperature in K (300)
- ***box_size*** - side-length of the simulation box in m (1e-6)
- ***charge*** - charge of the new particles in elementary charges (-1.0)
- ***density*** - the particle's density in kg/m3 (7874.0 Fe)
- ***dt*** - time step in s (1e-9)
- ***E*** - electric field in the chamber in V/m (300e3)
- ***n_substrate*** - refractive index of the substrate (1.4585 Si)
- ***n_gas*** - refractive index of the gas (1.0) 
- ***dynamic_viscosity*** - dynamic viscosity of the gas in kg/m/s (18.13e-6)
- ***mean_free_path*** - mean free path of the particles in m (66.5e-9)
- ***particle_susceptibility*** - the particles' magnetic susceptibility (-2.2e-5)
- ***interaction_length*** - max distance in m between two particles for interaction between them to be included (0.5e-6)
- ***v_x*** - gas velocity in m/s x-component (0.0)
- ***v_y*** - gas velocity in m/s y-component (0.0)
- ***v_z*** - gas velocity in m/s z-component (0.0)
- ***dielectric_substrate*** - relative dielectric constant of the substrate (3.9)
- ***dielectric_gas*** - relative dielectric constant of the gas (1.0)
- ***Bx*** - external magnetic field (B-field) x-component (0.0)
- ***By*** - external magnetic field (B-field) y-component (0.0)
- ***Bz*** - external magnetic field (B-field) z-component (0.0)
- ***m_saturation*** - saturation magnetization of the particles in A/m2 (1.707e6)
- ***alignment_field_strength*** - local H-field strength above which a particle's magnetization is aligned with the field (1e-5)
- ***diameter_std*** - standard deviation in m for the diameter of the particles in a log-norm distribution (0.0)
- ***diameter_std2*** - standard deviation in m for the diameter of doubly charged particles in a log-norm distribution (0.0)
- ***double_charge_fraction*** - fraction of doubly charged particles (0.0)
- ***print_trajectory*** - print trajectories for all particles? (false)
- ***remove_surface_charge*** - remove particle charge upon collision? (true)
- ***magnetic*** - include magnetic interactions? (true)
- ***magnetic_type*** - type of magnetic interaction to use in force calculations. Options are 'ferro' and 'para'. Default is 'ferro'.

The second kind of infile, the particle input file must have the same format as the particle output file described below.

### Output

The main result from running the AeroDep program is a file containing the positions of all the deposited particles in the simulation. This output is always generated and is written to the file `bin/particle_positions`. The file is updated for every 50th generated particle, as well as at the end of the program. 

Output file structure: each row corresponds to one particle. There are seven columns with data where column 1-3 is the particle's position, 4-6 its magnetization, and column 7 its diameter. 

The second kind of output files the program can generate are files containing the trajectories of an individual particles. If the `print_trajectory` option is set to `true` one trajectory file is generated for each deposited particle. These files are written to `bin/trajectories` (this folder must exist!). Each row in the file corresponds to the current position (x,y and z coordinate).

### Documentation

The documentation for the AeroDep project uses doxygen to generate nice documentation in html format. Thus, doxygen needs to be installed in order to generate the documentation. For unix-based systems (Linux ans macOS) the html output is generated using 
```bash
make docs
```
in `AeroDep/src`. For windows-based systems it is instead recommended to open and run `AeroDep/docs/doxygen_config` using doxywizard.

The generated documentation is accessed by opening `AeroDep/docs/html/index.html` in your browser of choice.
