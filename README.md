# AeroDep
Molecular dynamics-type simulations of the deposition process of aerosol nanoparticles.


### Overview
This program performs molecular dynamics-type calculations in order to simulate the final steps of the depositon process 
of (magnetic) nanoparticles in an aerosol phase onto an substrate. The particle concentration in the 
gas is assumed to be very low such that only one particle is in the aerosol phase in the simulation volume at any given
time. When the particle collides with the substrate, or another particle, its properties (such as position and magnetization)
are frozen, and a new particle is spawned.

The calculations are based on Euler's method for solving Newton's force equation by taking many small successive 
time steps. The forces included are of electrostratic, magnetic, and van der Waals nature. Interactions between the frozen 
particles and the incoming particle is taken into account using pair-wise interactions. In addition,
stochastic motion governed by Brownian motion is included as it has a significat effect on the particles' 
trajectories for small nanoparticles.

The core of the program builds upon the work in:
Krinke et al. "Microscopic aspects of the deposition of nanoparticles from the gas phase." Journal of Aerosol Science 33.10 (2002)
but with added functionality (such as magnetism).

### Installation and requirements 
The program is entirely run and install used the command line. The only requirement for compiling the parogram is to use the gcc compiler. This requirement can be relaxed if the flag `-ffast-math` is removed from `src/Makefile`, and as long as the other compiler supports OpenMP. 

For compilation `cd` to AeroDep/src and type `make`, which till put the compiled executable in `AeroDep/bin.`

### Usage
AeroDep is run from `AeroDep/bin` by typing `./AeroDep [infile] [particle_file]` where the commands in `[]` are optional. If no optional arguments are provided the program will execute using its default parameters (see below). The first optional argument is the filename of the input files containing the simulation parameters (see below). The second parameter is the filename for a file containing already deposited particles (i.e. an output file from the program). If the second option is provided, the particles specified in the files are added to the simulation volume as frozen particles before the simulation is run.

### Infiles
There are two kinds of input files to AeroDep, parameter input files and particle input files. The parameter input file contains the simumlation specific parameters. The syntax in the file have to be is `key = value`. Lines starting `#` are ignored, and all  spaces and tabs are always ignored. An example of the contents of a minimal infile is
```
nbr_of_particles = 10
```
which sets the total number of particles to 10. All other parameters take their default values. 

All available key are (default values within parenthesis):

- **particle_number** - total number of particles in the deposition (100)
- **diameter** - particle diameter in m (30e-9)
- **start_height** - distance between a newly generated particle and the highest frozen particle or substrate in m (150e-9)
- **temperature** - Chamber/gas temperature in K(300)
- **box_size** - side-length of the simulation box in m (1e-6)
- **charge** - charge of the new particles in elementary charges (-1)
- **density** - the particle's density in kg/m3 (7310)
- **dt** - time spep in s (1e-9)
- **E** - electric field in the chamber in V/m (300e3)
- **n_substrate** - refractive index of the substrate (1.4585)
- **n_gas** - refractive index of the gas (1.0) 
- **dynamic_viscocity** - dynamic viscocity of the gas in kg/m/s (18.13e-6)
- **mean_free_path** - mean free path of the particles iin m (66.5e-9)
- **particle_susceptibility** - the particles' magnetic susceptibility (-2.2e-5)
- **interaction_length** - max distance in m between two particles for interaction between them to be included (0.5e-6)
- **v_x** - gas velocity in m/s x-component (0.0)
- **v_y** - gas velocity in m/s y-component (0.0)
- **v_z** - gas velocity in m/s z-component (0.0)
- **dielectric_substrate** - relative dielectric constant of the substrate (3.9)
- **dielectric_gas** - relative dielectric constant of the gas (1.0)
- **Bx** - external magnetic field (B-field) x-component (0.0)
- **By** - external magnetic field (B-field) y-component (0.0)
- **Bz** - external magnetic field (B-field) z-component (0.0)
- **m_saturation** - saturation magnetization of the particles in A/m2 (1.713e6)
- **alignment_field_strength** - local H-field strength above which a particle's magnetization is aligned with the field (1e-5)
- **diameter_std** - standard deviation in m for the diameter of the particles in a log-norm distribution (0.0)
- **diameter_std2** - standard deviation in m for the diameter of doubly charged particles in a log-norm distribution (0.0)
- **double_charge_fraction** - fraction of particles with double charge (0.0)
- **print_trajectory** - print trajectories for all particles? (false)
- **remove_surface_charge** - remove particle charge upon collision? (true)
- **magnetic** - include magnetic interactions? (false)
- **magnetic_type** - type of magnetic interaction to use in force calculations. Options are 'ferro' and 'para' (ferro)

The particle input file takes has the same format as the particle output file discussed below.

### Output

The main result from running the AeroDep program is a file containing the positions of all the particles deposited in the simulation. This output is always generated and is written to the file `bin/particle_positions`. The file is updated for every 50 particles generated, as well as at the end of the program. 

Output file structure: each row corresponds to one particle. There are seven columns with data where column 1-3 is the particle's position, 4-6 its magnetization, and column 7 its diameter. 

The second kind of output file the program can generate is a file containing the trajectory of a single particle. If the `print_trajectory` option is set to true one trajectory file is generated for each deposited particle. These files are written to `bin/trajectories` (this folder must exist!). Each row in the file corresponds to the current position (x,y and z coordinate).
