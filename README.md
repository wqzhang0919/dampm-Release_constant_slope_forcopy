# DAMPM - Strucutre

An implementation of the Material Point Method (MPM) in the finite-strain elastoplastic framework. The code DAMPM is based om the code Larsie (developped by Lars Blatny). The main branches of the codes are 

+ Release_constant_slope
+ Release_curvilinear
+ Release_projected_gradient
+ Flow_complex_terrain
+ Flow_projected_gradient

The orthers branches are either obsoletes (meaned to be deleted) or in developpement. The branches used for the simulation in the paper "A Depth-Averaged Material Point Method for Shallow Landslide : Application to Snow Slab Avalanches Release" are Release_constant_slope and Release_curvileanar


# Release_contant_slope

## How to compile and Run

Make sure you have set the parameters of you choice in the file `mpm.cpp`. Open a terminal session in the `dampm` folder

+ To create the makefiles :

`mkdir build`

`cd build`

`cmake -DCMAKE_BUILD_TYPE=Release ..`

+ To compile

`make`

+ To run 

`./src/mpm`


### Dependencies

Dampm only relies only on CMake and the linear algebra library Eigen.
To use OpenMP, use the CMake option `-DUSE_OMP=True`, and make sure to select the appropriate parallelized functions you want to use in `simulation.cpp`. Please make sure you have cmake install as well as a recent enough compilator for C++.


## Main Classes and variables

Ths code is composed ot three classes, the class particles, the class grid and the class simulations. The first two classes are used to store the quantities of respectively the particles and the grid. The list attribute, their type and their shape can be found in the file `data_arrays.hpp`.
The attribute of the class particles are the machanical and numerical parameters used in the simulation. the methods (described in next section allow to run the simulation). the list of the method and attributes are in the file `simulation.hpp`. the vaklue of the attribute are set in the file `mpm.cpp`. 

## Structure of the code
This is in a nutshell the main methods of the class Simulation. The equation numbers refer to the preprint

```
- initialize() 
- initialize_particles() : Initialiaze particles position
- simulate()
    - advanceStep()
            - updateDt() : according to Eqs.(42) and (43)
            - remesh_on_complex_terrain() : recreate grid nodes
            - P2G() : according to Eq.(30)
            - compute_forces() : according to Eqs.(27) and (28)
            - update_grid_velocity() : according to Eq.(36)
            - G2P() : according to Eqs.(34) and (35)
            - elasticity() : according to Eq.(41)
            - plasticity() 
            - update_height() : according to Eq.(38)
	    - update_cohesion() : allow us to manually remove cohesion at a given time
        - saveParticleData()
```


## Output data

In the `mpm.cpp` file, the directory where the simulation data will be stored. *Please make sure that the directory exists before running the simulation*. The type of particle saving can also be change in the file, depemading of the type of post-process. 

## Post Process files

Some python `.py` and houdini `.hipnc` files for the post-process or the vizualisation of the simulated data are in the folder `post_process`. Please make sure that the format of the saved data and the path indicated in the post-process file are the right ones.

## Performance

Dampm is continously implemented according to the principle *premature optimization is the root of all evil*. Basic optimizations (precomputations, clever particle-grid loops, memory reads, etc...) are still being explored. Although still under development, Dampm also supports **parallelization** on shared memory with **OpenMP**.

## Contact or informations

Even if the code is not always well comented, its structure and the name of the variables allow a good readability. However, if you have any questions concerning the code or its application, feel free to contact the current of previous developper of the code at `louis.guillet@inria.fr` or `johan.gaume@epfl.ch`. 
