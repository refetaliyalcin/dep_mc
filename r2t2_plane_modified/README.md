# Radiative Transfer with Reciprocal Transactions
The Radiative Transfer with Reciprocal Transactions (R2T2) computes the light-scattering characteristics of dense discrete random media. All you need are incoherent T-matrices, corresponding albedos and incoherent mean free path. T-matrices for volume elements made of spherical scatterers can be computed with the [IVEGen](https://bitbucket.org/planetarysystemresearch/ivegen_pub)

R2T2 updates required packages from Bitbucket's servers, so you must add [public ssh key](https://confluence.atlassian.com/bitbucket/set-up-an-ssh-key-728138079.html#) to your Bitbucket account. 

Another option is to download the R2T2 from a [PLOSONE branch](https://bitbucket.org/planetarysystemresearch/r2t2_pub/src/PLOSONE/), which does not require any registration. 

## Requirements
* CMake (version >=3.0.2)
* Fortran compiler (supports F2003)
* MPI libraries
* hdf5
* LAPACK and BLAS / (or mkl)

## Compile
The program can be compiled with:

```console
mkdir build
cd build
cmake ..
make
```
which will create an executable **R2T2** into folder _/bin/_.

The R2T2 supports plane (not tested enough) and sphere geometry. The geometry is selected by adjusting the cmake flag GEOMETRY 


```console
cmake ..  -DGEOMETRY=SPHERE
```

```console
cmake ..  -DGEOMETRY=PLANE
```

The sphere geometry is used as a default geometry if the flag is not provided.

## Run
The R2T2 is run from the command line with
```console
R2T2 input.file TIME
```
in which _input.file_ contains the arguments. TIME defines how long the R2T2 should run. The algorithm will stop itself even if it has not finished tracing all the rays. If the argument is not given, the time is read from the input.file. 

## Input
Input file must have format
```console
PARAMETER = ARGUMENT
```
There are multiple parameters which you can tweak and all of them are listed in "definitions.f90" or "src/var_list_n.csv". Some of the frequently needed are:

* **wavelength**: define the wavelength in the same units as the **mean_free_path**
* **mean_free_path**: incoherent mean free path
* **cutoff_intensity**: stop tracing the ray after the intensity goes below this limit
* **tau_threshold**: distance exponential attenuation is killed
* **seed**: seed for the random number generator. 0 means random seed
* **medium_thickness_or_radius**: Radius or thickness of the medium in the same unit as the **mean_free_path**
* **volume_element_radius**: the radius of the volume element the medium in the same unit as the **mean_free_path**
* **default_run_time**: running time (same as the TIME given to the command line)
* **number of rays**: the number of traced rays
* **max_sca_procs**: how many consecutive scattering processes are allowed before new ray is generated
* **rt_the_angles**: how many theta angles are used to discretize the solid angle
* **rt_phi_angles**: how many phi angles are used to discretize the solid angle
* **pol_state_X_Y**: the initial Stokes configurations. X can be Q,U,V and Y MINUS or PLUS.
* **tmatrix_location**: path to the file containing T-matrices and albedos
* **output_details**: path to file to which details will be written
* **output_rt**: path to file to which rt output will be saved
* **output_cb**: path to file to which cb output will be saved

Notice that incoherent T-matrices are required, e.g., from the IVEGEN, in order to be able to run the R2T2.  

## Output
The output files defined by **output_rt** and **output_cb** need to be processed further with a python script "sph_output.py" (sphere geometry) located in _/scripts/_. This requires the polarization state Q-. Circular polarization can be added by computing the polarization state V- too.

The output from the plane geometry is processed with "pln_output.py". The plane geometry R2T2 should be run with all the six polarization states (Q+-,U+-,V+-).


## Contribute
Why even bother? Well, develop your Fortran, testing, and physics skills, while helping us to do reproducible and reliable science. As you see, the code is in lousy shape concerning the software development (e.g., testing is done by comparing the current results with the old results. See those long functions). 

[List of contributors](contributors).

So please, fork the project or create a branch where you can play with the code. Important changes will be merged to the master branch after they have been reviewed.


### To do list
* Combine RT and CB output
* improve python script
* improve cmake files
* Automatic testing
* Just trying to build the code (robust cmake confs) and run it
* Refactoring
* Better documentation
* Maybe something from [IVEGen](https://bitbucket.org/planetarysystemresearch/ivegen_pub)

## Projects included

[The CMake Fortran template](https://github.com/SethMMorton/cmake_fortran_template).

[Mersenne Twister and jumper](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/JUMP/dsfmt-jump.html)

[The interface for the Mersenne Twister ](https://github.com/jsspencer/dSFMT_F03_interface).
 



## License

See [LICENSE](LICENSE.txt)

## References 

1. [Timo Väisänen, Johannes Markkanen, Antti Penttilä, Karri Muinonen, Radiative transfer with reciprocal transactions:  Numerical method and its implementation, PLOS ONE, in review]()

2. [Karri Muinonen, Johannes Markkanen, Timo Väisänen, Jouni Peltoniemi, and Antti Penttilä, "Multiple scattering of light in discrete random media using incoherent interactions," Opt. Lett. 43, 683-686 (2018) ](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-43-4-683)

3. [J. Markkanen, T. Väisänen, A. Penttilä, and K. Muinonen, "Scattering and absorption in dense discrete random media of irregular particles," Opt. Lett. 43, 2925-2928 (2018)](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-43-12-2925)


