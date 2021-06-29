# Incoherent volume element generator (IVEGen)

The Incoherent volume element generator (IVEGen) supplies the [Radiative Transfer with Reciprocal Transactions](https://bitbucket.org/planetarysystemresearch/r2t2_pub) (R2T2) or [SIRIS4](https://bitbucket.org/planetarysystemresearch/siris4) with incoherent T-matrices and Mueller matrices. IVEGen uses the superposition T-matrix method (STMM) to compute the incoherent scattering properties of spherical volume elements made of spherical scatterers (normal or power law distributed). 

If you notice a bug or weird results, please create an [issue](https://bitbucket.org/planetarysystemresearch/stmm_t_matrix/issues?status=new&status=open) or contact the authors.


## Requirements
* CMake (version >=3.0.2)
* Fortran compiler (supports F2003)
* MPI libraries
* hdf5
* LAPACK and BLAS

## Compile
To compile the program type:

```console
mkdir build
cd build
cmake ..
make
```
After make, an executable IVEGen is built into folder /bin/.

The default setting is to build an optimized version of the program with normal distributed radiuses. If you wish to compile a debug build, or use a power law distribution, the user must type


**Power law**
```console
#DISCRETION REQUIRED: It is important to notice that the power law distribution does not work correctly in the current version.
cmake .. -DUSE_DISTRIBUTION:STRING=powerlaw
```

**Debug Build**
```console
cmake .. -DCMAKE_BUILD_TYPE=DEBUG 
```
instead. 

**NOTE:** 
NOTE: If you want to switch between builds (e.g., RELEASE to DEBUG, or powerlaw to normal), you must clear the content of the _/build/_ folder and start with the cmake command before compiling.


## Run

The IVEgen uses syntax
```console
./IVEgen -parameter1 value1 -parameter2 value2
```
The program stops running if invalid input is given. 

Keep the unit of the arguments consistent. (volrad, stdev, wavenumber, mean, upB, etc.)

Valid parameters are

* **n_samples**: How many incoherent T-matrices will be written to **T_out2**. Default: 128
* **k**: wavenumber, 2*pi/wavelength. efault: 1
* **vol_rad**: the radius of the volume element. efault: 2*pi
* **N_ave**: How many samples are generated for an ensemble average. efault: 128
* **read_lmk**: Read wavenumbers and permittivies from a **read_file**. Default: 0 (do not read)
* **read_file**: A file which contains wavenumbers and permittivies. Default: 'data.txt'
* **N_theta**: How many theta angles are generated for a S-matrix presentation: Default: 180
* **N_phi**: How many phi angles are generated for a S-matrix presentation: Default: 64
* **density**: Density of the volume element: Default: 0.15d0
* **eps_r**: Real permittivity of the material: Default: 1.7161
* **eps_i**: Imaginary permittivity of the material: Default: 0.01
* **T_out**: The output file for the coherent T_matrix: Default: T.h5
* **T_out2**: The output file for incoherent T-matrices: Default: 'T_multi.h5'
* **S_out_prefix**: The prefix of the S-matrix output file: Default: "mueller"
* **kappa_out_prefix**: The prefix of the general output file: Default: "kappa"

### Normal distribution specific


* **mean**: The mean of a normal distribution: Default: 2.0
* **stdev**: The standard deviation of a normal distribution: Default: 0.0 (constant radius)

### Power law distribution specific
* **lowB**: The lower bound of the distribution : Default: 1.0
* **upB**: The upper bound of the distribution : Default: 2.0
* **n**: power law index (must be positive): Default: 2.0


## Contribute

Why even bother? Well, develop your Fortran, testing, or physics skills, while helping us to do reproducible and reliable science. As you can see, the code is in lousy shape regarding the software development (e.g., testing is done by comparing the current results with the old results. See those long functions). [List of contributors](contributors).

So please, fork the project and play with the code. Important changes can be merged into the master branch after they have been reviewed.

### To do list
* Automatic testing
* Just trying to build the code (robust cmake confs) and run it
* Refactoring
* Better documentation
* Get rid of permittivity. It is just an unnecessary step
* Fix the power law distribution. (One branch may contain experimental version) 

## Projects included

[CMake Fortran Template](https://github.com/SethMMorton/cmake_fortran_template).


## License

See [LICENSE](LICENSE.txt)

## References 

1. [Timo Väisänen, Johannes Markkanen, Antti Penttilä, Karri Muinonen, Radiative transfer with reciprocal transactions:  Numerical method and its implementation, PLOS ONE, submitted]()

2. [Karri Muinonen, Johannes Markkanen, Timo Väisänen, Jouni Peltoniemi, and Antti Penttilä, "Multiple scattering of light in discrete random media using incoherent interactions," Opt. Lett. 43, 683-686 (2018) ](https://www.osapublishing.org/ol/abstract.cfm?uri=ol-43-4-683)

1. [Markkanen , J & Yuffa , A J 2017 , ' Fast superposition T-matrix solution for clusters with arbitrarily-shaped constituent particles ', JQSRT , vol 189 , pp. 181-188](https://helda.helsinki.fi/handle/10138/180562)



