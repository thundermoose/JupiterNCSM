# JupiterNCSM
This is JupiterNCSM, a no-core shell model code developed in C at Chalmers
University of Technology. 

## License

This software is released under the terms and conditions of GPL version 2.
As required by the license a full copy of it can be found in the LICENSE file.
This software is presented as is and absolutely no warranty is included.

## Dependencies

JupiterNCSM depends on the following libraries 
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [BLAS](http://www.netlib.org/blas/)
- [wigxjpf](http://fy.chalmers.se/subatom/wigxjpf/)
- [libconfig](https://hyperrealm.github.io/libconfig/)
- [OpenMP](https://www.openmp.org/)
- [GSL](https://www.gnu.org/software/gsl/)

At the moment the index-lists are computed with [Anicre](https://github.com/thundermoose/Anicre)
which is therefor an indirect dependence.

## Installation

To install JupiterNCSM first run
```
>> git clone git@git.chalmers.se:djarv/jupiterncsm.git
>> cd jupiterncsm
>> git submodule init
>> git submodule update
```

Before running make for the first time, you may need to direct it to where 
the dependencies are set, and also change compiler and linker. This can be done
in a makefile.local in which the following variables can be set:

``` makefile.local
compiler= # if other than gcc is desired
linker= # if other than gcc is desired
test_data_path=/path/to/jupiterncsm_test_data # (default is ./jupiterncsm_test_data)
wigxjpf_path=/path/to/wigxjpf
hdf5_comp_flags= -I/path/to/hdf5/include # if hdf5 is not installed under /usr/
hdf5_link_flags= -L/path/to/hdf5/lib -lhdf5 # if hdf5 is not installed under /usr/ 
blas_comp_flags= -I/path/to/blas/include
blas_link_flags= -L/path/to/blas/lib -lblas -llapack
libconfig_comp_flags= -I/path/to/libconfig
libconfig_link_flags=-lconfig
```

However, their default values works in many Linux distributions.

To run JupiterNCSM:s internal tests, run 
```
>> make test
>> ./test_scripts/run_all_tests.sh
```

To compile the release version
```
>> make
```

## Example

To illustrate how to use JupiterNCSM an example of how to compute the 
eigen spectrum of ⁴He for the N2LOsat interaction at ħΩ = 20 MeV up to Nmax = 6, 
is included in the subdirectory **example**. To try it out, download 
[jupiter_ncsm_example_4he_data.tar.gz](http://fy.chalmers.se/subatom/jupiterncsm_example_data/jupiter_ncsm_example_4he_data.tar.gz) 
and place it in **example** after which you can run the  script
```
>> cd example
>> ./example_4_helium.sh
Usage ./example_4_helium.sh <nmax> [<k-particle-forces>] [--max-loaded-memory <memstr>]
>> ./example_4_helium.sh 6 3 --max-loaded-memory 4GB
```
The tar package will be automatically unpacked if it is not already unpacked,
and then JupiterNCSM will be executed. See the script for details on how 
JupiterNCSM should be run.
The first argument of the script is the NCSM truncation parameter Nmax and is
required. The second argument, which is optional, can either be 2 or 3, if
it is 2 JupiterNCSM will only use the 2NF part of the interaction, and if it
is 3 both the 2NF and 3NF part will be used, the default value 2. The flag 
`--max-loaded-memory` limit how much RAM both Bacchus and Mercury may use,
default is 16GB.

When the script is ran, as in the example above, a new subdirectory is created 
named `example_run_2nf_and_3nf_nmax6`. This is the work directory of JupiterNCSM
where examples of the Bacchus configuration file, the evaluation order file can
be found as well as index lists, Krylow vectors and M-scheme matrix elements. 
The result of the calculation can be found in `bacchus_results` and in 
`eigen_vectors`. 

To view the resulting eigen spectrum you can run
```
>> tail -n 7 example_run_2nf_and_3nf_nmax6/bacchus_results
Diagonalization end after 935722 µs
eigensystem:
{
	num_eigenvalues = 11;
	eigenvalues = ( -25.90722535033397, 17.521187442926934, 23.599734543083045, 35.906134844316433, 53.329043345505333, 80.410359916207241, 104.87999820268072, 119.21925859073369, 129.67968655743007, 142.2911766330169, 167.74531258016918 );
}
Bacchus ends after 1.02472e+06 µs
```

## Authors

See [authors](AUTHORS.md).
