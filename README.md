# JupiterNCSM
This is JupiterNCSM, a no-core shell model code developed in C at Chalmers
University of Technology. It is licensed under GPL v2

## Installation

## Dependencies

## Example

To illustrate how to use JupiterNCSM an example of how to compute the 
eigen spectrum of ⁴He for the N2LOsat interaction at ħΩ = 20 MeV up to Nmax = 6, 
is included in the subdirectory **example**. To try it out, download 
[jupiter_ncsm_example_4he_data.tar.gz](https://unknownlink.nonexistent.se) 
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
