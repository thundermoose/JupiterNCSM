#!/usr/bin/bash

# This is script uses JupiterNCSM to compute the eigenspectrum of 4He.
# To use it you need also the jupiter_ncsm_example_4he_data.tar.gz
# found at

# Usage ./script.sh <nmax> [<k-particle-forces>] [--max-loaded-memory <memstr>]

if [ $# -lt 1 ]
then
	echo "Usage $0 <nmax> [<k-particle-forces>] [--max-loaded-memory <memstr>]"
fi

# Parsing script arguments

# Default values
nmax=0
k_particle_forces=2
Z=2
N=2
max_loaded_memory=16GB

variable_to_set="nmax"
for arg in $@
do
	if [ $arg = "--max-loaded-memory" ]
	then
		variable_to_set="max_loaded_memory"
	else
		eval "$variable_to_set=$arg"
		variable_to_set="k_particle_forces"
	fi
done

# External resources


# Extracting example data if not already extracted
if [ ! -d example_4he_data ]
then
	if [ ! -f jupiter_ncsm_example_4he_data.tar.gz ]
	then
		# In the future I might add wget to download the tar.gz package
		# but for now I will rely on the user to download it
		echo "Can't find jupiter_ncsm_example_4he_data.tar.gz in current directory $PWD"
		exit 1
	fi
	tar -xf jupiter_ncsm_example_4he_data.tar.gz
fi

# Setting up data location variables
example_data_path=$PWD/example_4he_data
index_lists_path=$example_data_path/anicr_runs/4He_nmax$nmax.jtd
comb_file_path=$index_lists_path/comb.txt
two_nucleon_force_path=$example_data_path/two_nucleon_forces/TBME_a4_bare_n2lo_sat_Nmax$((nmax))_0.hw20.0.bin
three_nucleon_force_path=$example_data_path/three_nucleon_forces/3NF_N2LO_450_N6_hw20_split.h5

# JupiterNCSM program
aurora=$PWD/../release/Aurora/aurora.x
mars=$PWD/../release/Mars/mars.x
mercury=$PWD/../release/Mercury/mercury_J_scheme_to_internal.x
bacchus=$PWD/../release/Bacchus/bacchus.x

# Setting up run environment
run_directory="./example_run_4_helium"
mkdir -p $run_directory
cd $run_directory

# Generating M-scheme elements, generating the evaluation-order formating the index-lists 

mkdir -p m_scheme_interaction
mkdir -p index_lists
if [ $k_particle_forces -eq 3 ]
then
	# 2NF and 3NF case
	# Creating the M-scheme forces
	$mercury $comb_file_path $index_lists_path m_scheme_interaction --num-protons $Z --num-neutrons $N --single-particle-energy $nmax --two-particle-energy $nmax --max-loaded-memory $max_loaded_memory --finished-blocks-file finished_blocks $two_nucleon_force_path $three_nucleon_force_path
	
	# Creating the evaluation order
	$mars $comb_file_path evaluation_order --num-protons $Z --num-neutrons $N

	# Formating the index lists
	$aurora $comb_file_path $index_lists_path index_lists $Z $N
else
	# 2NF only case
	# Creating the M-scheme forces
	$mercury $comb_file_path $index_lists_path m_scheme_interaction --num-protons $Z --num-neutrons $N --single-particle-energy $nmax --two-particle-energy $nmax --max-loaded-memory $max_loaded_memory $two_nucleon_force_path
	
	# Creating the evaluation order
	$mars $comb_file_path evaluation_order --num-protons $Z --num-neutrons $N --only-two-nucleon-forces

	# Formating the index lists
	$aurora $comb_file_path $index_lists_path index_lists $Z $N --no-3NF
fi

# Preparing for Lanczos
mkdir -p krylow_vectors
mkdir -p eigen_vectors

# Setting up bacchus configuration file
echo "interaction:" > bacchus.conf
echo '{' >> bacchus.conf
echo "combination_table_file = \"$comb_file_path\";" >> bacchus.conf
echo "evaluation_order_file = \"evaluation_order\";" >> bacchus.conf
echo "index_lists_base_directory = \"index_lists\";" >> bacchus.conf
echo "matrix_file_base_directory = \"m_scheme_interaction\";" >> bacchus.conf
echo "num_neutrons = $Z;" >> bacchus.conf
echo "num_protons = $N;" >> bacchus.conf
echo "}" >> bacchus.conf
echo "lanczos:" >> bacchus.conf
echo "{" >> bacchus.conf
echo "krylow_vector_directory = \"krylow_vectors\";" >> bacchus.conf
echo "max_num_lanczos_iterations = 20;" >> bacchus.conf
echo "convergence_tolerance = 1.0e-7;" >> bacchus.conf
echo "target_eigenvector = 1;" >> bacchus.conf
echo "converge_eigenvectors = false;" >> bacchus.conf
echo "eigenvector_directory  = \"eigen_vectors\";" >> bacchus.conf
echo "}" >> bacchus.conf
$bacchus --max-memory-load $max_loaded_memory | tee bacchus_results
