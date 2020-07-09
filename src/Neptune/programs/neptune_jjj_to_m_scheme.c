#include <stdlib.h>
#include <stdio.h>
// omp here
#include <input/read_3nf_file.h>
#include <bases/m_scheme_3p_basis.h>
#include <bases/jjj_coupled_3p.h>
#include <matrix_transform/matrix_transform.h>
#include <transform_scheduller/transform_scheduller_3p.h>
#include <output/out_file.h>
#include <thundertester/test.h>


void usage(char* prg)
{
	printf("Usage: %s <in_name> <output_name> [-f hdf5/unified-3p]"
	       "[-n <e-max-single-particle> <e-max-three-particle>] "
	       "[-w <kind> <value>]\n\n"
	       "The in_name is either a hdf5 file or a folder containing a .cnf and .mat file.\n",
	       prg);
}

int main(int argc, char** argv)
{
	if (argc<3)
	{
		usage(*argv);
		return 0;
	}

	// Instantiate the 3nf read system
	// Note for future: Should change 
	Data_File* data_3nf = open_data_file(argv[1]);

	// Get output filename

	char* output_filename = argv[2];

	// determine optional arguments

	quantum_number e_max1 = data_3nf->e_max;
	quantum_number e_max2 = data_3nf->e_max;
	file_type_t output_file_type = hdf5_file;
	size_t i;
	for (i = 3; i<argc; i++)
	{
		if (argv[i][0] == '-' && argv[i][1] == 'n')
		{
			e_max1 = atoi(argv[++i]);
			e_max2 = atoi(argv[++i]);
		}
		else if (argv[i][0] == '-' && argv[i][1] == 'w')
		{
			weight_t w = identify_weight(argv[++i]);
			set_weight(data_3nf,w,atof(argv[++i]));
		}
		else if (argv[i][0] == '-' && argv[i][1] == 'f')
			output_file_type = file_type_from_string(argv[++i]);
		else
		{
			fprintf(stderr,"Unknown option \"%s\"\n",argv[i]);
			return 1;
		}
	}
	// Setup the single particle bases
	Shells* shells = new_shells(e_max1);
	SP_States* sp_states = new_sp_states(shells);

	// Setup the m-scheme basis
	M_Scheme_3p_Basis *mp_basis = 
		new_m_scheme_3p_basis_no_m_rest(e_max2,
						sp_states);

	// Prepairing the outputfile for writing
	out_file_t outfile =
		create_new_out_file(output_filename,
				    mp_basis,
				    output_file_type);
	// Everything is now prepared
	// so we can just run the main
	// code
	transform_3p_data(data_3nf,
			  outfile,
			  mp_basis);
	// if we come this far
	// we can now finalize everything
	close_out_file(outfile);
	free_m_scheme_3p_basis(mp_basis);
	free_sp_states(sp_states);
	free_shells(shells);
	free_data_file(data_3nf);
	// Everything went as expected
	return 0;
}

new_test(test_of_test,
	 printf("test works\n");
	 assert_that(1 == 1);
	);
