#ifndef __READ_3NF_FILE__
#define __READ_3NF_FILE__
#include <stdlib.h>
#include <bases/jjj_coupled_3p.h>
#include <bases/shells.h>
#include <matrix_transform/matrix_transform.h>



typedef unsigned int format_signature;

#define HDF5 1
#define ASCII 2

typedef enum
{
	CE = 0,
	CD = 1,
	C1 = 2,
	C3 = 3,
	C4 = 4,
	U = 1000
} weight_t;

typedef struct _data_file_
{
	format_signature signature;
	void* data_pointer;
	quantum_number e_max;
} Data_File;

Data_File* open_data_file(const char* file_name);

void set_max_loaded_memory(Data_File* data_file,size_t max_loaded_memory);

weight_t identify_weight(const char *weight);

void set_weight(Data_File* file,
		weight_t weight,
		double value);
		


Dens_Matrix* get_matrix(Data_File* data_file,
		        void* m_basis,
			void* n_basis);


/* Check if j_scheme is a subset of the
 * j_scheme basis in datafile. If any
 * basis state in j_scheme do not exists in 
 * data file j_scheme is not a subset and this
 * method yields 0, otherwise 1
 */
int is_subset_of_basis(JJJ_Basis* j_scheme,
		       Data_File* datafile);


void free_data_file(Data_File* data_file);

#endif
