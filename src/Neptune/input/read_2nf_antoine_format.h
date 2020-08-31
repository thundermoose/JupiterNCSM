#ifndef __READ_2NF_ANTOINE_FORMAT__
#define __READ_2NF_ANTOINE_FORMAT__

#include <stdlib.h>
#include <bases/jt_basis.h>
#include <matrix_transform/matrix_transform.h>
struct _antoine_2nf_file_;
typedef struct _antoine_2nf_file_ *antoine_2nf_file_t;

antoine_2nf_file_t open_antoine_2nf_file(const char* file_name,
					 size_t num_particles,
					 quantum_number e_max1,
					 quantum_number e_max2);

jt_basis_t get_jt_basis(antoine_2nf_file_t data_file);

Dens_Matrix *get_antoine_matrix(antoine_2nf_file_t data_file,
		jt_basis_t bra_basis,
		jt_basis_t ket_basis,
		quantum_number Tz);

void free_antoine_2nf_file(antoine_2nf_file_t data_file);

#endif
