#ifndef __READ_3NF_ASCII_FORMAT__
#define __READ_3NF_ASCII_FORMAT__
#include "../bases/shells.h"
#include "../bases/jjj_coupled_3p.h"
#include "../matrix_transform/matrix_transform.h"

typedef int basis_type;
#define NO_BASIS 0
#define JJJ_BASIS 1
#define JT_BASIS 2

typedef struct _ascii_data_
{
  basis_type basis_format;
  Shells *shells;
  void *basis;

  size_t *ket_configs;
  size_t *bra_configs;
  double *elements;
  size_t num_elements;
} ASCII_Data;

ASCII_Data* open_ascii_data(const char* directory_name);

Dens_Matrix* get_matrix_ascii(ASCII_Data* data_file,
			      void* m_basis,
			      void* n_basis);



void free_ascii_data(ASCII_Data* data_file);

#endif
