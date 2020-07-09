#ifndef __FORTRAN_BLOCK__
#define __FORTRAN_BLOCK__

#include <stdio.h>

struct _fortran_block_;
typedef struct _fortran_block_ *fortran_block_t;

fortran_block_t read_fortran_block(FILE *file);

int get_num_bytes(fortran_block_t block);

const char *get_fortran_block_data(fortran_block_t block);

void free_fortran_block(fortran_block_t block);

void skip_fortran_block(FILE *file);

#endif
