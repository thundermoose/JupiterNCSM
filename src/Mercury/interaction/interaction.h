#ifndef __INTERACTION__
#define __INTERACTION__

#include <stdlib.h>
#include <interaction/header/header.h>

struct _interaction_;
typedef struct _interaction_ *interaction_t;

interaction_t new_interaction(const char *interaction_path);

const header_t get_header(const interaction_t interaction);

double get_matrix_element(interaction_t interaction,
			  int *bra_state,
			  int *ket_state,
			  size_t num_particles,
			  int E1,int E2,int Tz, int M);

void free_interaction(interaction_t interaction);

#endif
