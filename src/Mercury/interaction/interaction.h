#ifndef __INTERACTION__
#define __INTERACTION__

#include <stdlib.h>

struct _interaction_;
typedef struct _interaction_ *interaction_t;

interaction_t new_interaction(const char *interaction_path);

double get_matrix_element(interaction_t interaction,
			  int *bra_state,
			  int *ket_state,
			  size_t num_particles,
			  int E1,int E2,int Tz, int M);

void free_interaction(interaction_t interaction);

#endif
