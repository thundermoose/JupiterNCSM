#ifndef __SUBSPACE_OPERATOR__
#define __SUBSPACE_OPERATOR__

void create_subspace_operator(const char *subspace_operator_path,
			      const char *operator_path,
			      const char *workspace_path,
			      vector_t *training_vectors,
			      size_t num_training_vectors);

#endif
