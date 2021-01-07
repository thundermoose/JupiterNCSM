#ifndef __SUBSPACE_OPERATOR__
#define __SUBSPACE_OPERATOR__

#include <stdlib.h>
#include <vector/vector.h>
#include <evaluation_order/evaluation_order.h>

void create_subspace_operator(const char *subspace_operator_path,
			      const char *operator_path,
			      const char *workspace_path,
			      const char *index_list_base_directory,
			      combination_table_t combination_table,
			      evaluation_order_t evaluation_order,
			      vector_settings_t vector_setting,
			      vector_t *training_vectors,
			      size_t num_training_vectors,
			      size_t maximum_loaded_memory);

#endif
