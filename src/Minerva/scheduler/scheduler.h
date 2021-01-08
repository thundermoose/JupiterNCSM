#ifndef __SCHEDULER__
#define __SCHEDULER__

#include <evaluation_order/evaluation_order.h>
#include <combination_table/combination_table.h>

struct _scheduler_;
typedef struct _scheduler_ *scheduler_t;

scheduler_t new_scheduler(evaluation_order_t evaluation_order,
			  combination_table_t combination_table,
			  const char *index_lists_base_directory,
			  const char *matrix_file_base_directory,
			  size_t maximum_loaded_memory);

void run_matrix_vector_multiplication(const char *output_vector_base_directory,
				      const char *input_vector_base_directory,
				      scheduler_t scheduler);

void free_scheduler(scheduler_t scheduler);

#endif
