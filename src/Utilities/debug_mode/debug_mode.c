#define BACKEND_DEBUG_MODE
#include <debug_mode/debug_mode.h>
#include <log/log.h>
#ifndef NDEBUG
void *logging_malloc(size_t bytes_to_allocate,
		     const char *file_name,
		     const char *function_name,
		     const int line_number)
{
	void *allocated_pointer = malloc(bytes_to_allocate);
	static size_t num_mallocs = 0;
	num_mallocs++;
	log_log(file_name,
		function_name,
		line_number,
		num_mallocs,
		"%p = malloc(%lu) memory",
		allocated_pointer,
		bytes_to_allocate);
	return allocated_pointer;
}

void logging_free(void *pointer_to_free,
		  const char *file_name,
		  const char *function_name,
		  const int line_number)
{
	static size_t num_frees = 0;
	num_frees++;
	log_log(file_name,
		function_name,
		line_number,
		num_frees,
		"free(%p) memory",
		pointer_to_free);
	free(pointer_to_free);
}
#endif
