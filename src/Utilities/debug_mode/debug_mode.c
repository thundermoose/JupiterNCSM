#define BACKEND_DEBUG_MODE
#include <debug_mode/debug_mode.h>
#include <log/log.h>
#include <string.h>
#include <errno.h>
#ifdef DEBUG

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

void *logging_calloc(size_t num_elements,
		     size_t size,
		     const char *file_name,
		     const char *function_name,
		     const int line_number)
{
	void *allocated_pointer = calloc(num_elements,size);
	static size_t num_callocs = 0;
	num_callocs++;
	log_log(file_name,
		function_name,
		line_number,
		num_callocs,
		"%p = calloc(%lu,%lu) memory",
		allocated_pointer,
		num_elements,
		size);
	return allocated_pointer;
}

void *logging_realloc(void *origin,
		      size_t new_size,
		      const char *file_name,
		      const char *function_name,
		      const int line_number)
{
	void *allocated_pointer = realloc(origin,new_size);
	static size_t num_reallocs = 0;
	num_reallocs++;
	log_log(file_name,
		function_name,
		line_number,
		num_reallocs,
		"%p = realloc(%p,%lu) memory",
		allocated_pointer,
		origin,
		new_size);
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

FILE *logging_fopen(const char *path_name,
		    const char *mode,
		    const char *file_name,
		    const char *function_name,
		    const int line_number)
{
	FILE *file_handle = fopen(path_name,mode);
	static size_t num_fopen = 0;
	num_fopen++;
	log_log(file_name,
		function_name,
		line_number,
		num_fopen,
		"%p = fopen(%s,%s) files",
		file_handle,
		path_name,mode);
	if (file_handle == NULL)
		log_log(file_name,
			function_name,
			line_number,
			1,
			"%s: %s",
			path_name,
			strerror(errno));
	return file_handle;
}

void logging_fclose(FILE *file_handle,
		    const char *file_name,
		    const char *function_name,
		    const int line_number)
{
	static size_t num_fclose = 0;
	num_fclose++;
	log_log(file_name,
		function_name,
		line_number,
		num_fclose,
		"fclose(%p) files",
		file_handle);
	fclose(file_handle);
}

#endif
