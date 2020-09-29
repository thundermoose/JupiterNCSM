#ifndef __DEBUG_MODE_
#define __DEBUG_MODE_

#include <stdlib.h>
#include <stdio.h>
#ifdef DEBUG
#ifndef BACKEND_DEBUG_MODE
#define malloc(bytes_to_allocate) logging_malloc(bytes_to_allocate,\
						 __builtin_FILE(),\
						 __builtin_FUNCTION(),\
						 __builtin_LINE())
#define calloc(num_elements,size) logging_calloc(num_elements,size,\
						 __builtin_FILE(),\
						 __builtin_FUNCTION(),\
						 __builtin_LINE())
#define realloc(origin,new_size) logging_realloc(origin,\
						 new_size,\
						 __builtin_FILE(),\
						 __builtin_FUNCTION(),\
						 __builtin_LINE())
#define free(pointer_to_free) logging_free(pointer_to_free,\
					   __builtin_FILE(),\
					   __builtin_FUNCTION(),\
					   __builtin_LINE())


#define fopen(pathname,mode) logging_fopen(pathname,mode,\
					   __builtin_FILE(),\
					   __builtin_FUNCTION(),\
					   __builtin_LINE())

#define fclose(file_handle) logging_fclose(file_handle,\
					   __builtin_FILE(),\
					   __builtin_FUNCTION(),\
					   __builtin_LINE())
#endif

void *logging_malloc(size_t bytes_to_allocate,
		     const char *file_name,
		     const char *function_name,
		     const int line_number);

void *logging_calloc(size_t num_elements,
		     size_t size,
		     const char *file_name,
		     const char *function_name,
		     const int line_number);

void *logging_realloc(void *origin,
		      size_t new_size,
		      const char *file_name,
		      const char *function_name,
		      const int line_number);
void logging_free(void *pointer_to_free,
		  const char *file_name,
		  const char *function_name,
		  const int line_number);


FILE *logging_fopen(const char *path_name,
		    const char *mode,
		    const char *file_name,
		    const char *function_name,
		    const int line_number);

void logging_fclose(FILE *file_handle,
		    const char *file_name,
		    const char *function_name,
		    const int line_number);

#endif

#endif
