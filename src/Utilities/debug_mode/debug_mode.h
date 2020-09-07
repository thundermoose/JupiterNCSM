#ifndef __DEBUG_MODE_
#define __DEBUG_MODE_

#include <stdlib.h>
#ifndef NDEBUG
#ifndef BACKEND_DEBUG_MODE
#define malloc(bytes_to_allocate) logging_malloc(bytes_to_allocate,\
						 __builtin_FILE(),\
						 __builtin_FUNCTION(),\
						 __builtin_LINE())
#define free(pointer_to_free) logging_free(pointer_to_free,\
					   __builtin_FILE(),\
					   __builtin_FUNCTION(),\
					   __builtin_LINE())
#endif

void *logging_malloc(size_t bytes_to_allocate,
		     const char *file_name,
		     const char *function_name,
		     const int line_number);

void logging_free(void *pointer_to_free,
		  const char *file_name,
		  const char *function_name,
		  const int line_number);

#endif

#endif
