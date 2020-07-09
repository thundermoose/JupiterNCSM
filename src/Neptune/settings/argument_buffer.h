#ifndef __ARGUMENT_BUFFER__
#define __ARGUMENT_BUFFER__
#include <stdlib.h>

typedef struct
{
  size_t num_arguments;
  size_t next_argument_pos;
  char **arguments;
} Argument_Buffer;

Argument_Buffer *new_argument_buffer(size_t num_arguments,
				     char **arguments);
char *next_argument(Argument_Buffer* argument_buffer);

#endif
