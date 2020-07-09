#include "argument_buffer.h"

Argument_Buffer *empty_argument_buffer()
{
  return (Argument_Buffer*)calloc(1,sizeof(Argument_Buffer));
}

Argument_Buffer *new_argument_buffer(size_t num_arguments,
				     char **arguments)
{
  Argument_Buffer *argument_buffer = empty_argument_buffer();
  argument_buffer->num_arguments = num_arguments;
  argument_buffer->arguments = arguments;
  return argument_buffer;
}
   
char *next_argument(Argument_Buffer* argument_buffer)
{
  if (argument_buffer->next_argument_pos == argument_buffer->num_arguments)
    return NULL;
  return argument_buffer->arguments[argument_buffer->next_argument_pos++];
}
