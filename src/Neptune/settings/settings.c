#include "settings.h"
#include <stdlib.h>
#include "argument_buffer.h"
#include "argument_callback.h"
#include <debug_mode/debug_mode.h>

Settings *empty_settings()
{
  return (Settings*)calloc(1,sizeof(Settings));
}

void parse_arguments(Settings *settings,
		     size_t num_arguments,
		     char** arguments)
{
  Argument_Buffer *argument_buffer =
    new_argument_buffer(num_arguments,
			arguments);
  char* argument = NULL;
  while((argument = next_argument(argument_buffer)))
    {
      Argument_Callback *argument_callback =
	identify_argument(argument);
      (*argument_callback)(settings,
			   argument_buffer);
    }
  free(argument_buffer);
}

Settings *read_settings(size_t num_arguments,
			char** arguments)
{
  Settings *settings = empty_settings();
  parse_arguments(settings,
		  num_arguments,
		  arguments);
  return settings;
}

void free_settings(Settings *settings)
{
  free(settings);
}
