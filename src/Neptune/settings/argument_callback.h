#ifndef __ARGUMENT_CALLBACK__
#define __ARGUMENT_CALLBACK__

#include "settings.h"
#include "argument_buffer.h"

typedef void(*Argument_Callback)(Settings *,
				 Argument_Buffer *);

Argument_Callback *identify_argument(const char *argument);

#endif
