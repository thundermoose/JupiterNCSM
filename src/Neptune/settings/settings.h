#ifndef __SETTINGS__
#define __SETTINGS__

#include <stdlib.h>

typedef struct _settings_ {
  
} Settings;

Settings *read_settings(size_t num_arguments,
			char** arguments);

void free_settings(Settings *settings);

#endif
