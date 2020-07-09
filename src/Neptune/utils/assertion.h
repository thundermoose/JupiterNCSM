#ifndef __ASSERTION__
#define __ASSERTION__
#include <stdlib.h>
#include <stdio.h>

#include "terminal_colors.h"

#ifdef DEBUG
#define ASSERT(test,action,message...)			\
  if(!(test))						\
    {							\
      fprintf(stderr,RED "ASSERTION %s@%s:%d FAILED!\n",	\
	      __builtin_FUNCTION(),			\
	      __builtin_FILE(),				\
	      __builtin_LINE());			\
      fprintf(stderr,MAGENTA message);			\
      fprintf(stderr,NORMAL);				\
      action;						\
    }
#else
#define ASSERT(test,action,message...)
#endif
#endif
