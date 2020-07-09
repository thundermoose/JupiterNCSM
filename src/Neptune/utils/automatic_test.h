#ifndef __AUTOMATIC_TEST__
#define __AUTOMATIC_TEST__

#include "terminal_colors.h"

void test_failed();

#define TEST_FAILED(message...)			\
  fprintf(stderr,RED "TEST %s@%s:%d FAILED\n",	\
	  __builtin_FUNCTION(),			\
	  __builtin_FILE(),			\
	  __builtin_LINE());			\
  fprintf(stderr,MAGENTA message);		\
  fprintf(stderr,NORMAL);			\
  test_failed();				\
  return;



#ifdef DEBUG

#define BEGIN_TEST(name)			\
  __attribute__((constructor(101)))		\
  void name()					\
  {						\
  printf(BLUE "RUNS TEST %s:\n" NORMAL,#name);



#else

#define BEGIN_TEST(name)			\
  __attribute__((cold))				\
  void name()					\
  {						\
  printf("RUNS TEST %s:\n",#name);


#endif

#define END_TEST				\
  printf(GREEN"TEST %s@%s PASSED!\n",		\
	 __builtin_FUNCTION(),			\
	 __builtin_FILE());			\
  printf(NORMAL);				\
  }

#endif
