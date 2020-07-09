#include "automatic_test.h"
#include <stdlib.h>
#include <stdio.h>

int num_failed_tests = 0;

void test_failed()
{
#ifdef DEBUG
  num_failed_tests++;
#else
  exit(1);
#endif
}

#ifdef DEBUG
__attribute__((constructor(200)))
void test_passed_checker()
{
  if (num_failed_tests)
    {
      fprintf(stderr,RED"%d TESTS FAILED\n",
	      num_failed_tests);
      exit(1);
    }
  printf(GREEN "ALL TESTS PASSED\n" NORMAL);
}
#endif
