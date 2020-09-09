#include <stdlib.h>
#include <stdio.h>
#include <log/log.h>

__attribute__((constructor(101)))
void initialization()
{
	initiate_logging("ASIMPLE_LOGFILE",
			 "asimple_program.log");
}

int main()
{
  printf("This is just a simple program\n");
}
