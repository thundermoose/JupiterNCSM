#include <stdlib.h>
#include <stdio.h>
#include <log/log.h>
__attribute__((constructor(101)))
void initialization()
{
	initiate_logging("MINERVA_LOGFILE",
			 "minerva.log");
}
int main()
{
	printf("I am Minerva\n");
	return EXIT_SUCCESS;
}
