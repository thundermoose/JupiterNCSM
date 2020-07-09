#ifndef __ERROR__
#define __ERROR__

#include <stdlib.h>
#include <stdio.h>

#define error(message...)\
{\
	fprintf(stderr,"Error in %s @ %s:%d :",\
		__func__,__FILE__,__LINE__);\
	fprintf(stderr,message);\
	exit(EXIT_FAILURE);\
}

#endif
