#ifndef __DIRECTORY_TOOLS__
#define __DIRECTORY_TOOLS__


#include <stdlib.h>
#include <stdio.h>
#include <dirent.h>
#include <error/error.h>
#include <errno.h>
#include <string.h>
#include <sys/stat.h>

/*
#define error(message...)\
{\
	fprintf(stderr,"Error %s@%s:%d : ",\
		__func__,__FILE__,__LINE__);\
	fprintf(stderr,message);\
	exit(EXIT_FAILURE);\
}
*/

int directory_exists(const char *directory_name);

int create_directory(const char *directory_name);

void clear_directory(const char *directory_name);

#endif
