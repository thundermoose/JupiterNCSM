#include "directory_tools.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <string.h>
#include <errno.h>

int directory_exists(const char *directory_name)
{
	struct stat directory_status;
	return stat(directory_name,
			&directory_status) == 0 &&
		S_ISDIR(directory_status.st_mode);
}

void create_directory(const char *directory_name)
{
	mode_t directory_mode = S_IFDIR | 0777;
	if (mkdir(directory_name,directory_mode) == -1)
	{
		fprintf(stderr,"Could not create directory \"%s\", %s.\n",
				directory_name,
				strerror(errno));
		exit(EXIT_FAILURE);
	}
}
