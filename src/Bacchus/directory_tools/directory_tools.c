#include <directory_tools/directory_tools.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <string.h>
#include <error/error.h>

int directory_exists(const char *directory_name)
{
	DIR* directory = opendir(directory_name);
	if (directory)
	{
		closedir(directory);
		return 1;
	}
	else if (errno == ENOENT)
		return 0;
	else
		error("opendir(%s) failed. %s.\n",
		      directory_name,
		      strerror(errno));
}

int create_directory(const char *directory_name)
{
	return mkdir(directory_name,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}
