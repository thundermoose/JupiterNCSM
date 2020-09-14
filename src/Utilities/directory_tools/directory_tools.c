#include <directory_tools/directory_tools.h>

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
	if (mkdir(directory_name,
		  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
		return 1;
	else
		return 0;
}

void clear_directory(const char *directory_name)
{
	if (!directory_exists(directory_name))
		return;
	char *command = (char*)calloc(strlen(directory_name)+20,
				      sizeof(char));
	sprintf(command,
		"rm -rf %s*\n",directory_name);
	system(command);
	free(command);
}
