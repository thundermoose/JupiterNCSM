#include <settings/settings.h>
#include <error/error.h>
#include <string.h>

struct _settings_
{
	char *program_name;
     	char *combination_table_path;
	char *output_path;	
	size_t num_protons;
	size_t num_neutrons;
	int mode;
};

settings_t parse_settings(size_t num_arguments,
			  char **argument_list)
{
	settings_t settings = (settings_t)malloc(sizeof(struct _settings_));
	settings->program_name = *argument_list;
	if (num_arguments < 2)
	{
		settings->mode = 0;
	}
	else
	{
		settings->combination_table_path = argument_list[1];
		settings->output_path = argument_list[2];
		settings->num_protons = 0;
		settings->num_neutrons = 0;
		settings->mode = 3;
		for (size_t i = 3; i < num_arguments; i++)
		{
			if (strcmp(argument_list[i],"--num-protons") == 0)
				settings->num_protons =
				       	atoi(argument_list[++i]);
			else if (strcmp(argument_list[i],
					"--num-neutrons") == 0)
				settings->num_neutrons =
					atoi(argument_list[++i]);
			else if (strcmp(argument_list[i],
					"--only-two-nucleon-forces") == 0)
				settings->mode = 2;
			else
				error("Unknown argument \"%s\"\n",
				      argument_list[i]);
		}
	}
	return settings;
}

int should_display_usage(settings_t settings)
{
	return settings->mode == 0;
}

void display_usage(settings_t settings)
{
	printf("Usage: %s <combination table file> <output file> "
	       "[--num-protons <integer>] "
	       "[--num-neutrons <integer>] "
	       "[--only-two-nucleon-forces]\n",
	       settings->program_name);
}

const char *get_combination_table_path(settings_t settings)
{
	return settings->combination_table_path;
}

const char *get_output_path(settings_t settings)
{
	return settings->output_path;
}

size_t get_num_protons_argument(settings_t settings)
{
	return settings->num_protons;
}

size_t get_num_neutrons_argument(settings_t settings)
{
	return settings->num_neutrons;
}

int two_nucleon_force_only_mode(settings_t settings)
{
	return settings->mode == 2;
}

void free_settings(settings_t settings)
{
	free(settings);
}
