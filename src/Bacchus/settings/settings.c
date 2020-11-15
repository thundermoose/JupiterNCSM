#include <settings/settings.h>
#include <libconfig.h>
#include <string.h>
#include <string_tools/string_tools.h>
#include <log/log.h>
#include <error/error.h>

struct _settings_
{
	char *program_name;
	int show_help;
	char *combination_table_path;
	char *execution_order_path;
	char *index_lists_base_directory;
	char *matrix_file_base_directory;
	char *krylow_vector_directory;
	size_t num_neutrons;
	size_t num_protons;
	size_t max_num_lanczos_iterations;
	size_t maximum_loaded_memory;
	double tollerance;
};

settings_t parse_settings(size_t num_arguments,
			  char **argument_list)
{
	char *settings_file_name = "bacchus.conf";
	int show_help = 0;
	size_t maximum_loaded_memory = 0;
	for (size_t i = 1; i<num_arguments; i++)
	{
		if (strcmp(argument_list[i],"--settings-file") == 0)
		{
			settings_file_name = argument_list[++i];	
			if (!ends_with(settings_file_name, ".conf"))
			{
				error("\"%s\" is not a settings file,"
				      " the filename must end in .conf\n",
				      settings_file_name);
			}
		}
		else if (strcmp(argument_list[i],"--help") == 0 ||
			 strcmp(argument_list[i],"-h") == 0)
		{
			show_help = 1;	
		}
		else if (strcmp(argument_list[i],"--max-memory-load") == 0)
		{
			char *memory_string = argument_list[++i];
			if (!is_memory_string(memory_string))
				error("--max-memory-load followed by unknown "
				      " string \"%s\".\n",
				      memory_string);
			maximum_loaded_memory = 
				parse_memory_string(memory_string);		
		}
		else
		{
			error("Unknown argument \"%s\"\n",
			      argument_list[i]);
		}
	}	
	settings_t settings = (settings_t)malloc(sizeof(struct _settings_));
	settings->program_name = argument_list[0];
	settings->show_help = show_help;
	if (show_help)
	{
		return settings;
	}
	config_t config;
	config_init(&config);
	if (config_read_file(&config,settings_file_name) == CONFIG_FALSE)
		error("Could not read config \"%s\". %s\n",
		      settings_file_name,
		      config_error_text(&config));
	config_setting_t *interaction_setting =
	       	config_lookup(&config,"interaction");
	config_setting_t *lanczos_setting = 
		config_lookup(&config,"lanczos");
	char *string_buffer = NULL;
	if (config_setting_lookup_string(interaction_setting,
					 "combination_table_file",
					 (const char **)
					 &string_buffer)
	    == CONFIG_FALSE)
		error("No interaction.combination_table_file found in \"%s\"."
		      " %s\n",
		      settings_file_name,
		      config_error_text(&config));
	settings->combination_table_path = copy_string(string_buffer);
	log_entry("settings->combination_table_path = %s",
		  settings->combination_table_path);
	if (config_setting_lookup_string(interaction_setting,
					 "execution_order_file",
					 (const char **)
					 &string_buffer)
	    == CONFIG_FALSE)
		error("No interaction.execution_order_file found in \"%s\"."
		      " %s\n",
		      settings_file_name,
		      config_error_text(&config));
	settings->execution_order_path = copy_string(string_buffer);
	if (config_setting_lookup_string(interaction_setting,
					 "index_lists_base_directory",
					 (const char **)
					 &string_buffer)
	    == CONFIG_FALSE)
		error("No interaction.index_lists_base_directory found in "
		      "\"%s\". %s\n",
		      settings_file_name,
		      config_error_text(&config));
	settings->index_lists_base_directory = copy_string(string_buffer);
	if (config_setting_lookup_string(interaction_setting,
					 "matrix_file_base_directory",
					 (const char **)
					 &string_buffer)
	    == CONFIG_FALSE)
		error("No interaction.matrix_file_base_directory found in "
		      "\"%s\". %s\n",
		      settings_file_name,
		      config_error_text(&config));
	settings->matrix_file_base_directory = copy_string(string_buffer);
	if (config_setting_lookup_int64(interaction_setting,
					"num_neutrons",
					(long long*)
					&settings->num_neutrons)
	    == CONFIG_FALSE)
		error("No interaction.num_neutrons found in \"%s\"."
		      " %s",
		      settings_file_name,
		      config_error_text(&config));
	if (config_setting_lookup_int64(interaction_setting,
					"num_protons",
					(long long*)
					&settings->num_protons)
	    == CONFIG_FALSE)
		error("No interaction.num_protons found in \"%s\"."
		      " %s",
		      settings_file_name,
		      config_error_text(&config));
	if (config_setting_lookup_string(lanczos_setting,
					 "krylow_vector_directory",
					 (const char **)
					 &string_buffer)
	    == CONFIG_FALSE)
		error("No lanczos.krylow_vector_directory found in \"%s\"."
		      " %s\n",
		      settings_file_name,
		      config_error_text(&config));
	settings->krylow_vector_directory = copy_string(string_buffer);
	if (config_setting_lookup_int64(lanczos_setting,
					"max_num_lanczos_iterations",
					(long long*)
					&settings->max_num_lanczos_iterations)
	    == CONFIG_FALSE)
		error("No lanczos.max_num_lanczos_iterations found in \"%s\"."
		      " %s",
		      settings_file_name,
		      config_error_text(&config));
	if (config_setting_lookup_float(lanczos_setting,
				"convergence_tollerance",
				&settings->tollerance)
    		== CONFIG_FALSE)
    		error("No lanczos.convergence_tollerance found in \"%s\"."
		      " %s\n",
		      settings_file_name,
		      config_error_text(&config));
	if (config_setting_lookup_string(lanczos_setting,
					 "max_memory_load",
					 (const char **)
					 &string_buffer) == CONFIG_FALSE)
		settings->maximum_loaded_memory = (size_t)(16)<<30;
	else if (is_memory_string(string_buffer))
		settings->maximum_loaded_memory =
		       	parse_memory_string(string_buffer);
	else
		error("max_memory_load is not set to correct memory string\n");
	// Command argument has presidence over settings file
	if (maximum_loaded_memory > 0)
		settings->maximum_loaded_memory = maximum_loaded_memory;
	config_destroy(&config);
	return settings;	
}

int should_show_help_text(const settings_t settings)
{
	return settings->show_help;
}

void show_help_text(const settings_t settings)
{
	printf("Usage: %s [--settings-file <file-path>] [-h/--help] "
	       "[--max-memory-load <memory size>\n"
	       "Flags:\n"
	       "\t--settings-file <file-path>: To provide %s with an "
	       "alternative settings file than \"bacchus.conf\".\n"
	       "\t-h/--help: To display this message.\n"
	       "\t--max-memory-load <memory size>: To provide a different "
	       "limit on how much memory the matrix vector multiplication "
	       "uses\n"
	       "The settings file:\n"
	       "The settings file is read using the libconfig library."
	       "Therefore, the user is referred to the libconfig documentation"
	       "for details.\n"
	       "In the settings file there should be two groups, interaction"
	       " and lanczos. The interaction group should contain the"
	       " following fields:\n"
	       "\tcombination_table_file: A string containing the path to the"
	       " comb.txt file produced by anicr.\n"
	       "\texecution_order_file: A string containing the path to the"
	       " execution order file, e.g. greedy_2_16.txt.\n"
	       "\tindex_lists_base_directory: A string containing the path to"
	       " the base directory for the index lists generated by anicr.\n"
	       "\tmatrix_file_base_directory: A string containing the path to"
	       " the base direcotry for the matrix files generated by "
	       "mercury.\n"
	       "\tnum_neutrons: An integer giving the number of neutrons in"
	       " the nucleus of interest.\n"
	       "\tnum_protons: An integer giving the number of protons in"
	       " the nucleus of interest.\n"
	       "The lanczos group should contain the following fileds:\n"
	       "\tkrylow_vector_directory: A string containing the path to"
	       " where the Lanczos algorithm can store the krylow vectors on"
	       " disk."
	       "\tmax_num_lanczos_iterations: An integer limiting the maximum"
	       " number of iterations that the Lanczos algorithm should do\n"
	       "\tconvergence_tollerance: If the lowest eigenvalue differ with" 
	       " less than this number from the previous lowest eigenvalue,"
	       " the Lanczos algorithm is assumed to be converged\n",
		settings->program_name,
		settings->program_name);
}

const char *get_combination_table_path_setting(const settings_t settings)
{
	return settings->combination_table_path;
}

const char *get_execution_order_path_setting(const settings_t settings)
{
	return settings->execution_order_path;
}

const char *get_krylow_vector_directory_setting(const settings_t settings)
{
	return settings->krylow_vector_directory;
}

const char *get_index_lists_base_directory_setting(const settings_t settings)
{
	return settings->index_lists_base_directory;
}

const char *get_matrix_file_base_directory_setting(const settings_t settings)
{
	return settings->matrix_file_base_directory;
}

size_t get_num_neutrons_setting(const settings_t settings)
{
	return settings->num_neutrons;
}

size_t get_num_protons_setting(const settings_t settings)
{
	return settings->num_protons;
}

size_t get_max_num_lanczos_iterations_setting(const settings_t settings)
{
	return settings->max_num_lanczos_iterations;
}

size_t get_maximum_loaded_memory_setting(const settings_t settings)
{
	return settings->maximum_loaded_memory;
}

double get_tollerance_setting(const settings_t settings)
{
	return settings->tollerance;
}

void free_settings(settings_t settings)
{
	free(settings->combination_table_path);
	free(settings->execution_order_path);
	free(settings->index_lists_base_directory);
	free(settings->matrix_file_base_directory);
	free(settings->krylow_vector_directory);
	free(settings);
}
