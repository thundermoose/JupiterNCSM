#include <settings/settings.h>
#include <libconfig.h>
#include <string.h>
#include <error/error.h>
#include <string_tools/string_tools.h>

struct _settings_
{
	int should_show_help_text;
	char *program_name;
	char *combination_table_path;
	char **training_vector_paths;
	char *norm_matrix_path;
	char *workspace_path;
	char *index_list_path;
	char **evaluation_order_path;
	char **operator_paths;
	char **subspace_operator_paths;
	size_t num_protons;
	size_t num_neutrons;
	size_t num_training_vectors;
	size_t num_operators;
	size_t max_loaded_memory;
};

settings_t parse_settings(size_t num_arguments, char **argument_list)
{
	settings_t settings = (settings_t)calloc(1,sizeof(struct _settings_));
	settings->program_name = argument_list[0];
	if (num_arguments>1)
	{
		settings->should_show_help_text =  
			strcmp(argument_list[1],"-h") == 0 || 
			strcmp(argument_list[1],"--help") == 0;
		if (settings->should_show_help_text)
			return settings;
	}

	char settings_file[256] = "ceres.conf";
	if (num_arguments > 2 && 
	    (strcmp(argument_list[1],"-c") == 0 ||
	     strcmp(argument_list[1],"--config") == 0))
		memcpy(settings_file,
		       argument_list[2],
		       strlen(argument_list[2]));
	config_t config;
	config_init(&config);
	if (config_read_file(&config,settings_file) == CONFIG_FALSE)
		error("Could not open file %s. %s\n",
		      settings_file,
		      config_error_text(&config));
	char *string_buffer = NULL;
	if (config_lookup_string(&config,
				 "combination_table_path",
				 (const char**)&string_buffer) == CONFIG_FALSE)
		error("Could not get combination_table_path from %s. %s\n",
		      settings_file,
		      config_error_text(&config));
	settings->combination_table_path = copy_string(string_buffer);
	if (config_lookup_string(&config,
				 "norm_matrix_path",
				 (const char**)&string_buffer) == CONFIG_FALSE)
		error("Could not get norm_matrix_path from %s. %s\n",
		      settings_file,
		      config_error_text(&config));
	settings->norm_matrix_path = copy_string(string_buffer);
	if (config_lookup_string(&config,
				 "indexlists_path",
				 (const char**)&string_buffer) == CONFIG_FALSE)
		error("Could not get indexlists_path from %s. %s\n",
		      settings_file,
		      config_error_text(&config));
	settings->index_list_path = copy_string(string_buffer);
	if (config_lookup_string(&config,
				 "workspace_path",
				 (const char**)&string_buffer) == CONFIG_FALSE)
		error("Could not get workspace_path fromm %s. %s\n",
		      settings_file,
		      config_error_text(&config));
	settings->workspace_path = copy_string(string_buffer);
	if (config_lookup_int64(&config,
				"num_protons",
				(long long int*)
				&settings->num_protons) == CONFIG_FALSE)
		error("Could not get num_protons from %s. %s\n",
		      settings_file,
		      config_error_text(&config));
	if (config_lookup_int64(&config,
				"num_neutrons",
				(long long int*)
				&settings->num_neutrons) == CONFIG_FALSE)
		error("Could not get num_neutrons from %s. %s\n",
		      settings_file,
		      config_error_text(&config));
	if (config_lookup_string(&config,
				 "max_loaded_memory",
				 (const char**)&string_buffer) == CONFIG_FALSE)
		error("Could not get max_loaded_memory from %s. %s\n",
		      settings_file,
		      config_error_text(&config));			
	settings->max_loaded_memory = parse_memory_string(string_buffer);
	config_setting_t *training_vectors =
		config_lookup(&config,"training_vectors");
	if (training_vectors == NULL)
		error("Could not find training_vectors in %s. %s\n",
		      settings_file,
		      config_error_text(&config));
	if (config_setting_is_list(training_vectors) ||
	    config_setting_is_array(training_vectors))
	{
		settings->num_training_vectors = 
			config_setting_length(training_vectors);
		settings->training_vector_paths = 
			(char**)
			malloc(settings->num_training_vectors*sizeof(char*));
		for (size_t i = 0; i<settings->num_training_vectors; i++)
		{
			settings->training_vector_paths[i] =
				copy_string
				(config_setting_get_string_elem
				 (training_vectors,i));
		}
	}
	else if (config_setting_type(training_vectors) == CONFIG_TYPE_STRING)
	{
		settings->num_training_vectors = 1;
		settings->training_vector_paths = 
			(char**)malloc(sizeof(char*));
		settings->training_vector_paths[0] =
			copy_string
			(config_setting_get_string(training_vectors));
	}
	else
		error("training_vector in %s is not list, array or string"
		      " and is therefore invalid\n",
		      settings_file);
	config_setting_t *operators = config_lookup(&config,"operators");
	if (operators == NULL)
		error("Could not find operators in %s. %s\n",
		      settings_file,
		      config_error_text(&config));
	if (config_setting_is_list(operators))
	{
		settings->num_operators = config_setting_length(operators);
		settings->operator_paths = 
			(char**)malloc(settings->num_operators*sizeof(char*));
		settings->subspace_operator_paths = 
			(char**)malloc(settings->num_operators*sizeof(char*));
		settings->evaluation_order_path =
			(char**)malloc(settings->num_operators*sizeof(char*));
		for (size_t i = 0; i<settings->num_operators; i++)
		{
			config_setting_t *operator_setting =
				config_setting_get_elem(operators,i);
			if (!config_setting_is_group(operator_setting))
				error("Element %lu of operators in %s "
				      "is not a group\n", i, settings_file);
			if (config_setting_lookup_string
			    (operator_setting,
			     "operator_path",
			     (const char**)&string_buffer) == CONFIG_FALSE)
				error("Could not read operator path from "
				      "element %lu in operators in %s. %s\n",
				      i,settings_file,
				      config_error_text(&config));
			settings->operator_paths[i] =
			       	copy_string(string_buffer);
			if (config_setting_lookup_string
			    (operator_setting,
			     "subspace_operator_path",
			     (const char**)&string_buffer) == CONFIG_FALSE)
				error("Could not read operator path from "
				      "element %lu in operators in %s. %s\n",
				      i,settings_file,
				      config_error_text(&config));
			settings->subspace_operator_paths[i] =
				copy_string(string_buffer);
			if (config_setting_lookup_string
			    (operator_setting,
			     "evaluation_order_path",
			     (const char**)&string_buffer) == CONFIG_FALSE)
				error("Could not read evaluation_order_path "
				      "from element %lu in operators in "
				      "%s. %s\n",
				      i,settings_file,
				      config_error_text(&config));
			settings->evaluation_order_path[i] =
				copy_string(string_buffer);
		}
	}
	else if (config_setting_is_group(operators))
	{
		settings->num_operators = 1;
		settings->operator_paths = 
			(char**)malloc(sizeof(char*));
		settings->subspace_operator_paths = 
			(char**)malloc(sizeof(char*));
		if (config_setting_lookup_string
		    (operators,
		     "operator_path",
		     (const char**)&string_buffer) == CONFIG_FALSE)
			error("Could not read operator path from "
			      "in operators in %s. %s\n",
			      settings_file,
			      config_error_text(&config));
		settings->operator_paths[0] =
			copy_string(string_buffer);
		if (config_setting_lookup_string
		    (operators,
		     "subspace_operator_path",
		     (const char**)&string_buffer) == CONFIG_FALSE)
			error("Could not read operator path from "
			      "in operators in %s. %s\n",
			      settings_file,
			      config_error_text(&config));
		settings->subspace_operator_paths[0] =
			copy_string(string_buffer);

	}
	else
		error("operators in %s is neither an array or a group\n",
		      settings_file);
	config_destroy(&config);
	return settings;
}

int settings_should_show_help_text(const settings_t settings)
{
	return settings->should_show_help_text;
}

void settings_show_help_text(const settings_t settings)
{
	printf("Usage: %s [-h/--help] [-c/--config <configuration_file>]\n"
	       "This program transforms operators on to a subspace\n"
	       "defined by a set of given basis vectors.\n"
	       "Configurations are given either in a file named \"ceres.conf\""
	       " or a file provided with the flags \"-c\" or \"--config\".\n"
	       "A configuration file follow the format syntax provided by"
	       " libconfig, and should contain the following variables:\n"
	       "\tcombination_table_path: A string containing the path to \n"
	       "\t\tthe \"comb.txt\" file generated by anicr\n"
	       "\tnorm_matrix_path: A string containing the path to store\n"
	       "\t\tthe resulting norm-matrix for the given basis vectors\n"
	       "\tindexlists_path: The path where the index lists are stored\n"
	       "\tnum_protons: The proton number of the given nucleus\n"
	       "\tnum_neutrons: The neutron number of the given nucleus\n"
	       "\tmax_loaded_memory: A string setting a limit on how much\n"
	       "\t\tof RAM may be used for storing arrays during the matrix-\n"
	       "\t\tvector multiplications\n"
	       "\ttraining_vectors: Either a string, a list of strings or an\n"
	       "\t\tarray of strings listing the paths to the basis vectors\n"
	       "\t\tof the subspace\n"
	       "\toperators: Either a group or list of groups. The groups\n"
	       "\t\t must contain two strings:\n"
	       "\t\toperator_path: Where the matrix elements needed for\n"
	       "\t\t\tthe matrix-vector multiplications\n"
	       "\t\tsubspace_operator_path: Where to store the resulting\n"
	       "\t\t\tsubspace-operator\n"
	       "\t\tevaluation_order_path: Path to the file containing the\n"
	       "\t\t\torder in which the matrix vector multiplications should\n"
	       "\t\t\tdone\n"
	       "The norm-matrix and the subspace-operators are store as\n"
	       "numpy matrices and can be loaded directly in to python\n",
		settings->program_name);
}

const char *settings_combination_table_path(const settings_t settings)
{
	return settings->combination_table_path;
}

const char *settings_training_vector_path(const settings_t settings,
					  size_t index)
{
	return settings->training_vector_paths[index];
}

const char *settings_norm_matrix_path(const settings_t settings)
{
	return settings->norm_matrix_path;
}

const char *settings_workspace_path(const settings_t settings)
{
	return settings->workspace_path;
}

const char *settings_index_list_path(const settings_t settings)
{
	return settings->index_list_path;
}

const char *settings_operator_path(const settings_t settings, size_t index)
{
	return settings->operator_paths[index];
}

const char *settings_subspace_operator_path(const settings_t settings,
					    size_t index)
{
	return settings->subspace_operator_paths[index];
}

const char *settings_evaluation_order_path(const settings_t settings,
					   size_t index)
{
	return settings->evaluation_order_path[index];
}

size_t settings_num_protons(const settings_t settings)
{
	return settings->num_protons;
}

size_t settings_num_neutrons(const settings_t settings)
{
	return settings->num_neutrons;
}

size_t settings_num_training_vectors(const settings_t settings)
{
	return settings->num_training_vectors;
}

size_t settings_num_operators(const settings_t settings)
{
	return settings->num_operators;
}

size_t settings_max_loaded_memory(const settings_t settings)
{
	return settings->max_loaded_memory;
}

void free_settings(settings_t settings)
{
	free(settings->combination_table_path);
	free(settings->norm_matrix_path);
	free(settings->index_list_path);
	free(settings->evaluation_order_path);
	for (size_t i = 0; i<settings->num_training_vectors; i++)
		free(settings->training_vector_paths[i]);
	free(settings->training_vector_paths);
	for (size_t i = 0; i<settings->num_operators; i++)
	{
		free(settings->operator_paths[i]);
		free(settings->subspace_operator_paths[i]);
	}	
	free(settings->operator_paths);
	free(settings->subspace_operator_paths);
}
