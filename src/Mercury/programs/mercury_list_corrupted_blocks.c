#include <stdlib.h>
#include <stdio.h>
#include <combination_table/combination_table.h>
#include <iterator/iterator.h>
#include <matrix_block_setting/matrix_block_setting.h>
#include <string.h>
#include <errno.h>

static
void print_usage(const char *program_name);

int main(int num_arguments,
	 char **argument_list)
{
	if (num_arguments != 5)
	{
		print_usage(*argument_list);
		return EXIT_FAILURE;
	}
	const char *combination_table_path = argument_list[1];
	const char *interaction_path = argument_list[2];
	const size_t num_protons = atoll(argument_list[3]);
	const size_t num_neutrons = atoll(argument_list[4]);
	combination_table_t table =
		new_combination_table(combination_table_path,
				      num_protons,
				      num_neutrons);
	iterator_t matrix_blocks = new_matrix_block_settings_iterator(table);
	matrix_block_setting_t current_block;
	char *filename_buffer = 
		(char*)calloc(strlen(interaction_path)+128,sizeof(char));
	FILE *human_readable_output = fopen("corrupted_blocks.txt","w");
	FILE *output = fopen("corrupted_blocks","w");
	for (initialize(matrix_blocks,&current_block);
	     has_next_element(matrix_blocks);
	     next_element(matrix_blocks,&current_block))
	{
		sprintf(filename_buffer,
			"%s/%lu_matrix_elements",
			interaction_path,
			current_block.matrix_block_id);	
		FILE *file = fopen(filename_buffer,"r");
		if (file == NULL)
		{
			fprintf(human_readable_output,
				"Could not open matrix file %lu. %s\n",
				current_block.matrix_block_id,
				strerror(errno));	
			fwrite(&current_block.matrix_block_id,
			       sizeof(size_t),
			       1,
			       output);
			continue;
		}
		fseek(file,0,SEEK_END);
		size_t file_length = ftell(file);
		fseek(file,0,SEEK_SET);
		size_t expected_lenth = 
			get_matrix_block_length(current_block)*sizeof(double)+
			2*sizeof(size_t);
		if (file_length != expected_lenth)
		{
	
			fprintf(human_readable_output,
				"Block %lu is of length %lu B but should be %lu B.\n",
				current_block.matrix_block_id,
				file_length,
				expected_lenth);
			fwrite(&current_block.matrix_block_id,
			       sizeof(size_t),
			       1,
			       output);
		}
		fclose(file);
	}
	fclose(output);
	fclose(human_readable_output);
	free(filename_buffer);
	free_iterator(matrix_blocks);
	free_combination_table(table);
	return EXIT_SUCCESS;
}

static
void print_usage(const char *program_name)
{
	printf("%s <comb.txt> <interaction> <Z> <N>\n",
	       program_name);
}
