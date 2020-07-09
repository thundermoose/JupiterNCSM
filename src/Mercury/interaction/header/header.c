#include <interaction/header/header.h>
#include <energy_block_info/energy_block_info.h>
#include <log/log.h>
#include <array_builder/array_builder.h>
#include <global_constants/global_constants.h>
#include <error/error.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>


struct _header_
{
	size_t num_particles;
	int energy_max;
	size_t dimension;
	size_t num_energy_blocks;
	energy_block_info_t *energy_block_info;
};


static
void next_row(char **row, size_t *row_length,FILE *file);

header_t read_header(const char *interaction_path)
{
	log_entry("Reading interaction header file");
	size_t file_name_length = strlen(interaction_path)+8;
	char *file_name = (char*)malloc(file_name_length);
	sprintf(file_name,
		"%s/header",
		interaction_path);
	FILE *header_file = fopen(file_name,"r");
	if (header_file == NULL)
		error("Could not open header file %s. %s\n",
		      file_name,
		      strerror(errno));
	size_t row_length = 0;
	char *row = NULL;
	next_row(&row,&row_length,header_file);		
	if (strcmp(row,"=====Header=====\n") != 0)
		error("The file %s is not a header file\n",
		      file_name);
	log_entry("File \"%s\" is a header file",
		  file_name);
	next_row(&row,&row_length,header_file);
	size_t num_particles = 0;
	if (sscanf(row,
		   "num particles: %lu",
		   &num_particles) != 1)
		error("Could not read number of particles from header\n");
	log_entry("Correctly parsed the number of particles to %lu",
		  num_particles);
	next_row(&row,&row_length,header_file);
	//if (strncmp(row+5,interaction_path,strlen(interaction_path)) != 0)
	//	error("Header directory path is not the same as path given\n");
	log_entry("Header file contains the correct path");
	next_row(&row,&row_length,header_file);
	int energy_max = 0;
	if (sscanf(row,
		   "Nmax: %d",
		   &energy_max) != 1)
		error("Could not read energy max form header\n");
	log_entry("Correctly parsed the energy max to %d",energy_max);
	next_row(&row,&row_length,header_file);
	size_t dimension = 0;
	if (sscanf(row,
		   "dim: %lu",
		   &dimension) != 1)
		error("Could not read dimension from header\n");
	log_entry("Correctly parsed the dimension to %lu",dimension);
	next_row(&row,&row_length,header_file);
	energy_block_info_t *energy_block_info = NULL;
	size_t num_energy_blocks = 0;
	array_builder_t energy_block_builder =
		new_array_builder((void**)&energy_block_info,
				  &num_energy_blocks,
				  sizeof(energy_block_info_t));
	const char *new_table_header =
		"\tE1\tE2\tTz\tM\tconf_file\telem_file\n";
	const char *old_table_header =
		"\tdN\tTz\tM\tconf_file\telem_file\n";
	if (strcmp(row,new_table_header) == 0)
	{
		while (!feof(header_file))
		{
			if (getline(&row,&row_length,header_file) < 0)
				break;
			energy_block_info_t current_energy_block;
			current_energy_block.is_old = 0;
			if (sscanf(row,
				   "\t%d\t%d\t%d\t%d\t%s\t%s",
				   &current_energy_block.E1,
				   &current_energy_block.E2,
				   &current_energy_block.Tz,
				   &current_energy_block.M,
				   current_energy_block.configuration_file,
				   current_energy_block.element_file) != 6)
				break;
			current_energy_block.dE =
			       	current_energy_block.E2-
				current_energy_block.E1;
			append_array_element(energy_block_builder,
					     &current_energy_block);

		}
	}
	else if (strcmp(row,old_table_header) == 0)
	{
		while (!feof(header_file))
		{
			if (getline(&row,&row_length,header_file) < 0)
				break;
			energy_block_info_t current_energy_block;
			current_energy_block.is_old = 1;
			if (sscanf(row,
				   "\t%d\t%d\t%d\t%s\t%s",
				   &current_energy_block.dE,
				   &current_energy_block.Tz,
				   &current_energy_block.M,
				   current_energy_block.configuration_file,
				   current_energy_block.element_file) != 5)
				break;
			append_array_element(energy_block_builder,
					     &current_energy_block);
		}
	}
	else
	{
		error("Wrong table header\n");
	}
	log_entry("The table header is correct\n");
	free_array_builder(energy_block_builder);
	if (row)
		free(row);
	fclose(header_file);
	free(file_name);
	header_t header = (header_t)malloc(sizeof(struct _header_));
	header->num_particles = num_particles;
	header->energy_max = energy_max;
	header->dimension = dimension;
	header->num_energy_blocks = num_energy_blocks;
	header->energy_block_info = energy_block_info;
	return header;

}

size_t get_num_particles(const header_t header)
{
	return header->num_particles;
}

size_t find_block_index(const header_t header,
			int E1, int E2, int Tz, int M)
{
	if (E1 > E2)
	{
		int tmp = E1;
		E1 = E2;
		E2 = tmp;
	}
	log_entry("finding block %d %d %d %d",
		  E1,E2,Tz,M);
	int dE = abs(E1-E2);
	log_entry("dE = %d",dE);
	log_entry("header->num_energy_blocks = %lu",
		  header->num_energy_blocks);
	for (size_t i = 0; i < header->num_energy_blocks; i++)
	{
		energy_block_info_t block_info =
			header->energy_block_info[i];
		log_entry("block_info ={\n"
			  ".dE = %d\n"
			  ".E1 = %d\n"
			  ".E2 = %d\n"
			  ".Tz = %d\n"
			  ".M = %d\n"
			  ".configuration_file = \"%s\"\n"
			  ".element_file = \"%s\"\n"
			  ".is_old = %d\n"
			  "}",
			  block_info.dE,
			  block_info.E1,
			  block_info.E2,
			  block_info.Tz,
			  block_info.M,
			  block_info.configuration_file,
			  block_info.element_file,
			  block_info.is_old);
		if (((block_info.is_old == 0 &&
		      block_info.E1 == E1 &&
		      block_info.E2 == E2) ||
		     (block_info.is_old == 1 &&
		      block_info.dE == dE)) &&
		    block_info.Tz == Tz &&
		    block_info.M == M)
			return i;
	}
	return no_index;
}

energy_block_info_t get_energy_block_info(header_t header,
					  size_t index)
{
	return header->energy_block_info[index];
}

void free_header(header_t header)
{
	free(header->energy_block_info);
	free(header);
}

	static
void next_row(char **row, size_t *row_length,FILE *file)
{
	if (getline(row,row_length,file) < 0)
		error("Could not fetch row from file\n");
	log_entry("Fetched row \"%s\"",row);
}
