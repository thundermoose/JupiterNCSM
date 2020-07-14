#include <index_list/index_list.h>
#include <array_builder/array_builder.h>
#include <error/error.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

typedef struct
{
	int i,j,k;
} index_triple_t;

struct _index_list_
{
	index_triple_t *indices;
	size_t num_indices;	
};

index_list_t parse_human_readable_index_list(const char *file_name)
{
	index_list_t index_list =
	       	(index_list_t)calloc(1,sizeof(struct _index_list_));
	array_builder_t indices_builder =
		new_array_builder((void**)&index_list->indices,
				  &index_list->num_indices,
				  sizeof(index_triple_t));
	FILE *file = fopen(file_name,"r");
	if (file == NULL)
		error("Could not open file %s. %s\n",
		      file_name,
		      strerror(errno));
	char *row = NULL;
	size_t row_length = 0;
	size_t row_index = 0;
	while (!feof(file))
	{
		if (getline(&row,&row_length,file)<0)
			break;
		row_index++;
		index_triple_t current_triple;
		char sign_char;
		if (sscanf(row,"%d %d %c%d",
			   &current_triple.i,
			   &current_triple.j,
			   &sign_char,
			   &current_triple.k) != 4)
			error("Could not parse row %lu in file %s\n",
				row_index,file_name);
		assert(current_triple.k>=0);
		if (sign_char == '-')
			current_triple.k |= 0x80000000;
		append_array_element(indices_builder,
				     &current_triple);
	}
	fclose(file);
	free_array_builder(indices_builder);
	if (row != NULL)
		free(row);
	return index_list;	
}

void save_index_list(index_list_t index_list,
		     const char *file_name)
{
	FILE *file = fopen(file_name,"w");
	if (file == NULL)
		error("Could not open file %s. %s\n",
	      		file_name,
	  		strerror(errno));		
	if (fwrite(index_list->indices,
		   sizeof(index_triple_t),
		   index_list->num_indices,
		   file)<index_list->num_indices)
		error("Could not write indices to file %s\n",
		      file_name);
	fclose(file);
}

void free_index_list(index_list_t index_list)
{
	free(index_list->indices);
	free(index_list);
}
