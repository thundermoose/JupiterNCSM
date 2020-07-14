#include <index_list/index_list.h>
#include <array_builder/array_builder.h>
#include <log/log.h>
#include <error/error.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

struct _index_list_
{
	size_t num_elements;
	index_triple_t *elements;
};

index_list_t new_index_list(const char *base_directory,
			    const sub_basis_block_t in_block,
			    const sub_basis_block_t out_block,
			    const int sign)
{
	assert(in_block.depth == out_block.depth);
	index_list_t index_list =
		(index_list_t)calloc(1,sizeof(struct _index_list_));
	array_builder_t index_list_builder =
		new_array_builder((void**)&index_list->elements,
				  &index_list->num_elements,
				  sizeof(index_triple_t));
	char index_list_file_name[256];
	sprintf(index_list_file_name,
		"%s/index_list_E_in%d_E_out%d_M_in%d_M_out%d_dE%d_dM%d_depth%d_%s",
		base_directory,
		in_block.E,
		out_block.E,
		in_block.M,
		out_block.M,
		out_block.E-in_block.E,
		out_block.M-in_block.M,
		in_block.depth,
		sign > 0 ? "pos" : "neg");
	FILE *index_list_file = fopen(index_list_file_name,"r");
	if (index_list_file == NULL)
		error("Could not open file %s. %s\n",
		      index_list_file_name,
		      strerror(errno));
	char *row = NULL;
	size_t row_length = 0;
	while (!feof(index_list_file))
	{
		if (getline(&row,&row_length,index_list_file) < 0)
			break;
		index_triple_t current_index_triple;
		if (sscanf(row,
			   "%d %d %d",
			   &current_index_triple.out_index,
			   &current_index_triple.in_index,
			   &current_index_triple.matrix_index) != 3)
			error("Could not parse index list row \"%s\"\n",
			      row);
		append_array_element(index_list_builder,&current_index_triple);
	}
	free_array_builder(index_list_builder);
	fclose(index_list_file);
	if (row != NULL)
		free(row);
	return index_list;
}

index_list_t new_index_list_from_id(const char *base_directory,
				    const size_t id)
{
	char index_list_file_name[2048];
	sprintf(index_list_file_name,
		"%s/index_list_%lu",
		base_directory,id);
	FILE *index_list_file = fopen(index_list_file_name,"r");
	if (index_list_file == NULL)
		error("Could not open file %s. %s\n",
		      index_list_file_name,
		      strerror(errno));
	fseek(index_list_file,0,SEEK_END);
	size_t num_bytes_in_file = ftell(index_list_file);	
	fseek(index_list_file,0,SEEK_SET);
	assert(num_bytes_in_file % sizeof(index_triple_t) == 0);
	index_list_t index_list =
		(index_list_t)malloc(sizeof(struct _index_list_));
	index_list->num_elements = num_bytes_in_file / sizeof(index_triple_t);
	index_list->elements = (index_triple_t*)malloc(num_bytes_in_file);		
	if (fread(index_list->elements,
		  sizeof(index_triple_t),
		  index_list->num_elements,
		  index_list_file) < index_list->num_elements)
		error("Could not read the index_list elements from %s\n",
		      index_list_file_name);
	fclose(index_list_file);
	return index_list;
}

size_t length_index_list(const index_list_t index_list)
{
	return index_list->num_elements;
}

index_triple_t *get_index_list_elements(const index_list_t index_list)
{
	return index_list->elements;
}

void free_index_list(index_list_t index_list)
{
	log_entry("free_index_list(%p)",index_list);
	free(index_list->elements);
	free(index_list);
}
