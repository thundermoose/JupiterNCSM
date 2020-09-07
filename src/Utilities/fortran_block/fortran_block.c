#include <fortran_block/fortran_block.h>
#include <stdlib.h>
#include <debug_mode/debug_mode.h>

struct _fortran_block_
{
	int num_bytes;
	char *data;	
};

fortran_block_t read_fortran_block(FILE *file)
{
	int block_size = 0;
	if (fread(&block_size,
		  sizeof(int),
		  1,file) != 1)
	{
		fprintf(stderr, "Could not beginning of block\n");
		exit(EXIT_FAILURE);
	}
	char *data = (char*)malloc(block_size);
	if (fread(data,
		  sizeof(char),
		  block_size,
		  file) != block_size)
	{
		fprintf(stderr, "Could not read data of block of size %dB\n",
			block_size);
		exit(EXIT_FAILURE);
	}
	int block_end = 0;
	if (fread(&block_end,
		  sizeof(int),
		  1,file) != 1)
	{
		fprintf(stderr, "Could not read end of block of size %dB\n",
			block_size);
		exit(EXIT_FAILURE);
	}
	if (block_size != block_end)
	{
		fprintf(stderr, "This is not a fortran block\n");
		exit(EXIT_FAILURE);
	}
	fortran_block_t block =
	       	(fortran_block_t)malloc(sizeof(struct _fortran_block_));
	block->num_bytes = block_size;
	block->data = data;
	return block;
}

int get_num_bytes(fortran_block_t block)
{
	return block->num_bytes;
}

const char *get_fortran_block_data(fortran_block_t block)
{
	return block->data;
}

void free_fortran_block(fortran_block_t block)
{
	free(block->data);
	free(block);
}

void skip_fortran_block(FILE *file)
{
	int block_size = 0;
	if (fread(&block_size,
		  sizeof(int),
		  1,file) != 1)
	{
		fprintf(stderr, "Could not beginning of block\n");
		exit(EXIT_FAILURE);
	}
	long current_possition = ftell(file);
	fseek(file,current_possition+block_size,SEEK_SET);
	int block_end = 0;
	if (fread(&block_end,
		  sizeof(int),
		  1,file) != 1)
	{
		fprintf(stderr, "Could not read end of block of size %dB\n",
			block_size);
		exit(EXIT_FAILURE);
	}
	if (block_size != block_end)
	{
		fprintf(stderr, "This is not a fortran block\n");
		exit(EXIT_FAILURE);
	}
}
