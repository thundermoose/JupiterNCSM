#include "read_packed_states.h"
#include <error/error.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <debug_mode/debug_mode.h>

void read_packed_states(void **states,
			size_t *num_states,
			size_t state_size,
			const char *filename)
{
	FILE *basis_file = fopen(filename,"r");
	if (basis_file == NULL)
		error("Could not open file %s. %s\n",
		      filename,
		      strerror(errno));
	fseek(basis_file,0,SEEK_END);
	const size_t num_bytes = ftell(basis_file);
	fseek(basis_file,0,SEEK_SET);
	if (num_bytes % state_size != 0)
		error("Basis file %s do not include %lu bytes states.\n",
		      filename,
		      state_size);
	*num_states = num_bytes/state_size;
	*states = (unsigned int*)malloc(num_bytes);
	if (fread(*states,
		  state_size,
		  *num_states,
		  basis_file) != (*num_states))
		error("Could not read %lu bytes from %s.\n",
		      num_bytes,filename);
	fclose(basis_file);	
}

