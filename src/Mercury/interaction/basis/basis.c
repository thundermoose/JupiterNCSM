#include <interaction/basis/basis.h>
#include <log/log.h>
#include <error/error.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdint.h>

struct _basis_
{
	size_t num_states;
	int *states;
	size_t *hash_buckets;
	size_t num_hash_buckets;
};

static
size_t find_basis_hash_bucket(int *state,
			      int *states,
			      size_t *hash_buckets,
			      size_t num_hash_buckets,
			      size_t num_particles);

static
uint64_t basis_state_hash(int *state, size_t num_particles);

static
int compar_states(int *state_one,
		  int *state_two,
		  size_t num_particles);


basis_t read_basis(const char *interaction_path,
		   size_t num_particles)
{
	size_t file_name_length = strlen(interaction_path)+11;
	char *file_name = (char*)malloc(file_name_length);
	sprintf(file_name,
		"%s/basis.bin",
		interaction_path);
	FILE *basis_file = fopen(file_name,"r");
	if (basis_file == NULL)
		error("Could not open %s. %s\n",
		      file_name,
		      strerror(errno));
	fseek(basis_file,0,SEEK_END);
	size_t num_bytes = ftell(basis_file);
	fseek(basis_file,0,SEEK_SET);
	if (num_bytes % sizeof(int) != 0)
		error("%s do not contain integers\n",
		      file_name);
	size_t num_ints = num_bytes / sizeof(int);
	if (num_ints % num_particles != 0)
		error("%s do not conform to the set number"
		      " of particles in the header\n",
		      file_name);
	size_t num_states = num_ints/num_particles;
	size_t state_length = sizeof(int)*num_particles;	
	int *states = (int*)malloc(num_bytes);
	if (fread(states,state_length,num_states,basis_file) < num_states)
		error("Could not read states from %s\n",
		      file_name);
	fclose(basis_file);
	free(file_name);
	size_t num_hash_buckets = num_states*2;
	size_t *hash_buckets =
		(size_t*)calloc(num_hash_buckets,sizeof(size_t));
	for (size_t i = 0; i < num_states; i++)
	{
		int *state = states+i*num_particles;
		size_t index =find_basis_hash_bucket(state,
						     states,
						     hash_buckets,
						     num_hash_buckets,
						     num_particles) ;
		hash_buckets[index] = i+1;
		if (num_particles == 1)
		{
			log_entry("state = %d",state[0]);
		}
		else if (num_particles == 2)
		{
			log_entry("state = %d %d", state[0], state[1]);
		}
		else if (num_particles == 3)
		{
			log_entry("state = %d %d %d",
				  state[0], state[1], state[3]);
		}
		log_entry("index = %lu",index);
		log_entry("i = %lu",i);
	}
	basis_t basis = (basis_t)malloc(sizeof(struct _basis_));
	basis->num_states = num_states;
	basis->states = states;
	basis->hash_buckets = hash_buckets;
	basis->num_hash_buckets = num_hash_buckets;
	return basis;	
}

size_t find_basis_state(basis_t basis,
		  int *state,
		  size_t num_particles)
{
	if (num_particles == 1)
	{
		log_entry("state = %d",state[0]);
	}
	else if (num_particles == 2)
	{
		log_entry("state = %d %d", state[0], state[1]);
	}
	else if (num_particles == 3)
	{
		log_entry("state = %d %d %d",
			  state[0], state[1], state[3]);
	}
	size_t index = find_basis_hash_bucket(state,
					      basis->states,
					      basis->hash_buckets,
					      basis->num_hash_buckets,
					      num_particles);
	log_entry("index = %lu",index);
	log_entry("i = %lu",basis->hash_buckets[index]-1);
	return basis->hash_buckets[index]-1;
}

void free_basis(basis_t basis)
{
	free(basis->states);
	free(basis->hash_buckets);
	free(basis);
}

static
size_t find_basis_hash_bucket(int *state,
			      int *states,
			      size_t *hash_buckets,
			      size_t num_hash_buckets,
			      size_t num_particles)
{
	uint64_t hash = basis_state_hash(state,num_particles);
	if (num_particles == 1)
	{
		log_entry("state = %d",state[0]);
	}
	else if (num_particles == 2)
	{
		log_entry("state = %d %d", state[0], state[1]);
	}
	else if (num_particles == 3)
	{
		log_entry("state = %d %d %d",
			  state[0], state[1], state[3]);
	}
	log_entry("hash = %lx",hash);
	size_t index = hash % num_hash_buckets;
	while (hash_buckets[index] != 0 && 
	       compar_states(state,
			     states+
			     num_particles*(hash_buckets[index]-1),
			     num_particles) != 0)
	{
		hash++;
		index = hash % num_hash_buckets;
	}
	return index;
}
static
uint64_t basis_state_hash(int *state, size_t num_particles)
{
	uint64_t hash = 0;
	for (size_t i = 0; i<num_particles; i++)
	{
		uint64_t s = (uint64_t)(state[i]);
		hash^=s<<((13*i)%64)|(s>>(64-(13*i)%64));
	}
	return hash;
}

static
int compar_states(int *state_one,
		  int *state_two,
		  size_t num_particles)
{
	int diff = 0;
	for (size_t i = 0; i<num_particles; i++)
	{
		diff = state_one[i]-state_two[i];
		if (diff)
			return diff;
	}
	return 0;
}

