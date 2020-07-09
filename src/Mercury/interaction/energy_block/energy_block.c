#include <interaction/energy_block/energy_block.h>
#include <log/log.h>
#include <error/error.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>


typedef struct
{
	size_t bra_index;
	size_t ket_index;
} config_t;

struct _energy_block_
{
	config_t *configurations;
	double *elements; 
	size_t num_configurations;
	size_t *hash_buckets;
	size_t num_hash_buckets;
};

static
FILE *open_block_file(const char *interaction_path,
		      const char *file_name);

static
size_t compute_configuration_hash_bucket(config_t key,
					 config_t *keys,
					 size_t *hash_buckets,
					 size_t num_hash_buckets);

static
int compare_config(config_t a, config_t b);

static
uint64_t configuration_hash(config_t key);

energy_block_t open_energy_block(const char *interaction_path,
				 energy_block_info_t info)
{
	FILE* configuration_file =
		open_block_file(interaction_path,
				info.configuration_file);	   
	fseek(configuration_file,0,SEEK_END);
	size_t num_bytes = ftell(configuration_file);
	fseek(configuration_file,0,SEEK_SET);
	if (num_bytes % sizeof(config_t) != 0)
		error("The configuration file %s is corrupt\n",
		      info.configuration_file);
	size_t num_configurations = num_bytes/sizeof(config_t);
	config_t *configurations = (config_t*)malloc(num_bytes);
	if (fread(configurations,
		  sizeof(config_t),
		  num_configurations,
		  configuration_file) < num_configurations)
		error("Could not read configurations from %s\n",
		      info.configuration_file);
	fclose(configuration_file);
	FILE *element_file =
		open_block_file(interaction_path,
				info.element_file);
	fseek(element_file,0,SEEK_END);
	num_bytes = ftell(element_file);
	fseek(element_file,0,SEEK_SET);
	if (num_bytes % sizeof(double) != 0)
		error("The element file %s is corrupt\n",
		      info.element_file);
	size_t num_elements = num_bytes/sizeof(double);
	if (num_elements != num_configurations)
		error("The element file %s do not have as many"
		      " elements as %s has configurations\n",
		      info.element_file,info.element_file);
	double *elements = (double*)malloc(num_bytes);
	if (fread(elements,
		  sizeof(double),
		  num_elements,
		  element_file) < num_elements)
		error("Could not read elements for the %s file\n",
		      info.element_file);
	fclose(element_file);

	size_t num_hash_buckets = 2*num_configurations;
	size_t *hash_buckets = (size_t*)calloc(num_hash_buckets,
					       sizeof(size_t));
	for (size_t i = 0; i<num_configurations; i++)
	{
		config_t key = configurations[i];
		if (key.bra_index > key.ket_index)
		{
			size_t tmp = key.bra_index;
			key.bra_index = key.ket_index;
			key.ket_index = tmp;
		}
		log_entry("key[%lu] = (%lu %lu)",
			  i,
			  key.bra_index,
			  key.ket_index);
		size_t bucket_index =
			compute_configuration_hash_bucket(key,
							  configurations,
							  hash_buckets,
							  num_hash_buckets);
		log_entry("bucket_index = %lu",bucket_index);
		log_entry("elements[%lu] = %lg",i,elements[i]);
		hash_buckets[bucket_index] = i+1;
	}
	energy_block_t energy_block =
		(energy_block_t)malloc(sizeof(struct _energy_block_));
	energy_block->configurations = configurations;
	energy_block->elements = elements;
	energy_block->num_configurations = num_configurations;
	energy_block->hash_buckets = hash_buckets;
	energy_block->num_hash_buckets = num_hash_buckets;
	return energy_block;
}

double get_energy_block_element(energy_block_t energy_block,
				size_t bra_index,
				size_t ket_index)
{
	config_t key = {bra_index,ket_index};
	size_t index = 
		compute_configuration_hash_bucket
		(key,
		 energy_block->configurations,
		 energy_block->hash_buckets,
		 energy_block->num_hash_buckets);
	log_entry("index = %lu",index);
	if (energy_block->hash_buckets[index] == 0)
	{
		if (bra_index > ket_index)
		{
			key.bra_index = ket_index;
			key.ket_index = bra_index;
		}
		log_entry("key = (%lu %lu)",
			  key.bra_index,
			  key.ket_index);
		index = compute_configuration_hash_bucket
			(key,
			 energy_block->configurations,
			 energy_block->hash_buckets,
			 energy_block->num_hash_buckets);
		log_entry("index = %lu",index);
		if (energy_block->hash_buckets[index] == 0)
			return 0.0;
	}
	return energy_block->elements[energy_block->hash_buckets[index]-1];

}

void free_energy_block(energy_block_t energy_block)
{
	free(energy_block->configurations);
	free(energy_block->elements);
	free(energy_block->hash_buckets);
	free(energy_block);
}

static
FILE *open_block_file(const char *interaction_path,
		      const char *file_name)
{
	const size_t full_file_name_length = 
		strlen(interaction_path) + strlen(file_name)+2;
	char *full_file_name = (char*)malloc(full_file_name_length);
	sprintf(full_file_name,
		"%s/%s",
		interaction_path,
		file_name);
	FILE *file = fopen(full_file_name,"r");
	if (file == NULL)
		error("Could not open file %s. %s\n",
		      full_file_name,
		      strerror(errno));
	return file;
}

static
size_t compute_configuration_hash_bucket(config_t key,
					 config_t *keys,
					 size_t *hash_buckets,
					 size_t num_hash_buckets)
{
	log_entry("key = (%lu %lu)",key.bra_index,key.ket_index);
	uint64_t hash = configuration_hash(key);
	log_entry("hash = %lx",hash);
	size_t index = hash % num_hash_buckets;
	while (hash_buckets[index] != 0 &&
	       !compare_config(key,keys[hash_buckets[index]-1]))
	{
		log_entry("keys[%lu] = (%lu %lu)",
			  hash_buckets[index]-1,
			  keys[hash_buckets[index]-1].bra_index,
			  keys[hash_buckets[index]-1].ket_index);
		hash++;
		index = hash % num_hash_buckets;
	}
	return index;	
}

static
uint64_t configuration_hash(config_t key)
{
	if (key.bra_index < key.ket_index)
		return key.bra_index ^ __builtin_bswap64(key.ket_index);
	else
		return key.ket_index ^ __builtin_bswap64(key.bra_index);
}

static
int compare_config(config_t a, config_t b)
{
	log_entry("compares (%lu %lu) and (%lu %lu)",
		  a.bra_index,a.ket_index,
		  b.bra_index,b.ket_index);
	return (a.ket_index == b.ket_index &&
		a.bra_index == b.bra_index) ||
		(a.ket_index == b.bra_index &&
		 a.bra_index == b.ket_index);
}
