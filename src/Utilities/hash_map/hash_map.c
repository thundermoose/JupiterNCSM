#include <hash_map/hash_map.h>
#include <debug_mode/debug_mode.h>
#include <string.h>

struct _hash_map_
{
	size_t num_required_elements;
       	size_t margin_size;
	size_t element_size;	
	size_t key_size;
	__compar_fn_t key_compare_function;
	__hash_fn_t hash_function;
	size_t num_buckets;
	char *element_buckets;
	char *key_buckets;
};

static inline
size_t find_bucket(hash_map_t hash_map,
		   const void *key);

static
int is_zero_key(const char* key, size_t key_size);

hash_map_t new_hash_map(size_t num_required_elements,
			size_t margin_size,
			size_t element_size,
			size_t key_size,
			__compar_fn_t key_compare_function,
			__hash_fn_t hash_function)
{
	hash_map_t hash_map = (hash_map_t)malloc(sizeof(struct _hash_map_));
	hash_map->num_required_elements = num_required_elements;
	hash_map->margin_size = margin_size;
	hash_map->element_size = element_size;
	hash_map->key_size = key_size;
	hash_map->key_compare_function = key_compare_function;
	hash_map->hash_function = hash_function;
	hash_map->num_buckets = num_required_elements*margin_size;
	hash_map->element_buckets =
	       	(char*)malloc(hash_map->num_buckets*element_size);
	hash_map->key_buckets =
		(char*)calloc(hash_map->num_buckets,key_size);
	return hash_map;
}

void hash_map_insert(hash_map_t hash_map,
		     const void *key,
		     const void *element)
{
	size_t bucket_index = find_bucket(hash_map,key);
	memcpy(hash_map->key_buckets + bucket_index*hash_map->key_size,
	       key,
	       hash_map->key_size);
	memcpy(hash_map->element_buckets + bucket_index*hash_map->element_size,
	       element,
	       hash_map->element_size);
}

int hash_map_get(const hash_map_t hash_map,
		 const void *key,
		 void *element)
{
	size_t bucket_index = find_bucket(hash_map,key);
	if (hash_map->key_compare_function
	    (key, 
	     (void*)hash_map->key_buckets+
	     bucket_index*hash_map->key_size) == 0)
	{
		memcpy(element,
		       hash_map->element_buckets + 
		       bucket_index*hash_map->element_size,
		       hash_map->element_size);
		return 1;

	}
	else
	{
		return 0;
	}
}

void free_hash_map(hash_map_t hash_map)
{
	free(hash_map->element_buckets);
	free(hash_map->key_buckets);
	free(hash_map);
}

static inline
size_t find_bucket(hash_map_t hash_map,
		   const void *key)
{
	uint64_t hash = hash_map->hash_function(key);
	size_t index = hash % hash_map->num_buckets;
	__compar_fn_t compare = hash_map->key_compare_function;
	void *current_key = 
		(void*)(hash_map->key_buckets + index * hash_map->key_size);
	while (!is_zero_key(current_key,hash_map->key_size) &&
	       compare(key,current_key) != 0)
	{
		hash++;
		index = hash % hash_map->num_buckets;
		current_key =
			(void*)(hash_map->key_buckets + 
			index * hash_map->key_size);
	}
	return index;
}

static
int is_zero_key(const char* key, size_t key_size)
{
	for (size_t i = 0; i < key_size; i++)
		if (key[i] != 0)
			return 0;
	return 1;
}
