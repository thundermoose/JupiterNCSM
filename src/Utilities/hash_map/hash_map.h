#ifndef __HASH_MAP__
#define __HASH_MAP__

#include <stdlib.h>
#include <stdint.h>

struct _hash_map_;
typedef struct _hash_map_ *hash_map_t;

typedef uint64_t (*__hash_fn_t)(const void *key);

hash_map_t new_hash_map(size_t num_required_elements,
			size_t margin_size,
			size_t element_size,
			size_t key_size,
			__compar_fn_t key_compare_function,
			__hash_fn_t hash_function);

void hash_map_insert(hash_map_t hash_map,
		     const void *key,
		     const void *element);

int hash_map_get(const hash_map_t hash_map,
		 const void *key,
		 void *element);

void free_hash_map(hash_map_t hash_map);

#endif
