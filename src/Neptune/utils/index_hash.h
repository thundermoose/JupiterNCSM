#ifndef __INDEX_HASH__
#define __INDEX_HASH__

#include <stdlib.h>

#define no_index -1

struct _index_hash_;
typedef struct _index_hash_ *index_hash_t;

index_hash_t new_index_hash(size_t num_bins);

size_t get_index(const index_hash_t index_hash,
		size_t left_index,
		size_t right_index);

void set_index(index_hash_t index_hash,
		size_t left_index,
		size_t right_index,
		size_t combined_index);

void free_index_hash(index_hash_t index_hash);
#endif
