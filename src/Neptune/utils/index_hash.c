#include <utils/index_hash.h>
#include <utils/bit_operations.h>
#include <stdint.h>
#include <unit_testing/test.h>


typedef struct
{
	size_t left_index;
	size_t right_index;
} index_key_t;

typedef struct
{
	index_key_t key;
	size_t combined_index;
} bin_t;

struct _index_hash_
{
	size_t num_bins;
	bin_t *bins;
};

const bin_t empty_bin = 
{
	.key = 
	{
		.left_index = no_index,
		.right_index = no_index
	},
	.combined_index = no_index
};

static inline
uint64_t compute_hash(index_key_t key);

static inline
int has_right_key(bin_t bin, index_key_t key);

static inline
int is_empty(bin_t bin);

index_hash_t new_index_hash(size_t num_bins)
{
	index_hash_t index_hash =
		(index_hash_t)malloc(sizeof(struct _index_hash_));
	index_hash->num_bins = num_bins;
	index_hash->bins = (bin_t*)malloc(num_bins*sizeof(bin_t));
	for (size_t i = 0; i<num_bins; i++)
		index_hash->bins[i] = empty_bin;
	return index_hash;
}

size_t get_index(const index_hash_t index_hash,
		size_t left_index,
		size_t right_index)
{
	index_key_t key = {left_index,right_index};
	uint64_t hash = compute_hash(key);
	while (!has_right_key(index_hash->bins[hash % index_hash->num_bins],
				key))
		hash++;
	return index_hash->bins[hash % index_hash->num_bins].combined_index;
}

void set_index(index_hash_t index_hash,
		size_t left_index,
		size_t right_index,
		size_t combined_index)
{
	index_key_t key = {left_index,right_index};
	uint64_t hash = compute_hash(key);
	while (!is_empty(index_hash->bins[hash % index_hash->num_bins]))
		hash++;
	index_hash->bins[hash % index_hash->num_bins].key = key;
	index_hash->bins[hash % index_hash->num_bins].combined_index =
		combined_index;
}

void free_index_hash(index_hash_t index_hash)
{
	free(index_hash->bins);
	free(index_hash);
}

static inline
uint64_t compute_hash(index_key_t key)
{
	uint64_t h1 = key.left_index ^ key.right_index;
#define magic_number 0xFEDCBA9876543210
	uint64_t h2 = key.left_index ^ magic_number;
	uint64_t h3 = key.right_index ^ magic_number;
	return h1^(cyclic_rshift(h2,17))^(cyclic_lshift(h3,13));
}

static inline
int has_right_key(bin_t bin, index_key_t key)
{
	return (bin.key.left_index == key.left_index &&
		bin.key.right_index == key.right_index) ||
		bin.key.left_index == no_index ||
		bin.key.right_index == no_index ||
		bin.combined_index == no_index;
}

static inline
int is_empty(bin_t bin)
{
	return bin.combined_index == no_index &&
		bin.key.left_index == no_index &&
		bin.key.right_index == no_index;
}

int generate_and_retrive_works()
{
	size_t width = 100, height = 100;
	index_hash_t hash_map = new_index_hash(width*height*2);
	// just create an index that is chaotic but deterministic
	// to test against
#define compute_combined_index(i,j) (((i<<3)^j)%(width*height))
	for (size_t i = 0; i < width; i++)
		for (size_t j = 0; j<height; j++)
		{
			size_t combined_index = compute_combined_index(i,j);
			set_index(hash_map,i,j,combined_index);
		}
	for (size_t i = 0; i < width; i++)
		for (size_t j = 0; j<height; j++)
		{
			size_t combined_index = compute_combined_index(i,j);
			if (combined_index != get_index(hash_map,i,j))
				return 0;
		}
	return 1;
}

new_test(generate_and_retrived_test,
	 assert_that(generate_and_retrive_works());
	);
