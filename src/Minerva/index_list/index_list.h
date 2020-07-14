#ifndef __INDEX_LIST__
#define __INDEX_LIST__

#include <sub_basis_block/sub_basis_block.h>
#include <index_triple/index_triple.h>

struct _index_list_;
typedef struct _index_list_ *index_list_t;

index_list_t new_index_list(const char *base_directory,
			    const sub_basis_block_t in_block,
			    const sub_basis_block_t out_block,
			    const int sign);

index_list_t new_index_list_from_id(const char *base_directory,
				    const size_t id);

size_t length_index_list(const index_list_t index_list);

index_triple_t *get_index_list_elements(const index_list_t index_list);

void free_index_list(index_list_t index_list);

#endif
