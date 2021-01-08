#ifndef __INDEX_LIST__
#define __INDEX_LIST__

struct _index_list_;
typedef struct _index_list_ *index_list_t;

index_list_t parse_human_readable_index_list(const char *file_name);

index_list_t parse_binary_index_lists(const char *file_name);

void save_index_list(index_list_t index_list,
		     const char *file_name);

void free_index_list(index_list_t index_list);

#endif
