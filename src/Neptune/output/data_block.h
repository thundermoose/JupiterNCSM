#ifndef __DATA_BLOCK__
#define __DATA_BLOCK__
#include <bases/shells.h>
typedef enum
{
	unused_block,
	used_block
} block_status_t;

typedef struct _data_block_
{
  quantum_number Tz,M,E1,E2;
} Data_Block;

#endif
