#ifndef __GLOBAL_CONSTANTS__
#define __GLOBAL_CONSTANTS__

#if __GNUC__ <= 6

#define max_interaction_order 3

#define no_index (size_t)(-1)

#else

static const size_t max_interaction_order = 3;

static const size_t no_index = -1;

#endif

#endif
