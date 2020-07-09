#ifndef __PERMUTATION_TOOLS__
#define __PERMUTATION_TOOLS__

#define SWAP(a,b)				\
  {						\
    typeof(a) t = a;				\
    a = b;					\
    b = t;					\
  }

#define CYCLIC_SHIFT_R(a,b,c)			\
  {						\
    typeof(c) t = c;				\
    c = b;					\
    b = a;					\
    a = t;					\
  }

#define CYCLIC_SHIFT_L(a,b,c)			\
  {						\
    typeof(c) t = c;				\
    c = a;					\
    a = b;					\
    b = t;					\
  }

#endif
