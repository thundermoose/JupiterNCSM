#ifndef __DEBUG_MESSAGES__
#define __DEBUG_MESSAGES__

#include "terminal_colors.h"

#ifdef DEBUG
#ifdef COLOR

#define DEBUG_MESS(message...)			\
  fprintf(stderr,YELLOW"DM %s@%s:%d :",	\
	  __builtin_FUNCTION(),			\
	  __builtin_FILE(),			\
	  __builtin_LINE());			\
  fprintf(stderr,GREEN message);		\
  fprintf(stderr,NORMAL);

#define DEBUG_PRINT_MATRIX(matrix)\
{\
  fprintf(stderr,YELLOW"DM %s@%s:%d : matrix:\n",	\
	  __builtin_FUNCTION(),			\
	  __builtin_FILE(),			\
	  __builtin_LINE());			\
	for (size_t i = 0; i<matrix->m; i++)\
	{\
		for (size_t j = 0; j<matrix->n; j++)\
			fprintf(stderr,GREEN "%lg ",matrix->elements[i*matrix->n+j]);\
		fprintf(stderr,"\n");\
	}\
	fprintf(stderr,NORMAL "\n");\
}
#else

#define DEBUG_MESS(message...)			\
  fprintf(stderr,"DM %s@%s:%d :",	\
	  __builtin_FUNCTION(),			\
	  __builtin_FILE(),			\
	  __builtin_LINE());			\
  fprintf(stderr,message);		\
 
#define DEBUG_PRINT_MATRIX(matrix)\
{\
  fprintf(stderr,"DM %s@%s:%d : matrix:\n",	\
	  __builtin_FUNCTION(),			\
	  __builtin_FILE(),			\
	  __builtin_LINE());			\
	for (size_t i = 0; i<matrix->m; i++)\
	{\
		fprintf(stderr,"DM %s@%s:%d: ",\
			__builtin_FUNCTION(),\
			__builtin_FILE(),\
			__builtin_LINE());\
		for (size_t j = 0; j<matrix->n; j++)\
			fprintf(stderr,"%lg ",matrix->elements[i*matrix->n+j]);\
		fprintf(stderr,"\n");\
	}\
	fprintf(stderr,"\n");\
}
	

#endif

#define DEBUG_CALL(action...) action

#else

#define DEBUG_MESS(message...)

#define DEBUG_CALL(action...)

#define DEBUG_PRINT_MATRIX(matrix)\

#endif

#endif
