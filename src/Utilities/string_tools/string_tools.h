#ifndef __STRING_TOOLS__
#define __STRING_TOOLS__

#include <stdlib.h>

/* Allocates memory and copies "origin" to
 * that new memory
 */
char *copy_string(const char *origin);

/* Extract all words, ie, clusters of characters surounded by atleast one 
 * whitespace character ('\n', '\t' or ' '), from "string". The words are 
 * stored in an array that is returned through "words". The return value
 * is the number of extracted words.
 */
size_t extract_words(char ***words,
		     const char *string);

int begins_with(const char* str,
		const char* beg);

int ends_with(const char* str,
	      const char* end);

int is_integer(const char *string);

char *concatinate_strings(const char *string_1, const char *string_2);

#endif
