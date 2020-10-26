#include <string_tools/string_tools.h>
#include <unit_testing/test.h>
#include <stdlib.h>
#include <string.h>
#include <debug_mode/debug_mode.h>

static
size_t count_words(const char *string);

static
const char *find_next_word(const char *string);

static
size_t get_word_length(const char *string);

static
int white_character(char c);

char *copy_string(const char *origin)
{
	const size_t length = strlen(origin);
	char *target = (char*)calloc(length+1,sizeof(char));
	memcpy(target,origin,length);
	return target;
}

size_t extract_words(char ***words,
		     const char *string)
{
	size_t num_words = count_words(string);
	*words = (char**)malloc(num_words*sizeof(char*));
	for (size_t i = 0; i<num_words; i++)
	{
		string = find_next_word(string);
		size_t word_length = get_word_length(string);
		(*words)[i] = (char*)calloc(word_length+1,
					    sizeof(char));
		memcpy((*words)[i],
		       string,
		       word_length);
		string+=word_length;
	}
	return num_words;
}

int begins_with(const char* str,
		const char* beg)
{
  char *beg_end = (char*)beg;
  while (*(++beg_end) != 0);

  char* s = (char*)str;
  char* b = (char*)beg;
  while (b!=beg_end)
    {
      if (*(s++) != *(b++))
	return 0; 
    }
  return 1;
}

int ends_with(const char* str,
	      const char* end)
{
  char* str_end = (char*)str;
  while (*(++str_end) != 0);

  char* end_end = (char*)end;
  while (*(++end_end) != 0);
  
  while (end_end>end)
    {
      if (*(end_end--) != *(str_end--))
	{
	  return 0;
	}
    }
  return 1;
}

int is_integer(const char *string)
{
	if (*string == '-')
		string++;
	do
	{
		char current = *(string++);
		if (current<'0' || current>'9')
			return 0;
	}
	while (*string != '\0');
	return 1;
}

char *concatinate_strings(const char *string_1,
			  const char *string_2)
{
	const size_t length_string_1 = strlen(string_1);
	const size_t length_string_2 = strlen(string_2);
	const size_t length_string = length_string_1 + length_string_2;
	char *string =
	       	(char*)malloc(length_string + 1);
	memcpy(string,
	       string_1,
	       length_string_1);
	memcpy(string+length_string_1,
	       string_2,
	       length_string_2);
	string[length_string] = 0;
	return string;
}

char *jump_back_word(char *string, char *current_word)
{
	current_word--;
	while (current_word > string &&
	       white_character(*current_word))
		current_word--;	
	while (current_word > string &&
	       !white_character(*current_word))
		current_word--;	
	return current_word++;
}

static
size_t count_words(const char *string)
{
	size_t num_words = !white_character(*string) ? 1 : 0;
	while (*(++string) != 0)
	{
		if (!white_character(*string) && white_character(*(string-1)))
			num_words++;
	}
	return num_words;
}

static
const char *find_next_word(const char *string)
{
	while (white_character(*string))
	{
		string++;
	}
	return string;
}

static
size_t get_word_length(const char *string)
{
	size_t word_length = 0;
	while (!white_character(*string))
	{
		word_length++;
		string++;
	}
	return word_length;
}

static
int white_character(char c)
{
	return c == ' ' || c == '\n' || c == '\t';
}

new_test(extracting_words_from_simple_sentence,
	 const char *message = "This is a test";
	 const char *expected_words[4] =
	 {
	 	"This",
		"is",
		"a",
		"test"
	 };
	 char **words = NULL;
	 size_t num_words = extract_words(&words,message);
	 printf("There are %lu words and expects 4\n",
		num_words);
	 for (size_t i = 0; i<num_words; i++)
	 	printf("(%lu): %s\n",
		       i,words[i]);
	 assert_that(num_words == 4);
	for (size_t i = 0; i<4; i++)
	{
		assert_that(strcmp(words[i],expected_words[i]) == 0);
	}	
	for (size_t i = 0; i<4; i++)
		free(words[i]);
	free(words);
	);

new_test(extracting_words_special_case,
	 const char *message = "n      0   0    0   0    0           2  # ARRAY:     3=        24\n";
	 const char *expected_words[11] =
	 {
	 "n",
	 "0",
	 "0",
	 "0",
	 "0",
	 "0",
	 "2",
	 "#",
	 "ARRAY:",
	 "3=",
	 "24"
	 };
	 char **words = NULL;
	 size_t num_words = extract_words(&words,message);
	 printf("There are %lu words and expected 11\n",
		num_words);
	 for (size_t i = 0; i<num_words; i++)
		printf("(%lu): %s\n",
		       i,words[i]);
	assert_that(num_words == 11);
	for (size_t i = 0; i<num_words; i++)
		assert_that(strcmp(words[i],expected_words[i]) == 0);
	for (size_t i = 0; i<num_words; i++)
		free(words[i]);
	free(words);
	);
