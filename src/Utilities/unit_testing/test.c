/*
*************** License ***************
thunder-tester is a simple system for in code unit testsing in C using gcc
Copyright (C) 2019  Tor Djärv email: tordjarv@gmail.com
This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published b
    the Free Software Foundation, either version 3 of the License, o
    (at your option) any later version

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty o
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
GNU General Public License for more details
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "test.h"
#include "memory.h"
#include <directory_tools/directory_tools.h>
#include <time_tools/time_tools.h>
#include <string.h>

#define MULTIPLE_CALL_GUARD				\
	{						\
		static int has_already_been_called = 0;	\
		if (has_already_been_called)		\
			return;				\
		has_already_been_called = 1;		\
	}
typedef enum
{
	NORMAL_MODE,
	LIST_TESTS,
	RUN_SPECIFIED_TESTS,
	RUN_ALL_TESTS
} test_mode_t;

typedef struct
{
	char name[256];
	char filename[256];
} test_name_t;

test_mode_t test_mode = NORMAL_MODE;

test_name_t *tests_to_run = NULL;
size_t num_tests_to_run = 0;

char current_test_directory[2048];

char **test_file_paths = NULL;
size_t num_test_file_paths = 0;
size_t allocated_num_test_file_paths = 0;

static
void push_file_path(char *file_path)
{
	if (num_test_file_paths == allocated_num_test_file_paths)
	{
		allocated_num_test_file_paths+=allocated_num_test_file_paths+1;
		test_file_paths = REALLOC(test_file_paths,
					  allocated_num_test_file_paths);
	}
	test_file_paths[num_test_file_paths++] = file_path;
}

const char* get_test_file_path(const char* file_name)
{
	size_t buffer_length =
	       	strlen(current_test_directory) + strlen(file_name) + 1;
	char *buffer = (char*)calloc(buffer_length,sizeof(char));
	sprintf(buffer,"%s%s",current_test_directory,file_name);
	push_file_path(buffer);
	return buffer;
}

static
test_mode_t determine_test_mode(const size_t num_arguments,
				char **arguments)
{
	size_t i;
	for (i = 0; i<num_arguments; i++)
	{
		if (strcmp(arguments[i],"--run-test") == 0||
		    strcmp(arguments[i],"--run-tests-in") == 0)
			return RUN_SPECIFIED_TESTS;
		else if (strcmp(arguments[i],"--run-all-tests") == 0)
			return RUN_ALL_TESTS;
		else if (strcmp(arguments[i],"--list-tests") == 0) 
			return LIST_TESTS;
	}
	return NORMAL_MODE;
}

static
void determine_tests_to_run(const size_t num_arguments,
			    char **arguments)
{
	size_t i = 0;
	for (i = 1; i<num_arguments; i++)
		if (strcmp(arguments[i],"--run-test") == 0 ||
		    strcmp(arguments[i],"--run-tests-in") == 0)
			num_tests_to_run++;
	tests_to_run = CALLOC(num_tests_to_run,test_name_t);
	size_t current_test = 0;
	for (i = 1; i<num_arguments; i++)
	{
		if (strcmp(arguments[i],"--run-test") == 0)
		{
			i++;
			memcpy(tests_to_run[current_test++].name,
			       arguments[i],
			       strlen(arguments[i])+1);
		}
		else if(strcmp(arguments[i],"--run-tests-in") == 0)
		{
			i++;
			memcpy(tests_to_run[current_test++].filename,
			       arguments[i],
			       strlen(arguments[i])+1);
		}
	}
}

void test_header(int num_arguments,
		 char **arguments)
{
	MULTIPLE_CALL_GUARD;
	test_mode = determine_test_mode(num_arguments,
					arguments);
	switch (test_mode)
	{
	case NORMAL_MODE:
	case RUN_ALL_TESTS:
		printf(TC_BLUE "Running all tests\n" TC_NORMAL);
		break;
	case RUN_SPECIFIED_TESTS:
		printf(TC_BLUE "Running specified tests\n" TC_NORMAL);
		determine_tests_to_run(num_arguments,
				       arguments);
		break;
	case LIST_TESTS:
		printf(TC_BLUE "List tests\n" TC_NORMAL);
		break;
	}
}

int run_this_test(const char *name,
		  const char *filename)
{
	if (test_mode != RUN_SPECIFIED_TESTS)
		return 1;
	size_t i;
	for (i = 0; i<num_tests_to_run; i++)
		if (strcmp(name,tests_to_run[i].name) == 0 ||
		    strcmp(filename,tests_to_run[i].filename) == 0)
			return 1;
	return 0;
}

int print_name_only()
{
	return test_mode == LIST_TESTS;
}

test_name_t *failed_tests = NULL;
size_t num_failed_tests = 0;
void test_failed(const char *name,
		 const char *filename)
{
	num_failed_tests++;
	failed_tests = REALLOC(failed_tests,
			       num_failed_tests);
	memcpy(failed_tests[num_failed_tests-1].name,
	       name,
	       strlen(name)+1);
	memcpy(failed_tests[num_failed_tests-1].filename,
	       filename,
	       strlen(filename)+1);
}

static
const char *remove_directory_part(const char *program_name)
{
	const char *head = program_name;
	while (*head != 0)
	{
		if (*head == '/')
			program_name = head+1;
		head++;
	}
	return program_name;
}

void prepare_test_environment(const char *test_name,
			      const char *program_name)
{
	char directory_name_buffer[2048]={0};
	const char *name = remove_directory_part(program_name);
	sprintf(directory_name_buffer,
		"/tmp/%s/",name);
	if (!directory_exists(directory_name_buffer))
		create_directory(directory_name_buffer);
	sprintf(directory_name_buffer,
		"/tmp/%s/%s/",
		name,test_name);
	if (!directory_exists(directory_name_buffer))
		create_directory(directory_name_buffer);
	else
		clear_directory(directory_name_buffer);
	memcpy(current_test_directory,
	       directory_name_buffer,
	       strlen(directory_name_buffer));
	current_test_directory[strlen(directory_name_buffer)] = 0;
}

void test_summary()
{
	MULTIPLE_CALL_GUARD;
	if (tests_to_run != NULL)
		free(tests_to_run);
	for (size_t i = 0; i<num_test_file_paths; i++)
		free(test_file_paths[i]);
	if (test_file_paths != NULL)
		free(test_file_paths);
	if (num_failed_tests>0)
	{
		printf(TC_RED "Listing all failed tests:\n");
		for (size_t i = 0; i<num_failed_tests; i++)
			printf(TC_RED "test %s@%s failed\n" TC_NORMAL,
			       failed_tests[i].name,
			       failed_tests[i].filename);
		exit(EXIT_FAILURE);
	}
	if (test_mode != NORMAL_MODE)
		exit(EXIT_SUCCESS);
	printf(TC_GREEN "all tests passed\n" TC_NORMAL);
}
