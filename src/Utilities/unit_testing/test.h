/*
 *************** License ***************
 thunder-tester is a simple system for in code unit testing in C using gcc
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

#ifndef __TEST__
#define __TEST__

#include <stdlib.h>
#include <stdio.h>
#include "tc_colors.h"

/* A function that returns a file path
 * that depends on the test name
 */
const char* get_test_file_path(const char* file_name);

#ifdef TEST
/* The following type is used
 * internally to indicate the
 * status of the test 
 */
typedef enum
{
	TEST_PASSED,
	TEST_FAILED
} test_status_t;

typedef struct
{
	test_status_t status;
	char file_name[1024];
	size_t line_number;
} test_result_t;

/* Displays the test header and handles
 * potential test arguments
 */
void test_header(int num_arguments,
		 char **arguments);

/* Used internally to register a failed test
 */
void test_failed(const char *name,
		 const char *filename);

/* Setup a directory tree under /tmp
 */
void prepare_test_environment(const char *test_name,
			      const char *program_name);
/* Used internally to determine if the 
 * current test should be executed
 */
int run_this_test(const char *name,
		  const char *filename);
int print_name_only();
void test_summary();
#define assert_that(code...)				\
	if (!(code))					\
{						\
	test_result_t result =			\
	{					\
		.status = TEST_FAILED,          \
		.file_name = __FILE__,		\
		.line_number = __LINE__		\
	};					\
	return result;				\
}

#define new_test(name,code...)				\
	test_result_t test_##name()			\
{						\
	code;					\
	{					\
		test_result_t result = 		\
		{				\
			.status = TEST_PASSED,  \
			.file_name = __FILE__,	\
			.line_number = __LINE__	\
		};				\
		return result;			\
	}					\
}								\
__attribute__((constructor(500)))				\
void runtest_##name(int num_arguments,char **argument_list)	\
{								\
	if (!run_this_test(#name,__FILE__))			\
	return;						\
	if (print_name_only())					\
	{							\
		printf(#name " @ %s\n",__FILE__);		\
		return;						\
	}							\
	prepare_test_environment(#name,*argument_list);		\
	printf(TC_BLUE "running test " #name "\n" TC_NORMAL);	\
	test_result_t result = test_##name();			\
	switch (result.status)					\
	{							\
		case TEST_PASSED:					\
			printf(TC_GREEN "test " #name " passed\n" TC_NORMAL); \
		break;						\
		case TEST_FAILED:					\
			printf(TC_RED "test " #name 			\
			       " failed at %s@%lu\n" TC_NORMAL,	\
			       result.file_name,			\
			       result.line_number); 			\
		test_failed(#name, __FILE__);			\
		break;						\
	}							\
}

#define new_test_silent(name,code...)				

	__attribute__((constructor(450)))
static void run_test_header(
			    int num_arguments,
			    char **arguments)
{
	test_header(
		    num_arguments,
		    arguments);
}

	__attribute__((constructor(550)))
static void run_test_summary()
{
	test_summary();
}
#else

#define new_test(name,code...)

#define new_test_silent(name,code...)				
#endif

#endif
