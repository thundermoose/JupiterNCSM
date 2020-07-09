#include <stdlib.h>
#include <stdio.h>
#include <thundertester/test.h>
#include <log/log.h>
#include <matrix/matrix.h>

__attribute__((constructor(101)))
void setting_up_log()
{
	initiate_logging("BACCHUS_LOG_FILE",
			 "bacchus.log");
}

int main(int num_arguments, char **argument_list)
{
	printf("Just a dummy\n");
}

__attribute__((destructor(900)))
void finalize()
{
	finalize_logging();
}

new_test(a_simple_test,
	 printf("The test works\n");
	 assert_that(1);
	);
