#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "questions.h"
#include <string_tools/string_tools.h>

int yes_no_question(const char* question)
{
    size_t trial;
    for (trial = 0; trial<3; trial++)
	{
	    if (trial > 0)
		{
		    printf("The answer you provided is not valid, try again");
		    fflush(stdin); // To remove any
		}
	    char answer[64];
	    printf("%s [yes(y)/no(n)/quit(q)]: ",
		   question);
	    if (scanf("%s",answer) != 1)
		continue;
	    if (strcmp(answer,"yes")== 0 ||
		strcmp(answer,"y")==0)
		{
		    return 1;
		}
	    if (strcmp(answer,"no")== 0 ||
		strcmp(answer,"n")==0)
		{
		    return 0;
		}
	    if (strcmp(answer,"quit")== 0 ||
		strcmp(answer,"q")==0)
		{
		    return -1;
		}
	
	}
    printf("The answer you provided is not valid, interpreting it as quit\n");
    return -1;
}


ssize_t ask_for_index(const char* question)
{
    size_t trial;
    for (trial = 0; trial<3; trial++)
	{
	    if (trial >0)
		{
		    printf("The answer you provided is not valid, try again");
		    fflush(stdin); // To remove any
		}
	    char answer[64];
	    printf("%s [<positiv integer>/quit(q)]: ",
		   question);
	    if (scanf("%s",answer) != 1)
		continue;
	    if (is_integer(answer))
		{
		    return atoi(answer);
		}
	    if (strcmp(answer,"quit")== 0 ||
		strcmp(answer,"q")==0)
		{
		    return -1;
		}
	
	}
    printf("The answer you provided is not valid, interpreting it as quit\n");
    return -1;
}
