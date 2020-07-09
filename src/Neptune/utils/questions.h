#ifndef __QUESTIONS__
#define __QUESTIONS__
#include <stdlib.h>

/* This method asks a provided
 * yes no question to the use,
 * and wait for  yes [y], no [n] or
 * quit (q) (it is not case sensitive).
 * If some invalid string is entered,
 * it complains and asks the question
 * again up to 3 times. 
 * It returns 1 for yes, 0 for no and
 * -1 for quit or unanswered question
 */
int yes_no_question(const char* question);


/* This method asks a provided
 * question for which the answere is
 * either a positive integer or quit.
 * If some invalid string is entered,
 * it complains and asks the again,
 * up to 3 times.
 * It returns either a positive integer,
 * that corresponds to the integer answer
 * provided by the user or -1 corresponding
 * to quit or threefold invalid answer
 */
ssize_t ask_for_index(const char* question);

#endif
