#ifndef __LOG__
#define __LOG__

#include <stdlib.h>
#ifndef NLOGING
#define log_entry(message...)						\
{									\
	static size_t num_calls = 0;					\
	num_calls++; 							\
 	log_log(__FILE__,__func__,__LINE__,num_calls,message);		\
}
#else
#define log_entry(message...)
#endif

void initiate_logging(const char *log_file_name_variable,
		      const char *default_log_file_name);

void log_log(const char *file_name,
		const char *function_name,
		size_t line_number,
		size_t num_calls,
		const char *format,...);

void finalize_logging();
#endif
