#include <log/log.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>


FILE *log_file;

void initiate_logging(const char *log_file_name_variable,
		      const char *default_log_file_name)
{
	const char *file_name = getenv(log_file_name_variable);
	if (file_name == NULL)
		file_name = default_log_file_name;
	log_file = fopen(file_name,"w");
}

void time_stamp()
{
	time_t t= time(NULL);
	struct tm local_time = *localtime(&t);
	int year = local_time.tm_year+1900;
	int month = local_time.tm_mon+1;
	int day = local_time.tm_mday;
	int hour = local_time.tm_hour;
	int min = local_time.tm_min;
	int sec = local_time.tm_sec;
	fprintf(log_file,"Y%04d-M%02d-D%02d %02d:%02d:%02d:\n",
		year,
		month,
		day,
		hour,
		min,
		sec);
}

void new_log_entry(const char *file_name,
		const char *function_name,
		size_t line_number,
		size_t num_calls,
		const char *user_message)
{
	time_stamp();
	fprintf(log_file,"(%lu) In %s@%s:%lu: \"%s\"\n",
		num_calls,
		function_name,
		file_name,
		line_number,
		user_message);
	fflush(log_file);
}

void log_log(const char *file_name,
		const char *function_name,
		size_t line_number,
		size_t num_calls,
		const char *format,...)
{
	va_list last_arguments;
	va_start(last_arguments,format);
	char user_message[1024];
	vsprintf(user_message,format,last_arguments);
	new_log_entry(file_name,
			function_name,
			line_number,
			num_calls,
			user_message);	
}


void finalize_logging()
{
	time_stamp();
	fprintf(log_file,
		"End of Log\n");
	fflush(log_file);
	fclose(log_file);
}
