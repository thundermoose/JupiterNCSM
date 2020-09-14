#include <time_tools/time_tools.h>
#include <time.h>

void get_time_stamp(int *year,
		    int *month,
		    int *day,
		    int *hour,
		    int *minute,
		    int *second)
{
	time_t now = time(NULL);
	struct tm *current_time = localtime(&now);
	*year = current_time->tm_year+1900;
	*month = current_time->tm_mon + 1;
	*day = current_time->tm_mday;
	*hour = current_time->tm_hour;
	*minute = current_time->tm_min;
	*second = current_time->tm_sec;
}
