#include <stdlib.h>
#include <stdio.h>

int main(int num_arguments,
	 char **argument_list)
{
	arguments_t arguments = parse_arguments(num_arguments,
						argument_list);
	combination_table_t combination_table = 
		new_combination_table(get_combination_table_argument(arguments),
				      get_num_protons_argument(arguments),
				      get_num_neutron_argument(arguments));
}
