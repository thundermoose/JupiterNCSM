#include <interaction/interaction.h>
#include <thundertester/test.h>
#include <log/log.h>
#include <math.h>

new_test(finding_pp_diagonal_state,
	 const char *interaction_path = 
	 TEST_DATA "n2lo_sat_2nf_hw20_Nmax2_noho";
	 interaction_t interaction = new_interaction(interaction_path);
	 int bra_state[2] = {2,14};
	 int ket_state[2] = {2,14};
	 double expected_value = 10.274936;
	 double value = get_matrix_element(interaction,
					   bra_state,
					   ket_state,
					   2,
					   1,1,-2,-4);
	 log_entry("value = %lg, expected %lg",
		   value,expected_value);
	assert_that(fabs(value-expected_value)<1e-5);
	free_interaction(interaction);
	);

new_test(finding_np_diagonal_state,
	 const char *interaction_path = 
	 TEST_DATA "n2lo_sat_2nf_hw20_Nmax4_noho";
	 interaction_t interaction = new_interaction(interaction_path);
	 int bra_state[2] = {1,6};
	 int ket_state[2] = {23,12};
	 double expected_value = 0.346070;
	 double value = get_matrix_element(interaction,
					   bra_state,
					   ket_state,
					   2,
					   1,3,0,0);
	 log_entry("value = %lg, expected %lg",
		   value,expected_value);
	assert_that(fabs(value-expected_value)<1e-5);
	free_interaction(interaction);
	);
