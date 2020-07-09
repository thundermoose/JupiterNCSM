#include <matrix_vector_multiplication/matrix_vector_multiplication.h>
#include <log/log.h>
#include <stdlib.h>
#include <assert.h>

void multiplication_neutrons(vector_block_t out_block,
				 const vector_block_t in_block,
				 const matrix_block_t block,
				 const index_list_t neutron_list,
				 const int sign)
{
	const size_t num_neutron_indices =
		length_index_list(neutron_list);
	assert(get_proton_dimension(out_block) ==
	       get_proton_dimension(in_block));
	const size_t num_proton_states =
		get_proton_dimension(in_block);
	const size_t num_in_neutron_states =
		get_neutron_dimension(in_block);
	const size_t num_out_neutron_states =
		get_neutron_dimension(out_block);
	double *out_vector_elements =
		get_vector_block_elements(out_block);
	double *in_vector_elements =
		get_vector_block_elements(in_block);
	index_triple_t *neutron_indices =
		get_index_list_elements(neutron_list);
	double *matrix_elements = get_matrix_block_elements(block);
	log_entry("num_neutron_indices = %lu",
		  num_neutron_indices);
	for (size_t i = 0; i < num_neutron_indices; i++)
	{
		const size_t neutron_out_index = neutron_indices[i].out_index;
		const size_t neutron_in_index = neutron_indices[i].in_index;
		const double matrix_element =
			matrix_elements[neutron_indices[i].matrix_index];
		for (size_t proton_state = 0;
		     proton_state < num_proton_states;
		     proton_state++)
		{
			size_t out_index =
				neutron_out_index +
				num_out_neutron_states*proton_state;
			size_t in_index =
				neutron_in_index +
				num_in_neutron_states*proton_state;
			log_entry("%lg(%lu) %c= %lg(%lu) * %lg(%lu)",
				  out_vector_elements[out_index],
				  out_index,
				  sign > 0 ? '+' : '-',
				  matrix_element,
				  neutron_indices[i].matrix_index,
				  in_vector_elements[in_index],
				  in_index);
			out_vector_elements[out_index] += 
				sign * matrix_element * 
				in_vector_elements[in_index];	
		}
	}
}

void multiplication_protons(vector_block_t out_block,
				const vector_block_t in_block,
				const matrix_block_t block,
				const index_list_t proton_list,
				const int sign)
{
	const size_t num_proton_indices =
		length_index_list(proton_list);
	assert(get_neutron_dimension(out_block) ==
	       get_neutron_dimension(in_block));
	const size_t num_neutron_states =
		get_neutron_dimension(in_block);
	double *out_vector_elements =
		get_vector_block_elements(out_block);
	double *in_vector_elements =
		get_vector_block_elements(in_block);
	index_triple_t *proton_indices =
		get_index_list_elements(proton_list);
	double *matrix_elements = get_matrix_block_elements(block);
	log_entry("num_proton_indices = %lu",
		  num_proton_indices);
	for (size_t i = 0; i < num_proton_indices; i++)
	{
		const size_t proton_out_index = proton_indices[i].out_index;
		const size_t proton_in_index = proton_indices[i].in_index;
		const double matrix_element =
			matrix_elements[proton_indices[i].matrix_index];
		for (size_t neutron_state = 0;
		     neutron_state < num_neutron_states;
		     neutron_state++)
		{
			size_t out_index =
				num_neutron_states*proton_out_index +
				neutron_state;
			size_t in_index =
				num_neutron_states*proton_in_index +
				neutron_state;
			log_entry("%lg(%lu) %c= %lg(%lu) * %lg(%lu)",
				  out_vector_elements[out_index],
				  out_index,
				  sign > 0 ? '+' : '-',
				  matrix_element,
				  proton_indices[i].matrix_index,
				  in_vector_elements[in_index],
				  in_index);
			out_vector_elements[out_index] += 
				sign * matrix_element * 
				in_vector_elements[in_index];	
		}
	}
}

void multiplication_neutrons_protons(vector_block_t out_block,
					 const vector_block_t in_block,
					 const matrix_block_t block,
					 const index_list_t neutron_list,
					 const index_list_t proton_list,
					 const int sign)
{
	const size_t num_neutron_indices =
		length_index_list(neutron_list);
	const size_t num_proton_indices =
		length_index_list(proton_list);
	const size_t num_in_neutron_states =
		get_neutron_dimension(in_block);
	const size_t num_out_neutron_states =
		get_neutron_dimension(out_block);
	const size_t neutron_matrix_dimension =
		get_neutron_matrix_dimension(block);
	double *in_vector_elements =
		get_vector_block_elements(in_block);
	double *out_vector_elements =
		get_vector_block_elements(out_block);
	double *matrix_elements =
		get_matrix_block_elements(block);
	index_triple_t *neutron_indices =
		get_index_list_elements(neutron_list);
	index_triple_t *proton_indices =
		get_index_list_elements(proton_list);
	log_entry("num_neutron_indices = %lu",
		  num_neutron_indices);
	log_entry("num_proton_indices = %lu",
		  num_proton_indices);
	for (size_t i = 0; i < num_neutron_indices; i++)
	{
		const size_t neutron_in_index =
			neutron_indices[i].in_index;
		const size_t neutron_out_index =
			neutron_indices[i].out_index;
		const size_t neutron_matrix_index =
			neutron_indices[i].matrix_index;
		for (size_t j = 0; j<num_proton_indices; j++)
		{
			const size_t proton_in_index =
				proton_indices[j].in_index;
			const size_t proton_out_index =
				proton_indices[j].out_index;
			const size_t proton_matrix_index =
				proton_indices[j].matrix_index;
			const size_t in_index = 
				neutron_in_index +
				num_in_neutron_states * proton_in_index;
			const size_t out_index = 
				neutron_out_index +
				num_out_neutron_states * proton_out_index;
			const size_t matrix_index = 
				neutron_matrix_index +
				neutron_matrix_dimension * proton_matrix_index;
			const double matrix_element =
				matrix_elements[matrix_index];
			log_entry("%lg(%lu) %c= %lg(%lu,%lu/%lu,%lu) * %lg(%lu)",
				  out_vector_elements[out_index],
				  out_index,
				  sign > 0 ? '+' : '-',
				  matrix_elements[matrix_index],
				  matrix_index,
				  neutron_matrix_index,
				  neutron_matrix_dimension,
				  proton_matrix_index,
				  in_vector_elements[in_index],
				  in_index);
			out_vector_elements[out_index] +=
				sign * matrix_element *
				in_vector_elements[in_index];
		}
	}
}

void multiplication_neutrons_off_diag(vector_block_t out_block_left,
				      vector_block_t out_block_right,
				      const vector_block_t in_block_left,
				      const vector_block_t in_block_right,
				      const matrix_block_t block,
				      const index_list_t neutron_list,
				      const int sign)
{
	const size_t num_neutron_indices =
		length_index_list(neutron_list);
	assert(get_proton_dimension(out_block_left) ==
	       get_proton_dimension(in_block_left));
	assert(get_proton_dimension(out_block_right) ==
	       get_proton_dimension(in_block_right));
	const size_t num_proton_states =
		get_proton_dimension(in_block_left);
	const size_t num_in_neutron_states =
		get_neutron_dimension(in_block_left);
	const size_t num_out_neutron_states =
		get_neutron_dimension(out_block_left);
	double *out_vector_elements_left =
		get_vector_block_elements(out_block_left);
	double *out_vector_elements_right =
		get_vector_block_elements(out_block_right);
	double *in_vector_elements_left =
		get_vector_block_elements(in_block_left);
	double *in_vector_elements_right =
		get_vector_block_elements(in_block_right);
	index_triple_t *neutron_indices =
		get_index_list_elements(neutron_list);
	double *matrix_elements = get_matrix_block_elements(block);
	log_entry("num_neutron_indices = %lu",
		  num_neutron_indices);
	for (size_t i = 0; i < num_neutron_indices; i++)
	{
		const size_t neutron_out_index = neutron_indices[i].out_index;
		const size_t neutron_in_index = neutron_indices[i].in_index;
		const double matrix_element =
			matrix_elements[neutron_indices[i].matrix_index];
		for (size_t proton_state = 0;
		     proton_state < num_proton_states;
		     proton_state++)
		{
			size_t out_index =
				neutron_out_index +
				num_out_neutron_states*proton_state;
			size_t in_index =
				neutron_in_index +
				num_in_neutron_states*proton_state;
			log_entry("%lg(%lu) %c= %lg(%lu) * %lg(%lu)",
				  out_vector_elements_left[out_index],
				  out_index,
				  sign > 0 ? '+' : '-',
				  matrix_element,
				  neutron_indices[i].matrix_index,
				  in_vector_elements_left[in_index],
				  in_index);
			log_entry("%lg(%lu) %c= %lg(%lu) * %lg(%lu)",
				  out_vector_elements_right[in_index],
				  in_index,
				  sign > 0 ? '+' : '-',
				  matrix_element,
				  neutron_indices[i].matrix_index,
				  in_vector_elements_right[out_index],
				  out_index);
			out_vector_elements_left[out_index] += 
				sign * matrix_element * 
				in_vector_elements_left[in_index];	
			out_vector_elements_right[in_index] += 
				sign * matrix_element * 
				in_vector_elements_right[out_index];	
		}
	}
}

void multiplication_protons_off_diag(vector_block_t out_block_left,
				     vector_block_t out_block_right,
				     const vector_block_t in_block_left,
				     const vector_block_t in_block_right,
				     const matrix_block_t block,
				     const index_list_t proton_list,
				     const int sign)
{
	const size_t num_proton_indices =
		length_index_list(proton_list);
	assert(get_neutron_dimension(out_block_left) ==
	       get_neutron_dimension(in_block_left));
	assert(get_neutron_dimension(out_block_right) ==
	       get_neutron_dimension(in_block_right));
	const size_t num_neutron_states =
		get_neutron_dimension(in_block_left);
	double *out_vector_elements_left =
		get_vector_block_elements(out_block_left);
	double *in_vector_elements_left =
		get_vector_block_elements(in_block_left);
	double *out_vector_elements_right =
		get_vector_block_elements(out_block_right);
	double *in_vector_elements_right =
		get_vector_block_elements(in_block_right);
	index_triple_t *proton_indices =
		get_index_list_elements(proton_list);
	double *matrix_elements = get_matrix_block_elements(block);
	log_entry("num_proton_indices = %lu",
		  num_proton_indices);
	for (size_t i = 0; i < num_proton_indices; i++)
	{
		const size_t proton_out_index = proton_indices[i].out_index;
		const size_t proton_in_index = proton_indices[i].in_index;
		const double matrix_element =
			matrix_elements[proton_indices[i].matrix_index];
		for (size_t neutron_state = 0;
		     neutron_state < num_neutron_states;
		     neutron_state++)
		{
			size_t out_index =
				num_neutron_states*proton_out_index +
				neutron_state;
			size_t in_index =
				num_neutron_states*proton_in_index +
				neutron_state;
			log_entry("%lg(%lu) %c= %lg(%lu) * %lg(%lu)",
				  out_vector_elements_left[out_index],
				  out_index,
				  sign > 0 ? '+' : '-',
				  matrix_element,
				  proton_indices[i].matrix_index,
				  in_vector_elements_right[in_index],
				  in_index);
			log_entry("%lg(%lu) %c= %lg(%lu) * %lg(%lu)",
				  out_vector_elements_left[in_index],
				  in_index,
				  sign > 0 ? '+' : '-',
				  matrix_element,
				  proton_indices[i].matrix_index,
				  in_vector_elements_right[out_index],
				  out_index);
			out_vector_elements_left[out_index] += 
				sign * matrix_element * 
				in_vector_elements_left[in_index];	
			out_vector_elements_right[in_index] += 
				sign * matrix_element * 
				in_vector_elements_right[out_index];	
		}
	}
}

void multiplication_neutrons_protons_off_diag(vector_block_t out_block_left,
					      vector_block_t out_block_right,
					      const vector_block_t 
					      in_block_left,
					      const vector_block_t 
					      in_block_right,
					      const matrix_block_t block,
					      const index_list_t neutron_list,
					      const index_list_t proton_list,
					      const int sign)
{
	const size_t num_neutron_indices =
		length_index_list(neutron_list);
	const size_t num_proton_indices =
		length_index_list(proton_list);
	const size_t num_in_neutron_states =
		get_neutron_dimension(in_block_left);
	const size_t num_out_neutron_states =
		get_neutron_dimension(out_block_left);
	const size_t neutron_matrix_dimension =
		get_neutron_matrix_dimension(block);
	double *in_vector_elements_left =
		get_vector_block_elements(in_block_left);
	double *out_vector_elements_left =
		get_vector_block_elements(out_block_left);
	double *in_vector_elements_right =
		get_vector_block_elements(in_block_right);
	double *out_vector_elements_right =
		get_vector_block_elements(out_block_right);
	double *matrix_elements =
		get_matrix_block_elements(block);
	index_triple_t *neutron_indices =
		get_index_list_elements(neutron_list);
	index_triple_t *proton_indices =
		get_index_list_elements(proton_list);
	log_entry("num_neutron_indices = %lu",
		  num_neutron_indices);
	log_entry("num_proton_indices = %lu",
		  num_proton_indices);
	for (size_t i = 0; i < num_neutron_indices; i++)
	{
		const size_t neutron_in_index =
			neutron_indices[i].in_index;
		const size_t neutron_out_index =
			neutron_indices[i].out_index;
		const size_t neutron_matrix_index =
			neutron_indices[i].matrix_index;
		for (size_t j = 0; j<num_proton_indices; j++)
		{
			const size_t proton_in_index =
				proton_indices[j].in_index;
			const size_t proton_out_index =
				proton_indices[j].out_index;
			const size_t proton_matrix_index =
				proton_indices[j].matrix_index;
			const size_t in_index = 
				neutron_in_index +
				num_in_neutron_states * proton_in_index;
			const size_t out_index = 
				neutron_out_index +
				num_out_neutron_states * proton_out_index;
			const size_t matrix_index = 
				neutron_matrix_index +
				neutron_matrix_dimension * proton_matrix_index;
			const double matrix_element =
				matrix_elements[matrix_index];
			log_entry("%lg(%lu) %c= %lg(%lu,%lu/%lu,%lu) * %lg(%lu)",
				  out_vector_elements_left[out_index],
				  out_index,
				  sign > 0 ? '+' : '-',
				  matrix_elements[matrix_index],
				  matrix_index,
				  neutron_matrix_index,
				  neutron_matrix_dimension,
				  proton_matrix_index,
				  in_vector_elements_left[in_index],
				  in_index);
			log_entry("%lg(%lu) %c= %lg(%lu,%lu/%lu,%lu) * %lg(%lu)",
				  out_vector_elements_right[in_index],
				  in_index,
				  sign > 0 ? '+' : '-',
				  matrix_elements[matrix_index],
				  matrix_index,
				  neutron_matrix_index,
				  neutron_matrix_dimension,
				  proton_matrix_index,
				  in_vector_elements_right[out_index],
				  out_index);
			out_vector_elements_left[out_index] +=
				sign * matrix_element *
				in_vector_elements_left[in_index];
			out_vector_elements_right[in_index] +=
				sign * matrix_element *
				in_vector_elements_right[out_index];
		}
	}
}
