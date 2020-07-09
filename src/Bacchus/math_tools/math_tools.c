#include <math_tools/math_tools.h>
#include <math.h>
#include <thundertester/test.h>

size_t sum_sizes(const size_t *sizes,size_t num_sizes)
{
	size_t accumulator = 0;
	for (size_t i = 0; i<num_sizes; i++)
		accumulator+=sizes[i];
	return accumulator;
}

double array_scalar_product(const double *first_array,
			    const double *second_array,
			    size_t num_elements)
{
	double accumulator = 0;
	for (size_t i = 0; i<num_elements; i++)
		accumulator += first_array[i]*second_array[i];
	return accumulator;
}

void subtract_array_projection(double *target_array,
			       double projection,
			       const double *direction_array,
			       size_t num_elements)
{
	for (size_t i = 0; i<num_elements; i++)
		target_array[i]-=projection*direction_array[i];
}

double array_square_norm(const double *array,
			 size_t num_elements)
{
	double accumulator = 0.0;
	for (size_t i = 0; i<num_elements; i++)
		accumulator+=array[i]*array[i];
	return accumulator;
}

void scale_array(double *array,
		 double scaling,
		 size_t num_elements)
{
	for (size_t i = 0; i<num_elements; i++)
		array[i]*=scaling;
}

new_test(summing_1_2_3_expects_6,
	 const size_t array[3] = {1,2,3};
	 assert_that(sum_sizes(array,3) == 6);
	);

new_test(scalar_product_of_1_1_1_and_1_1_1_expects_3,
	 const double first_array[3] = {1,1,1};
	 const double second_array[3] = {1,1,1};
	 assert_that(fabs(array_scalar_product(first_array,second_array,3)-3)<1e-10)
	);

