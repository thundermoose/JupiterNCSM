#ifndef __MATH_TOOLS__
#define __MATH_TOOLS__

#include <stdlib.h>

#define square(a) (a)*(a)

size_t sum_sizes(const size_t *sizes,size_t num_sizes);

double array_scalar_product(const double *first_array,
			    const double *second_array,
			    size_t num_elements);

void subtract_array_projection(double *target_array,
			       double projection,
			       const double *direction_array,
			       size_t num_elements);

double array_square_norm(const double *array,
			 size_t num_elements);

void scale_array(double *array,
		 double scaling,
		 size_t num_elements);

void array_add_scaled(double *target_array,
		      double scaling_factor,
		      const double *term,
		      size_t num_elements);
#endif
