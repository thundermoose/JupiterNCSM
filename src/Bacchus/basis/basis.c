#include <basis/basis.h>
#include <assert.h>

struct _basis_
{
	vector_t *vectors;
	size_t num_vectors;
};

basis_t new_basis_from_vectors(
		vector_t *vectors,
		size_t num_vectors)
{
	basis_t basis = (basis_t)malloc(sizeof(struct _basis_));
	basis->vectors = vectors;
	basis->num_vectors = num_vectors;
	return basis;
}

vector_t construct_vector(basis_t basis,
		double *amplitudes,
		size_t num_amplitudes)
{
	assert(basis != NULL);
	assert(amplitudes != NULL);
	assert(basis->vectors != NULL);
	assert(num_amplitudes < basis->num_vectors);	
	return NULL;
}

void free_basis(basis_t basis)
{
}
