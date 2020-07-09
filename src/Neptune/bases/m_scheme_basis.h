#ifndef __M_SCHEME_BASIS__
#define __M_SCHEME_BASIS__

struct _m_scheme_basis_;
typedef struct _m_scheme_basis_ *m_scheme_basis_t;

m_scheme_basis_t new_m_scheme_basis(quantum_number e_max1,
				    quantum_number e_max2);
#endif
