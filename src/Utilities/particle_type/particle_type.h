#ifndef __PARTICLE_TYPE__
#define __PARTICLE_TYPE__

typedef enum
{
	NEUTRON,PROTON
} particle_type_t;

particle_type_t parse_particle_type(const char *string);

#endif
