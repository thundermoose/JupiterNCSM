#include <particle_type/particle_type.h>
#include <log/log.h>
#include <error/error.h>

particle_type_t parse_particle_type(const char *string)
{
	log_entry("parse_particle_type(%s)",string);
	switch (*string)
	{
		case 'n':
			return NEUTRON;
		case 'p':
			return PROTON;
		default:
			error("Not a particle type %s\n",
			      string);
	}	
}
