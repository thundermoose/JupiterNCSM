#ifndef __SP_STATES__
#define __SP_STATES__
#include <stdlib.h>
#include "shells.h"

// To indicate intention
// to index a single particle
// list
typedef int sp_state_index;

// A single particle state,
// corresponding to shell shell
typedef struct _sp_state_{
	shell_index shell;
	quantum_number m;
} SP_State;

// Contains all single particle states
typedef struct _sp_states_{
	SP_State* sp_states;
	sp_state_index dimension;
	Shells* shells;
} SP_States;


// To generate a single particle list
// from a list of all shells
SP_States* new_sp_states(Shells* shells);

// Lists the sp_states
void list_sp_states(SP_States* s);

// free all memory used by sp_states
void free_sp_states(SP_States* sp_states);
#endif
