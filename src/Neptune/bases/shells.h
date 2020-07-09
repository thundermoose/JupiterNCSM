#ifndef __SHELLS__
#define __SHELLS__

#include <stdlib.h>
#include <stdio.h>

// To indicate intention
// to reperesent some integer
// quantum number
typedef int quantum_number;

// To indicate intention
// to index the shells
typedef size_t true_shell_index;
typedef ssize_t strue_shell_index;
typedef size_t shell_index;
typedef ssize_t sshell_index;

typedef struct _true_shell_
{
  quantum_number e,n,l,j;
} True_Shell;

// We use harmonic oscilator
// states as our single particle
// basis
typedef struct _shell_
{
  true_shell_index tse;
  quantum_number e,n,l,j,tz;
} Shell;

typedef struct _shells_
{
  Shell* shells;
  shell_index num_of_shells;
  True_Shell* true_shells;
  true_shell_index num_of_true_shells;
  quantum_number e_max;
} Shells;

// To generate all shells
// with energy below e_max GJ style
Shells* new_shells(quantum_number e_max);

// To generate all shells
// with energy below e_max Antoine style
Shells* new_antoine_shells(quantum_number e_max);

// To read all shells given
// in a ascii filestream
Shells* new_shells_from_ascii_file(FILE* f);

// To find the lower limit of an energy range
true_shell_index find_minimum_true_shell(Shells* shells,
					 quantum_number e);
// To find the upper limit of an energy range
true_shell_index find_maximum_true_shell(Shells* shells,
					 quantum_number e);

// Lists all shells in shells
// and prints them to stdout
void list_shells(Shells* shells);


// compare two shells
// returns 0 if equal
// 1 if a is a higher shell
// -1 if b is a higher shell
int cmp_shells(Shell *a,
	       Shell *b);

// for every shell in shells_a
// find the corresponding shell in shells_b
// if any, and store the corresponding
// index, if no corresponding shell exists,
// -1 is used instead
sshell_index* matching_shells(Shells* shells_a,
			      Shells* shells_b);

strue_shell_index* matching_true_shells(Shells* shells_a,
					Shells* shells_b);

int triangle_in_equality(True_Shell shell_a,
		True_Shell shell_b,
		quantum_number J);

quantum_number shell_parity(True_Shell shell);

// To free all memory used
// by the Shells object
void free_shells(Shells* shells);

#endif
