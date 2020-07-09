#include <string_tools/string_tools.h>
#include <bases/shells.h>
#include <search.h>
#include <string.h>

/* The following function generate a 
 * Shells object with all single particle
 * shells up to e_max
 */
Shells* new_shells(quantum_number e_max)
{
	Shells* shells =
		(Shells*)malloc(sizeof(Shells));
	shells->e_max = e_max;
	// It is easy to show that the following
	// expression is the exact number of
	// shells below e_max
	shells->num_of_true_shells =
		(1+e_max+(e_max*(e_max+1))/2); 
	shells->num_of_shells =
		2*shells->num_of_true_shells;
	shells->shells =
		(Shell*)malloc(sizeof(Shell)*
			       shells->num_of_shells);
	shells->true_shells =
		(True_Shell*)malloc(sizeof(True_Shell)*
				    shells->num_of_true_shells);
	Shell current;
	shell_index i = 0;
	true_shell_index j = 0;
	for (current.e=0;
	     current.e<=e_max;
	     current.e++)
	{
		// For constant harmonic oscillator energy
		// l must change with 2
		for (current.l = current.e;
		     current.l>=(current.e&1);
		     current.l-=2)
		{
			current.n = (current.e-current.l)/2;
			// By the triangle inequality for
			// coupling the orbital angular momentum
			// and fermionic spin total j, here
			// repersented as 2*j, must be (l+1/2)
			// and |l-1/2|
			for (current.j=current.l*2+1;
			     current.j>=abs(current.l*2-1);
			     current.j-=2)
			{
				current.tse = j;
				True_Shell ts ={current.e,
					current.n,
					current.l,
					current.j};
				shells->true_shells[j++] = ts;
				// protons have tz = -1
				// neutrons have tz = 1
				for(current.tz=-1;
				    current.tz<=1;
				    current.tz+=2)
				{
					shells->shells[i++] = current;
				}
			}
		}
	}  
	return shells;
}

Shells* new_antoine_shells(quantum_number e_max)
{
	Shells* shells =
		(Shells*)malloc(sizeof(Shells));
	shells->e_max = e_max;
	// It is easy to show that the following
	// expression is the exact number of
	// shells below e_max
	shells->num_of_true_shells =
		(1+e_max+(e_max*(e_max+1))/2); 
	shells->num_of_shells =
		2*shells->num_of_true_shells;
	shells->shells =
		(Shell*)malloc(sizeof(Shell)*
			       shells->num_of_shells);
	shells->true_shells =
		(True_Shell*)malloc(sizeof(True_Shell)*
				    shells->num_of_true_shells);
	Shell current;
	shell_index i = 0;
	true_shell_index j = 0;
	for (current.e=0;
	     current.e<=e_max;
	     current.e++)
	{
		// For constant harmonic oscillator energy
		// l must change with 2
		for (current.l = (current.e&1);
		     current.l<=current.e;
		     current.l+=2)
		{
			current.n = (current.e-current.l)/2;
			// By the triangle inequality for
			// coupling the orbital angular momentum
			// and fermionic spin total j, here
			// repersented as 2*j, must be (l+1/2)
			// and |l-1/2|
			for (current.j=abs(current.l*2-1);
			     current.j<=current.l*2+1;
			     current.j+=2)
			{
				current.tse = j;
				True_Shell ts ={current.e,
					current.n,
					current.l,
					current.j};
				shells->true_shells[j++] = ts;
				// protons have tz = -1
				// neutrons have tz = 1
				for(current.tz=-1;
				    current.tz<=1;
				    current.tz+=2)
				{
					shells->shells[i++] = current;
				}
			}
		}
	}  
	return shells;
}
strue_shell_index find_true_shell(True_Shell* true_shells,
				  true_shell_index num_true_shells,
				  True_Shell lookfor)
{
	true_shell_index i;
	for (i = 0; i<num_true_shells; i++)
	{
		if (true_shells[i].e == lookfor.e &&
		    true_shells[i].n == lookfor.n &&
		    true_shells[i].l == lookfor.l &&
		    true_shells[i].j == lookfor.j)
		{
			return i;
		}
	}
	return -1;
}


/* The following function reads a given file  
 * and generate a Shells object.
 * It expects that the first non empty row
 * contians the text
 * "---single-particle-states---"
 * then it skips empty lines or lines begining 
 * with "---".
 * then it takes a list of integers
 * formated as 
 * i n l j tz
 * then when either the file ends
 * or a new line begining with --- comes
 * it assumes that the shell block is over
 */
Shells* new_shells_from_ascii_file(FILE* f)
{
	size_t n = 0;
	char* row = NULL;
	Shells* out = NULL;
	size_t allocated_size = 0;
	// We allow for any number of blank
	// rows in the begining of the file
	do
	{
		if (getline(&row,&n,f)<0)
		{
			if (row != NULL)
				free(row);

			return NULL;
		}
	}
	while (strcmp(row,"\n")==0);

	// We need to verifying that the row it reads is
	// the header of a single particle block
	if (strcmp(row,
		   "---single-particle-states---\n") == 0)
	{
		out = (Shells*)malloc(sizeof(Shells));
		// At the moment what out->e_max should be
		// is unknown, however we know for a fact
		// that it atleast is >= 0.
		out->e_max = 0;
		// We do not know the final size
		// of the array at this point
		// so we assume that it is empty
		out->num_of_shells = 0;
		out->shells = NULL;
	}
	else
	{
		free(row);
		return NULL;
	}
	// As always skip any following blank lines
	// or comments
	do
	{
		if (getline(&row,&n,f)<0)
		{
			free(row);
			free(out);
			return NULL;
		}
	}
	while (strcmp(row,"\n")==0 ||
	       begins_with(row,"---"));
	size_t position;
	// Now we are ready to read
	// the actual data of the ascii file
	do
	{

		printf("row: \"%s\"\n",
		       row);
		Shell shell;
		if (sscanf(row,
			   "%*d  %d  %d  %d  %d",
			   &shell.n,
			   &shell.l,
			   &shell.j,
			   &shell.tz)<4 &&
		    strcmp(row,"\n")!=0)
		{
			fprintf(stderr,
				"wrong line format in \"%s\"\n",
				row);
			free(row);
			exit(1);
		}
		else if(strcmp(row,"\n")!=0)
		{
			shell.e = 2*shell.n+shell.l;
			// Our current best estimate of the final size
			// of out->shells could be wrong, and if we
			// assumed a to small size we need to increase it
			if (allocated_size==out->num_of_shells)
			{
				// We double the size of the array
				// this will only make us have to
				// do roughly log(n) allocations where
				// n is the acctual final size
				allocated_size=2*allocated_size+1;
				out->shells =
					(Shell*)realloc(out->shells,
							sizeof(Shell)*
							allocated_size);
			}
			out->shells[out->num_of_shells++]=shell;	  
			if (shell.e>out->e_max)
			{
				out->e_max = shell.e;
			}
		}
		position =ftell(f);
		if (getline(&row,&n,f)<0)
		{
			break;
		}
	}
	while (!begins_with(row,"---"));

	printf("out->num_of_shells = %ld\n",
	       out->num_of_shells);

	// any comments that follows the
	// data are assumed to belong to
	// the next block
	if (begins_with(row,"---"))
	{
		fseek(f,position,SEEK_SET);
	}
	free(row);
	// We do now know the actual final size
	// and can therefore rescale our out->shells
	// array to fit, since it most certainly is
	// to large. In fact only if n=2^k-1 the array
	// has the correct size
	out->shells =
		(Shell*)realloc(out->shells,
				sizeof(Shell)*
				out->num_of_shells);

	// We can from the shells now construct
	// the true_shells array
	// we do not know how many they should be
	// but know that the must be fewer than
	// the shells
	true_shell_index allocated_true_shells =
		out->num_of_shells;
	out->true_shells =
		(True_Shell*)malloc(sizeof(True_Shell)*
				    allocated_true_shells);
	out->num_of_true_shells = 0;
	shell_index i;
	for (i = 0; i<out->num_of_shells; i++)
	{
		// construct corresponding true_shell
		True_Shell current =
		{out->shells[i].e,
			out->shells[i].n,
			out->shells[i].l,
			out->shells[i].j};
		// check if true_shell already exists,
		// if so add true shell index to the shell
		strue_shell_index cur =
			find_true_shell(out->true_shells,
					out->num_of_true_shells,
					current);
		if (cur!=-1)
		{
			out->shells[i].tse = cur;

		}
		else // if it does not exists, append it last
		{
			out->true_shells[out->num_of_true_shells++] = current;
		}

	}
	// The program might have overestimated the needed size of
	// the true_shells array
	if (allocated_true_shells>out->num_of_true_shells)
	{
		allocated_true_shells = out->num_of_true_shells;
		out->true_shells =
			(True_Shell*)realloc(out->true_shells,
					     sizeof(True_Shell)*
					     allocated_true_shells);
	}				   
	return out;
}


void list_shells(Shells* shells)
{
	shell_index shell_i;
	printf("Shells:\n");
	printf("e\tn\tl\tj\ttz\n");
	for (shell_i = 0;
	     shell_i<shells->num_of_shells;
	     shell_i++)
	{
		printf("(%ld): %d\t%d\t%d\t%d\t%d\n",
		       shell_i,
		       shells->shells[shell_i].e,
		       shells->shells[shell_i].n,
		       shells->shells[shell_i].l,
		       shells->shells[shell_i].j,
		       shells->shells[shell_i].tz);
	}
	printf("\n\n\n");
	printf("True Shells:\n");
	printf("e\tn\tl\tj\ttz\n");
	for (shell_i = 0;
	     shell_i<shells->num_of_true_shells;
	     shell_i++)
	{
		printf("(%ld): %d\t%d\t%d\t%d\n",
		       shell_i,
		       shells->true_shells[shell_i].e,
		       shells->true_shells[shell_i].n,
		       shells->true_shells[shell_i].l,
		       shells->true_shells[shell_i].j);
	}
	printf("\n\n\n");
}

int cmp_shells(Shell *a,
	       Shell *b)
{
	int diff = a->e-b->e;
	if (diff)
		return diff;
	diff = a->l-b->l;
	if (diff)
		return diff;
	diff = a->j-b->j;
	if (diff)
		return diff;
	return a->tz-b->tz;
}

sshell_index find_shell(Shells* shells,
			Shell look_for)
{
	Shell* candidate = (Shell*)
		lfind((void*)&look_for,
		      (void*)shells->shells,
		      &shells->num_of_shells,
		      sizeof(Shell),
		      (__compar_fn_t)cmp_shells);
	if (candidate == NULL)
		return -1;
	return (sshell_index)(candidate-shells->shells);

}

sshell_index* matching_shells(Shells* shells_a,
			      Shells* shells_b)
{
	sshell_index *conversion =
		(sshell_index*)malloc(sizeof(sshell_index)*
				      shells_a->num_of_shells);
	shell_index i;
	for (i = 0;
	     i<shells_a->num_of_shells;
	     i++)
	{
		conversion[i] =
			find_shell(shells_b,
				   shells_a->shells[i]);
	}
	return conversion;
}
/*
   strue_shell_index find_true_shell(Shells *shells,
   True_Shell ts)
   {
   true_shell_index i;
   for (i = 0; i<shells->num_of_true_shells; i++)
   {
   if (shells->true_shells[i].e = ts.e &&
   shells->true_shells[i].n = ts.n &&
   shells->true_shells[i].l = ts.l &&
   shells->true_shells[i].j = ts.j)
   return i;
   }
   return -1;
   }
 */
strue_shell_index* matching_true_shells(Shells* shells_a,
					Shells* shells_b)
{
	strue_shell_index* conversion =
		(strue_shell_index*)malloc(sizeof(strue_shell_index)*
					   shells_a->num_of_true_shells);

	true_shell_index i;
	for (i = 0;
	     i < shells_a->num_of_true_shells;
	     i++)
	{
		conversion[i] =
			find_true_shell(shells_b->true_shells,
					shells_b->num_of_true_shells,
					shells_a->true_shells[i]);
	}
	return conversion;
}


// To find the lower limit of an energy range
true_shell_index find_minimum_true_shell(Shells* shells,
					 quantum_number e)
{
	// we assume that the shells are ordered in
	// increasing energy, so we can use binary search
	true_shell_index min_i = 0;
	true_shell_index max_i = shells->num_of_true_shells;
	true_shell_index i = 0;
	while (max_i-min_i>1)
	{
		i = (max_i+min_i)>>1;
		quantum_number e_i = shells->true_shells[i].e;
		if (e_i<e)
		{
			min_i = i;
		}
		else if (e_i>e)
		{
			max_i = i;
		}
		else
		{
			break;
		}
	}
	while (i-1>0 &&
	       shells->true_shells[i-1].e==e)
	{
		i--;
	}
	return i;
}
// To find the upper limit of an energy range
true_shell_index find_maximum_true_shell(Shells* shells,
					 quantum_number e)
{
	// we assume that the shells are ordered in
	// increasing energy, so we can use binary search
	true_shell_index min_i = 0;
	true_shell_index max_i = shells->num_of_true_shells;
	true_shell_index i = 0;
	while (max_i-min_i>1)
	{
		i = (max_i+min_i)>>1;
		quantum_number e_i = shells->true_shells[i].e;
		if (e_i<e)
		{
			min_i = i;
		}
		else if (e_i>e)
		{
			max_i = i;
		}
		else
		{
			break;
		}
	}
	while (i+1<shells->num_of_true_shells &&
	       shells->true_shells[i+1].e==e)
	{
		i++;
	}
	return i;
}

int triangle_in_equality(True_Shell shell_a,
			 True_Shell shell_b,
			 quantum_number J)
{
	return shell_a.j+shell_b.j >= 2*J && abs(shell_a.j-shell_b.j)<=2*J;
}

quantum_number shell_parity(True_Shell shell)
{
	return shell.l&1;
}


/* The following method frees all
 * memory allocated by either of the
 * two functions above
 */
void free_shells(Shells* shells){
	free(shells->shells);
	free(shells->true_shells);
	free(shells);
}
