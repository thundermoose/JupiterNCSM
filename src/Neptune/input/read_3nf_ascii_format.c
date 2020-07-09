#include <stdlib.h>
#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <errno.h>
#include <input/read_3nf_ascii_format.h>
#include <bases/jt_basis_3p.h>
#include <string_tools/string_tools.h>

#define HAS_CONFIG_FILE 1
#define HAS_MATRIX_FILE 2
#define IS_VALID 3

int scan_directory(const char* directory_name,
		   DIR* simple_directory,
		   char config_file[],
		   char matrix_file[])
{
  struct dirent* child;
  int is_valid = 0;
  while ((child=readdir(simple_directory)))
    {
      // A config filename ends in this
      // context with .cnf
      if (ends_with(child->d_name,
		    ".cnf"))
	{
	  sprintf(config_file,
		  "%s/%s",
		  directory_name,
		  child->d_name);
	  is_valid |=HAS_CONFIG_FILE;
	}
      // A matrix filename ends in this
      // context with .mat
      if (ends_with(child->d_name,
		    ".mat"))
	{
	  sprintf(matrix_file,
		  "%s/%s",
		  directory_name,
		  child->d_name);
	  is_valid |= HAS_MATRIX_FILE;
	}
    }
  // if the directory is missing either a
  // config file or a matrix file this is
  // not a valid directory
  if (!(is_valid & HAS_CONFIG_FILE) ||
      !(is_valid & HAS_MATRIX_FILE))
    {
      if (!(is_valid & HAS_CONFIG_FILE))
	{
	  fprintf(stderr,
		  "The input directory is missing a config file\n");
	}
      if (!(is_valid & HAS_MATRIX_FILE))
	{
	  fprintf(stderr,
		  "The input directory is missing a matrix file\n");
	}
      
      return 1;
    }
  return 0;
}




int read_matrix_file(const char *matrix_filename,
		     ASCII_Data *data_file)
{
  size_t allocation_size = 0;
  data_file->ket_configs = NULL;
  data_file->bra_configs = NULL;
  data_file->elements = NULL;
  data_file->num_elements = 0;
  FILE *matrix_file =
    fopen(matrix_filename,"r");

  if (matrix_file == NULL)
    {
      fprintf(stderr,
	      "Could not open file \"%s\". %s\n",
	      matrix_filename,
	      strerror(errno));

      return 1;
    }

  // read line by line of the file
  // and determine if it is a comment
  // or a matrix element
  size_t n = 0;
  char *row = NULL;
  while (!feof(matrix_file))
    {
      if (getline(&row,&n,matrix_file)<0)
	{
	  break;
	}
      if (begins_with(row,"---"))
	{
	  continue;
	}
      size_t b,k;
      double e;
      if (sscanf(row," %ld %ld %lg ",&b,&k,&e))
	{
	  // Since we do not know the necessary allocation size
	  // we double the array sizes when ever we run out
	  // of memory, this means that we only need to do
	  // O(log(n)) memory allocations for n matrix elements
	  if (data_file->num_elements == allocation_size)
	    {
	      allocation_size += allocation_size+1;
	      data_file->ket_configs =
		(size_t*)realloc(data_file->ket_configs,
				 sizeof(size_t)*allocation_size);
	      data_file->bra_configs =
		(size_t*)realloc(data_file->bra_configs,
				 sizeof(size_t)*allocation_size);
	      data_file->elements = 
		(double*)realloc(data_file->elements,
				 sizeof(double)*allocation_size);
	    }
	  data_file->bra_configs[data_file->num_elements] = b;
	  data_file->ket_configs[data_file->num_elements] = k;
	  data_file->elements[data_file->num_elements] = e;
	  data_file->num_elements++;
	}
    }
  // now we know what sizes the arrays need to have
  // but we might have allocated more memory so we
  // trim the arrays so that the data fit exactly
  if (allocation_size > data_file->num_elements)
    {
      
      data_file->ket_configs =
	(size_t*)realloc(data_file->ket_configs,
				 sizeof(size_t)*data_file->num_elements);
      data_file->bra_configs =
	(size_t*)realloc(data_file->bra_configs,
				 sizeof(size_t)*data_file->num_elements);
      data_file->elements = 
	(double*)realloc(data_file->elements,
			 sizeof(double)*data_file->num_elements);

    }
  free(row);
  return 0;
}

// This function takes a peek forward in the file
// to determin if the basis format, if it is jjj
// of JT, after the format is determined, it will
// rewind the file to position at entry
basis_type determine_basis_format(FILE* config)
{
  // Since we want to rewind the file
  // to the starting point at the end
  // we need to get the current position
  unsigned long entyrpos = ftell(config);

  // we are going to read the file row by row
  // using getline, therefore we need to
  // initiate getline's arguments
  char* row = NULL; // getline will allocate the necessary memory
  size_t row_len = 0;

  // skip ahead to the three particle states
  do
    {
      if (getline(&row,&row_len,config) < 0)
	{
	  fprintf(stderr,"No three particle states found\n");
	  exit(1);
	}
    }
  while (row != NULL &&
	 strcmp(row,"---three-particle-states---\n") != 0);
  // The data we are interested in is the row directly after,
  // since it contatains of what type of basis it is
  if (getline(&row,&row_len,config) < 0)
    {
      fprintf(stderr,"File ended directly after three particle states heading\n");
      exit(1);
    }

  // Determining the output format
  basis_type format = NO_BASIS;
  if (strcmp(row,"---i---a---b---c---Jab-Jabc---\n") == 0)
    {
      format =  JJJ_BASIS;
    }
  else if (strcmp(row,"---i---a---b---c---Jab-Jabc-Tab-Tabc---\n") == 0)
    {
      format =  JT_BASIS;
    }
  // Rewind the file to the entry position
  fseek(config,entyrpos,SEEK_SET);

  if (row != NULL)
    {
      free(row);
    }
  
  return format;
}


ASCII_Data* open_ascii_data(const char* directory_name)
{
  // Since the given file_name is not a hdf5 file,
  // we want it to be a directory
  DIR* simple_directory = (DIR*)opendir(directory_name);
  if (simple_directory==NULL)
    {
      return NULL;
    }

  

  // It is however not enough for it to just be a directory
  // the directory must contain a .cnf file and a .mat file
  char config_file[256];
  char matrix_file[256];

  if (scan_directory(directory_name,
		     simple_directory,
		     config_file,
		     matrix_file))
    {
      return NULL;
    }

  // Reading the configuration file:
  // We need to extract first the shells
  // and then the jjj-coupled configurations
  FILE* config = fopen(config_file,
		       "r");
  ASCII_Data* data_file =
    (ASCII_Data*)malloc(sizeof(ASCII_Data));
  
  // Here we read the shells
  data_file->shells =
    new_shells_from_ascii_file(config);

  // Determine if we are dealing with a jjj basis or a jt basis
  data_file->basis_format =
    determine_basis_format(config);
  
  // Here we read the basis,
  // note that these should come after the shells
  if (data_file->basis_format == JJJ_BASIS)
    {
      data_file->basis =
	(void*)new_jjj_basis_from_ascii(config,
					data_file->shells);
    }
  else if (data_file->basis_format == JT_BASIS)
    {
      data_file->basis =
	(void*)new_jt_basis_from_ascii_file(config,
					    data_file->shells);
    }
  
  fclose(config);

  
  if (read_matrix_file(matrix_file,
		       data_file))
    {
      return NULL;
    }

  return data_file;
}


int is_diagonal_jjj(JJJ_Basis* m_basis,
		    JJJ_Basis* n_basis)
{
  if (m_basis->dimension != n_basis->dimension)
    return 0;
  size_t i;
  for (i = 0;
       i<m_basis->dimension;
       i++)
    {
      if (m_basis->states[i].a != n_basis->states[i].a ||
	  m_basis->states[i].b != n_basis->states[i].b ||
	  m_basis->states[i].c != n_basis->states[i].c ||
	  m_basis->states[i].j_ab != n_basis->states[i].j_ab ||
	  m_basis->states[i].j_abc != n_basis->states[i].j_abc ||
	  m_basis->states[i].tz != n_basis->states[i].tz ||
	  m_basis->states[i].parity != n_basis->states[i].parity)
	return 0;
    }
  return 1;
}

int is_diagonal_jt(JT_Basis* m_basis,
		   JT_Basis* n_basis)
{
  if (m_basis->dimension != n_basis->dimension)
    return 0;
  size_t i;
  for (i = 0;
       i<m_basis->dimension;
       i++)
    {
      if (m_basis->states[i].a != n_basis->states[i].a ||
	  m_basis->states[i].b != n_basis->states[i].b ||
	  m_basis->states[i].c != n_basis->states[i].c ||
	  m_basis->states[i].jab != n_basis->states[i].jab ||
	  m_basis->states[i].jabc != n_basis->states[i].jabc ||
	  m_basis->states[i].tab != n_basis->states[i].tab ||
	  m_basis->states[i].tabc != n_basis->states[i].tabc)
	return 0;
    }
  return 1;
}

int sgn(int v)
{
  if (v == 0)
    return 0;
  return v/abs(v);
}


Dens_Matrix* get_matrix_ascii(ASCII_Data* data_file,
			      void* m_basis_p,
			      void* n_basis_p)
{
  if (data_file->basis_format == JJJ_BASIS)
    {
      JJJ_Basis* m_basis = (JJJ_Basis*)m_basis_p;
      JJJ_Basis* n_basis = (JJJ_Basis*)n_basis_p;
      JJJ_Basis* file_basis = data_file->basis;
      // prepare output matrix
      Dens_Matrix *output =
	new_zero_matrix(m_basis->dimension,
			n_basis->dimension);
  
      // match shells
      ssize_t *m_shell_matches =
	matching_shells(data_file->shells,
			m_basis->shells);
      ssize_t *n_shell_matches = m_shell_matches;
      if (m_basis->shells != n_basis->shells)
	{
	  n_shell_matches =
	    matching_shells(data_file->shells,
			    n_basis->shells);
	}
  
      // match m_basis with file_basis,
      // and compute translation array
      ssize_t *m_state_matches =
	matching_jjj_states(file_basis,
			    m_basis,
			    m_shell_matches);
						
      // match n_basis in the exact same way
      ssize_t *n_state_matches = m_state_matches;
  
      if (m_basis != n_basis)
	{
	  n_state_matches =
	    matching_jjj_states(file_basis,
				n_basis,
				n_shell_matches);
	}
      int diagonal_block =
	is_diagonal_jjj(m_basis,
			n_basis);
  
      // loop over all data_file->elements
      // and if they have a corresponding matrix
      // element in the output matrix put it there
      // utilizing that we might be dealing with symmetric blocks
      size_t i;
      for (i = 0;
	   i<data_file->num_elements;
	   i++)
	{
	  if (m_state_matches[data_file->bra_configs[i]] == 0 ||
	      n_state_matches[data_file->ket_configs[i]] == 0)
	    continue;
	  ssize_t index_m =
	    abs(m_state_matches[data_file->bra_configs[i]])-1;
	  ssize_t index_n =
	    abs(n_state_matches[data_file->ket_configs[i]])-1;
	  int phase = sgn(n_state_matches[data_file->ket_configs[i]])*
	    sgn(m_state_matches[data_file->bra_configs[i]]);
	  double element = data_file->elements[i]*phase;
	  
      
	  ELEM(output,index_m,index_n) = element;
	  if (diagonal_block)
	    ELEM(output,index_n,index_m) = element;
	}
      free(m_shell_matches);
      if (m_shell_matches != n_shell_matches)
	{
	  free(n_shell_matches);
	}
      free(m_state_matches);
      if(m_basis != n_basis){
	free(n_state_matches);
      }
  
      return output;
    }
  else if (data_file->basis_format == JT_BASIS)
    {
      JT_Basis* m_basis = m_basis_p;
      JT_Basis* n_basis = n_basis_p;
      JT_Basis* data_basis = data_file->basis;

      // prepare output matrix
      Dens_Matrix *output = 
	new_zero_matrix(m_basis->dimension,
			n_basis->dimension);
      // match m_shells to data_basis->shells
      // the resulting array has one element
      // for each state in data_file->shells
      // that are an index in m_basis->shells
      strue_shell_index *m_shell_matches=
	matching_true_shells(data_file->shells,
			     m_basis->shells);
			     
      // match n_shells to data_basis->shells
      strue_shell_index *n_shell_matches=
	m_shell_matches;
      if (n_basis->shells != m_basis->shells)
	{
	  n_shell_matches = matching_true_shells(data_file->shells,
						 m_basis->shells);
	}
      
      // match m_basis to data_basis
      ssize_t *m_basis_matches =
	matching_jt_states(m_basis,
			   data_basis,
			   m_shell_matches);
      // match n_basis to data_basis
      ssize_t *n_basis_matches =
	m_basis_matches;
      if (n_basis != m_basis)
	{
	  n_basis_matches =
	    matching_jt_states(n_basis,
			       data_basis,
			       m_shell_matches);
	}
      free(m_shell_matches);
      if (m_basis->shells != n_basis->shells)
	free(n_shell_matches);
      int diagonal_block =
	is_diagonal_jt(m_basis,
		       n_basis);
      // loop through the matrix elements in
      // data_file, and use the index maps
      // to get correct indices
      size_t i;
      for (i = 0; i<data_file->num_elements; i++)
	{
	   ssize_t index_m =
	    m_basis_matches[data_file->bra_configs[i]];
	   ssize_t index_n =
	     n_basis_matches[data_file->ket_configs[i]];
	   double element = data_file->elements[i];
	   if (index_m<0 || index_n<0)
	     continue;
	   
	   ELEM(output,index_m,index_n) = element;
	   if (diagonal_block)
	     ELEM(output,index_n,index_m) = element;
	}
      free(m_basis_matches);
      if (m_basis != n_basis)
	free(n_basis_matches);
      return output;
    }
  return NULL;
}


void free_ascii_data(ASCII_Data* data_file)
{
  free_shells(data_file->shells);
  free_jjj_basis(data_file->basis);
  free(data_file->ket_configs);
  free(data_file->bra_configs);
  free(data_file->elements);
  free(data_file);
}
