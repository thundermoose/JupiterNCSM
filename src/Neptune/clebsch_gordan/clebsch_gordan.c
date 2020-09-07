#include "clebsch_gordan.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <debug_mode/debug_mode.h>

Clebsch_Gordan_Data* initiate_clebsch_gordan(int max_j)
{
  Clebsch_Gordan_Data* out;
#pragma omp single copyprivate(out)
  {
    out = (Clebsch_Gordan_Data*)malloc(sizeof(Clebsch_Gordan_Data));
    out->max_num_threads = omp_get_num_threads();
    out->wig_data =
      (struct wigxjpf_temp**)malloc(sizeof(struct wigxjpf_temp*)*out->max_num_threads);
    wigxjpf_fill_factors(3*max_j+1);
  }
#pragma omp barrier
  out->wig_data[omp_get_thread_num()] = wigxjpf_temp_alloc(max_j+1);
#pragma omp barrier
  if (out->wig_data[omp_get_thread_num()] == NULL)
    {
      fprintf(stderr,"thread: %d could not initiate wigxjpf\n",omp_get_thread_num());
      exit(1);
    }
  
  return out;
}

double clebsch_gordan(int j1,int j2,int j3,
		      int m1,int m2,int m3,
		      Clebsch_Gordan_Data* cgd)
{
  double result;
  int num_thread = omp_get_thread_num();
  calc_3j_double(&result,
		 j1,j2,j3,
		 m1,m2,-m3,
		 cgd->wig_data[num_thread]);
  return result*sqrt(j3+1)*(1-((j1-j2+m3)&2)); 
}

void free_clebsch_gordan(Clebsch_Gordan_Data* cgd)
{
  wigxjpf_temp_free(cgd->wig_data[omp_get_thread_num()]);
#pragma omp barrier
#pragma omp single copyprivate(cgd)
  {
    free(cgd->wig_data);
    free(cgd);
  }
}
