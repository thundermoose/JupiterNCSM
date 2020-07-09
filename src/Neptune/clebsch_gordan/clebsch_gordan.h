#ifndef __CLEBSCH_GORDAN__
#define __CLEBSCH_GORDAN__
#include <stdlib.h>
#include <wigxjpf.h>

typedef struct _clebsch_gordan_data_{
  struct wigxjpf_temp** wig_data;
  size_t max_num_threads;
} Clebsch_Gordan_Data; // Should be shared between threads

Clebsch_Gordan_Data* initiate_clebsch_gordan(int max_j); // Should be runed once, inside omp

double clebsch_gordan(int j1,int j2,int j3,
		      int m1,int m2,int m3,
		      Clebsch_Gordan_Data* cgd);
void free_clebsch_gordan(Clebsch_Gordan_Data* cgd);
#endif
