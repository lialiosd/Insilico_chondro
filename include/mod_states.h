#ifndef MOD_STATES_H
#define MOD_STATES_H

  #include<iostream>
  #include "mkl.h"

  // int check_states(float* state1,float* state2, const int max_id);
  float* all_close(MKL_INT n, float* v1, float* v2, float rtol = 1e-5, float atol=1e-8);

#endif
