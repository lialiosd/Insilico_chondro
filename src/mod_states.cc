#include "mod_states.h"

// int check_states(float* state1,float* state2, const int max_id){
  
//   return 0;
// }

float* all_close(MKL_INT n, float* v1, float* v2, double rtol, double atol){
  /*
    Returns TRUE if the 2 arrays are close within tolerance.
    The following expression is evaluated:
    abs(v1 - v2) <= (atol + rtol * abs(v2))
         p1                          p2
  */

  // Initialize pointers
  float *p1, *p2;
  // Allocate the memory for the intermediate steps
  p1 = (float*)mkl_calloc(n, sizeof(float), 64);
  p2 = (float*)mkl_calloc(n, sizeof(float), 64);

  vsSub(n, v1, v2, p1);
  vsAbs(n, p1, p1);
  for (int i; i<n; i++){
    printf("%f", p1[i]);
  }
  vsAbs(n, v2, p2);

  // Free memory
  mkl_free(p1);
  mkl_free(p2);

  return 0;

}
