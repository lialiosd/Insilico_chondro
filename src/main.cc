#include <stdio.h>

#include "mkl.h"

#include "mod_states.h"
// #include "test.h"
// #include "z_update.h"

int sandbox() {
  MKL_INT n=4;
  //float Z[60] = {0.0};
  //float res = z_update(Z, 13);
  // Initialization
  float a[4] = {1.0, 2.0, 3.0, 4.0};
  float b[4] = {1.0001, 2.0001, 3.00003, 4.00004};
  const float *ap;
  const float *bp;
  float *c;
  float *res;

  c = (float*)mkl_calloc(n*n, sizeof(float), 32);
  res = (float*)mkl_calloc(n*n, sizeof(float), 4);

  ap = a;
  bp = b;

  printf("pre\n");
  printf("%f\n", all_close(n, a, b, 1e-3, 1e-5));
  printf("post\n");
  vsAdd(n, a, b, res);

  mkl_free(c);
  mkl_free(res);
  return 0;
}


int main(){
  MKL_INT n = 60;
  float *global_state;
  float *fast_state;
  float *slow_state;

  global_state = (float *) mkl_calloc(n, sizeof(float), 64);
  fast_state   = (float *) mkl_calloc(n, sizeof(float), 64);
  slow_state   = (float *) mkl_calloc(n, sizeof(float), 64);
  

  VSLStreamStatePtr stream;
  vslNewStream(&stream, VSL_BRNG_SFMT19937, 111);
  vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n, global_state, 0.0, 1.0);

  for (unsigned i=0; i<n; i++){
    printf("%f ", global_state[i]);
  }

  printf("\n \n");

  // Clean memory / streams
  mkl_free(global_state);
  mkl_free(fast_state);
  mkl_free(slow_state);
  // vslDeleteStream(&stream);

  float ran[n];

  vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n, ran, 0.0, 1.0);

  for (unsigned i=0; i<n; i++){
    printf("%f ", global_state[i]);
  }
  // Clean memory / streams
  vslDeleteStream(&stream);
  return 0;
}
