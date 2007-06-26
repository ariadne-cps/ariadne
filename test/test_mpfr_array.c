#include <stdio.h>
#include <assert.h>
#include <mpfr.h>
#include "mpfr_array.h"

int main() {
  printf("test_mpfr_array\n");

  mpfr_ptr x;
  char xstr[100];


  x = mpfr_array_alloc_init(1);
  mpfr_array_clear_free(x,1);

  x = mpfr_array_alloc_init2(1,2);
  mpfr_array_clear_free(x,1);

  x = mpfr_array_alloc_init2(1,53);
  mpfr_array_clear_free(x,1);

  x = mpfr_array_alloc_init2(1,65);
  mpfr_array_clear_free(x,1);

  printf("\n");

  x = mpfr_array_alloc_init(2);
  printf("%f\n",mpfr_get_d(x+0,GMP_RNDN));
  printf("%f\n",mpfr_get_d(x+1,GMP_RNDN));
  mpfr_array_clear_free(x,2);
  
  printf("\n");

  x = mpfr_array_alloc_init(3);
  printf("%f\n",mpfr_get_d(x+2,GMP_RNDN));
  mpfr_set_ui(x,3,GMP_RNDN);
  printf("%f\n",mpfr_get_d(x,GMP_RNDN));
  mpfr_set_d(x+1,2.3,GMP_RNDN);
  printf("%f\n",mpfr_get_d(x+1,GMP_RNDN));
  mpfr_set_d(x+2,-4.2,GMP_RNDN);
  printf("%f\n",mpfr_get_d(x+2,GMP_RNDN));
  mpfr_array_clear_free(x,3);

  printf("\n");
 
  x = mpfr_array_alloc_init(2);
  mpfr_set_ui(x,0,GMP_RNDN);
  printf("%f\n",mpfr_get_d(x+0,GMP_RNDN));
  printf("%f\n",mpfr_get_d(x+1,GMP_RNDN));
  mpfr_array_clear_free(x,2);

  x = mpfr_array_alloc_init2(2,53);
  mpfr_set_ui(x,5,GMP_RNDN);
  printf("%f\n",mpfr_get_d(x+0,GMP_RNDN));
  printf("%f\n",mpfr_get_d(x+1,GMP_RNDN));
  mpfr_array_clear_free(x,2);
  
  printf("\n");

  x = mpfr_array_alloc_init2(2,53);
  mpfr_set_ui(x+1,1,GMP_RNDN);
  printf("%f\n",mpfr_get_d(x+0,GMP_RNDN));
  printf("%f\n",mpfr_get_d(x+1,GMP_RNDN));
  mpfr_array_clear_free(x,2);

  printf("\n");

  x = mpfr_array_alloc_init2(2,65);
  mpfr_set_ui(x+0,1,GMP_RNDN);
  mpfr_set_ui(x+1,2,GMP_RNDN);
  printf("%f\n",mpfr_get_d(x+0,GMP_RNDN));
  printf("%f\n",mpfr_get_d(x+1,GMP_RNDN));
  mpfr_array_clear_free(x,2);

  printf("\n\n");

  x = mpfr_array_alloc(7);
  mpfr_array_init(x,2);
  mpfr_set_ui(x+0,0,GMP_RNDN);
  printf("%f\n",mpfr_get_d(x+0,GMP_RNDN));
  printf("%f\n",mpfr_get_d(x+1,GMP_RNDN));
  
  printf("\n");
 
  mpfr_array_set_prec(x,2,128);
  printf("%f\n",mpfr_get_d(x+0,GMP_RNDN));
  mpfr_set_ui(x+0,1,GMP_RNDN);
  printf("%f\n",mpfr_get_d(x+0,GMP_RNDN));
  mpfr_array_clear(x,2);
  mpfr_array_free(x,7);
 


}
