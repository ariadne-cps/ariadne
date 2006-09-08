#include <iostream>
#include <cassert>

#include <mpfr.h>

using namespace std;

int main() {

  double one;
  double one_approx;
  mpfr_t mpfr_one;
  mpfr_t mpfr_exp_one_lower;
  mpfr_t mpfr_exp_one_upper;
  mpf_t mpf_one;
  mpf_t mpf_one_lower;
  mpf_t mpf_one_upper;
  mpf_t mpf_exp_one_lower;
  mpf_t mpf_exp_one_upper;

  one=1;
  mpfr_init_set_d(mpfr_one,one,GMP_RNDN);
  one_approx=mpfr_get_d(mpfr_one,GMP_RNDN);
  assert(one=one_approx);

  mpf_init_set_d(mpf_one,one);
  mpf_init(mpf_one_lower);
  mpf_init(mpf_one_upper);
  mpfr_get_f(mpf_one_lower,mpfr_one,GMP_RNDD);
  mpfr_get_f(mpf_one_upper,mpfr_one,GMP_RNDD);

  assert(mpf_cmp(mpf_one_lower,mpf_one)<=0);
  assert(mpf_cmp(mpf_one,mpf_one_upper)<=0);

  mpfr_init(mpfr_exp_one_lower);
  mpfr_init(mpfr_exp_one_upper);
  mpfr_exp(mpfr_exp_one_lower,mpfr_one,GMP_RNDD);
  mpfr_exp(mpfr_exp_one_upper,mpfr_one,GMP_RNDU);
  assert(mpfr_cmp(mpfr_exp_one_lower,mpfr_exp_one_upper)<0);

  mpf_init(mpf_exp_one_lower);
  mpf_init(mpf_exp_one_upper);
  mpfr_get_f(mpf_exp_one_lower,mpfr_exp_one_lower,GMP_RNDD);
  mpfr_get_f(mpf_exp_one_upper,mpfr_exp_one_upper,GMP_RNDD);
  
  if(mpf_cmp(mpf_exp_one_lower,mpf_exp_one_upper)>=0) {
    std::cout << "ERROR: Bug in mpfr_get_f\n  " << std::flush;
  }
  assert(mpf_cmp(mpf_exp_one_lower,mpf_exp_one_upper)<0);
  

  return 0;
}
