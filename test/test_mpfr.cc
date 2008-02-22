#include <iostream>
#include <cassert>

#include <mpfr.h>
#include "test/test.h"

using namespace std;

class TestMpfr {
  double one;
  mpf_t mpf_one;
  mpfr_t mpfr_one;
 
 public:
  TestMpfr() {
		one = 1.0;
		mpfr_init_set_d(mpfr_one,one,GMP_RNDN);
		mpf_init_set_d(mpf_one,one);		
	}
	
	void test_approx() {
		double one_approx;
		
		one_approx=mpfr_get_d(mpfr_one,GMP_RNDN);
		ARIADNE_TEST_ASSERT(one==one_approx);	
	}
	
	void test_rounding() {
		mpf_t mpf_one_lower, mpf_one_upper;
		
		mpf_init(mpf_one_lower);
		mpf_init(mpf_one_upper);
		mpfr_get_f(mpf_one_lower,mpfr_one,GMP_RNDD);
		mpfr_get_f(mpf_one_upper,mpfr_one,GMP_RNDD);

		ARIADNE_TEST_ASSERT(mpf_cmp(mpf_one_lower,mpf_one)<=0);
		ARIADNE_TEST_ASSERT(mpf_cmp(mpf_one,mpf_one_upper)<=0);
	}
	
	void test_exp() {
		mpfr_t mpfr_exp_one_lower, mpfr_exp_one_upper;
		
		mpfr_init(mpfr_exp_one_lower);
		mpfr_init(mpfr_exp_one_upper);
		mpfr_exp(mpfr_exp_one_lower,mpfr_one,GMP_RNDD);
		mpfr_exp(mpfr_exp_one_upper,mpfr_one,GMP_RNDU);
		ARIADNE_TEST_ASSERT(mpfr_cmp(mpfr_exp_one_lower,mpfr_exp_one_upper)<0);	
	}
	
  // REGRESSION TEST: Test for bug in mpfr_get_f
	void test_mpfr_get_f() {
		mpfr_t mpfr_exp_one_lower, mpfr_exp_one_upper;
		mpf_t mpf_exp_one_lower, mpf_exp_one_upper;

		mpfr_init(mpfr_exp_one_lower);
		mpfr_init(mpfr_exp_one_upper);
		mpfr_exp(mpfr_exp_one_lower,mpfr_one,GMP_RNDD);
		mpfr_exp(mpfr_exp_one_upper,mpfr_one,GMP_RNDU);
				
		mpf_init(mpf_exp_one_lower);
		mpf_init(mpf_exp_one_upper);
		mpfr_get_f(mpf_exp_one_lower,mpfr_exp_one_lower,GMP_RNDD);
		mpfr_get_f(mpf_exp_one_upper,mpfr_exp_one_upper,GMP_RNDD);
		
		if(mpf_cmp(mpf_exp_one_lower,mpf_exp_one_upper)>=0) {
			std::cout << "ERROR: Bug in mpfr_get_f\n  " << std::flush;
		}
		ARIADNE_TEST_ASSERT(mpf_cmp(mpf_exp_one_lower,mpf_exp_one_upper)<0);	
	}
	
	void test() {
		ARIADNE_TEST_CALL(test_approx());
		ARIADNE_TEST_CALL(test_rounding());
		ARIADNE_TEST_CALL(test_exp());
		ARIADNE_TEST_CALL(test_mpfr_get_f());
	}
	
};

int main() {
  TestMpfr().test();
  return ARIADNE_TEST_FAILURES;
}
