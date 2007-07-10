/* mpfr_init -- initialize a floating-point number

Copyright 1999, 2001, 2002, 2004, 2006 Free Software Foundation, Inc.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdlib.h>
#include <stdio.h>
#include "mpfr.h"
#include "mpfr_array.h"
#include "mpfr-impl.h"

#define MPFR_ARRAY_MALLOC_SIZE(n,s)                                     \
  ( n * (sizeof(mpfr_t) + sizeof(mpfr_size_limb_t) + BYTES_PER_MP_LIMB * ((size_t) s) ) )


mpfr_ptr
mpfr_array_alloc_init (size_t n)
{
  return mpfr_array_alloc_init2 (n, __gmpfr_default_fp_bit_precision);
}

mpfr_ptr
mpfr_array_alloc_init2 (size_t n, mp_prec_t p)
{
  mpfr_ptr a;
  mpfr_ptr x;
  /* mpfr_t x; */
  mp_size_t limbsize;
  mp_ptr limbs;
  mp_ptr limb;
  size_t i;

  /* Check if we can represent the number of limbs
   * associated to the maximum of mpfr_prec_t*/
  MPFR_ASSERTN( MP_SIZE_T_MAX >= (MPFR_PREC_MAX/BYTES_PER_MP_LIMB) );

  /* Check for correct BITS_PER_MP_LIMB and BYTES_PER_MP_LIMB */
  MPFR_ASSERTN( BITS_PER_MP_LIMB == BYTES_PER_MP_LIMB * CHAR_BIT
                && sizeof(mp_limb_t) == BYTES_PER_MP_LIMB );

  /* Check for correct EXP NAN, ZERO & INF in both mpfr.h and in mpfr-impl.h */
  MPFR_ASSERTN( __MPFR_EXP_NAN  == MPFR_EXP_NAN  );
  MPFR_ASSERTN( __MPFR_EXP_ZERO == MPFR_EXP_ZERO );
  MPFR_ASSERTN( __MPFR_EXP_INF  == MPFR_EXP_INF  );

  MPFR_ASSERTN( MPFR_EMAX_MAX <= (MPFR_EXP_MAX >> 1)  );
  MPFR_ASSERTN( MPFR_EMIN_MIN >= -(MPFR_EXP_MAX >> 1) );

  /* p=1 is not allowed since the rounding to nearest even rule requires at
     least two bits of mantissa: the neighbours of 3/2 are 1*2^0 and 1*2^1,
     which both have an odd mantissa */
  MPFR_ASSERTN(p >= MPFR_PREC_MIN && p <= MPFR_PREC_MAX);

  limbsize = (mp_size_t) ((p - 1) / BITS_PER_MP_LIMB) + 1;
  a   = (mpfr_ptr) (*__gmp_allocate_func)(MPFR_ARRAY_MALLOC_SIZE(n, limbsize)); 
  limbs = (mp_ptr) (a + n);
  
  /*
    printf("a=%p, l=%p, d=%i, n=%i, p=%i, ls=%i, s=%i\n",a,limbs,(char*)limbs-(char*)a,n,p,limbsize,MPFR_ARRAY_MALLOC_SIZE(n, limbsize));
  */

  for(i=0; i!=n; ++i) {
    x=a+i;
    limb=limbs+i*limbsize;
    MPFR_PREC(x) = p;                   /* Set prec */
    MPFR_EXP (x) = MPFR_EXP_INVALID;    /* make sure that the exp field has a
                                           valid value in the C point of view */
    MPFR_SET_POS(x);                    /* Set a sign */
    MPFR_SET_MANT_PTR(x, limb);         /* Set Mantissa ptr */
    MPFR_SET_ALLOC_SIZE(x, limbsize);   /* Fix alloc size of Mantissa */
    MPFR_SET_NAN(x);                    /* initializes to NaN */
  }

  return a;
}



mpfr_ptr
mpfr_array_realloc_set_prec (mpfr_ptr a, size_t n, mpfr_prec_t p)
{
  mpfr_array_clear_free(a,n);
  return mpfr_array_alloc_init2(n,p);
}


void
mpfr_array_clear_free (mpfr_ptr a, size_t n)
{
  /*
    printf("f=%p\n",a);
  */

  free(a);
  /*
  (*__gmp_free_func) (a,
                      MPFR_ARRAY_MALLOC_SIZE (n,MPFR_GET_ALLOC_SIZE (a)));
  */
}




mpfr_ptr
mpfr_array_alloc (size_t n)
{
  return (mpfr_ptr) (*__gmp_allocate_func) (n*sizeof(mpfr_t));
}

void
mpfr_array_init (mpfr_ptr a, size_t n)
{
  mpfr_array_init2 (a, n, __gmpfr_default_fp_bit_precision);
}

void
mpfr_array_init2 (mpfr_ptr a, size_t n, mp_prec_t p)
{
  mpfr_ptr x;
  mp_size_t limbsize;
  mp_ptr limbs;
  mp_ptr limb;
  size_t i;

  /* Check if we can represent the number of limbs
   * associated to the maximum of mpfr_prec_t*/
  MPFR_ASSERTN( MP_SIZE_T_MAX >= (MPFR_PREC_MAX/BYTES_PER_MP_LIMB) );

  /* Check for correct BITS_PER_MP_LIMB and BYTES_PER_MP_LIMB */
  MPFR_ASSERTN( BITS_PER_MP_LIMB == BYTES_PER_MP_LIMB * CHAR_BIT
                && sizeof(mp_limb_t) == BYTES_PER_MP_LIMB );

  /* Check for correct EXP NAN, ZERO & INF in both mpfr.h and in mpfr-impl.h */
  MPFR_ASSERTN( __MPFR_EXP_NAN  == MPFR_EXP_NAN  );
  MPFR_ASSERTN( __MPFR_EXP_ZERO == MPFR_EXP_ZERO );
  MPFR_ASSERTN( __MPFR_EXP_INF  == MPFR_EXP_INF  );

  MPFR_ASSERTN( MPFR_EMAX_MAX <= (MPFR_EXP_MAX >> 1)  );
  MPFR_ASSERTN( MPFR_EMIN_MIN >= -(MPFR_EXP_MAX >> 1) );

  /* p=1 is not allowed since the rounding to nearest even rule requires at
     least two bits of mantissa: the neighbours of 3/2 are 1*2^0 and 1*2^1,
     which both have an odd mantissa */
  MPFR_ASSERTN(p >= MPFR_PREC_MIN && p <= MPFR_PREC_MAX);

  limbsize = (mp_size_t) ((p - 1) / BITS_PER_MP_LIMB) + 1;
  limbs   = (mp_ptr) (*__gmp_allocate_func)(n*MPFR_MALLOC_SIZE(limbsize));

  for(i=0; i!=n; ++i) {
    x=a+i;
    limb=limbs+i*limbsize;
    MPFR_PREC(x) = p;                   /* Set prec */
    MPFR_EXP (x) = MPFR_EXP_INVALID;    /* make sure that the exp field has a
                                           valid value in the C point of view */
    MPFR_SET_POS(x);                    /* Set a sign */
    MPFR_SET_MANT_PTR(x, limb);         /* Set Mantissa ptr */
    MPFR_SET_ALLOC_SIZE(x, limbsize);   /* Fix alloc size of Mantissa */
    MPFR_SET_NAN(x);                    /* initializes to NaN */
  }
}

void
mpfr_array_set_prec (mpfr_ptr a, size_t n, mpfr_prec_t p)
{
  mp_size_t limbsize, limboldsize;
  mp_ptr limbs;
  mpfr_ptr x;
  size_t i;

  /* first, check if p is correct */
  MPFR_ASSERTN (p >= MPFR_PREC_MIN && p <= MPFR_PREC_MAX);

  /* Calculate the new number of limbs */
  limbsize = (p - 1) / BITS_PER_MP_LIMB + 1;

  /* Realloc only if the new size is not equal to the old */
  /* We realloc even if new size is smaller to save memory of large arrays */
  limboldsize = MPFR_GET_ALLOC_SIZE (a);
  if (limbsize != limboldsize)
    {
      limbs = (mp_ptr) (*__gmp_reallocate_func)
        (MPFR_GET_REAL_PTR(a), n*MPFR_MALLOC_SIZE(limboldsize), n*MPFR_MALLOC_SIZE(limbsize));
      for(i=0; i!=n; ++i) 
        {
          x=a+i;
          MPFR_SET_MANT_PTR(x, limbs+i*limbsize);
          MPFR_SET_ALLOC_SIZE(x, limbsize);
          MPFR_PREC (x) = p;
          MPFR_SET_NAN (x); /* initializes to NaN */
        }
    }
}

mpfr_prec_t
mpfr_array_get_prec (mpfr_srcptr a, size_t n)
{
  return MPFR_PREC(a);
}

void
mpfr_array_clear (mpfr_ptr m, size_t n)
{
  size_t i;
  (*__gmp_free_func) (MPFR_GET_REAL_PTR (m),
                      n*MPFR_MALLOC_SIZE (MPFR_GET_ALLOC_SIZE (m)));
  for(i=0; i!=n; ++i) {
    MPFR_MANT (m+i) = (mp_limb_t *) 0;
  }
}

void
mpfr_array_free (mpfr_ptr a, size_t n)
{
  (*__gmp_free_func) (a,n*sizeof(mpfr_t));
}



void
mpfr_array_init_zero (mpfr_ptr a, size_t n)
{
  mpfr_array_init(a,n);
  size_t i;
  for(i=0; i!=n; ++i) 
    {
      mpfr_set_si(a+i,0,GMP_RNDN);
    }
}


void
mpfr_array_init_zero2 (mpfr_ptr a, size_t n, mpfr_prec_t p)
{
  mpfr_array_init2(a,n,p);
  size_t i;
  for(i=0; i!=n; ++i) 
    {
      mpfr_set_si(a+i,0,GMP_RNDN);
    }
}



void
mpfr_array_copy (mpfr_ptr a, size_t n, mpfr_srcptr x, mp_rnd_t rnd_mode)
{
  size_t i;
  for(i=0; i!=n; ++i) 
    {
      mpfr_set(a+i,x+i,rnd_mode);
    }
}


void
mpfr_array_set (mpfr_ptr a, size_t n, mpfr_t x, mp_rnd_t rnd_mode)
{
  size_t i;
  for(i=0; i!=n; ++i) 
    {
      mpfr_set(a+i,x,rnd_mode);
    }
}



void 
mpfr_array_set_ui (mpfr_ptr rop, size_t n, unsigned long int op, mp_rnd_t rnd)
{
  size_t i;
  for(i=0; i!=n; ++i) 
    {
      mpfr_set_ui(rop+i,op,rnd);
    }
}

void 
mpfr_array_set_si (mpfr_ptr rop, size_t n, long int op, mp_rnd_t rnd)
{
  size_t i;
  for(i=0; i!=n; ++i) 
    {
      mpfr_set_si(rop+i,op,rnd);
    }
}


#ifdef _MPFR_H_HAVE_INTMAX_T

void mpfr_array_set_uj (mpfr_ptr rop, size_t n, uintmax_t op, mp_rnd_t rnd)
{
  size_t i;
  for(i=0; i!=n; ++i) 
    {
      mpfr_set_uj(rop+i,op,rnd);
    }
}

void mpfr_array_set_sj (mpfr_ptr rop, size_t n, intmax_t op, mp_rnd_t rnd)
{
  size_t i;
  for(i=0; i!=n; ++i) 
    {
      mpfr_set_sj(rop+i,op,rnd);
    }
}

#endif


void mpfr_array_set_d (mpfr_ptr rop, size_t n, double op, mp_rnd_t rnd)
{
  size_t i;
  for(i=0; i!=n; ++i) 
    {
      mpfr_set_d(rop+i,op,rnd);
    }
}


void mpfr_array_set_ld (mpfr_t rop, size_t n, long double op, mp_rnd_t rnd)
{
  size_t i;
  for(i=0; i!=n; ++i) 
    {
      mpfr_set_ld(rop+i,op,rnd);
    }
}


void mpfr_array_set_z (mpfr_t rop, size_t n, mpz_t op, mp_rnd_t rnd)
{
  size_t i;
  for(i=0; i!=n; ++i) 
    {
      mpfr_set_z(rop+i,op,rnd);
    }
}


void mpfr_array_set_q (mpfr_t rop, size_t n, mpq_t op, mp_rnd_t rnd)
{
  size_t i;
  for(i=0; i!=n; ++i) 
    {
      mpfr_set_q(rop+i,op,rnd);
    }
}


void mpfr_array_set_f (mpfr_t rop, size_t n, mpf_t op, mp_rnd_t rnd)
{
  size_t i;
  for(i=0; i!=n; ++i) 
    {
      mpfr_set_f(rop+i,op,rnd);
    }
}




