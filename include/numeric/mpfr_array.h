/* mpfr_array.h -- Include file for mpfr arrays.

Copyright 2007 Pieter Collins, CWI

This file is intended to become part of the MPFR Library.

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

#ifndef _MPFR_ARRAY_H
#define _MPFR_ARRAY_H

#include "gmp.h"
#include "mpfr.h"

_MPFR_DECLSPEC mpfr_ptr mpfr_array_alloc_init _MPFR_PROTO ((size_t));
_MPFR_DECLSPEC mpfr_ptr mpfr_array_alloc_init2 _MPFR_PROTO ((size_t, mpfr_prec_t));
_MPFR_DECLSPEC mpfr_ptr mpfr_array_realloc_set_prec _MPFR_PROTO ((mpfr_ptr, size_t, mpfr_prec_t));
_MPFR_DECLSPEC void mpfr_array_clear_free _MPFR_PROTO ((mpfr_ptr, size_t));

_MPFR_DECLSPEC mpfr_ptr mpfr_array_alloc _MPFR_PROTO ((size_t));
_MPFR_DECLSPEC void mpfr_array_init _MPFR_PROTO ((mpfr_ptr, size_t));
_MPFR_DECLSPEC void mpfr_array_init2 _MPFR_PROTO ((mpfr_ptr, size_t, mpfr_prec_t));
_MPFR_DECLSPEC void mpfr_array_set_prec _MPFR_PROTO ((mpfr_ptr, size_t, mpfr_prec_t));
_MPFR_DECLSPEC mpfr_prec_t mpfr_array_get_prec _MPFR_PROTO ((mpfr_srcptr, size_t));
_MPFR_DECLSPEC void mpfr_array_clear _MPFR_PROTO ((mpfr_ptr, size_t));
_MPFR_DECLSPEC void mpfr_array_free _MPFR_PROTO ((mpfr_ptr, size_t));

_MPFR_DECLSPEC void mpfr_array_init_zero2 _MPFR_PROTO ((mpfr_ptr, size_t, mpfr_prec_t));

_MPFR_DECLSPEC void mpfr_array_copy _MPFR_PROTO ((mpfr_ptr, size_t, mpfr_srcptr, mp_rnd_t));
_MPFR_DECLSPEC void mpfr_array_set _MPFR_PROTO ((mpfr_ptr, size_t, mpfr_t, mp_rnd_t));

_MPFR_DECLSPEC void mpfr_array_set_ui _MPFR_PROTO ((mpfr_ptr, size_t, unsigned long int, mp_rnd_t));
_MPFR_DECLSPEC void mpfr_array_set_si _MPFR_PROTO ((mpfr_ptr, size_t, long int, mp_rnd_t));
#ifdef _MPFR_H_HAVE_INTMAX_T
_MPFR_DECLSPEC void mpfr_array_set_uj _MPFR_PROTO ((mpfr_ptr, size_t, uintmax_t, mp_rnd_t));
_MPFR_DECLSPEC void mpfr_array_set_sj _MPFR_PROTO ((mpfr_ptr, size_t, intmax_t, mp_rnd_t));
#endif
_MPFR_DECLSPEC void mpfr_array_set_d _MPFR_PROTO ((mpfr_ptr, size_t, double, mp_rnd_t));
_MPFR_DECLSPEC void mpfr_array_set_ld _MPFR_PROTO ((mpfr_ptr, size_t, long double, mp_rnd_t));
_MPFR_DECLSPEC void mpfr_array_set_z _MPFR_PROTO ((mpfr_ptr, size_t, mpz_t, mp_rnd_t));
_MPFR_DECLSPEC void mpfr_array_set_q _MPFR_PROTO ((mpfr_ptr, size_t, mpq_t, mp_rnd_t));
_MPFR_DECLSPEC void mpfr_array_set_f _MPFR_PROTO ((mpfr_ptr, size_t, mpf_t, mp_rnd_t));

#endif /* _MPFR_ARRAY_H */
