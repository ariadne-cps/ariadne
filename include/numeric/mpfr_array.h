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

#ifndef __MPFR_ARRAY_H
#define __MPFR_ARRAY_H

#include "mpfr.h"

__MPFR_DECLSPEC mpfr_ptr mpfr_array_alloc_init _MPFR_PROTO ((size_t));
__MPFR_DECLSPEC mpfr_ptr mpfr_array_alloc_init2 _MPFR_PROTO ((size_t, mpfr_prec_t));
__MPFR_DECLSPEC void mpfr_array_clear_free _MPFR_PROTO ((mpfr_ptr, size_t));

__MPFR_DECLSPEC mpfr_ptr mpfr_array_alloc _MPFR_PROTO ((size_t));
__MPFR_DECLSPEC void mpfr_array_init _MPFR_PROTO ((mpfr_ptr, size_t));
__MPFR_DECLSPEC void mpfr_array_init2 _MPFR_PROTO ((mpfr_ptr, size_t, mpfr_prec_t));
__MPFR_DECLSPEC void mpfr_array_set_prec _MPFR_PROTO ((mpfr_ptr, size_t, mpfr_prec_t));
__MPFR_DECLSPEC mpfr_prec_t mpfr_array_get_prec _MPFR_PROTO ((mpfr_srcptr, size_t));
__MPFR_DECLSPEC void mpfr_array_set _MPFR_PROTO ((mpfr_ptr, size_t, mpfr_t, mp_rnd_t));
__MPFR_DECLSPEC void mpfr_array_set_si _MPFR_PROTO ((mpfr_ptr, size_t, long, mp_rnd_t));
__MPFR_DECLSPEC void mpfr_array_set_d _MPFR_PROTO ((mpfr_ptr, size_t, double, mp_rnd_t));
__MPFR_DECLSPEC void mpfr_array_clear _MPFR_PROTO ((mpfr_ptr, size_t));
__MPFR_DECLSPEC void mpfr_array_free _MPFR_PROTO ((mpfr_ptr, size_t));

#endif /* __MPFR_ARRAY_H */
