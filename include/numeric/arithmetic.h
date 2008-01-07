/***************************************************************************
 *            arithmetic.h
 *
 *  Copyright 2007  Alberto Casagrande, Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
/*! \file numeric/arithmetic.h
 *  \brief Arithmetic operators and rounded arithmetic.
 */
 
#ifndef ARIADNE_NUMERIC_ARITHMETIC_H
#define ARIADNE_NUMERIC_ARITHMETIC_H

#include "numeric/rounding.h"

namespace Ariadne {
namespace Numeric {

template<class T> class Float;

template<class T>
inline int floor(const Float<T>& x) { 
  int r; floor_(r,x); return r; }

template<class T>
inline int ceil(const Float<T>& x) { 
  int r; ceil_(r,x); return r; }


template<class T> inline
void
decr(Float<T>& x) {
  next_(x,x,round_up);
}

template<class T> inline
void
incr(Float<T>& x) {
  next_(x,x,round_up);
}

template<class T> inline
Float<T>
prec(const Float<T>& x) {
  return next(x,round_down);
}

template<class T> inline
Float<T>
succ(const Float<T>& x) {
  return next(x,round_up);
}

template<class T>
inline Float<T> next_up(const Float<T>& x) { 
  Float<T> r; next_(r,x,round_up); return r; }
template<class T>
inline Float<T> next_down(const Float<T>& x) { 
  Float<T> r; next_(r,x,round_down); return r; }


template<class T> inline
Float<T>
abs(const Float<T>& x) {
  Float<T> r; abs_(r,x); return r;
}



template<class T> 
inline Float<T> add_up(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; add_(r,x,y,round_up); return r; }
template<class T> 
inline Float<T> add_down(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; add_(r,x,y,round_down); return r; }
template<class T> 
inline Float<T> add_approx(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; add_(r,x,y,round_approx); return r; }

template<class T> 
inline Float<T> sub_up(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; sub_(r,x,y,round_up); return r; }
template<class T> 
inline Float<T> sub_down(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; sub_(r,x,y,round_down); return r; }
template<class T> 
inline Float<T> sub_approx(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; sub_(r,x,y,round_approx); return r; }

template<class T> 
inline Float<T> mul_up(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; mul_(r,x,y,round_up); return r; }
template<class T> 
inline Float<T> mul_down(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; mul_(r,x,y,round_down); return r; }
template<class T> 
inline Float<T> mul_approx(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; mul_(r,x,y,round_approx); return r; }

template<class T> 
inline Float<T> mul_up(const Float<T>& x, const int& y) { 
  Float<T> r; mul_(r,x,y,round_up); return r; }
template<class T> 
inline Float<T> mul_down(const Float<T>& x, const int& y) { 
  Float<T> r; mul_(r,x,y,round_down); return r; }
template<class T> 
inline Float<T> mul_approx(const Float<T>& x, const int& y) { 
  Float<T> r; mul_(r,x,y,round_approx); return r; }

template<class T> 
inline Float<T> mul_approx(const Float<T>& x, const double& y) { 
  Float<T> r; mul_(r,x,y,round_approx); return r; }
template<class T> 
inline Float<T> mul_approx(const double& x, const Float<T>& y) { 
  Float<T> r; mul_(r,x,y,round_approx); return r; }


template<class T> 
inline Float<T> div_up(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; div_(r,x,y,round_up); return r; }
template<class T> 
inline Float<T> div_down(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; div_(r,x,y,round_down); return r; }
template<class T> 
inline Float<T> div_approx(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; div_(r,x,y,round_approx); return r; }

template<class T> 
inline Float<T> div_up(const Float<T>& x, const int& y) { 
  Float<T> r; div_(r,x,y,round_up); return r; }
template<class T> 
inline Float<T> div_down(const Float<T>& x, const int& y) { 
  Float<T> r; div_(r,x,y,round_down); return r; }
template<class T> 
inline Float<T> div_approx(const Float<T>& x, const int& y) { 
  Float<T> r; div_(r,x,y,round_approx); return r; }
template<class T> 
inline Float<T> div_approx(const Float<T>& x, const double& y) { 
  Float<T> r; div_(r,x,y,round_approx); return r; }

template<class T> 
inline Float<T> med_approx(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; med_(r,x,y,round_approx); return r; }
template<class T> 
inline Float<T> rad_up(const Float<T>& x, const Float<T>& y) { 
  Float<T> r; rad_(r,x,y,round_up); return r; }

template<class T, class N> 
inline Float<T> pow_up(const Float<T>& x, const N& n) { 
  Float<T> r; pow_(r,x,n,round_up); return r; }
template<class T, class N> 
inline Float<T> pow_down(const Float<T>& x, const N& n) { 
  Float<T> r; pow_(r,x,n,round_down); return r; }
template<class T, class N> 
inline Float<T> pow_approx(const Float<T>& x, const N& n) { 
  Float<T> r; pow_(r,x,n,round_approx); return r; }




template<class T> 
inline Float<T> sqrt_up(const Float<T>& x) { 
  Float<T> r; sqrt_(r,x,round_up); return r; }
template<class T> 
inline Float<T> sqrt_down(const Float<T>& x) { 
  Float<T> r; sqrt_(r,x,round_down); return r; }
template<class T> 
inline Float<T> sqrt_approx(const Float<T>& x) { 
  Float<T> r; sqrt_(r,x,round_approx); return r; }

template<class T> 
inline Float<T> exp_up(const Float<T>& x) { 
  Float<T> r; exp_(r,x,round_up); return r; }
template<class T> 
inline Float<T> exp_down(const Float<T>& x) { 
  Float<T> r; exp_(r,x,round_down); return r; }
template<class T> 
inline Float<T> exp_approx(const Float<T>& x) { 
  Float<T> r; exp_(r,x,round_approx); return r; }

template<class T> 
inline Float<T> log_up(const Float<T>& x) { 
  Float<T> r; log_(r,x,round_up); return r; }
template<class T> 
inline Float<T> log_down(const Float<T>& x) { 
  Float<T> r; log_(r,x,round_down); return r; }
template<class T> 
inline Float<T> log_approx(const Float<T>& x) { 
  Float<T> r; log_(r,x,round_approx); return r; }


} // namespace Ariadne
} // namespace Numeric



#endif /* ARIADNE_ARITHMETIC_H */
