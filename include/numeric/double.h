/***************************************************************************
 *            numeric/double.h
 *
 *  Copyright  2004-7  Pieter Collins
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
 
/*! \file numeric/double.h
 *  \brief Rounded arithmetic for double precision
 */

#ifndef ARIADNE_NUMERIC_DOUBLE_H
#define ARIADNE_NUMERIC_DOUBLE_H

namespace Ariadne {
namespace Numeric {

typedef unsigned short rounding_mode;

void incr_(double& r);
void decr_(double& r);

void add_(double& r, const double& x, const double& y, rounding_mode rnd);
void sub_(double& r, const double& x, const double& y, rounding_mode rnd);
void mul_(double& r, const double& x, const double& y, rounding_mode rnd);
void div_(double& r, const double& x, const double& y, rounding_mode rnd);

void add_(double& r, const int& x, const double& y, rounding_mode rnd);
void sub_(double& r, const int& x, const double& y, rounding_mode rnd);
void mul_(double& r, const int& x, const double& y, rounding_mode rnd);
void div_(double& r, const int& x, const double& y, rounding_mode rnd);

void add_(double& r, const double& x, const int& y, rounding_mode rnd);
void sub_(double& r, const double& x, const int& y, rounding_mode rnd);
void mul_(double& r, const double& x, const int& y, rounding_mode rnd);
void div_(double& r, const double& x, const int& y, rounding_mode rnd);

void add_(double& r, const unsigned int& x, const double& y, rounding_mode rnd);
void sub_(double& r, const unsigned int& x, const double& y, rounding_mode rnd);
void mul_(double& r, const unsigned int& x, const double& y, rounding_mode rnd);
void div_(double& r, const unsigned int& x, const double& y, rounding_mode rnd);

void add_(double& r, const double& x, const unsigned int& y, rounding_mode rnd);
void sub_(double& r, const double& x, const unsigned int& y, rounding_mode rnd);
void mul_(double& r, const double& x, const unsigned int& y, rounding_mode rnd);
void div_(double& r, const double& x, const unsigned int& y, rounding_mode rnd);

void med_(double& r, const double& x, const double& y, rounding_mode rnd);
void rad_(double& r, const double& x, const double& y, rounding_mode rnd);

void pow_(double& r, const double& x, const unsigned int& n, rounding_mode rnd);
void pow_(double& r, const double& x, const int& n, rounding_mode rnd);

void sqrt_(double& r, const double& x, rounding_mode rnd);
void hypot_(double& r, const double& x, const double& x, rounding_mode rnd);

void exp_(double& r, const double& x, rounding_mode rnd);
void log_(double& r, const double& x, rounding_mode rnd);
void sin_(double& r, const double& x, rounding_mode rnd);
void cos_(double& r, const double& x, rounding_mode rnd);
void tan_(double& r, const double& x, rounding_mode rnd);
void asin_(double& r, const double& x, rounding_mode rnd);
void acos_(double& r, const double& x, rounding_mode rnd);
void atan_(double& r, const double& x, rounding_mode rnd);
void sinh_(double& r, const double& x, rounding_mode rnd);
void cosh_(double& r, const double& x, rounding_mode rnd);
void tanh_(double& r, const double& x, rounding_mode rnd);




/*
template<class Rnd> inline void add_(double& r, const double& x, const double& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x+y; }
template<class Rnd> inline void sub_(double& r, const double& x, const double& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x-y; }
template<class Rnd> inline void mul_(double& r, const double& x, const double& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x*y; }
template<class Rnd> inline void div_(double& r, const double& x, const double& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x/y; }

template<class Rnd> inline void add_(double& r, const int& x, const double& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x+y; }
template<class Rnd> inline void sub_(double& r, const int& x, const double& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x-y; }
template<class Rnd> inline void mul_(double& r, const int& x, const double& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x*y; }
template<class Rnd> inline void div_(double& r, const int& x, const double& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x/y; }

template<class Rnd> inline void add_(double& r, const double& x, const int& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x+y; }
template<class Rnd> inline void sub_(double& r, const double& x, const int& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x-y; }
template<class Rnd> inline void mul_(double& r, const double& x, const int& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x*y; }
template<class Rnd> inline void div_(double& r, const double& x, const int& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x/y; }

template<class Rnd> inline void add_(double& r, const unsigned int& x, const double& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x+y; }
template<class Rnd> inline void sub_(double& r, const unsigned int& x, const double& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x-y; }
template<class Rnd> inline void mul_(double& r, const unsigned int& x, const double& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x*y; }
template<class Rnd> inline void div_(double& r, const unsigned int& x, const double& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x/y; }

template<class Rnd> inline void add_(double& r, const double& x, const unsigned int& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x+y; }
template<class Rnd> inline void sub_(double& r, const double& x, const unsigned int& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x-y; }
template<class Rnd> inline void mul_(double& r, const double& x, const unsigned int& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x*y; }
template<class Rnd> inline void div_(double& r, const double& x, const unsigned int& y, Rnd) {
  set_rounding_mode<Rnd>(); r=x/y; }

inline void med_(double& r, const double& x, const double& y, RoundApprox) {
  set_rounding_mode<RoundApprox>(); r=(x+y)/2; }
inline void rad_(double& r, const double& x, const double& y, RoundUp) {
  set_rounding_mode<RoundUp>(); r=(y-x)/2; }

template<class Rnd> inline void pow_(double& r, const double& x, const unsigned int& n, Rnd rnd) {
  set_rounding_mode<Rnd>(); 
  double p(x); unsigned int m=n; r=1.0; 
  while(m) { if(m%2) { r*=p; } p*=p; m/=2; }
}

template<class Rnd> inline void pow_(double& r, const double& x, const int& n, Rnd rnd) {
  set_rounding_mode<Rnd>(); 
  unsigned int m=(n>=0) ? n : -n;
  double p=(n>=0) ? x : 1/x;
  r=1.0; while(m) { if(m%2) { r*=p; } p*=p; m/=2; }
}

*/

}
}

#include "double.inline.h"

#endif /* ARIADNE_NUMERIC_DOUBLE_H */
