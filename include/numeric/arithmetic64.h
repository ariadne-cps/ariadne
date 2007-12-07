/***************************************************************************
 *            arithmetic64.h
 *
 *  Copyright 2005-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file arithmetic64.h
 *  \brief Simple arithmetic functions using expression templates.
 */
 
#ifndef ARIADNE_ARITHMETIC64_H
#define ARIADNE_ARITHMETIC64_H

#include "numeric/expression.h"
#include "numeric/rounding.h"

namespace Ariadne {
  namespace Numeric {

    inline void floor_(Float64& r, const Float64& x) { r=std::floor(x._value); }
    inline void floor_(Integer& r, const Float64& x) { r=int(std::floor(x._value)); }
    inline void floor_(int& r, const Float64& x) { r=int(std::floor(x._value)); }
	  
    inline void ceil_(Float64& r, const Float64& x) { r=std::ceil(x._value); }
    inline void ceil_(Integer& r, const Float64& x) { r=int(std::ceil(x._value)); }
    inline void ceil_(int& r, const Float64& x) { r=int(std::ceil(x._value)); }
	  
    template<class Rnd> 
    inline void add_(Float64& r, const Float64& x1,const Float64& x2, Rnd) {
      set_rounding_mode<Rnd>(); r._value=x1._value+x2._value; }

    inline void add_(Float64& r, const Float64& x1,const Float64& x2, RoundDown) {
      r._value=Float64::rounding().add_down(x1._value,x2._value); }
    inline void add_(Float64& r, const Float64& x1,const Float64& x2, RoundUp) {
      r._value=Float64::rounding().add_up(x1._value,x2._value); }
    inline void add_(Float64& r, const Float64& x1,const Float64& x2, RoundApprox) {
      r._value=x1._value+x2._value; }

    
    inline void sub_(Float64& r, const Float64& x1,const Float64& x2, RoundDown) {
      r._value=Float64::rounding().sub_down(x1._value,x2._value); }
    inline void sub_(Float64& r, const Float64& x1,const Float64& x2, RoundUp) {
      r._value=Float64::rounding().sub_up(x1._value,x2._value); }
    inline void sub_(Float64& r, const Float64& x1,const Float64& x2, RoundApprox) {
      r._value=x1._value-x2._value; }
    
    inline void mul_(Float64& r, const Float64& x1,const Float64& x2, RoundDown) {
      r._value=Float64::rounding().mul_down(x1._value,x2._value); }
    inline void mul_(Float64& r, const Float64& x1,const Float64& x2, RoundUp) {
      r._value=Float64::rounding().mul_up(x1._value,x2._value); }
    inline void mul_(Float64& r, const Float64& x1,const Float64& x2, RoundApprox) {
      r._value=x1._value*x2._value; }
     
    inline void div_(Float64& r, const Float64& x1,const Float64& x2, RoundDown) {
      r._value=Float64::rounding().div_down(x1._value,x2._value); }
    inline void div_(Float64& r, const Float64& x1,const Float64& x2, RoundUp) {
      r._value=Float64::rounding().div_up(x1._value,x2._value); }
    inline void div_(Float64& r, const Float64& x1,const Float64& x2, RoundApprox) {
      r._value=x1._value/x2._value; }

  /*
    inline void div_(Float64& r, const Float64& x1,const Float64& x2, RoundDown) {
      r._value=Float64::rounding().div_down(x1._value,x2._value); }
    inline void div_(Float64& r, const Float64& x1,const Float64& x2, RoundUp) {
      r._value=Float64::rounding().div_up(x1._value,x2._value); }
    inline void div_(Float64& r, const Float64& x1,const Float64& x2, RoundApprox) {
      r._value=x1._value/x2._value; }
  */
    
    inline void med_(Float64& r, const Float64& x1,const Float64& x2, RoundApprox) {
      r._value=(x1._value+x2._value)/2; }
    

	  
	  
    inline void mul_(Float64& r, const Float64& x, const int& n, RoundDown) {
      r._value=Float64::rounding().mul_down(x._value,n); }
    inline void mul_(Float64& r, const Float64& x,const int& n, RoundUp) {
      r._value=Float64::rounding().mul_up(x._value,n); }
    inline void mul_(Float64& r, const Float64& x, const int& n, RoundApprox) {
      r._value=x._value*n; }
     
    inline void mul_(Float64& r, const int& n, const Float64& x, RoundDown) {
      r._value=Float64::rounding().mul_down(x._value,n); }
    inline void mul_(Float64& r, const int& n, const Float64& x, RoundUp) {
      r._value=Float64::rounding().mul_up(x._value,n); }
    inline void mul_(Float64& r, const int& n, const Float64& x, RoundApprox) {
      r._value=x._value*n; }

    inline void mul_(Float64& r, const Float64& x, const double& d, RoundDown) {
      r._value=Float64::rounding().mul_down(x._value,d); }
    inline void mul_(Float64& r, const Float64& x,const double& d, RoundUp) {
      r._value=Float64::rounding().mul_up(x._value,d); }
    inline void mul_(Float64& r, const Float64& x, const double& d, RoundApprox) {
      r._value=x._value*d; }
     
    inline void mul_(Float64& r, const uint& n, const Float64& x, RoundApprox) {
      r._value=x._value*n; }
    inline void mul_(Float64& r, const Float64& x, const uint& n, RoundApprox) {
      r._value=x._value*n; }
      
    inline void mul_(Float64& r, const double& d, const Float64& x, RoundDown) {
      r._value=Float64::rounding().mul_down(x._value,d); }
    inline void mul_(Float64& r, const double& d, const Float64& x, RoundUp) {
      r._value=Float64::rounding().mul_up(x._value,d); }
    inline void mul_(Float64& r, const double& d, const Float64& x, RoundApprox) {
      r._value=x._value*d; }


    inline void div_(Float64& r, const Float64& x,const int& n, RoundDown) {
      r._value=Float64::rounding().div_down(x._value,n); }
    inline void div_(Float64& r, const Float64& x,const int& n, RoundUp) {
      r._value=Float64::rounding().div_up(x._value,n); }
    inline void div_(Float64& r, const Float64& x,const int& n, RoundApprox) {
      r._value=x._value/n; }
    

    template<class T, class Rnd>
    inline void pow_(Float<T>& r, const Float<T>& x, const uint& n, Rnd rnd) {
      Float<T> p=x; r=1; uint m=n;
      while(m) { if(m%2) { mul_(r,r,p,rnd); } mul_(p,p,p,rnd); m/=2; }
    }      

    template<class T, class Rnd>
    inline void pow_(Float<T>& r, const Float<T>& x, const int& n, Rnd rnd) {
      if(n>=0) { pow_(r,x,uint(n),rnd); }
      else { Float<T> t; div_(t,1,x,rnd); pow_(r,t,uint(-n),rnd); }
    }

    inline void next_(Float64& r, const Float64& x, RoundDown) { 
      Float64 min=std::numeric_limits<double>::min();
      sub_(r,x,min,round_down); }
    inline void next_(Float64& r, const Float64& x, RoundUp) { 
      Float64 min=std::numeric_limits<double>::min();
      add_(r,x,min,round_up); }
  
    inline void incr(Float64& x) { 
      Float64 m=std::numeric_limits<double>::min(); add_(x,x,m,round_up); }
    inline void decr(Float64& x) { 
      Float64 m=std::numeric_limits<double>::min(); sub_(x,x,m,round_down); }

    inline Float64 next(Float64& x, RoundDown) { 
      Float64 r; next_(r,x,round_down); return r; }
    inline Float64 next(Float64& x, RoundUp) { 
      Float64 r; next_(r,x,round_down); return r; }

    inline void set_(Rational& q, const Float64& x) { mpq_set_d(q._value,x._value); mpq_canonicalize(q._value); }
    template<> inline Rational& Rational::operator=(const Float64& x) { set_(*this,x); return *this; }

    inline void set_(Float64& r, const Rational& q, RoundDown) { r._value=mpq_get_d(q._value); next_(r,r,round_down); }
    inline void set_(Float64& r, const Rational& q, RoundUp) { r._value=mpq_get_d(q._value); next_(r,r,round_up); }
    inline void set_(Float64& r, const Rational& q, RoundApprox) { r._value=mpq_get_d(q._value); }



  }
}

#endif /* ARIADNE_ARITHMETIC64_H */
