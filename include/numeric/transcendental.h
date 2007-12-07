/***************************************************************************
 *            transcendental.h
 *
 *  Copyright  2007  Pieter Collins
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
 
/*! \file numeric/transcendental.h
 *  \brief Templated algebraic and transcendental functions for floating point types
 */

#ifndef ARIADNE_TRANSCENDENTAL_H
#define ARIADNE_TRANSCENDENTAL_H


#include "numeric/interval.class.h"
#include "numeric/rounding.h"

namespace Ariadne {
  namespace Numeric {

    template<class R, class X> void sqrt_(Interval<R>& r, const Interval<X>& x);
    template<class R, class X, class Y> void hypot_(Interval<R>& r, const Interval<X>& x, const Interval<Y>& y);

    template<class R, class X> void exp_(Interval<R>& r, const Interval<X>& x);
    template<class R, class X> void log_(Interval<R>& r, const Interval<X>& x);

    template<class R> void pi_(Interval<R>& r);
    template<class R, class X> void sin_(Interval<R>& r, const Interval<X>& x);
    template<class R, class X> void cos_(Interval<R>& r, const Interval<X>& x);
    template<class R, class X> void tan_(Interval<R>& r, const Interval<X>& x);
    template<class R, class X> void asin_(Interval<R>& r, const Interval<X>& x);
    template<class R, class X> void acos_(Interval<R>& r, const Interval<X>& x);
    template<class R, class X> void atan_(Interval<R>& r, const Interval<X>& x);







        
    template<class T> inline void sqrt_(Interval< Float<T> >& r, const Float<T>& x) {
      sqrt_(r._lower,x,round_down); sqrt_(r._upper,x,round_up); } 
    template<class R> inline void sqrt_(Interval<R>& r, const Interval<R>& x) {
      sqrt_(r._lower,x._lower,round_down); sqrt_(r._upper,x._upper,round_up); } 



    template<class T, class Rnd> inline 
    void hypot_(Float<T>& r, const Float<T>& x, const Float<T>& y, Rnd rnd) {
      if(abs_(x)>abs_(y)) { hypot_(r,y,x,Rnd()); }
      Float<T> t; div_(t,x,y,rnd); mul_(t,t,t,rnd); add_(t,1,t,rnd); sqrt_(t,t,rnd); mul_(r,t,y,rnd);
    }

    template<class T> inline 
    void hypot_(Interval< Float<T> >& r, const Float<T>& x, const Float<T>& y) {
      hypot_(r._lower,x,y,round_down); hypot_(r._upper,x,y,round_up); } 

    template<class R, class X, class Y> inline 
    void hypot_(Interval<R> r, const Interval<X>& x, const Interval<Y>& y) {
      hypot_(r._lower,x._lower,y._lower,round_down); hypot_(r._upper,x._lower,y._lower,round_up); } 
  


    template<class T> inline
    void exp_(Interval< Float<T> >& r, const Float<T>& x) {
      exp_(r._lower,x,round_down); exp_(r._upper,x,round_up); } 
    
    template<class R, class X> inline 
    void exp_(Interval<R>& r, const Interval<X>& x) {
      exp_(r._lower,x._lower,round_down); exp_(r._upper,x._lower,round_up); } 


    template<class T> inline
    void log_(Interval< Float<T> >& r, const Float<T>& x) {
      log_(r._lower,x,round_down); log_(r._upper,x,round_up); } 

    template<class R, class X> inline 
    void log_(Interval<R>& r, const Interval<X>& x) {
      log_(r._lower,x._lower,round_down); log_(r._upper,x._lower,round_up); } 






    template<class R, class X> inline
    void sin_(Interval<R>& r, const X& x) {
      sin_(r,Interval<R>(x)); }

    template<class R, class X> inline
    void cos_(Interval<R>& r, const X& x) {
      cos_(r,Interval<R>(x)); }

    template<class R, class X> inline
    void tan_(Interval<R>& r, const X& x) {
      tan_(r,Interval<R>(x)); }

    template<class R, class X> inline
    void asin_(Interval<R>& r, const X& x) {
      asin_(r,Interval<R>(x)); }

    template<class R, class X> inline
    void acos_(Interval<R>& r, const X& x) {
      acos_(r,Interval<R>(x)); }

    template<class R, class X> inline
    void atan_(Interval<R>& r, const X& x) {
      atan_(r,Interval<R>(x)); }


 
  }

}

#endif /* ARIADNE_TRANSCENDENTAL_H */
