/***************************************************************************
 *            transcendental64.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file numeric/transcendental64.h
 *  \brief Transcendental functions for 64 bit floating-point type.
 */

#ifndef ARIADNE_TRANSCENDENTAL64_H
#define ARIADNE_TRANSCENDENTAL64_H

#include <cmath>
#include <limits>
#include <mpfr.h>
#include <boost/numeric/interval/rounded_arith.hpp>
#include <boost/numeric/interval/rounded_transc.hpp>
#include <boost/numeric/interval/hw_rounding.hpp>

#include "numeric/interval.class.h"
#include "numeric/rounding.h"

namespace Ariadne {
  namespace Numeric {
        
    inline void sqrt_(Float64& r, const Float64& x, RoundDown) { 
      r._value=Float64::rounding().sqrt_down(x._value); }
    inline void sqrt_(Float64& r, const Float64& x, RoundUp) { 
      r._value=Float64::rounding().sqrt_up(x._value); }
    inline void sqrt_(Float64& r, const Float64& x, RoundApprox) { 
      r._value=std::sqrt(x._value); }


    inline void exp_(Float64& r, const Float64& x, RoundApprox) {
      r._value=std::exp(x._value); }
    inline void exp_(Float64& r, const Float64& x, RoundDown) {
      r._value=Float64::rounding().exp_down(x._value); }
    inline void exp_(Float64& r, const Float64& x, RoundUp) {
      r._value=Float64::rounding().exp_up(x._value); }
    
    inline void log_(Float64& r, const Float64& x, RoundApprox) {
      r._value=std::log(x._value); }
    inline void log_(Float64& r, const Float64& x, RoundDown) {
      r._value=Float64::rounding().log_down(x._value); }
    inline void log_(Float64& r, const Float64& x, RoundUp) {
      r._value=Float64::rounding().log_up(x._value); }


    inline void pi_(Float64& r, RoundApprox) {
      r._value=std::acos(2)*2; }
    inline void pi_(Interval64& r) { 
      r._lower=3.14; r._upper=3.15; }

    inline void sin_(Float64& r, Float64& x, RoundApprox) {
      r._value=std::sin(x._value); }
    void sin_(Interval64& r, Interval64& x);

    inline void cos_(Float64& r, Float64& x, RoundApprox) {
      r._value=std::cos(x._value); }
    void cos_(Interval64& r, Interval64& x);

    inline void tan_(Float64& r, Float64& x, RoundApprox) {
      r._value=std::tan(x._value); }
    void tan_(Interval64& r, Interval64& x);

    inline void asin_(Float64& r, Float64& x, RoundApprox) {
      r._value=std::asin(x._value); }
    void asin_(Interval64& r, Interval64& x);

    inline void acos_(Float64& r, Float64& x, RoundApprox) {
      r._value=std::acos(x._value); }
    void acos_(Interval64& r, Interval64& x);

    inline void atan_(Float64& r, Float64& x, RoundApprox) {
      r._value=std::atan(x._value); }
    void atan_(Interval64& r, Interval64& x);
   


  }
}

#include "numeric/transcendental.h"

#endif /* ARIADNE_TRANSCENDENTAL64_H */
