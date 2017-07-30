/***************************************************************************
 *            float-raw.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file float-raw.hpp
 *  \brief Temporary header for raw floating-point type.
 */

#ifndef ARIADNE_FLOAT_RAW_HPP
#define ARIADNE_FLOAT_RAW_HPP

#include "float.decl.hpp"
#include "float64.hpp"
#include "floatmp.hpp"
#include "numeric/operators.hpp"

namespace Ariadne {

/*
inline Float64 nul(exact,Float64 x) { return nul(x); }
inline Float64 pos(exact,Float64 x) { return pos(x); }
inline Float64 neg(exact,Float64 x) { return neg(x); }
inline Float64 hlf(exact,Float64 x) { return hlf(x); }

inline Float64 max(exact,Float64 x1, Float64 x2) { return max(x1,x2); }
inline Float64 min(exact,Float64 x1, Float64 x2) { return min(x1,x2); }
inline Float64 abs(exact,Float64 x) { return abs(x); }


inline Float64 add(approx,Float64 x1, Float64 x2) { return add(near,x1,x2); }
inline Float64 add(near,Float64 x1, Float64 x2) { return add(near,x1,x2); }
inline Float64 add(up,Float64 x1, Float64 x2) { return add(up,x1,x2); }
inline Float64 add(down,Float64 x1, Float64 x2) { return add(down,x1,x2); }

inline Float64 sub(approx,Float64 x1, Float64 x2) { return sub(near,x1,x2); }
inline Float64 sub(near,Float64 x1, Float64 x2) { return sub(near,x1,x2); }
inline Float64 sub(up,Float64 x1, Float64 x2) { return sub(up,x1,x2); }
inline Float64 sub(down,Float64 x1, Float64 x2) { return sub(down,x1,x2); }

inline Float64 mul(approx,Float64 x1, Float64 x2) { return mul(near,x1,x2); }
inline Float64 mul(near,Float64 x1, Float64 x2) { return mul(near,x1,x2); }
inline Float64 mul(up,Float64 x1, Float64 x2) { return mul(up,x1,x2); }
inline Float64 mul(down,Float64 x1, Float64 x2) { return mul(down,x1,x2); }

inline Float64 div(approx,Float64 x1, Float64 x2) { return div(near,x1,x2); }
inline Float64 div(near,Float64 x1, Float64 x2) { return div(near,x1,x2); }
inline Float64 div(up,Float64 x1, Float64 x2) { return div(up,x1,x2); }
inline Float64 div(down,Float64 x1, Float64 x2) { return div(down,x1,x2); }

inline Float64 pow(approx,Float64 x, Int n) { return pow(near,x,n); }
inline Float64 pow(up,Float64 x, Int n) { return pow(up,x,n); }
inline Float64 pow(down,Float64 x, Int n) { return pow(down,x,n); }

inline Float64 next(down,Float64 x) { return sub(down,x,Float64::min(Precision64())); }
inline Float64 next(up,Float64 x) { return add(up,x,Float64::min(Precision64())); }

inline Float64 rec(approx,Float64 x) { return rec(near,x); }
inline Float64 rec(near,Float64 x) { return rec(near,x); }
inline Float64 rec(up,Float64 x) { return rec(up,x); }
inline Float64 rec(down,Float64 x) { return rec(down,x); }

inline Float64 sqrt(up,Float64 const& x) { return sqrt(up,x); }
inline Float64 sqrt(down,Float64 const& x) { return sqrt(down,x); }

inline Float64 exp(up,Float64 const& x) { return exp(up,x); }
inline Float64 exp(down,Float64 const& x) { return exp(down,x); }

inline Float64 log(up,Float64 const& x) { return log(up,x); }
inline Float64 log(down,Float64 const& x) { return log(down,x); }

inline Float64 sin(up,Float64 const& x) { return sin(up,x); }
inline Float64 sin(down,Float64 const& x) { return sin(down,x); }

inline Float64 cos(up,Float64 const& x) { return cos(up,x); }
inline Float64 cos(down,Float64 const& x) { return cos(down,x); }

inline Float64 atan(up,Float64 const& x) { return atan(up,x); }
inline Float64 atan(down,Float64 const& x) { return atan(down,x); }

inline Float64 pi_up(Precision64 const& pr) { return Float64::pi_up(pr); }
inline Float64 pi_approx(Precision64 const& pr) { return Float64::pi_near(pr); }
inline Float64 pi_down(Precision64 const& pr) { return Float64::pi_down(pr); }


inline FloatMP nul(exact,FloatMP x) { return nul(x); }
inline FloatMP pos(exact,FloatMP x) { return pos(x); }
inline FloatMP neg(exact,FloatMP x) { return neg(x); }
inline FloatMP hlf(exact,FloatMP x) { return hlf(x); }

inline FloatMP max(exact,FloatMP x1, FloatMP x2) { return max(x1,x2); }
inline FloatMP min(exact,FloatMP x1, FloatMP x2) { return min(x1,x2); }
inline FloatMP abs(exact,FloatMP x) { return abs(x); }


inline FloatMP add(approx,FloatMP const& x1, FloatMP const& x2) { return add(near,x1,x2); }
inline FloatMP add(near,FloatMP const& x1, FloatMP const& x2) { return add(near,x1,x2); }
inline FloatMP add(up,FloatMP const& x1, FloatMP const& x2) { return add(up,x1,x2); }
inline FloatMP add(down,FloatMP const& x1, FloatMP const& x2) { return add(down,x1,x2); }

inline FloatMP sub(approx,FloatMP const& x1, FloatMP const& x2) { return sub(near,x1,x2); }
inline FloatMP sub(near,FloatMP const& x1, FloatMP const& x2) { return sub(near,x1,x2); }
inline FloatMP sub(up,FloatMP const& x1, FloatMP const& x2) { return sub(up,x1,x2); }
inline FloatMP sub(down,FloatMP const& x1, FloatMP const& x2) { return sub(down,x1,x2); }

inline FloatMP mul(approx,FloatMP const& x1, FloatMP const& x2) { return mul(near,x1,x2); }
inline FloatMP mul(near,FloatMP const& x1, FloatMP const& x2) { return mul(near,x1,x2); }
inline FloatMP mul(up,FloatMP const& x1, FloatMP const& x2) { return mul(up,x1,x2); }
inline FloatMP mul(down,FloatMP const& x1, FloatMP const& x2) { return mul(down,x1,x2); }

inline FloatMP div(approx,FloatMP const& x1, FloatMP const& x2) { return div(near,x1,x2); }
inline FloatMP div(near,FloatMP const& x1, FloatMP const& x2) { return div(near,x1,x2); }
inline FloatMP div(up,FloatMP const& x1, FloatMP const& x2) { return div(up,x1,x2); }
inline FloatMP div(down,FloatMP const& x1, FloatMP const& x2) { return div(down,x1,x2); }

inline FloatMP pow(approx,FloatMP const& x, Int n) { return pow(near,x,n); }
inline FloatMP pow(up,FloatMP const& x, Int n) { return pow(up,x,n); }
inline FloatMP pow(down,FloatMP const& x, Int n) { return pow(down,x,n); }

inline FloatMP next(down,FloatMP const& x) { return sub(down,x,FloatMP::min(x.precision())); }
inline FloatMP next(up,FloatMP const& x) { return add(up,x,FloatMP::min(x.precision())); }

inline FloatMP rec(approx,FloatMP const& x) { return rec(near,x); }
inline FloatMP rec(near,FloatMP const& x) { return rec(near,x); }
inline FloatMP rec(up,FloatMP const& x) { return rec(up,x); }
inline FloatMP rec(down,FloatMP const& x) { return rec(down,x); }

inline FloatMP sqrt(up,FloatMP const& x) { return sqrt(up,x); }
inline FloatMP sqrt(down,FloatMP const& x) { return sqrt(down,x); }

inline FloatMP exp(up,FloatMP const& x) { return exp(up,x); }
inline FloatMP exp(down,FloatMP const& x) { return exp(down,x); }

inline FloatMP log(up,FloatMP const& x) { return log(up,x); }
inline FloatMP log(down,FloatMP const& x) { return log(down,x); }

inline FloatMP sin(up,FloatMP const& x) { return sin(up,x); }
inline FloatMP sin(down,FloatMP const& x) { return sin(down,x); }

inline FloatMP cos(up,FloatMP const& x) { return cos(up,x); }
inline FloatMP cos(down,FloatMP const& x) { return cos(down,x); }

inline FloatMP atan(up,FloatMP const& x) { return atan(up,x); }
inline FloatMP atan(down,FloatMP const& x) { return atan(down,x); }

inline FloatMP pi_up(PrecisionMP const& pr) { return FloatMP::pi_up(pr); }
inline FloatMP pi_approx(PrecisionMP const& pr) { return FloatMP::pi_near(pr); }
inline FloatMP pi_down(PrecisionMP const& pr) { return FloatMP::pi_down(pr); }

inline FloatMP med(near,FloatMP x1, FloatMP x2) { return hlf(add(near,x1,x2)); }
inline FloatMP rad(up,FloatMP x1, FloatMP x2) { return hlf(sub(up,x2,x1)); }

FloatMP add(up,FloatMP x1, Float64 x2);
FloatMP sub(down,FloatMP x1, Float64 x2);
*/
/*
inline Float64 max(exact,Float64 const& x1, Float64 const& x2) { return max(x1,x2); }
inline Float64 min(exact,Float64 const& x1, Float64 const& x2) { return min(x1,x2); }
inline Float64 abs(exact,Float64 const& x) { return abs(x); }

inline Float64 add(Float64 const& x1, Float64 const& x2, Float64::RoundingModeType rnd) { return apply(Add(),x1,x2,rnd); }
inline Float64 sub(Float64 const& x1, Float64 const& x2, Float64::RoundingModeType rnd) { return apply(Sub(),x1,x2,rnd); }
inline Float64 mul(Float64 const& x1, Float64 const& x2, Float64::RoundingModeType rnd) { return apply(Mul(),x1,x2,rnd); }
inline Float64 div(Float64 const& x1, Float64 const& x2, Float64::RoundingModeType rnd) { return apply(Div(),x1,x2,rnd); }
inline Float64 pow(Float64 const& x, Int n, Float64::RoundingModeType rnd) { return apply(Pow(),x,n,rnd); }

inline Float64 nul(exact,Float64 const& x) { return nul(x); }
inline Float64 pos(exact,Float64 const& x) { return pos(x); }
inline Float64 neg(exact,Float64 const& x) { return neg(x); }
inline Float64 hlf(exact,Float64 const& x) { return hlf(x); }

inline Float64 add(up,Float64 const& x1, Float64 const& x2) { return add(x1,x2,Float64::upward); }
inline Float64 add(down,Float64 const& x1, Float64 const& x2) { return add(x1,x2,Float64::downward); }
inline Float64 add(near,Float64 const& x1, Float64 const& x2) { return add(x1,x2,Float64::to_nearest); }
inline Float64 add(approx,Float64 const& x1, Float64 const& x2) { return add(x1,x2,Float64::to_nearest); }

inline Float64 sub(up,Float64 const& x1, Float64 const& x2) { return sub(x1,x2,Float64::upward); }
inline Float64 sub(down,Float64 const& x1, Float64 const& x2) { return sub(x1,x2,Float64::downward); }
inline Float64 sub(near,Float64 const& x1, Float64 const& x2) { return sub(x1,x2,Float64::to_nearest); }
inline Float64 sub(approx,Float64 const& x1, Float64 const& x2) { return sub(x1,x2,Float64::to_nearest); }

inline Float64 mul(up,Float64 const& x1, Float64 const& x2) { return mul(x1,x2,Float64::upward); }
inline Float64 mul(down,Float64 const& x1, Float64 const& x2) { return mul(x1,x2,Float64::downward); }
inline Float64 mul(near,Float64 const& x1, Float64 const& x2) { return mul(x1,x2,Float64::to_nearest); }
inline Float64 mul(approx,Float64 const& x1, Float64 const& x2) { return mul(x1,x2,Float64::to_nearest); }

inline Float64 div(up,Float64 const& x1, Float64 const& x2) { return div(x1,x2,Float64::upward); }
inline Float64 div(down,Float64 const& x1, Float64 const& x2) { return div(x1,x2,Float64::downward); }
inline Float64 div(near,Float64 const& x1, Float64 const& x2) { return div(x1,x2,Float64::to_nearest); }
inline Float64 div(approx,Float64 const& x1, Float64 const& x2) { return div(x1,x2,Float64::to_nearest); }

inline Float64 fma(up,Float64 const& x1, Float64 const& x2, Float64 const& y) { return add(mul(x1,x2,Float64::upward),y,upward); }
inline Float64 fma(down,Float64 const& x1, Float64 const& x2, Float64 const& y) { return add(mul(x1,x2,Float64::downward),y,downward); }
inline Float64 fma(approx,Float64 const& x1, Float64 const& x2, Float64 const& y) { return add(mul(x1,x2,Float64::to_nearest),y,to_nearest); }

inline Float64 next(up,Float64 const& x) { return add(up,x,Float64::min()); }
inline Float64 next(down,Float64 const& x) { return sub(down,x,Float64::min()); }

inline Float64 rad(up,Float64 const& x1, Float64 const& x2) { return hlf(sub(up,x2,x1)); }
inline Float64 med(near,Float64 const& x1, Float64 const& x2) { return hlf(add(near,x2,x1)); }

inline Float64 pow(up,Float64 const& x, Int n) { return pow(x,n,Float64::upward); }
inline Float64 pow(down,Float64 const& x, Int n) { return pow(x,n,Float64::downward); }
inline Float64 pow(approx,Float64 const& x, Int n) { return pow(x,n,Float64::to_nearest); }

Float64 pi_up();
Float64 pi_down();
Float64 pi_approx();

inline Float64 sqrt(approx,Float64 x) { return std::sqrt(x.dbl); }
inline Float64 exp(approx,Float64 x) { return std::exp(x.dbl); }
inline Float64 log(approx,Float64 x) { return std::log(x.dbl); }
inline Float64 sin(approx,Float64 x) { return std::sin(x.dbl); }
inline Float64 cos(approx,Float64 x) { return std::cos(x.dbl); }
inline Float64 tan(approx,Float64 x) { return std::tan(x.dbl); }
inline Float64 asin(approx,Float64 x) { return std::asin(x.dbl); }
inline Float64 acos(approx,Float64 x) { return std::acos(x.dbl); }
inline Float64 atan(approx,Float64 x) { return std::atan(x.dbl); }

inline FloatMP sqrt(approx,FloatMP const& x) { return sqrt(near,x); }
inline FloatMP exp(approx,FloatMP const& x) { return exp(near,x); }
inline FloatMP log(approx,FloatMP const& x) { return log(near,x); }
inline FloatMP sin(approx,FloatMP const& x) { return sin(near,x); }
inline FloatMP cos(approx,FloatMP const& x) { return cos(near,x); }
inline FloatMP tan(approx,FloatMP const& x) { return tan(near,x); }
inline FloatMP asin(approx,FloatMP const& x) { return asin(near,x); }
inline FloatMP acos(approx,FloatMP const& x) { return acos(near,x); }
inline FloatMP atan(approx,FloatMP const& x) { return atan(near,x); }
*/

} // namespace Ariadne

#endif /* ARIADNE_FLOAT_RAW_HPP */
