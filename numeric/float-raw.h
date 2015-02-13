/***************************************************************************
 *            float-raw.h
 *
 *  Copyright 2008-15  Pieter Collins
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

/*! \file float-raw.h
 *  \brief Temporary header for raw floating-point type.
 */

#ifndef ARIADNE_FLOAT_RAW_H
#define ARIADNE_FLOAT_RAW_H

#include "float.decl.h"
#include "float64.h"
#include "floatmp.h"
#include "expression/operators.h"

namespace Ariadne {

static const Float64 inf = std::numeric_limits<double>::infinity();

inline Float64 nul_exact(Float64 x) { return nul(x); }
inline Float64 pos_exact(Float64 x) { return pos(x); }
inline Float64 neg_exact(Float64 x) { return neg(x); }
inline Float64 half_exact(Float64 x) { return half(x); }

inline Float64 max_exact(Float64 x1, Float64 x2) { return max(x1,x2); }
inline Float64 min_exact(Float64 x1, Float64 x2) { return min(x1,x2); }
inline Float64 abs_exact(Float64 x) { return abs(x); }

inline Float64 add_approx(Float64 x1, Float64 x2) { return add(x1,x2,to_nearest); }
inline Float64 add_near(Float64 x1, Float64 x2) { return add(x1,x2,to_nearest); }
inline Float64 add_up(Float64 x1, Float64 x2) { return add(x1,x2,upward); }
inline Float64 add_down(Float64 x1, Float64 x2) { return add(x1,x2,downward); }

inline Float64 sub_approx(Float64 x1, Float64 x2) { return sub(x1,x2,to_nearest); }
inline Float64 sub_near(Float64 x1, Float64 x2) { return sub(x1,x2,to_nearest); }
inline Float64 sub_up(Float64 x1, Float64 x2) { return sub(x1,x2,upward); }
inline Float64 sub_down(Float64 x1, Float64 x2) { return sub(x1,x2,downward); }

inline Float64 mul_approx(Float64 x1, Float64 x2) { return mul(x1,x2,to_nearest); }
inline Float64 mul_near(Float64 x1, Float64 x2) { return mul(x1,x2,to_nearest); }
inline Float64 mul_up(Float64 x1, Float64 x2) { return mul(x1,x2,upward); }
inline Float64 mul_down(Float64 x1, Float64 x2) { return mul(x1,x2,downward); }

inline Float64 div_approx(Float64 x1, Float64 x2) { return div(x1,x2,to_nearest); }
inline Float64 div_near(Float64 x1, Float64 x2) { return div(x1,x2,to_nearest); }
inline Float64 div_up(Float64 x1, Float64 x2) { return div(x1,x2,upward); }
inline Float64 div_down(Float64 x1, Float64 x2) { return div(x1,x2,downward); }

inline Float64 pow_approx(Float64 x, Int n) { return pow(x,n,to_nearest); }
inline Float64 pow_up(Float64 x, Int n) { return pow(x,n,upward); }
inline Float64 pow_down(Float64 x, Int n) { return pow(x,n,downward); }

inline Float64 next_down(Float64 x) { return sub_down(x,Float64::min()); }
inline Float64 next_up(Float64 x) { return add_up(x,Float64::min()); }

inline Float64 rec_approx(Float64 x) { return rec(x,to_nearest); }
inline Float64 rec_up(Float64 x) { return rec(x,upward); }
inline Float64 rec_down(Float64 x) { return rec(x,downward); }

inline Float64 sqrt_up(Float64 const& x) { return sqrt(x,upward); }
inline Float64 sqrt_down(Float64 const& x) { return sqrt(x,downward); }

inline Float64 exp_up(Float64 const& x) { return exp(x,upward); }
inline Float64 exp_down(Float64 const& x) { return exp(x,downward); }

inline Float64 log_up(Float64 const& x) { return log(x,upward); }
inline Float64 log_down(Float64 const& x) { return log(x,downward); }

inline Float64 sin_up(Float64 const& x) { return sin(x,upward); }
inline Float64 sin_down(Float64 const& x) { return sin(x,downward); }

inline Float64 cos_up(Float64 const& x) { return cos(x,upward); }
inline Float64 cos_down(Float64 const& x) { return cos(x,downward); }

inline Float64 atan_up(Float64 const& x) { return atan(x,upward); }
inline Float64 atan_down(Float64 const& x) { return atan(x,downward); }


/*
inline Float64 max_exact(Float64 const& x1, Float64 const& x2) { return max(x1,x2); }
inline Float64 min_exact(Float64 const& x1, Float64 const& x2) { return min(x1,x2); }
inline Float64 abs_exact(Float64 const& x) { return abs(x); }

inline Float64 add(Float64 const& x1, Float64 const& x2, Float64::RoundingModeType rnd) { return apply(Add(),x1,x2,rnd); }
inline Float64 sub(Float64 const& x1, Float64 const& x2, Float64::RoundingModeType rnd) { return apply(Sub(),x1,x2,rnd); }
inline Float64 mul(Float64 const& x1, Float64 const& x2, Float64::RoundingModeType rnd) { return apply(Mul(),x1,x2,rnd); }
inline Float64 div(Float64 const& x1, Float64 const& x2, Float64::RoundingModeType rnd) { return apply(Div(),x1,x2,rnd); }
inline Float64 pow(Float64 const& x, Int n, Float64::RoundingModeType rnd) { return apply(Pow(),x,n,rnd); }

inline Float64 nul_exact(Float64 const& x) { return nul(x); }
inline Float64 pos_exact(Float64 const& x) { return pos(x); }
inline Float64 neg_exact(Float64 const& x) { return neg(x); }
inline Float64 half_exact(Float64 const& x) { return half(x); }

inline Float64 add_up(Float64 const& x1, Float64 const& x2) { return add(x1,x2,Float64::upward); }
inline Float64 add_down(Float64 const& x1, Float64 const& x2) { return add(x1,x2,Float64::downward); }
inline Float64 add_near(Float64 const& x1, Float64 const& x2) { return add(x1,x2,Float64::to_nearest); }
inline Float64 add_approx(Float64 const& x1, Float64 const& x2) { return add(x1,x2,Float64::to_nearest); }

inline Float64 sub_up(Float64 const& x1, Float64 const& x2) { return sub(x1,x2,Float64::upward); }
inline Float64 sub_down(Float64 const& x1, Float64 const& x2) { return sub(x1,x2,Float64::downward); }
inline Float64 sub_near(Float64 const& x1, Float64 const& x2) { return sub(x1,x2,Float64::to_nearest); }
inline Float64 sub_approx(Float64 const& x1, Float64 const& x2) { return sub(x1,x2,Float64::to_nearest); }

inline Float64 mul_up(Float64 const& x1, Float64 const& x2) { return mul(x1,x2,Float64::upward); }
inline Float64 mul_down(Float64 const& x1, Float64 const& x2) { return mul(x1,x2,Float64::downward); }
inline Float64 mul_near(Float64 const& x1, Float64 const& x2) { return mul(x1,x2,Float64::to_nearest); }
inline Float64 mul_approx(Float64 const& x1, Float64 const& x2) { return mul(x1,x2,Float64::to_nearest); }

inline Float64 div_up(Float64 const& x1, Float64 const& x2) { return div(x1,x2,Float64::upward); }
inline Float64 div_down(Float64 const& x1, Float64 const& x2) { return div(x1,x2,Float64::downward); }
inline Float64 div_near(Float64 const& x1, Float64 const& x2) { return div(x1,x2,Float64::to_nearest); }
inline Float64 div_approx(Float64 const& x1, Float64 const& x2) { return div(x1,x2,Float64::to_nearest); }

inline Float64 fma_up(Float64 const& x1, Float64 const& x2, Float64 const& y) { return add(mul(x1,x2,Float64::upward),y,upward); }
inline Float64 fma_down(Float64 const& x1, Float64 const& x2, Float64 const& y) { return add(mul(x1,x2,Float64::downward),y,downward); }
inline Float64 fma_approx(Float64 const& x1, Float64 const& x2, Float64 const& y) { return add(mul(x1,x2,Float64::to_nearest),y,to_nearest); }

inline Float64 next_up(Float64 const& x) { return add_up(x,Float64::min()); }
inline Float64 next_down(Float64 const& x) { return sub_down(x,Float64::min()); }

inline Float64 rad_up(Float64 const& x1, Float64 const& x2) { return half(sub_up(x2,x1)); }
inline Float64 med_near(Float64 const& x1, Float64 const& x2) { return half(add_near(x2,x1)); }

inline Float64 pow_up(Float64 const& x, Int n) { return pow(x,n,Float64::upward); }
inline Float64 pow_down(Float64 const& x, Int n) { return pow(x,n,Float64::downward); }
inline Float64 pow_approx(Float64 const& x, Int n) { return pow(x,n,Float64::to_nearest); }

Float64 pi_up();
Float64 pi_down();
Float64 pi_approx();

inline Float64 sqrt_approx(Float64 x) { return std::sqrt(x.dbl); }
inline Float64 exp_approx(Float64 x) { return std::exp(x.dbl); }
inline Float64 log_approx(Float64 x) { return std::log(x.dbl); }
inline Float64 sin_approx(Float64 x) { return std::sin(x.dbl); }
inline Float64 cos_approx(Float64 x) { return std::cos(x.dbl); }
inline Float64 tan_approx(Float64 x) { return std::tan(x.dbl); }
inline Float64 asin_approx(Float64 x) { return std::asin(x.dbl); }
inline Float64 acos_approx(Float64 x) { return std::acos(x.dbl); }
inline Float64 atan_approx(Float64 x) { return std::atan(x.dbl); }
*/
inline FloatMP sqrt_approx(FloatMP const& x) { return sqrt(x,FloatMP::to_nearest); }
inline FloatMP exp_approx(FloatMP const& x) { return exp(x,FloatMP::to_nearest); }
inline FloatMP log_approx(FloatMP const& x) { return log(x,FloatMP::to_nearest); }
inline FloatMP sin_approx(FloatMP const& x) { return sin(x,FloatMP::to_nearest); }
inline FloatMP cos_approx(FloatMP const& x) { return cos(x,FloatMP::to_nearest); }
inline FloatMP tan_approx(FloatMP const& x) { return tan(x,FloatMP::to_nearest); }
inline FloatMP asin_approx(FloatMP const& x) { return asin(x,FloatMP::to_nearest); }
inline FloatMP acos_approx(FloatMP const& x) { return acos(x,FloatMP::to_nearest); }
inline FloatMP atan_approx(FloatMP const& x) { return atan(x,FloatMP::to_nearest); }

} // namespace Ariadne

#endif /* ARIADNE_FLOAT_RAW_H */
