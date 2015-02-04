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

static const Float inf = std::numeric_limits<double>::infinity();

inline Float nul_exact(Float x) { return nul(x); }
inline Float pos_exact(Float x) { return pos(x); }
inline Float neg_exact(Float x) { return neg(x); }
inline Float half_exact(Float x) { return half(x); }

inline Float max_exact(Float x1, Float x2) { return max(x1,x2); }
inline Float min_exact(Float x1, Float x2) { return min(x1,x2); }
inline Float abs_exact(Float x) { return abs(x); }

inline Float add_approx(Float x1, Float x2) { return add(x1,x2,to_nearest); }
inline Float add_near(Float x1, Float x2) { return add(x1,x2,to_nearest); }
inline Float add_up(Float x1, Float x2) { return add(x1,x2,upward); }
inline Float add_down(Float x1, Float x2) { return add(x1,x2,downward); }

inline Float sub_approx(Float x1, Float x2) { return sub(x1,x2,to_nearest); }
inline Float sub_near(Float x1, Float x2) { return sub(x1,x2,to_nearest); }
inline Float sub_up(Float x1, Float x2) { return sub(x1,x2,upward); }
inline Float sub_down(Float x1, Float x2) { return sub(x1,x2,downward); }

inline Float mul_approx(Float x1, Float x2) { return mul(x1,x2,to_nearest); }
inline Float mul_near(Float x1, Float x2) { return mul(x1,x2,to_nearest); }
inline Float mul_up(Float x1, Float x2) { return mul(x1,x2,upward); }
inline Float mul_down(Float x1, Float x2) { return mul(x1,x2,downward); }

inline Float div_approx(Float x1, Float x2) { return div(x1,x2,to_nearest); }
inline Float div_near(Float x1, Float x2) { return div(x1,x2,to_nearest); }
inline Float div_up(Float x1, Float x2) { return div(x1,x2,upward); }
inline Float div_down(Float x1, Float x2) { return div(x1,x2,downward); }

inline Float pow_approx(Float x, Int n) { return pow(x,n,to_nearest); }
inline Float pow_up(Float x, Int n) { return pow(x,n,upward); }
inline Float pow_down(Float x, Int n) { return pow(x,n,downward); }

inline Float next_down(Float x) { return sub_down(x,Float::min()); }
inline Float next_up(Float x) { return add_up(x,Float::min()); }

inline Float rec_approx(Float x) { return rec(x,to_nearest); }
inline Float rec_up(Float x) { return rec(x,upward); }
inline Float rec_down(Float x) { return rec(x,downward); }

inline Float sqrt_up(Float const& x) { return sqrt(x,upward); }
inline Float sqrt_down(Float const& x) { return sqrt(x,downward); }

inline Float exp_up(Float const& x) { return exp(x,upward); }
inline Float exp_down(Float const& x) { return exp(x,downward); }

inline Float log_up(Float const& x) { return log(x,upward); }
inline Float log_down(Float const& x) { return log(x,downward); }

inline Float sin_up(Float const& x) { return sin(x,upward); }
inline Float sin_down(Float const& x) { return sin(x,downward); }

inline Float cos_up(Float const& x) { return cos(x,upward); }
inline Float cos_down(Float const& x) { return cos(x,downward); }


/*
inline Float max_exact(Float const& x1, Float const& x2) { return max(x1,x2); }
inline Float min_exact(Float const& x1, Float const& x2) { return min(x1,x2); }
inline Float abs_exact(Float const& x) { return abs(x); }

inline Float add(Float const& x1, Float const& x2, Float::RoundingModeType rnd) { return apply(Add(),x1,x2,rnd); }
inline Float sub(Float const& x1, Float const& x2, Float::RoundingModeType rnd) { return apply(Sub(),x1,x2,rnd); }
inline Float mul(Float const& x1, Float const& x2, Float::RoundingModeType rnd) { return apply(Mul(),x1,x2,rnd); }
inline Float div(Float const& x1, Float const& x2, Float::RoundingModeType rnd) { return apply(Div(),x1,x2,rnd); }
inline Float pow(Float const& x, Int n, Float::RoundingModeType rnd) { return apply(Pow(),x,n,rnd); }

inline Float nul_exact(Float const& x) { return nul(x); }
inline Float pos_exact(Float const& x) { return pos(x); }
inline Float neg_exact(Float const& x) { return neg(x); }
inline Float half_exact(Float const& x) { return half(x); }

inline Float add_up(Float const& x1, Float const& x2) { return add(x1,x2,Float::upward); }
inline Float add_down(Float const& x1, Float const& x2) { return add(x1,x2,Float::downward); }
inline Float add_near(Float const& x1, Float const& x2) { return add(x1,x2,Float::to_nearest); }
inline Float add_approx(Float const& x1, Float const& x2) { return add(x1,x2,Float::to_nearest); }

inline Float sub_up(Float const& x1, Float const& x2) { return sub(x1,x2,Float::upward); }
inline Float sub_down(Float const& x1, Float const& x2) { return sub(x1,x2,Float::downward); }
inline Float sub_near(Float const& x1, Float const& x2) { return sub(x1,x2,Float::to_nearest); }
inline Float sub_approx(Float const& x1, Float const& x2) { return sub(x1,x2,Float::to_nearest); }

inline Float mul_up(Float const& x1, Float const& x2) { return mul(x1,x2,Float::upward); }
inline Float mul_down(Float const& x1, Float const& x2) { return mul(x1,x2,Float::downward); }
inline Float mul_near(Float const& x1, Float const& x2) { return mul(x1,x2,Float::to_nearest); }
inline Float mul_approx(Float const& x1, Float const& x2) { return mul(x1,x2,Float::to_nearest); }

inline Float div_up(Float const& x1, Float const& x2) { return div(x1,x2,Float::upward); }
inline Float div_down(Float const& x1, Float const& x2) { return div(x1,x2,Float::downward); }
inline Float div_near(Float const& x1, Float const& x2) { return div(x1,x2,Float::to_nearest); }
inline Float div_approx(Float const& x1, Float const& x2) { return div(x1,x2,Float::to_nearest); }

inline Float fma_up(Float const& x1, Float const& x2, Float const& y) { return add(mul(x1,x2,Float::upward),y,upward); }
inline Float fma_down(Float const& x1, Float const& x2, Float const& y) { return add(mul(x1,x2,Float::downward),y,downward); }
inline Float fma_approx(Float const& x1, Float const& x2, Float const& y) { return add(mul(x1,x2,Float::to_nearest),y,to_nearest); }

inline Float next_up(Float const& x) { return add_up(x,Float::min()); }
inline Float next_down(Float const& x) { return sub_down(x,Float::min()); }

inline Float rad_up(Float const& x1, Float const& x2) { return half(sub_up(x2,x1)); }
inline Float med_near(Float const& x1, Float const& x2) { return half(add_near(x2,x1)); }

inline Float pow_up(Float const& x, Int n) { return pow(x,n,Float::upward); }
inline Float pow_down(Float const& x, Int n) { return pow(x,n,Float::downward); }
inline Float pow_approx(Float const& x, Int n) { return pow(x,n,Float::to_nearest); }

Float pi_up();
Float pi_down();
Float pi_approx();

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
