/***************************************************************************
 *            rounding.hpp
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

/*! \file rounding.hpp
 *  \brief Functions to set and retrieve the processor rounding mode.
 *  May be platform-dependent.
 */

#ifndef ARIADNE_ROUNDING_HPP
#define ARIADNE_ROUNDING_HPP


#if defined __GNUC__ && ( defined __i386__ || defined __x86_64 || defined _M_IX86 || defined _M_X86 )
    #if ( defined __SSE_MATH__ &&  defined __SSE2__ )
        #define ARIADNE_SSE_ROUNDING
    #elif __GNUC__ >= 5 || ( __GNUC__ == 4 && __GNUC_MINOR__ >= 3 )
        #define ARIADNE_GCC_ROUNDING
    #else
        #define ARIADNE_C99_ROUNDING
    #endif
#else
    #define ARIADNE_BOOST_ROUNDING
#endif


//#undef ARIADNE_SSE_ROUNDING
//#undef ARIADNE_GCC_ROUNDING
//#undef ARIADNE_C99_ROUNDING

//#define ARIADNE_GCC_ROUNDING


#if defined ARIADNE_SSE_ROUNDING

#include <xmmintrin.h>

namespace Ariadne {

typedef unsigned short rounding_mode_t;

const rounding_mode_t to_nearest = _MM_ROUND_NEAREST;
const rounding_mode_t downward = _MM_ROUND_DOWN;
const rounding_mode_t upward = _MM_ROUND_UP;
const rounding_mode_t toward_zero = _MM_ROUND_TOWARD_ZERO;

const rounding_mode_t ROUND_NEAR = to_nearest;
const rounding_mode_t ROUND_DOWN = downward;
const rounding_mode_t ROUND_UP   = upward;
const rounding_mode_t ROUND_ZERO = toward_zero;

inline void _set_rounding_to_nearest() { _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);  }
inline void _set_rounding_downward() { _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);  }
inline void _set_rounding_upward() { _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);  }
inline void _set_rounding_toward_zero() { _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);  }

inline void _set_rounding_mode(rounding_mode_t rnd) { _MM_SET_ROUNDING_MODE(rnd); }
inline rounding_mode_t _get_rounding_mode() { return _MM_GET_ROUNDING_MODE(); }

enum class RoundingMode : rounding_mode_t {
    ROUND_NEAR = to_nearest, ROUND_DOWN = downward, ROUND_UP   = upward
};

} // namespace Ariadne



#elif defined ARIADNE_C99_ROUNDING

#include <fenv.h>

namespace Ariadne {

typedef unsigned short rounding_mode_t;

const rounding_mode_t ROUND_NEAR = FE_TONEAREST;
const rounding_mode_t ROUND_DOWN = FE_DOWNWARD;
const rounding_mode_t ROUND_UP = FE_UPWARD;
const rounding_mode_t ROUND_ZERO = FE_TOWARDZERO;

const rounding_mode_t to_nearest = FE_TONEAREST;
const rounding_mode_t downward = FE_DOWNWARD;
const rounding_mode_t upward = FE_UPWARD;
const rounding_mode_t toward_zero = FE_TOWARDZERO;

inline void _set_rounding_to_nearest() { fesetround(FE_TONEAREST);  }
inline void _set_rounding_downward() { fesetround(FE_DOWNWARD);  }
inline void _set_rounding_upward() { fesetround(FE_UPWARD);  }
inline void _set_rounding_toward_zero() { fesetround(FE_TOWARDZERO);  }

inline void _set_rounding_mode(rounding_mode_t rnd) { fesetround(rnd); }
inline rounding_mode_t _get_rounding_mode() { return fegetround(); }

} // namespace Ariadne



#elif defined ARIADNE_BOOST_ROUNDING

#include <boost/numeric/interval/hw_rounding.hpp>

namespace Ariadne {

typedef boost::numeric::interval_lib::rounding_control<double>::rounding_mode rounding_mode_t;

class RoundToNearest { };
class RoundDownward { };
class RoundUpward { };
class RoundTowardZero { };

const RoundToNearest  to_nearest=RoundToNearest();
const RoundDownward   downward=RoundDownward();
const RoundUpward     upward=RoundUpward();
const RoundTowardZero toward_zero=RoundTowardZero();

inline void Float64::set_rounding_to_nearest() { boost::numeric::interval_lib::rounding_control<double>::to_nearest(); }
inline void Float64::set_rounding_downward() { boost::numeric::interval_lib::rounding_control<double>::downward(); }
inline void Float64::set_rounding_upward() { boost::numeric::interval_lib::rounding_control<double>::upward(); }
inline void Float64::set_rounding_toward_zero() { boost::numeric::interval_lib::rounding_control<double>::toward_zero(); }

inline void Float64::set_rounding_mode(RoundToNearest) { boost::numeric::interval_lib::rounding_control<double>::to_nearest(); }
inline void Float64::set_rounding_mode(RoundDownward) { boost::numeric::interval_lib::rounding_control<double>::downward(); }
inline void Float64::set_rounding_mode(RoundUpward) { boost::numeric::interval_lib::rounding_control<double>::upward(); }
inline void Float64::set_rounding_mode(RoundTowardZero) { boost::numeric::interval_lib::rounding_control<double>::toward_zero(); }

inline void Float64::set_rounding_mode(rounding_mode_t rnd) { boost::numeric::interval_lib::rounding_control<double>::Float64::set_rounding_mode(rnd); }
inline rounding_mode_t Float64::get_rounding_mode() { rounding_mode_t rnd; boost::numeric::interval_lib::rounding_control<double>::Float64::get_rounding_mode(rnd); return rnd; }

} // namespace Ariadne



#elif defined ARIADNE_GCC_ROUNDING

namespace Ariadne {

typedef unsigned short rounding_mode_t;

const rounding_mode_t ROUND_NEAR = 895;
const rounding_mode_t ROUND_DOWN = 895+1024;
const rounding_mode_t ROUND_UP = 895+2048;
const rounding_mode_t ROUND_ZERO = 895+3072;

const rounding_mode_t to_nearest   = ROUND_NEAR;
const rounding_mode_t downward     = ROUND_DOWN;
const rounding_mode_t upward       = ROUND_UP;
const rounding_mode_t toward_zero  = ROUND_ZERO;

inline rounding_mode_t _get_rounding_mode() { rounding_mode_t rnd; asm volatile ("fstcw %0" : "=m" (rnd) ); return rnd; }
inline void _get_rounding_mode(rounding_mode_t& rnd) { asm volatile ("fstcw %0" : "=m" (rnd) ); }
inline void _set_rounding_mode(rounding_mode_t rnd) { asm volatile ("fldcw %0" : : "m" (rnd) ); }
inline void _set_rounding_to_nearest() { asm volatile ("fldcw %0" : : "m" (ROUND_NEAR) ); }
inline void _set_rounding_downward() { asm volatile ("fldcw %0" : : "m" (ROUND_DOWN) ); }
inline void _set_rounding_upward() { asm volatile ("fldcw %0" : : "m" (ROUND_UP) ); }
inline void _set_rounding_toward_zero() { asm volatile ("fldcw %0" : : "m" (ROUND_ZERO) ); }

} // namespace Ariadne



#elif defined ARIADNE_MSVC_ROUNDING

static const unsigned short ARIADNE_FENV_BASE = 895;
static unsigned short ARIADNE_ROUND_TMP = ARIADNE_FENV_BASE;

namespace Ariadne {

typedef unsigned short rounding_mode_t;

const rounding_mode_t to_nearest   = ARIADNE_FENV_BASE;
const rounding_mode_t downward     = ARIADNE_FENV_BASE + 1024;
const rounding_mode_t upward       = ARIADNE_FENV_BASE + 2048;
const rounding_mode_t toward_zero  = ARIADNE_FENV_BASE + 3072;

inline void Float64::set_rounding_to_nearest() { __asm fldcw to_nearest; }
inline void Float64::set_rounding_downward() { __asm fldcw downward; }
inline void Float64::set_rounding_upward() { __asm fldcw upward; }
inline void Float64::set_rounding_toward_zero() { __asm fldcw toward_zero; }

inline void Float64::set_rounding_mode(rounding_mode_t rnd) { ARIADNE_RND_TMP=rnd; __asm fldcw ARIADNE_ROUND_TMP; }
inline rounding_mode_t Float64::get_rounding_mode() { __asm fstcw ARIADNE_ROUND_TMP; ARIADNE_RND_TMP; }

} // namespace Ariadne



#else // No rounding

namespace Ariadne {

typedef unsigned short rounding_mode_t;

const rounding_mode_t to_nearest   = 0000;
const rounding_mode_t downward     = 1024;
const rounding_mode_t upward       = 2048;
const rounding_mode_t toward_zero  = 3072;

inline void Float64::set_rounding_to_nearest() { }
inline void Float64::set_rounding_downward() { }
inline void Float64::set_rounding_upward() { }
inline void Float64::set_rounding_toward_zero() { }

inline void Float64::set_rounding_mode(rounding_mode_t rnd) { }
inline rounding_mode_t Float64::get_rounding_mode() { return 0 }

} // namespace Ariadne

#endif

/************  Publicly-accessible rounding-mode changing *******************/

namespace Ariadne {

//! \ingroup NumericModule
//! \brief The integral type used to represent the rounding mode.
typedef unsigned short RoundingModeType;

//! \brief The floating-point environment value for rounding arithmetic to the nearest exactly-representable value.
extern const RoundingModeType to_nearest;
extern const RoundingModeType ROUND_NEAR;
//! \brief The floating-point environment value for upwards-rounded arithmetic.
extern const RoundingModeType downward;
extern const RoundingModeType ROUND_DOWN;
//! \brief The floating-point environment value for downwards-rounded arithmetic.
extern const RoundingModeType downward;
extern const RoundingModeType ROUND_UP;
//! \brief The floating-point environment value for rounding arithmetic to zero.
extern const RoundingModeType toward_zero;
extern const RoundingModeType ROUND_ZERO;

//! \brief Set the builtin rounding mode. \ingroup NumericModule
void set_rounding_mode(RoundingModeType rnd);
//! \brief Get the current rounding mode \ingroup NumericModule
RoundingModeType get_rounding_mode();

//! \brief Set the rounding mode to nearest. \ingroup NumericModule
void set_rounding_to_nearest();
//! \brief Set the rounding mode to downwards rounding. \ingroup NumericModule
void set_rounding_downward();
//! \brief Set the rounding mode to upwards rounding. \ingroup NumericModule
void set_rounding_upward();
//! \brief Set the rounding mode to towards-zero rounding. \ingroup NumericModule
void set_rounding_toward_zero();

//! \brief Set the rounding mode to the expected default rounding mode.
void set_default_rounding();



} // namespace Ariadne

#endif // ARIADNE_ROUNDING_HPP

