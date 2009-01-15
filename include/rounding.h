/***************************************************************************
 *            rounding.h
 *
 *  Copyright 2008  Pieter Collins
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
 
/*! \file rounding.h
 *  \brief Functions to set and retrieve the processor rounding mode. 
 *  May be platform-dependent.
 */

#ifndef ARIADNE_ROUNDING_H
#define ARIADNE_ROUNDING_H

#ifdef DOXYGEN
namespace Ariadne {
//! \brief The unsigned integral type used to represent the rounding mode.
typedef unsigned short rounding_mode_t;

//! \brief The floating-point environment value for rounding arithmetic to the nearest exactly-representable value.
const rounding_mode_t to_nearest;
//! \brief The floating-point environment value for upwards-rounded arithmetic.
const rounding_mode_t downward;
//! \brief The floating-point environment value for downwards-rounded arithmetic.
const rounding_mode_t upward;
//! \brief The floating-point environment value for rounding arithmetic to zero.
const rounding_mode_t toward_zero;

//! \brief Set the rounding mode to nearest.
inline void set_rounding_to_nearest();
//! \brief Set the rounding mode to downwards rounding.
inline void set_rounding_downward();
//! \brief Set the rounding mode to upwards rounding.
inline void set_rounding_upward();
//! \brief Set the rounding mode to towards-zero rounding.
inline void set_rounding_toward_zero();

//! \brief Set the rounding mode.
inline void set_rounding_mode(const rounding_mode_t& rnd);
//! \brief Get the current rounding mode.
inline rounding_mode_t get_rounding_mode();
}
#endif

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



#if defined ARIADNE_SSE_ROUNDING

#include <xmmintrin.h>

namespace Ariadne {

typedef unsigned int rounding_mode_t;

const rounding_mode_t to_nearest = _MM_ROUND_NEAREST;
const rounding_mode_t downward = _MM_ROUND_DOWN;
const rounding_mode_t upward = _MM_ROUND_UP;
const rounding_mode_t toward_zero = _MM_ROUND_TOWARD_ZERO;

inline void set_rounding_to_nearest() { _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);  }
inline void set_rounding_downward() { _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);  }
inline void set_rounding_upward() { _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);  }
inline void set_rounding_toward_zero() { _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);  }

inline void set_rounding_mode(rounding_mode_t rnd) { _MM_SET_ROUNDING_MODE(rnd); }
inline rounding_mode_t get_rounding_mode() { return _MM_GET_ROUNDING_MODE(); }

} // namespace Ariadne



#elif defined ARIADNE_C99_ROUNDING

#include <fenv.h>

namespace Ariadne {

typedef unsigned short rounding_mode_t;

const rounding_mode_t to_nearest = FE_TONEAREST;
const rounding_mode_t downward = FE_DOWNWARD;
const rounding_mode_t upward = FE_UPWARD;
const rounding_mode_t toward_zero = FE_TOWARDZERO;

inline void set_rounding_to_nearest() { fesetround(FE_TONEAREST);  }
inline void set_rounding_downward() { fesetround(FE_DOWNWARD);  }
inline void set_rounding_upward() { fesetround(FE_UPWARD);  }
inline void set_rounding_toward_zero() { fesetround(FE_TOWARDZERO);  }

inline void set_rounding_mode(rounding_mode_t rnd) { fesetround(rnd); }
inline rounding_mode_t get_rounding_mode() { return fegetround(); }

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

inline void set_rounding_to_nearest() { boost::numeric::interval_lib::rounding_control<double>::to_nearest(); }
inline void set_rounding_downward() { boost::numeric::interval_lib::rounding_control<double>::downward(); }
inline void set_rounding_upward() { boost::numeric::interval_lib::rounding_control<double>::upward(); }
inline void set_rounding_toward_zero() { boost::numeric::interval_lib::rounding_control<double>::toward_zero(); }

inline void set_rounding_mode(RoundToNearest) { boost::numeric::interval_lib::rounding_control<double>::to_nearest(); }
inline void set_rounding_mode(RoundDownward) { boost::numeric::interval_lib::rounding_control<double>::downward(); }
inline void set_rounding_mode(RoundUpward) { boost::numeric::interval_lib::rounding_control<double>::upward(); }
inline void set_rounding_mode(RoundTowardZero) { boost::numeric::interval_lib::rounding_control<double>::toward_zero(); }

inline void set_rounding_mode(rounding_mode_t rnd) { boost::numeric::interval_lib::rounding_control<double>::set_rounding_mode(rnd); }
inline rounding_mode_t get_rounding_mode() { rounding_mode_t rnd; boost::numeric::interval_lib::rounding_control<double>::get_rounding_mode(rnd); return rnd; }

} // namespace Ariadne



#elif defined ARIADNE_GCC_ROUNDING

const unsigned short ARIADNE_FENV_BASE = 895;
const unsigned short ARIADNE_ROUND_TO_NEAREST = ARIADNE_FENV_BASE + 0000;
const unsigned short ARIADNE_ROUND_DOWNWARD = ARIADNE_FENV_BASE + 1024;
const unsigned short ARIADNE_ROUND_UPWARD = ARIADNE_FENV_BASE + 2048;
const unsigned short ARIADNE_ROUND_TOWARD_ZERO = ARIADNE_FENV_BASE + 3072;

namespace Ariadne {

typedef unsigned short rounding_mode_t;

const rounding_mode_t to_nearest = ARIADNE_ROUND_TO_NEAREST;
const rounding_mode_t downward = ARIADNE_ROUND_DOWNWARD;
const rounding_mode_t upward = ARIADNE_ROUND_UPWARD;
const rounding_mode_t toward_zero = ARIADNE_ROUND_TOWARD_ZERO;

//inline void set_rounding_mode(rounding_mode_t rnd) { asm volatile ("fldcw %0" : : "m" (rnd) ); }
inline void set_rounding_mode(const rounding_mode_t& rnd) { asm volatile ("fldcw %0" : : "m" (rnd) ); }
inline void get_rounding_mode(rounding_mode_t& rnd) { asm volatile ("fstcw %0" : "=m" (rnd) ); }
inline rounding_mode_t get_rounding_mode() { rounding_mode_t rnd; get_rounding_mode(rnd); return rnd; }

//inline void set_nearest() { asm volatile ("fldcw to_nearest"); }
inline void set_rounding_to_nearest() { set_rounding_mode(to_nearest); }
inline void set_rounding_downward() { set_rounding_mode(downward); }
inline void set_rounding_upward() { set_rounding_mode(upward); }
inline void set_rounding_toward_zero() { set_rounding_mode(toward_zero); }
} // namespace Ariadne



#elif defined ARIADNE_MSVC_ROUNDING

const unsigned short ARIADNE_FENV_BASE = 895;
unsigned short ARIADNE_ROUND_TMP = ARIADNE_FENV_BASE;

namespace Ariadne {

typedef unsigned short rounding_mode_t;

const rounding_mode_t to_nearest   = ARIADNE_FENV_BASE;
const rounding_mode_t downward     = ARIADNE_FENV_BASE + 1024;
const rounding_mode_t upward       = ARIADNE_FENV_BASE + 2048;
const rounding_mode_t toward_zero  = ARIADNE_FENV_BASE + 3072;

inline void set_rounding_to_nearest() { __asm fldcw to_nearest; }
inline void set_rounding_downward() { __asm fldcw downward; }
inline void set_rounding_upward() { __asm fldcw upward; }
inline void set_rounding_toward_zero() { __asm fldcw toward_zero; }

inline void set_rounding_mode(rounding_mode_t rnd) { ARIADNE_RND_TMP=rnd; __asm fldcw ARIADNE_ROUND_TMP; }
inline rounding_mode_t get_rounding_mode() { __asm fstcw ARIADNE_ROUND_TMP; ARIADNE_RND_TMP; }

} // namespace Ariadne



#else // No rounding

namespace Ariadne {

typedef unsigned short rounding_mode_t;

const rounding_mode_t to_nearest   = 0000;
const rounding_mode_t downward     = 1024;
const rounding_mode_t upward       = 2048;
const rounding_mode_t toward_zero  = 3072;

inline void set_rounding_to_nearest() { }
inline void set_rounding_downward() { }
inline void set_rounding_upward() { }
inline void set_rounding_toward_zero() { }

inline void set_rounding_mode(rounding_mode_t rnd) { }
inline rounding_mode_t get_rounding_mode() { return 0 }

}

#endif

#endif // ARIADNE_ROUNDING_H

