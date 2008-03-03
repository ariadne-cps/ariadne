/***************************************************************************
 *            rounding.h
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
 
/*! \file numeric/rounding.h
 *  \brief Rounding mode classes.
 */

#ifndef ARIADNE_ROUNDING_H
#define ARIADNE_ROUNDING_H

#include <fenv.h>
#include <mpfr.h>

namespace Ariadne {
namespace Numeric {
        
template<class RM> class Round { };
class Near; class Up; class Down; class Chop; class Approx; class Exact;

typedef Round<Near> RoundNear;
typedef Round<Down> RoundDown;
typedef Round<Up> RoundUp;
typedef Round<Chop> RoundChop;
typedef Round<Approx> RoundApprox;
typedef Round<Exact> RoundExact;

#ifdef DOXYGEN
//! \name Rounding classes.
//@{
//! \brief Round downwards
class RoundNear { };
//! \brief Round downwards
class RoundDown { };
//! \brief Round upwards
class RoundUp { };
//! \brief Round toward zero
class RoundChop { };
//! \brief Round without control on error
class RoundApprox { };
//! \brief No rounding allowed
class RoundExact { };
//@}
#endif

static const RoundNear round_near=RoundNear();    
static const RoundUp round_up=RoundUp();    
static const RoundDown round_down=RoundDown();
static const RoundChop round_chop=RoundChop();
static const RoundApprox round_approx=RoundApprox();; 
static const RoundExact round_exact=RoundExact(); 

typedef unsigned short rounding_mode_type;
const rounding_mode_type hardware_round_near = FE_TONEAREST;
const rounding_mode_type hardware_round_down = FE_DOWNWARD;
const rounding_mode_type hardware_round_up = FE_UPWARD;
const rounding_mode_type hardware_round_chop = FE_TOWARDZERO;
const rounding_mode_type hardware_round_approx = FE_TONEAREST;

inline rounding_mode_type hardware_rounding_mode(RoundNear) { return hardware_round_near; }
inline rounding_mode_type hardware_rounding_mode(RoundDown) { return hardware_round_down; }
inline rounding_mode_type hardware_rounding_mode(RoundUp) { return hardware_round_up; }
inline rounding_mode_type hardware_rounding_mode(RoundChop) { return hardware_round_chop; }
inline rounding_mode_type hardware_rounding_mode(RoundApprox) { return hardware_round_near; }

inline mpfr_rnd_t mpfr_rounding_mode(RoundNear) { return GMP_RNDN; }
inline mpfr_rnd_t mpfr_rounding_mode(RoundDown) { return GMP_RNDD; }
inline mpfr_rnd_t mpfr_rounding_mode(RoundUp) { return GMP_RNDU; }
inline mpfr_rnd_t mpfr_rounding_mode(RoundChop) { return GMP_RNDZ; }
inline mpfr_rnd_t mpfr_rounding_mode(RoundApprox) { return GMP_RNDN; }


inline rounding_mode_type get_rounding_mode() { return fegetround(); }
inline void set_rounding_mode(rounding_mode_type rnd) { fesetround(rnd); }
template<class RM> inline void set_rounding_mode(Round<RM> rnd) { fesetround(hardware_rounding_mode(rnd)); }

}
}
 
#endif /* ARIADNE_ROUNDING_H */
