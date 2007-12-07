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

namespace Ariadne {
  namespace Numeric {
        
    template<class RM> class Round { };
    class Up; class Down; class Approx; class Exact;

    typedef Round<Down> RoundDown;
    typedef Round<Up> RoundUp;
    typedef Round<Approx> RoundApprox;
    typedef Round<Exact> RoundExact;

#ifdef DOXYGEN
    //! \name Rounding classes.
    //@{
    //! \brief Round downwards
    class RoundDown { };
    //! \brief Round upwards
    class RoundUp { };
    //! \brief Round without control on error
    class RoundApprox { };
    //! \brief No rounding allowed
    class RoundExact { };
    //@}
#endif

    static const RoundUp round_up=RoundUp();    
    static const RoundDown round_down=RoundDown();
    static const RoundApprox round_approx=RoundApprox();; 
    static const RoundExact round_exact=RoundExact(); 

  }
}
 
#endif /* ARIADNE_ROUNDING_H */
