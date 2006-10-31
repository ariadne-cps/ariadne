/***************************************************************************
 *            conversion.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file conversion.h
 *  \brief Approximation of numerical types.
 */

#ifndef _ARIADNE_APPROXIMATION_H
#define _ARIADNE_APPROXIMATION_H

namespace Ariadne {
  namespace Numeric {
    
    //! \name Conversion operations. 
    //@{ 
    //! \ingroup Numeric
    /*! \brief Convert \a x to an element of R. */
    template<class R,class A> inline R conv_exact(const A& x) { return static_cast<R>(x); }
    /*! \brief Approximate \a x by an element of R. */
    template<class R,class A> inline R conv_approx(const A& x) { return conv_exact<R>(x); };
    /*! \brief Approximate \a x by an element of R, rounding down. */
    template<class R,class A> inline R conv_down(const A& x) { return conv_exact<R>(x); };
    /*! \brief Approximate \a x by an element of R, rounding up. */
    template<class R,class A> inline R conv_up(const A& x) { return conv_exact<R>(x); };

    /*! \brief Approximate \a x by an integer of type N, rounding down. */
    template<class N,class A> inline N int_down(const A& x);
    /*! \brief Approximate \a x by an integer of type N, rounding up. */
    template<class N,class A> inline N int_up(const A& x);
    //@}

  }    
}

#endif /* _ARIADNE_APPROXIMATION_H */
