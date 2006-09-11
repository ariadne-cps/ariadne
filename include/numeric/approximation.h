/***************************************************************************
 *            approximation.h
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
 
/*! \file approximation.h
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
    template<typename R,typename A> R conv_exact(const A& x);
    /*! \brief Approximate \a x by an element of R. */
    template<typename R,typename A> R conv_approx(const A& x);
    /*! \brief Approximate \a x by an element of R, rounding down. */
    template<typename R,typename A> R conv_down(const A& x);
    /*! \brief Approximate \a x by an element of R, rounding up. */
    template<typename R,typename A> R conv_up(const A& x);
    //@}

  }    
}

#endif /* _ARIADNE_APPROXIMATION_H */
