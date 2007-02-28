/***************************************************************************
 *            comparison.h
 *
 *  Copyright 2006  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file comparison.h
 *  \brief Comparison operators for numerical types.
 */
 
#ifndef ARIADNE_COMPARISON_H
#define ARIADNE_COMPARISON_H

namespace Ariadne {
  namespace Numeric {

    //! \name Comparison operators. 
    //@{
    //! \ingroup Numeric
    /*! \brief Equal to. */
    template<class R1, class R2> inline bool operator==(const R1& x1, const R2& x2);

    /*! \brief Not equal to. */
    template<class R1, class R2> inline bool operator!=(const R1& x1, const R2& x2);
      
    /*! \brief Strictly less than. */
    template<class R1, class R2> inline bool operator< (const R1& x1, const R2& x2);
      
    /*! \brief Strictly greater than. */
    template<class R1, class R2> inline bool operator> (const R1& x1, const R2& x2);
      
    /*! \brief Less than or equal to. */
    template<class R1, class R2> inline bool operator<=(const R1& x1, const R2& x2);
      
    /*! \brief Greater than or equal to. */
    template<class R1, class R2> inline bool operator>=(const R1& x1, const R2& x2);
    
    //!}

  }
}



#endif /* ARIADNE_COMPARISON_H */
