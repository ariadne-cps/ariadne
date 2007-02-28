/***************************************************************************
 *            rectangle.h
 *
 *  Mon 2 May 2005
 *  Copyright 2005  Alberto Casagrande, Pieter Collins
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
 
/*! \file rectangle_expression.h
 *  \brief Base class for rectangle expression templates.
 */

#ifndef ARIADNE_RECTANGLE_EXPRESSION_H
#define ARIADNE_RECTANGLE_EXPRESSION_H

#include <iosfwd>

#include "../declarations.h"

namespace Ariadne {
  namespace Geometry {

    /*! \brief %Base class tag for rectangle expressions
     *
     * Note that this class is \em not itself a model of a
     * %RectangleExpression, but must be a base class of any 
     * %RectangleExpression. */
    template<class E>
    class RectangleExpression 
    {
     public:
      /*!\brief Convert \a *this to a reference to E. */
      E& operator() () { return static_cast<E&>(*this); }
      /*!\brief Convert \a *this to a constant reference to E. */
      const E& operator() () const { return static_cast<const E&>(*this); }
    };

    
  }
}

#endif /* ARIADNE_RECTANGLE_EXPRESSION_H */
