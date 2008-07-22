/***************************************************************************
 *            geometry/box_expression.h
 *
 *  
 *  Copyright 2005-6  Pieter Collins
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
 
/*! \file geometry/box_expression.h
 *  \brief Base class for box expression templates.
 */

#ifndef ARIADNE_GEOMETRY_BOX_EXPRESSION_H
#define ARIADNE_GEOMETRY_BOX_EXPRESSION_H

namespace Ariadne {
  

    /*! \brief %Base class tag for box expressions
     *
     * Note that this class is \em not itself a model of a
     * %BoxExpression, but must be a base class of any 
     * %BoxExpression. */
    template<class E>
    class BoxExpression 
    {
     public:
      /*!\brief Convert \a *this to a reference to E. */
      E& operator() () { return static_cast<E&>(*this); }
      /*!\brief Convert \a *this to a constant reference to E. */
      const E& operator() () const { return static_cast<const E&>(*this); }
    };


} // namespace Ariadne

#endif /* ARIADNE_GEOMETRY_BOX_EXPRESSION_H */
