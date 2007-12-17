/***************************************************************************
 *            geometry/expression.h
 *
 *  Copyright 2007  Pieter Collins
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
 
/*! \file geometry/expression.h
 *  \brief Base class for rectangle expression templates.
 */

#ifndef ARIADNE_GEOMETRY_EXPRESSION_H
#define ARIADNE_GEOMETRY_EXPRESSION_H

namespace Ariadne {
  namespace Geometry {

    class OverApproximate 
    {
     public:
      template<class R, class A> void operator() (R& r, const A& a) {
        over_approximate_(r,a);
      }
    };

  }
}

#endif /* ARIADNE_GEOMETRY_EXPRESSION_H */
