/***************************************************************************
 *            geometry/interfaces.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file geometry/interfaces.h
 *  \brief Interfaces for function objects implementing standard geometrical operations.
 */

#ifndef ARIADNE_GEOMETRY_INTERFACES_H
#define ARIADNE_GEOMETRY_INTERFACES_H

#include "geometry/declarations.h"

namespace Ariadne { 
  
    
    template<class R, class BS> 
    RadiusInterface
    {
      RadiusInterface<R,BS>* clone() const = 0;
      typename R operator() (const BS&) const = 0;
    };

    template<class R, class BS> 
    BoundingBoxInterface
    {
      BoundingBoxInterface<R,BS>* clone() const = 0;
      typename Box<R> operator() (const BS&) const = 0;
    };

    template<class BS1, class BS2> 
    OverApproximationInterface
    {
      OverApproximateInterface<BS1,BS2>* clone() const = 0;
      BS1 operator() (const BS2&) const = 0;
    };

    template<class BS> 
    SubdivideInterface
    {
      SubdivideInterface<BS>* clone() const = 0;
      ListSet<BS> operator() (const BS2&) const = 0;
    };

    template<class R, class BS> 
    InnerApproximationInterface
    {
      InnerApproximationInterface<R>* clone() const = 0;
      GridCellListSet<R> operator() (const BS&, const Grid<R>&) const = 0;
    };

    template<class R, class BS> 
    OuterApproximationInterface
    {
      OuterApproximationInterface<R>* clone() const = 0;
      GridCellListSet<R> operator() (const BS&, const Grid<R>&) const = 0;
    };

    template<class BS, class S> 
    SubsetInterface
    {
      SubsetInterface<R>* clone() const = 0;
      tribool operator() (const BS&, const S&) const = 0;
    };

    template<class BS, class S> 
    IntersectInterface
    {
      IntersectInterface<R>* clone() const = 0;
      tribool operator() (const BS&, const S&) const = 0;
    };

    template<class BS, class S> 
    DisjointInterface
    {
      DisjointInterface<R>* clone() const = 0;
      tribool operator() (const BS&, const S&) const = 0;
    };


  
} // namespace Ariadne

#endif /* ARIADNE_GEOMETRY_INTERFACES_H */
