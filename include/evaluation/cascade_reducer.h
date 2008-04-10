/***************************************************************************
 *            cascade_reducer.h
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
 
/*! \file cascade_reducer.h
 *  \brief Methods for reducing the description of basic sets
 */

#ifndef ARIADNE_CASCADE_REDUCER_H
#define ARIADNE_CASCADE_REDUCER_H

#include "reducer_interface.h"
#include "geometry/zonotope.h"

namespace Ariadne {
  namespace Evaluation {

    template<class ES> class CascadeReducer;

    /*! \ingroup Approximators
     *  \brief Class for over-approximating a zonotope using Kuhn's cascade method.
     */ 
    template<class R>
    class CascadeReducer< Geometry::Zonotope<R> >
      : public ReducerInterface< Geometry::Zonotope<R> > 
    {
      typedef Geometry::Zonotope<R> ES;
     public:
      CascadeReducer<ES>(size_type maximum_number_of_blocks) : _max_blocks(maximum_number_of_blocks) { }
      virtual CascadeReducer<ES>* clone() const { return new CascadeReducer<ES>(*this); }
      virtual ES over_approximate(const ES& es) const { return Geometry::cascade_over_approximation(es,this->_max_blocks); }
     private:
      size_type _max_blocks;
    };



  }
}

#endif /* ARIADNE_CASCADE_REDUCER_H */
