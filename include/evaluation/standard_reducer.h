/***************************************************************************
 *            standard_reducer.h
 *
 *  Copyright  2007-8  Pieter Collins
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
 
/*! \file orthogonal_reducer.h
 *  \brief Methods for reducing the description of basic sets using orthogonal over-approximation
 */

#ifndef ARIADNE_STANDARD_REDUCER_H
#define ARIADNE_STANDARD_REDUCER_H

#include "geometry/zonotope.h"
#include "evaluation/reducer_interface.h"
#include "evaluation/orthogonal_reducer.h"

namespace Ariadne {
  
template<class ES> class StandardReducer;

/*! \ingroup Approximators
 *  \brief Default class for over-approximating a zonotope. Dispactches to orthogonal reducer.
 *
 */ 
template<class R>
class StandardReducer< Zonotope<R> >
  : public OrthogonalReducer< Zonotope<R> > 
{
  typedef Zonotope<R> ES;
 public:
  virtual StandardReducer<ES>* clone() const { return new StandardReducer<ES>(*this); }
};

  
} // namespace Ariadne

#endif /* ARIADNE_STANDARD_REDUCER_H */
