/***************************************************************************
 *            identity_reducer.h
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
 
/*! \file identity_reducer.h
 *  \brief Methods for reducing the description of basic sets
 */

#ifndef ARIADNE_IDENTITY_REDUCER_H
#define ARIADNE_IDENTITY_REDUCER_H

#include "reducer_interface.h"

namespace Ariadne {
  

    /*! \ingroup Approximators
     *  \brief Class for over-approximating an enclosure set by a simpler set. Returns the original set.
     */ 
    template<class ES>
    class IdentityReducer
      : public ReducerInterface<ES> {
     public:
      virtual IdentityReducer<ES>* clone() const { return new IdentityReducer<ES>(*this); }
      virtual ES over_approximate(const ES& es) const { return es; }
    };

 
  
} // namespace Ariadne

#endif /* ARIADNE_IDENTITY_REDUCER_H */
