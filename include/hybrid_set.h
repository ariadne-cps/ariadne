/***************************************************************************
 *            hybrid_set.h
 *
 *  Copyright 2008  Pieter Collins
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
 
/*! \file hybrid_set.h
 *  \brief Sets in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_SET_H
#define ARIADNE_HYBRID_SET_H

#include <map>
#include "macros.h"
#include "grid.h"

namespace Ariadne {

typedef uint DiscreteState;

class GridSet { 
  template<class S> void adjoin(const S&); 
  template<class S> void adjoin_outer_approximation(const S&); 
};

class HybridGridSet 
  : public std::map<DiscreteState,GridSet>
{
  bool has_location(DiscreteState q) const {
    return this->find(q)!=this->end(); }

  template<class S> void adjoin(DiscreteState q, const S& s) {
    this->operator[](q).adjoin(s); }

  template<class S> void adjoin_outer_approximation(DiscreteState q, const S& s) {
    this->operator[](q).adjoin_outer_approximation(s); }

  const GridSet& operator[](DiscreteState q) const {
    ARIADNE_ASSERT(this->find(q)!=this->end());
    return const_cast<HybridGridSet*>(this)->operator[](q);
  }

};

} // namespace Ariadne

#endif
