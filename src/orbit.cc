/***************************************************************************
 *            orbit.cc
 *
 *  Copyright 2007-8  Pieter Collins
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
 
#include "orbit.h"

#include "approximate_taylor_model.h"
#include "list_set.h"
#include "hybrid_set.h"

namespace Ariadne {

struct Orbit<HybridTaylorSetType>::Data {
  Data(const HybridTaylorSetType& initial_set) 
    : initial(initial_set) { }
  HybridTaylorSetType initial;
  HybridTaylorSetListType reach;
  HybridTaylorSetListType intermediate;
  HybridTaylorSetListType final;
};

Orbit<HybridTaylorSetType>::
Orbit(const HybridTaylorSetType& initial_set)
  : _data(new Data(initial_set))
{
}

void
Orbit<HybridTaylorSetType>::
adjoin_reach(const HybridTaylorSetType& set)
{
  this->_data->reach.adjoin(set);
}

void
Orbit<HybridTaylorSetType>::
adjoin_intermediate(const HybridTaylorSetType& set)
{
  this->_data->intermediate.adjoin(set);
}

void
Orbit<HybridTaylorSetType>::
adjoin_final(const HybridTaylorSetType& set)
{
  this->_data->final.adjoin(set);
}


void
Orbit<HybridTaylorSetType>::
adjoin_reach(const HybridTaylorSetListType& list_set)
{
  this->_data->reach.adjoin(list_set);
}

void
Orbit<HybridTaylorSetType>::
adjoin_intermediate(const HybridTaylorSetListType& list_set)
{
  this->_data->intermediate.adjoin(list_set);
}

void
Orbit<HybridTaylorSetType>::
adjoin_final(const HybridTaylorSetListType& list_set)
{
  this->_data->final.adjoin(list_set);
}


HybridTaylorSetType const&
Orbit<HybridTaylorSetType>::
initial() const
{
  return this->_data->initial;
}

HybridTaylorSetListType const&
Orbit<HybridTaylorSetType>::
reach() const
{
  return this->_data->reach;
}

HybridTaylorSetListType const&
Orbit<HybridTaylorSetType>::
intermediate() const
{
  return this->_data->intermediate;
}

HybridTaylorSetListType const&
Orbit<HybridTaylorSetType>::
final() const
{
  return this->_data->final;
}





struct Orbit<HybridGridCell>::Data {
  Data(const HybridGridCell& initial_set) 
    : initial(initial_set) { }
  HybridGridCell initial;
  HybridGridTreeSet reach;
  HybridGridTreeSet intermediate;
  HybridGridTreeSet final;
};

Orbit<HybridGridCell>::
Orbit(const HybridGridCell& initial_set)
  : _data(new Data(initial_set))
{
}

Orbit<HybridGridCell>::
Orbit(const HybridGridCell& initial_set,
      const HybridGridTreeSet& reach_set,
      const HybridGridTreeSet& intermediate_set,
      const HybridGridTreeSet& final_set)
  : _data(new Data(initial_set))
{
  this->_data->reach=reach_set;
  this->_data->intermediate=intermediate_set;
  this->_data->final=final_set;
}

HybridGridCell const&
Orbit<HybridGridCell>::
initial() const
{
  return this->_data->initial;
}

HybridGridTreeSet const&
Orbit<HybridGridCell>::
reach() const
{
  return this->_data->reach;
}

HybridGridTreeSet const&
Orbit<HybridGridCell>::
intermediate() const
{
  return this->_data->intermediate;
}

HybridGridTreeSet const&
Orbit<HybridGridCell>::
final() const
{
  return this->_data->final;
}


} // namespace Ariadne
