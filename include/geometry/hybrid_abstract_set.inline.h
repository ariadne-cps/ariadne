/***************************************************************************
 *            hybrid_abstract_set.inline.h
 *
 *  Copyright  2006-7  Alberto Casagrande,  Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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

#include "base/stlio.h"
#include "geometry/empty_set.h"
#include "geometry/rectangular_set.h"
#include "geometry/polyhedral_set.h"

namespace Ariadne { 




template<class R> inline
HybridSet<R>::~HybridSet() 
{
}
  
template<class R> inline
HybridSet<R>::HybridSet() 
  : _component_sets()
{
}
  
  
 
template<class R> inline
HybridSet<R>::HybridSet(const HybridSpace& locations)
  : _component_sets() 
{
  for(HybridSpace::const_iterator loc_iter=locations.begin();
      loc_iter!=locations.end(); ++loc_iter)
  {
    discrete_state_type q=loc_iter->id();
    dimension_type d=loc_iter->dimension();
    this->new_location(q,d);
  }
}


template<class R> inline
HybridSet<R>::HybridSet(const HybridSet<R>& hs)
  : _component_sets() 
{
  for(typename HybridSet<R>::locations_const_iterator loc_iter=hs.locations_begin();
      loc_iter!=hs.locations_end(); ++loc_iter)
  {
    discrete_state_type q=loc_iter->first;
    const S& s=*loc_iter->second;
    this->new_location(q,s);
  }
}


template<class R> inline
HybridSet<R>&
HybridSet<R>::operator=(const HybridSet<R>& hs)
{
  if(this!=&hs) {
    this->_component_sets.clear();
    for(typename HybridSet<R>::locations_const_iterator loc_iter=hs.locations_begin();
        loc_iter!=hs.locations_end(); ++loc_iter)
    {
      discrete_state_type q=loc_iter->first;
      const S& s=*loc_iter->second;
      this->new_location(q,s);
    }
  }
  return *this;
}


template<class R> template<class DS> inline
HybridSet<R>::HybridSet(const HybridDenotableSet<DS>& hds)
  : _component_sets() 
{
  for(typename HybridDenotableSet<DS>::locations_const_iterator loc_iter=hds.locations_begin();
      loc_iter!=hds.locations_end(); ++loc_iter)
  {
    discrete_state_type q=loc_iter->first;
    const DS& ds=loc_iter->second;
    this->new_location(q,ds);
  }
}




template<class R> inline
typename HybridSet<R>::S&
HybridSet<R>::new_location(discrete_state_type q, const S& s)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"S& HybridSet<R>::new_location(discrete_state_type q, S t)");
  this->_component_sets.insert(std::make_pair(q,shared_ptr<S>(s.clone())));
  return (*this)[q];
}

template<class R> inline
typename HybridSet<R>::S&
HybridSet<R>::new_location(discrete_state_type q, dimension_type d)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"S& HybridSet<R>::new_location(discrete_state_type q, dimension_type d)");
  this->_component_sets.insert(std::make_pair(q,shared_ptr<S>(new EmptySet<R>(d))));
  return (*this)[q];
}

template<class R> template<class T> inline
typename HybridSet<R>::S&
HybridSet<R>::new_location(discrete_state_type q, const T& s)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"S& HybridSet<R>::new_location(discrete_state_type q, SS s)");
  this->_component_sets.insert(std::make_pair(q,shared_ptr<S>(s.clone())));
  return (*this)[q];
}







template<class R> inline
HybridSpace
HybridSet<R>::locations() const 
{ 
  HybridSpace result;
  for(locations_const_iterator loc_iter=this->locations_begin(); 
      loc_iter!=this->locations_end(); ++loc_iter) 
  {
    result.new_location(loc_iter->first,loc_iter->second->space().dimension());
  } 
  return result;
}




template<class R> inline
size_type
HybridSet<R>::number_of_locations() const 
{ 
  return _component_sets.size(); 
}


template<class R> inline
bool 
HybridSet<R>::has_location(discrete_state_type q) const
{ 
  return this->_component_sets.find(q)!=this->_component_sets.end();
}


template<class R> inline
typename HybridSet<R>::S& 
HybridSet<R>::operator[](discrete_state_type q)
{ 
  if(!this->has_location(q)) {
    ARIADNE_THROW(InvalidLocation,"S& HybridSet::operator[](discrete_state_type q)","this->locations()="<<this->locations()<<", q="<<q);
  }
  return *this->_component_sets.find(q)->second;
}


template<class R> inline  
const typename HybridSet<R>::S& 
HybridSet<R>::operator[](discrete_state_type q) const 
{ 
  if(!this->has_location(q)) {
    ARIADNE_THROW(InvalidLocation,"S HybridSet::operator[](discrete_state_type q) const","this->locations()="<<this->locations()<<", q="<<q);
  }
  return *this->_component_sets.find(q)->second;
}


 




template<class R> inline
typename HybridSet<R>::locations_iterator 
HybridSet<R>::locations_begin()
{ 
  return this->_component_sets.begin();
}


template<class R> inline
typename HybridSet<R>::locations_const_iterator 
HybridSet<R>::locations_begin() const
{ 
  return this->_component_sets.begin();
}


template<class R> inline
typename HybridSet<R>::locations_iterator 
HybridSet<R>::locations_end()
{ 
  return this->_component_sets.end();
}


template<class R> inline
typename HybridSet<R>::locations_const_iterator 
HybridSet<R>::locations_end() const
{ 
  return this->_component_sets.end();
}




template<class R1, class R2> inline
tribool
subset(const HybridSet<R1>& hs1, const HybridSet<R2>& hs2)
{
  ARIADNE_CHECK_SAME_LOCATIONS(hs1,hs2,"tribool subset(HybridSet<S1> hs1,HybridSet<S2> hs2)"); 
  
  tribool result=true;
  for(typename HybridSet<R1>::locations_const_iterator hs1_iter=hs1.begin();
      hs1_iter!=hs1.end(); ++hs1_iter)
  {
    DiscreteState q=hs1_iter->first;
    result=result && subset(hs1[q],hs2[q]);
    if(!result) {
      break;
    }
  }
  return result;
}


template<class R> inline
std::ostream& 
operator<<(std::ostream& os, const HybridSet<R>& hs)
{ 
  return hs.write(os);
}










template<class R> inline
typename HybridSet<R>::S&
HybridSet<R>::new_location(discrete_state_type q, const Rectangle<real_type>& r)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"HybridSet<R>::new_location(discrete_state_type q, Box<R> r)");
  RectangularSet<real_type> rs(r);
  HybridSet<R>::new_location(q,static_cast<const S&>(rs));

  return (*this)[q];
}


template<class R> inline
typename HybridSet<R>::S&
HybridSet<R>::new_location(discrete_state_type q, const Polyhedron<real_type>& p)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"HybridSet<R>::new_location(discrete_state_type q, Polyhedron<R> p)");
  PolyhedralSet<real_type> ps(p);
  HybridSet<R>::new_location(q,static_cast<const S&>(ps));

  return (*this)[q];
}




}
