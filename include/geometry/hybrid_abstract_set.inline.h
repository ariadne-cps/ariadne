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

#include "../base/stlio.h"
#include "../geometry/rectangular_set.h"
#include "../geometry/polyhedral_set.h"

namespace Ariadne { 




template<class S> inline
Geometry::HybridAbstractSet<S>::~HybridAbstractSet() 
{
}
  
template<class S> inline
Geometry::HybridAbstractSet<S>::HybridAbstractSet() 
  : _component_sets()
{
}
  
  
 
template<class S> inline
Geometry::HybridAbstractSet<S>::HybridAbstractSet(const HybridSpace& locations)
  : _component_sets() 
{
  for(HybridSpace::const_iterator loc_iter=locations.begin();
      loc_iter!=locations.end(); ++loc_iter)
  {
    location_type q=loc_iter->id();
    dimension_type d=loc_iter->dimension();
    this->new_location(q,d);
  }
}


template<class S> inline
Geometry::HybridAbstractSet<S>::HybridAbstractSet(const HybridAbstractSet<S>& hs)
  : _component_sets() 
{
  for(typename HybridAbstractSet<S>::locations_const_iterator loc_iter=hs.locations_begin();
      loc_iter!=hs.locations_end(); ++loc_iter)
  {
    location_type q=loc_iter->first;
    const S& s=*loc_iter->second;
    this->new_location(q,s);
  }
}


template<class S> inline
Geometry::HybridAbstractSet<S>&
Geometry::HybridAbstractSet<S>::operator=(const HybridAbstractSet<S>& hs)
{
  if(this!=&hs) {
    this->_component_sets.clear();
    for(typename HybridAbstractSet<S>::locations_const_iterator loc_iter=hs.locations_begin();
        loc_iter!=hs.locations_end(); ++loc_iter)
    {
      location_type q=loc_iter->first;
      const S& s=*loc_iter->second;
      this->new_location(q,s);
    }
  }
  return *this;
}


template<class S> template<class DS> inline
Geometry::HybridAbstractSet<S>::HybridAbstractSet(const HybridDenotableSet<DS>& hds)
  : _component_sets() 
{
  for(typename HybridDenotableSet<DS>::locations_const_iterator loc_iter=hds.locations_begin();
      loc_iter!=hds.locations_end(); ++loc_iter)
  {
    location_type q=loc_iter->first;
    const DS& ds=loc_iter->second;
    this->new_location(q,ds);
  }
}

template<class S> template<class S1> inline
Geometry::HybridAbstractSet<S>::HybridAbstractSet(const HybridAbstractSet<S1>& hs)
  : _component_sets() 
{
  for(typename HybridAbstractSet<S1>::locations_const_iterator loc_iter=hs.locations_begin();
      loc_iter!=hs.locations_end(); ++loc_iter)
  {
    location_type q=loc_iter->first;
    const S1& s=*loc_iter->second;
    this->new_location(q,s);
  }
}




template<class S> inline
S&
Geometry::HybridAbstractSet<S>::new_location(location_type q, const S& s)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"S& HybridAbstractSet<S>::new_location(location_type q, S t)");
  this->_component_sets.insert(std::make_pair(q,shared_ptr<S>(s.clone())));
  return (*this)[q];
}

template<class S> template<class SS> inline
S&
Geometry::HybridAbstractSet<S>::new_location(location_type q, const SS& s)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"S& HybridAbstractSet<S>::new_location(location_type q, SS s)");
  this->_component_sets.insert(std::make_pair(q,shared_ptr<S>(s.clone())));
  return (*this)[q];
}







template<class S> inline
Geometry::HybridSpace
Geometry::HybridAbstractSet<S>::locations() const 
{ 
  HybridSpace result;
  for(locations_const_iterator loc_iter=this->locations_begin(); 
      loc_iter!=this->locations_end(); ++loc_iter) 
  {
    result.new_location(loc_iter->first,loc_iter->second->dimension());
  } 
  return result;
}




template<class S> inline
Geometry::location_type 
Geometry::HybridAbstractSet<S>::number_of_locations() const 
{ 
  return _component_sets.size(); 
}


template<class S> inline
bool 
Geometry::HybridAbstractSet<S>::has_location(location_type q) const
{ 
  return this->_component_sets.find(q)!=this->_component_sets.end();
}


template<class S> inline
S& 
Geometry::HybridAbstractSet<S>::operator[](location_type q)
{ 
  if(!this->has_location(q)) {
    ARIADNE_THROW(InvalidLocation,"S& HybridAbstractSet::operator[](location_type q)","this->locations()="<<this->locations()<<", q="<<q);
  }
  return *this->_component_sets.find(q)->second;
}


template<class S> inline  
const S& 
Geometry::HybridAbstractSet<S>::operator[](location_type q) const 
{ 
  if(!this->has_location(q)) {
    ARIADNE_THROW(InvalidLocation,"S HybridAbstractSet::operator[](location_type q) const","this->locations()="<<this->locations()<<", q="<<q);
  }
  return *this->_component_sets.find(q)->second;
}


 

template<class S> inline  
tribool 
Geometry::HybridAbstractSet<S>::bounded() const { 
  tribool result=true; 
  for(locations_const_iterator loc_iter=this->locations_begin(); loc_iter!=this->locations_end(); ++loc_iter) {
    result=result && loc_iter->second->bounded();
    if(!result) { 
      return result; 
    }
  }
  return result;
}
  



template<class S> inline
typename Geometry::HybridAbstractSet<S>::locations_iterator 
Geometry::HybridAbstractSet<S>::locations_begin()
{ 
  return this->_component_sets.begin();
}


template<class S> inline
typename Geometry::HybridAbstractSet<S>::locations_const_iterator 
Geometry::HybridAbstractSet<S>::locations_begin() const
{ 
  return this->_component_sets.begin();
}


template<class S> inline
typename Geometry::HybridAbstractSet<S>::locations_iterator 
Geometry::HybridAbstractSet<S>::locations_end()
{ 
  return this->_component_sets.end();
}


template<class S> inline
typename Geometry::HybridAbstractSet<S>::locations_const_iterator 
Geometry::HybridAbstractSet<S>::locations_end() const
{ 
  return this->_component_sets.end();
}




template<class S1, class S2> inline
tribool
Geometry::subset(const HybridAbstractSet<S1>& hs1, const HybridAbstractSet<S2>& hs2)
{
  ARIADNE_CHECK_SAME_LOCATIONS(hs1,hs2,"tribool subset(HybridAbstractSet<S1> hs1,HybridAbstractSet<S2> hs2)"); 
  
  tribool result=true;
  for(typename HybridAbstractSet<S1>::locations_const_iterator hs1_iter=hs1.begin();
      hs1_iter!=hs1.end(); ++hs1_iter)
  {
    location_type q=hs1_iter->first;
    result=result && subset(hs1[q],hs2[q]);
    if(!result) {
      break;
    }
  }
  return result;
}


template<class S> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const HybridAbstractSet<S>& hs)
{ 
  return hs.write(os);
}










template<class S> inline
S&
Geometry::HybridAbstractSet<S>::new_location(location_type q, const Geometry::Rectangle<real_type>& r)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"HybridAbstractSet<S>::new_location(location_type q, Rectangle<R> r)");
  RectangularSet<real_type> rs(r);
  HybridAbstractSet<S>::new_location(q,static_cast<const S&>(rs));

  return (*this)[q];
}


template<class S> inline
S&
Geometry::HybridAbstractSet<S>::new_location(location_type q, const Geometry::Polyhedron<real_type>& p)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"HybridAbstractSet<S>::new_location(location_type q, Polyhedron<R> p)");
  PolyhedralSet<real_type> ps(p);
  HybridAbstractSet<S>::new_location(q,static_cast<const S&>(ps));

  return (*this)[q];
}




}
