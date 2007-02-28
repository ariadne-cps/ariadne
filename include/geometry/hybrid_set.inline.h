/***************************************************************************
 *            hybrid_set.inline.h
 *
 *  Copyright  2006  Alberto Casagrande,  Pieter Collins
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

namespace Ariadne { namespace Geometry {


template<class S> inline
HybridSet<S>::HybridSet() 
  : _component_sets()
{
}
  
  
 
template<class S> inline
HybridSet<S>::HybridSet(const std::map<location_type,dimension_type>& locations)
  : _component_sets() 
{
  for(typename std::map<location_type,dimension_type>::const_iterator loc_iter=locations.begin();
      loc_iter!=locations.end(); ++loc_iter)
  {
    location_type q=loc_iter->first;
    dimension_type d=loc_iter->second;
    this->_component_sets.insert(std::make_pair(q,S(d)));
  }
}


template<class S> inline
HybridSet<S>::HybridSet(const HybridSet<S>& hs)
  : _component_sets(hs._component_sets) 
{
}


template<class S> template<class S1> inline
HybridSet<S>::HybridSet(const HybridSet<S1>& hs)
  : _component_sets() 
{
  for(typename HybridSet<S1>::const_iterator loc_iter=hs.begin();
      loc_iter!=hs.end(); ++loc_iter)
  {
    location_type q=loc_iter->first;
    const S1& s=loc_iter->second;
    this->_component_sets.insert(std::make_pair(q,s));
  }
}



template<class S> inline
S&
HybridSet<S>::new_location(location_type q, dimension_type d)
{
  if(this->has_location(q)) {
    std::ostringstream msg;
    msg << "The hybrid set already has location " << q << ".";
    throw HybridSystemError(msg.str());
  }
  this->_component_sets.insert(std::make_pair(q,S(d)));
  return (*this)[q];
}


template<class S> inline
S&
HybridSet<S>::new_location(location_type q, const S& s)
{
  if(this->has_location(q)) {
    std::ostringstream msg;
    msg << "The hybrid set already has location " << q << ".";
    throw HybridSystemError(msg.str());
  }
  this->_component_sets.insert(std::make_pair(q,s));
  return (*this)[q];
}


template<class S> template<class T> inline
S&
HybridSet<S>::new_location(location_type q, const T& t)
{
  if(this->has_location(q)) {
    std::ostringstream msg;
    msg << "The hybrid set already has location " << q << ".";
    throw HybridSystemError(msg.str());
  }
  this->_component_sets.insert(std::make_pair(q,S(t)));
  return (*this)[q];
}


template<class S> inline
std::map<location_type,dimension_type>
HybridSet<S>::locations() const 
{ 
  std::map<location_type,dimension_type> result;
  for(const_iterator loc_iter=this->_component_sets.begin(); 
      loc_iter!=this->_component_sets.end(); ++loc_iter) 
  {
    result.insert(std::make_pair(loc_iter->first,loc_iter->second.dimension()));
  } 
  return result;
}



template<class S> inline
location_type 
HybridSet<S>::number_of_locations() const 
{ 
  return _component_sets.size(); 
}


template<class S> inline
bool 
HybridSet<S>::has_location(location_type q) const
{ 
  return this->_component_sets.find(q)!=this->_component_sets.end();
}


template<class S> inline
S& 
HybridSet<S>::operator[](location_type q)
{ 
  this->check_location(q,"HybridSet<S>::operator[](location_type q)");
  return this->_component_sets.find(q)->second;
}


template<class S> inline  
const S& 
HybridSet<S>::operator[](location_type q) const 
{ 
  this->check_location(q,"HybridSet<S>::operator[](location_type q) const");
  return this->_component_sets.find(q)->second;
}


template<class S> inline  
void 
HybridSet<S>::clear()
{ 
  for(iterator loc_iter=this->_component_sets.begin(); 
      loc_iter!=this->_component_sets.end(); ++loc_iter) 
  {
    loc_iter->second.clear();
  } 
}


template<class S> inline  
tribool 
HybridSet<S>::empty() const { 
  tribool result=true; 
  for(const_iterator loc_iter=this->begin(); loc_iter!=this->end(); ++loc_iter) {
    result=result && loc_iter->second.empty();
    if(!result) { 
      return result; 
    }
  }
  return result;
}
  

template<class S> template<class S1> inline  
void 
HybridSet<S>::adjoin(location_type q, const S1& s) {
  check_location(q,"HybridSet<S>::adjoin(location_type q, const S1& s)");
  (*this)[q].adjoin(s);
}


template<class S> template<class S1> inline  
void 
HybridSet<S>::adjoin(const HybridSet<S1>& hs) {
  for(typename HybridSet<S1>::const_iterator loc_iter=hs.begin();
      loc_iter!=hs.end(); ++loc_iter)
  {
    location_type q=loc_iter->first;
    if(!this->has_location(q)) {
      throw std::runtime_error("Cannot adjoin a hybrid set to another with different discrete locations");
    }
  }

  for(typename HybridSet<S1>::const_iterator loc_iter=hs.begin();
      loc_iter!=hs.end(); ++loc_iter)
  {
    location_type q=loc_iter->first;
    const S1& s=loc_iter->second;
    (*this)[q].adjoin(s);
  }
}


template<class S> inline
typename HybridSet<S>::iterator 
HybridSet<S>::begin()
{ 
  return this->_component_sets.begin();
}


template<class S> inline
typename HybridSet<S>::const_iterator 
HybridSet<S>::begin() const
{ 
  return this->_component_sets.begin();
}


template<class S> inline
typename HybridSet<S>::iterator 
HybridSet<S>::end()
{ 
  return this->_component_sets.end();
}


template<class S> inline
typename HybridSet<S>::const_iterator 
HybridSet<S>::end() const
{ 
  return this->_component_sets.end();
}



template<class S> inline
void 
HybridSet<S>::check_location(location_type q, const char* where) const
{ 
  if (!this->has_location(q)) {
    std::ostringstream o;
    o << where << ": The hybrid set with locations " << this->locations() << " does not have a location with id " << q << ".";
    throw HybridSystemError(o.str());
  }
}

template<class S1, class S2> inline
tribool
subset(const HybridSet<S1>& hs1, const HybridSet<S2>& hs2)
{
  if(hs1.locations()!=hs2.locations()) {
    throw HybridSystemError("Comparing sets with different discrete locations");
  }
  
  tribool result=true;
  for(typename HybridSet<S1>::const_iterator hs1_iter=hs1.begin();
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

template<class S> 
std::ostream& 
HybridSet<S>::write(std::ostream& os) const
{ 
  return os << _component_sets;
}


template<class S> 
std::ostream& 
operator<<(std::ostream& os, const HybridSet<S>& hs)
{ 
  return hs.write(os);
}





template<class R> inline
size_type
HybridGridMaskSet<R>::capacity() const
{
  size_type result=0;
  for(typename HybridGridMaskSet<R>::const_iterator loc_iter=this->begin(); 
      loc_iter!=this->end(); ++loc_iter) 
  {
    result+=loc_iter->second.capacity();
  }
  return result;
}


template<class R> inline
size_type
HybridGridMaskSet<R>::size() const
{
  size_type result=0;
  for(typename HybridGridMaskSet<R>::const_iterator loc_iter=this->begin(); 
      loc_iter!=this->end(); ++loc_iter) 
  {
    result+=loc_iter->second.size();
  }
  return result;
}


template<class R> inline
size_type
HybridGridCellListSet<R>::size() const
{
  size_type result=0;
  for(typename HybridGridCellListSet<R>::const_iterator loc_iter=this->begin(); 
      loc_iter!=this->end(); ++loc_iter) 
  {
    result+=loc_iter->second.size();
  }
  return result;
}


template<class R> inline
void
HybridGridCellListSet<R>::unique_sort() 
{
  for(typename HybridGridCellListSet<R>::iterator loc_iter=this->begin(); 
      loc_iter!=this->end(); ++loc_iter) 
  {
    loc_iter->second.unique_sort();
  }
}


template<class R> inline
HybridGridCellListSet<R>
difference(const HybridGridCellListSet<R>& hgcl, const HybridGridMaskSet<R>& hgms) 
{
  if(hgcl.locations()!=hgms.locations()) {
    throw HybridSystemError("Intersection of sets with different discrete locations");
  }
  
  HybridGridCellListSet<R> result=hgcl;
  result.clear();
  for(typename HybridGridCellListSet<R>::const_iterator loc_iter=result.begin(); 
      loc_iter!=result.end(); ++loc_iter) 
  {
    location_type q=loc_iter->first;
    result[q]=difference(hgcl[q],hgms[q]);
  }
  return result;
}


template<class R> inline
HybridGridCellListSet<R>
regular_intersection(const HybridGridCellListSet<R>& hgcl, const HybridGridMaskSet<R>& hgms) 
{
  if(hgcl.locations()!=hgms.locations()) {
    throw HybridSystemError("Intersection of sets with different discrete locations");
  }
  
  HybridGridCellListSet<R> result=hgcl;
  result.clear();
  for(typename HybridGridCellListSet<R>::const_iterator loc_iter=result.begin(); 
      loc_iter!=result.end(); ++loc_iter) 
  {
    location_type q=loc_iter->first;
    result[q]=regular_intersection(hgcl[q],hgms[q]);
  }
  return result;
}


template<class R> inline
HybridGridCellListSet<R>
regular_intersection(const HybridGridMaskSet<R>& hgms, const HybridGridCellListSet<R>& hgcl) 
{
  return regular_intersection(hgcl,hgms);
}



}}
