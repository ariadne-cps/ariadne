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
#include "../geometry/rectangular_set.h"
#include "../geometry/polyhedral_set.h"

namespace Ariadne { 



template<class S> inline
Geometry::HybridSetBase<S>::~HybridSetBase() 
{
}
  
template<class S> inline
Geometry::HybridSetBase<S>::HybridSetBase() 
  : _component_sets()
{
}
  
  
 
template<class S> inline
Geometry::HybridSetBase<S>::HybridSetBase(const HybridSpace& locations)
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
Geometry::HybridSetBase<S>::HybridSetBase(const HybridSetBase<S>& hs)
  : _component_sets() 
{
  for(typename HybridSetBase<S>::const_iterator loc_iter=hs.begin();
      loc_iter!=hs.end(); ++loc_iter)
  {
    location_type q=loc_iter->first;
    const S& s=*loc_iter->second;
    this->new_location(q,s);
  }
}


template<class S> inline
Geometry::HybridSetBase<S>&
Geometry::HybridSetBase<S>::operator=(const HybridSetBase<S>& hs)
{
  if(this!=&hs) {
    this->_component_sets.clear();
    for(typename HybridSetBase<S>::const_iterator loc_iter=hs.begin();
        loc_iter!=hs.end(); ++loc_iter)
    {
      location_type q=loc_iter->first;
      const S& s=*loc_iter->second;
      this->new_location(q,s);
    }
  }
  return *this;
}


template<class S> template<class S1> inline
Geometry::HybridSetBase<S>::HybridSetBase(const HybridSetBase<S1>& hs)
  : _component_sets() 
{
  for(typename HybridSetBase<S1>::const_iterator loc_iter=hs.begin();
      loc_iter!=hs.end(); ++loc_iter)
  {
    location_type q=loc_iter->first;
    const S1& s=*loc_iter->second;
    this->new_location(q,s);
  }
}



template<class S> inline
S&
Geometry::HybridSetBase<S>::new_location(location_type q, dimension_type d)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"S& HybridSetBase<S>::new_location(location_type q, dimension_type d)");
  this->_component_sets.insert(std::make_pair(q,new S(d)));
  return (*this)[q];
}


template<class S> inline
S&
Geometry::HybridSetBase<S>::new_location(location_type q, const S& s)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"S& HybridSetBase<S>::new_location(location_type q, S t)");
  this->_component_sets.insert(std::make_pair(q,shared_ptr<S>(s.clone())));
  return (*this)[q];
}



template<class S> template<class T> inline
S&
Geometry::HybridSetBase<S>::new_location(location_type q, const T& t)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"SetInterface<R>& HybridSet<R>::new_location<T>(location_type q, T t)");
  typedef typename S::real_type R;
  this->_component_sets.insert(std::make_pair(q,shared_ptr<S>(new S(t))));

  return (*this)[q];
}


template<class S> template<class T1, class T2> inline
S&
Geometry::HybridSetBase<S>::new_location(location_type q, const T1& t1, const T2& t2)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"SetInterface<R>& HybridSet<R>::new_location<T>(location_type q, T1 t1, T2 t2)");
  typedef typename S::real_type R;
  this->_component_sets.insert(std::make_pair(q,shared_ptr<S>(new S(t1,t2))));

  return (*this)[q];
}



template<class S> inline
Geometry::HybridSpace
Geometry::HybridSetBase<S>::locations() const 
{ 
  HybridSpace result;
  for(const_iterator loc_iter=this->_component_sets.begin(); 
      loc_iter!=this->_component_sets.end(); ++loc_iter) 
  {
    result.new_location(loc_iter->first,loc_iter->second->dimension());
  } 
  return result;
}




template<class S> inline
Geometry::location_type 
Geometry::HybridSetBase<S>::number_of_locations() const 
{ 
  return _component_sets.size(); 
}


template<class S> inline
bool 
Geometry::HybridSetBase<S>::has_location(location_type q) const
{ 
  return this->_component_sets.find(q)!=this->_component_sets.end();
}


template<class S> inline
S& 
Geometry::HybridSetBase<S>::operator[](location_type q)
{ 
  if(!this->has_location(q)) {
    ARIADNE_THROW(InvalidLocation,"S& HybridSetBase::operator[](location_type q)","this->locations()="<<this->locations()<<", q="<<q);
  }
  return *this->_component_sets.find(q)->second;
}


template<class S> inline  
const S& 
Geometry::HybridSetBase<S>::operator[](location_type q) const 
{ 
  if(!this->has_location(q)) {
    ARIADNE_THROW(InvalidLocation,"S HybridSetBase::operator[](location_type q) const","this->locations()="<<this->locations()<<", q="<<q);
  }
  return *this->_component_sets.find(q)->second;
}


template<class S> inline  
void 
Geometry::HybridSetBase<S>::clear()
{ 
  for(iterator loc_iter=this->_component_sets.begin(); 
      loc_iter!=this->_component_sets.end(); ++loc_iter) 
  {
    loc_iter->second->clear();
  } 
}


template<class S> inline  
tribool 
Geometry::HybridSetBase<S>::empty() const { 
  tribool result=true; 
  for(const_iterator loc_iter=this->begin(); loc_iter!=this->end(); ++loc_iter) {
    result=result && loc_iter->second->empty();
    if(!result) { 
      return result; 
    }
  }
  return result;
}
  

template<class S> inline  
tribool 
Geometry::HybridSetBase<S>::bounded() const { 
  tribool result=true; 
  for(const_iterator loc_iter=this->begin(); loc_iter!=this->end(); ++loc_iter) {
    result=result && loc_iter->second->bounded();
    if(!result) { 
      return result; 
    }
  }
  return result;
}
  

template<class S> template<class T> inline  
void 
Geometry::HybridSetBase<S>::adjoin(location_type q, const T& t) {
  ARIADNE_CHECK_LOCATION(*this,q,"HybridSet<S>::adjoin<T>(location_type q, T t)");
  (*this)[q].adjoin(t);
}


template<class S> template<class S1> inline  
void 
Geometry::HybridSetBase<S>::adjoin(const HybridSetBase<S1>& hs) {
  ARIADNE_CHECK_SAME_LOCATIONS(*this,hs,"void HybridSet<S>::adjoin<SET>(HybridSetBase<SET>)");

  for(typename HybridSetBase<S1>::const_iterator loc_iter=hs.begin();
      loc_iter!=hs.end(); ++loc_iter)
  {
    location_type q=loc_iter->first;
    const S1& s=*loc_iter->second;
    (*this)[q].adjoin(s);
  }
}


template<class S> template<class S1> inline  
void 
Geometry::HybridSetBase<S>::restrict(const HybridSetBase<S1>& hs) {
  ARIADNE_CHECK_SAME_LOCATIONS(*this,hs,"void HybridSet<S>::restrict<SET>(HybridSetBase<SET>)");

  for(typename HybridSetBase<S1>::const_iterator loc_iter=hs.begin();
      loc_iter!=hs.end(); ++loc_iter)
  {
    location_type q=loc_iter->first;
    const S1& s=*loc_iter->second;
    (*this)[q].restrict(s);
  }
}


template<class S> template<class S1> inline  
void 
Geometry::HybridSetBase<S>::remove(const HybridSetBase<S1>& hs) {
  ARIADNE_CHECK_SAME_LOCATIONS(*this,hs,"void HybridSet<S>::remove<SET>(HybridSetBase<SET>)");

  for(typename HybridSetBase<S1>::const_iterator loc_iter=hs.begin();
      loc_iter!=hs.end(); ++loc_iter)
  {
    location_type q=loc_iter->first;
    const S1& s=*loc_iter->second;
    (*this)[q].remove(s);
  }
}


template<class S> inline
typename Geometry::HybridSetBase<S>::iterator 
Geometry::HybridSetBase<S>::begin()
{ 
  return this->_component_sets.begin();
}


template<class S> inline
typename Geometry::HybridSetBase<S>::const_iterator 
Geometry::HybridSetBase<S>::begin() const
{ 
  return this->_component_sets.begin();
}


template<class S> inline
typename Geometry::HybridSetBase<S>::iterator 
Geometry::HybridSetBase<S>::end()
{ 
  return this->_component_sets.end();
}


template<class S> inline
typename Geometry::HybridSetBase<S>::const_iterator 
Geometry::HybridSetBase<S>::end() const
{ 
  return this->_component_sets.end();
}




template<class S1, class S2> inline
tribool
Geometry::subset(const HybridSetBase<S1>& hs1, const HybridSetBase<S2>& hs2)
{
  ARIADNE_CHECK_SAME_LOCATIONS(hs1,hs2,"tribool subset(HybridSetBase<S1> hs1,HybridSetBase<S2> hs2)"); 
  
  tribool result=true;
  for(typename HybridSetBase<S1>::const_iterator hs1_iter=hs1.begin();
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
Geometry::operator<<(std::ostream& os, const HybridSetBase<S>& hs)
{ 
  return hs.write(os);
}






template<class R> inline
size_type
Geometry::HybridGridMaskSet<R>::capacity() const
{
  size_type result=0;
  for(typename HybridGridMaskSet<R>::const_iterator loc_iter=this->begin(); 
      loc_iter!=this->end(); ++loc_iter) 
  {
    result+=loc_iter->second->capacity();
  }
  return result;
}


template<class R> inline
size_type
Geometry::HybridGridMaskSet<R>::size() const
{
  size_type result=0;
  for(typename HybridGridMaskSet<R>::const_iterator loc_iter=this->begin(); 
      loc_iter!=this->end(); ++loc_iter) 
  {
    result+=loc_iter->second->size();
  }
  return result;
}


template<class R> inline
size_type
Geometry::HybridGridCellListSet<R>::size() const
{
  size_type result=0;
  for(typename HybridGridCellListSet<R>::const_iterator loc_iter=this->begin(); 
      loc_iter!=this->end(); ++loc_iter) 
  {
    result+=loc_iter->second->size();
  }
  return result;
}


template<class R> inline
void
Geometry::HybridGridCellListSet<R>::unique_sort() 
{
  for(typename HybridGridCellListSet<R>::iterator loc_iter=this->begin(); 
      loc_iter!=this->end(); ++loc_iter) 
  {
    loc_iter->second->unique_sort();
  }
}


template<class R> inline
Geometry::HybridGridCellListSet<R>
Geometry::difference(const HybridGridCellListSet<R>& hgcl, const HybridGridMaskSet<R>& hgms) 
{
  ARIADNE_CHECK_SAME_LOCATIONS(hgcl,hgms,"HybridGridCellListSet difference(HybridGridCellListSet hgcl, HybridGridMaskSet hgms)");
  
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
Geometry::HybridGridMaskSet<R>
Geometry::difference(const HybridGridMaskSet<R>& hgms1, const HybridGridMaskSet<R>& hgms2) 
{
  ARIADNE_CHECK_SAME_LOCATIONS(hgms1,hgms2,"HybridGridCellListSet difference(HybridGridCellListSet hgms1, HybridGridMaskSet hgms2)");
  
  HybridGridMaskSet<R> result;
  for(typename HybridGridMaskSet<R>::const_iterator loc_iter=result.begin(); 
      loc_iter!=result.end(); ++loc_iter) 
  {
    location_type q=loc_iter->first;
    result.new_location(q,difference(hgms1[q],hgms2[q]));
  }
  return result;
}


template<class R> inline
Geometry::HybridGridCellListSet<R>
Geometry::regular_intersection(const HybridGridCellListSet<R>& hgcl, const HybridGridMaskSet<R>& hgms) 
{
  ARIADNE_CHECK_SAME_LOCATIONS(hgcl,hgms,"HybridGridCellListSet regular_intersection(HybridGridCellListSet hgcl, HybridGridMaskSet hgms)");
  
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
Geometry::HybridGridCellListSet<R>
Geometry::regular_intersection(const HybridGridMaskSet<R>& hgms, const HybridGridCellListSet<R>& hgcl) 
{
  return regular_intersection(hgcl,hgms);
}


template<class R> inline
Geometry::HybridGridMaskSet<R>
Geometry::regular_intersection(const HybridGridMaskSet<R>& hgms1, const HybridGridMaskSet<R>& hgms2) 
{
  ARIADNE_CHECK_SAME_LOCATIONS(hgms1,hgms2,"HybridGridMaskSet regular_intersection(HybridGridMaskSet hgms1, HybridGridMaskSet hgms2)");

  HybridGridMaskSet<R> result;
  for(typename HybridGridMaskSet<R>::const_iterator loc_iter=result.begin(); 
      loc_iter!=result.end(); ++loc_iter) 
  {
    location_type q=loc_iter->first;
    result.new_location(q,regular_intersection(hgms1[q],hgms2[q]));
  }
  return result;
}





template<class R> inline
Geometry::HybridSet<R>::HybridSet()
  : HybridSetBase< SetInterface<R> >() 
{   
}


template<class R> inline
Geometry::HybridSet<R>::HybridSet(const HybridSet<R>& hs)
  : HybridSetBase< SetInterface<R> >(hs) 
{   
}


template<class R> template<class S> inline
Geometry::HybridSet<R>::HybridSet(const HybridSetBase<S>& hs)
  : HybridSetBase< SetInterface<R> >() 
{   
  for(typename HybridSetBase<S>::const_iterator loc_iter=hs.begin();
      loc_iter!=hs.end(); ++loc_iter)
  {
    location_type q=loc_iter->first;
    const S& s=*loc_iter->second;
    this->new_location(q,s);
  }
}


template<class R> inline
Geometry::SetInterface<R>&
Geometry::HybridSet<R>::new_location(location_type q, const Geometry::SetInterface<R>& s)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"HybridSet<R>::new_location(location_type q, SetInterface<R> s)");
  HybridSetBase< SetInterface<R> >::new_location(q,s);

  return (*this)[q];
}


template<class R> inline
Geometry::SetInterface<R>&
Geometry::HybridSet<R>::new_location(location_type q, const Geometry::Rectangle<R>& r)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"HybridSet<R>::new_location(location_type q, Rectangle<R> r)");
  RectangularSet<R> rs(r);
  HybridSetBase< SetInterface<R> >::new_location(q,static_cast<const SetInterface<R>&>(rs));

  return (*this)[q];
}


template<class R> inline
Geometry::SetInterface<R>&
Geometry::HybridSet<R>::new_location(location_type q, const Geometry::Polyhedron<R>& p)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"HybridSet<R>::new_location(location_type q, Polyhedron<R> p)");
  PolyhedralSet<R> ps(p);
  HybridSetBase< SetInterface<R> >::new_location(q,static_cast<const SetInterface<R>&>(ps));

  return (*this)[q];
}




}
