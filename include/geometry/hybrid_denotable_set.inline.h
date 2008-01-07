/***************************************************************************
 *            hybrid_denotable_set.inline.h
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

#include "base/stlio.h"
#include "base/types.h"
#include "geometry/rectangular_set.h"
#include "geometry/polyhedral_set.h"
#include "geometry/hybrid_space.h"

namespace Ariadne { 




template<class R> inline
Geometry::HybridSpace
Geometry::HybridGrid<R>::locations() const 
{ 
  HybridSpace result;
  for(locations_const_iterator loc_iter=this->locations_begin(); 
      loc_iter!=this->locations_end(); ++loc_iter) 
  {
    result.new_location(loc_iter->first,loc_iter->second.dimension());
  } 
  return result;
}




template<class S> inline
Geometry::HybridDenotableSet<S>::~HybridDenotableSet() 
{
}
  
template<class S> inline
Geometry::HybridDenotableSet<S>::HybridDenotableSet() 
  : _component_sets()
{
}
  
  
 
template<class S> inline
Geometry::HybridDenotableSet<S>::HybridDenotableSet(const HybridSpace& locations)
  : _component_sets() 
{
  for(HybridSpace::const_iterator loc_iter=locations.begin();
      loc_iter!=locations.end(); ++loc_iter)
  {
    DiscreteState q=loc_iter->id();
    dimension_type d=loc_iter->dimension();
    this->new_location(q,d);
  }
}


template<class S> inline
Geometry::HybridDenotableSet<S>::HybridDenotableSet(const HybridDenotableSet<S>& hs)
  : _component_sets() 
{
  for(typename HybridDenotableSet<S>::locations_const_iterator loc_iter=hs.locations_begin();
      loc_iter!=hs.locations_end(); ++loc_iter)
  {
    DiscreteState q=loc_iter->first;
    const S& s=loc_iter->second;
    this->new_location(q,s);
  }
}


template<class S> inline
Geometry::HybridDenotableSet<S>&
Geometry::HybridDenotableSet<S>::operator=(const HybridDenotableSet<S>& hs)
{
  if(this!=&hs) {
    this->_component_sets.clear();
    for(locations_const_iterator loc_iter=hs.locations_begin();
        loc_iter!=hs.locations_end(); ++loc_iter)
    {
      DiscreteState q=loc_iter->first;
      const S& s=loc_iter->second;
      this->new_location(q,s);
    }
  }
  return *this;
}


template<class S> template<class S1> inline
Geometry::HybridDenotableSet<S>::HybridDenotableSet(const HybridDenotableSet<S1>& hs)
  : _component_sets() 
{
  for(typename HybridDenotableSet<S1>::locations_const_iterator loc_iter=hs.locations_begin();
      loc_iter!=hs.locations_end(); ++loc_iter)
  {
    DiscreteState q=loc_iter->first;
    const S1& s=loc_iter->second;
    this->new_location(q,s);
  }
}



template<class DS> inline
DS&
Geometry::HybridDenotableSet<DS>::new_location(DiscreteState q, dimension_type d)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"DS& HybridDenotableSet<DS>::new_location(DiscreteState q, dimension_type d)");
  this->_component_sets.insert(std::make_pair(q,DS(d)));
  return (*this)[q];
}


template<class DS> inline
DS&
Geometry::HybridDenotableSet<DS>::new_location(DiscreteState q, const DS& s)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"DS& HybridDenotableSet<DS>::new_location(DiscreteState q, S t)");
  this->_component_sets.insert(std::make_pair(q,s));
  return (*this)[q];
}



template<class DS> template<class T> inline
DS&
Geometry::HybridDenotableSet<DS>::new_location(DiscreteState q, const T& t)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"DS& HybridDenotableSet<DS>::new_location<T>(DiscreteState q, T t)");
  this->_component_sets.insert(std::make_pair(q,static_cast<DS>(t)));

  return (*this)[q];
}


template<class DS> template<class T1, class T2> inline
DS&
Geometry::HybridDenotableSet<DS>::new_location(DiscreteState q, const T1& t1, const T2& t2)
{
  ARIADNE_CHECK_NEW_LOCATION(*this,q,"DS& HybridDenotableSet<DS>::new_location<T>(DiscreteState q, T1 t1, T2 t2)");
  this->_component_sets.insert(std::make_pair(q,static_cast<DS>(t1,t2)));

  return (*this)[q];
}



template<class DS> inline
Geometry::HybridSpace
Geometry::HybridDenotableSet<DS>::locations() const 
{ 
  HybridSpace result;
  for(locations_const_iterator loc_iter=this->locations_begin(); 
      loc_iter!=this->locations_end(); ++loc_iter) 
  {
    result.new_location(loc_iter->first,loc_iter->second.dimension());
  } 
  return result;
}


template<class DS> inline
Geometry::HybridSpace
Geometry::HybridDenotableSet<DS>::space() const 
{ 
  return this->locations();
}




template<class DS> inline
size_type 
Geometry::HybridDenotableSet<DS>::number_of_locations() const 
{ 
  return _component_sets.size(); 
}


template<class DS> inline
bool 
Geometry::HybridDenotableSet<DS>::has_location(DiscreteState q) const
{ 
  return this->_component_sets.find(q)!=this->_component_sets.end();
}


template<class DS> inline
DS& 
Geometry::HybridDenotableSet<DS>::operator[](DiscreteState q)
{ 
  if(!this->has_location(q)) {
    ARIADNE_THROW(InvalidLocation,"DS& HybridDenotableSet::operator[](DiscreteState q)","this->locations()="<<this->locations()<<", q="<<q);
  }
  return this->_component_sets.find(q)->second;
}


template<class DS> inline  
const DS& 
Geometry::HybridDenotableSet<DS>::operator[](DiscreteState q) const 
{ 
  if(!this->has_location(q)) {
    ARIADNE_THROW(InvalidLocation,"DS HybridDenotableSet::operator[](DiscreteState q) const","this->locations()="<<this->locations()<<", q="<<q);
  }
  return this->_component_sets.find(q)->second;
}


template<class DS> inline
size_type
Geometry::HybridDenotableSet<DS>::size() const
{
  size_type result=0;
  for(typename HybridDenotableSet<DS>::locations_const_iterator loc_iter=this->locations_begin(); 
      loc_iter!=this->locations_end(); ++loc_iter)
  {
    result+=loc_iter->second.size();
  }
  return result;
}


template<class DS> inline  
void 
Geometry::HybridDenotableSet<DS>::clear()
{ 
  for(locations_iterator loc_iter=this->locations_begin(); 
      loc_iter!=this->locations_end(); ++loc_iter) 
  {
    loc_iter->second.clear();
  } 
}


template<class DS> inline  
tribool 
Geometry::HybridDenotableSet<DS>::empty() const 
{ 
  tribool result=true; 
  for(locations_const_iterator loc_iter=this->locations_begin(); loc_iter!=this->locations_end(); ++loc_iter) {
    result=result && loc_iter->second.empty();
    if(!result) { 
      return result; 
    }
  }
  return result;
}
  

template<class DS> inline  
tribool 
Geometry::HybridDenotableSet<DS>::bounded() const 
{ 
  tribool result=true; 
  for(locations_const_iterator loc_iter=this->locations_begin(); loc_iter!=this->locations_end(); ++loc_iter) {
    result=result && loc_iter->second.bounded();
    if(!result) { 
      return result; 
    }
  }
  return result;
}
  

template<class DS> template<class T> inline  
void 
Geometry::HybridDenotableSet<DS>::adjoin(DiscreteState q, const T& t) 
{
  ARIADNE_CHECK_LOCATION(*this,q,"HybridDenotableSet<DS>::adjoin<T>(DiscreteState q, T t)");
  (*this)[q].adjoin(t);
}

template<class DS> template<class T> inline  
void 
Geometry::HybridDenotableSet<DS>::adjoin(const std::pair<DiscreteState,T>& pair)
{
  ARIADNE_CHECK_LOCATION(*this,pair.first,"HybridDenotableSet<DS>::adjoin<T>(pair<DS,T>)")
  (*this)[pair.first].adjoin(pair.second);
}



template<class DS>inline  
typename Geometry:: HybridDenotableSet<DS>::const_iterator
Geometry::HybridDenotableSet<DS>::begin() const
{
  return const_iterator(this->_component_sets,false);
}

template<class DS>inline  
typename Geometry:: HybridDenotableSet<DS>::const_iterator
Geometry::HybridDenotableSet<DS>::end() const
{
  return const_iterator(this->_component_sets,true);
}


template<class DS> template<class T> inline  
void 
Geometry::HybridDenotableSet<DS>::adjoin(const HybridBasicSet<T>& hbs) 
{
  ARIADNE_CHECK_LOCATION(*this,hbs.state(),"HybridDenotableSet<S>::adjoin<T>(DiscreteState q, T t)");
  (*this)[hbs.state()].adjoin(hbs.set());
}


template<class DS> template<class DSA> inline  
void 
Geometry::HybridDenotableSet<DS>::adjoin(const HybridDenotableSet<DSA>& hs) 
{
  ARIADNE_CHECK_SAME_LOCATIONS(*this,hs,"void HybridDenotableSet<DS>::adjoin(HybridDenotableSet<DSA>)");

  for(typename HybridDenotableSet<DSA>::locations_const_iterator loc_iter=hs.locations_begin();
      loc_iter!=hs.locations_end(); ++loc_iter)
  {
    DiscreteState q=loc_iter->first;
    const DSA& s=loc_iter->second;
    (*this)[q].adjoin(s);
  }
}


template<class DS> template<class DSA> inline  
void 
Geometry::HybridDenotableSet<DS>::restrict(const HybridDenotableSet<DSA>& hs) 
{
  ARIADNE_CHECK_SAME_LOCATIONS(*this,hs,"void HybridDenotableSet<DS>::restrict(HybridDenotableSet<DSA>)");

  for(typename HybridDenotableSet<DSA>::locations_const_iterator loc_iter=hs.locations_begin();
      loc_iter!=hs.locations_end(); ++loc_iter)
  {
    DiscreteState q=loc_iter->first;
    const DSA s=loc_iter->second;
    (*this)[q].restrict(s);
  }
}


template<class DS> template<class DSA> inline  
void 
Geometry::HybridDenotableSet<DS>::remove(const HybridDenotableSet<DSA>& hs) 
{
  ARIADNE_CHECK_SAME_LOCATIONS(*this,hs,"void HybridDenotableSet<DS>::remove(HybridDenotableSet<DSA>)");

  for(typename HybridDenotableSet<DSA>::locations_const_iterator loc_iter=hs.locations_begin();
      loc_iter!=hs.locations_end(); ++loc_iter)
  {
    DiscreteState q=loc_iter->first;
    const DSA& s=loc_iter->second;
    (*this)[q].remove(s);
  }
}


template<class DS> inline
typename Geometry::HybridDenotableSet<DS>::locations_iterator 
Geometry::HybridDenotableSet<DS>::locations_begin()
{ 
  return this->_component_sets.begin();
}


template<class DS> inline
typename Geometry::HybridDenotableSet<DS>::locations_const_iterator 
Geometry::HybridDenotableSet<DS>::locations_begin() const
{ 
  return this->_component_sets.begin();
}


template<class DS> inline
typename Geometry::HybridDenotableSet<DS>::locations_iterator 
Geometry::HybridDenotableSet<DS>::locations_end()
{ 
  return this->_component_sets.end();
}


template<class DS> inline
typename Geometry::HybridDenotableSet<DS>::locations_const_iterator 
Geometry::HybridDenotableSet<DS>::locations_end() const
{ 
  return this->_component_sets.end();
}




template<class DS1, class DS2> inline
tribool
Geometry::subset(const HybridDenotableSet<DS1>& hs1, const HybridDenotableSet<DS2>& hs2)
{
  ARIADNE_CHECK_SAME_LOCATIONS(hs1,hs2,"tribool subset(HybridDenotableSet<DS1> hs1,HybridDenotableSet<DS2> hs2)"); 
  
  tribool result=true;
  for(typename HybridDenotableSet<DS1>::locations_const_iterator hs1_iter=hs1.locations_begin();
      hs1_iter!=hs1.locations_end(); ++hs1_iter)
  {
    DiscreteState q=hs1_iter->first;
    result=result && subset(hs1[q],hs2[q]);
    if(!result) {
      break;
    }
  }
  return result;
}


template<class DS> inline
std::ostream& 
Geometry::HybridDenotableSet<DS>::write(std::ostream& os) const
{ 
  os << "HybridDenotableSet( { \n";
  for(locations_const_iterator iter=this->locations_begin(); iter!=this->locations_end(); ++iter)
  {
    DiscreteState loc=iter->first;
    const DS& set=iter->second;
    os << "  "<<loc<<": " << set << ",\n";
  }
  os << "} )";
  return os;
}


template<class R> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const HybridGrid<R>& hs)
{ 
  return os << "HybridGrid(...)";
}

template<class DS> inline
std::ostream& 
Geometry::operator<<(std::ostream& os, const HybridDenotableSet<DS>& hs)
{ 
  return hs.write(os);
}







template<class R> inline
Geometry::HybridGridCellListSet<R>::HybridGridCellListSet(const HybridGrid<R>& hg) 
{
  for(typename HybridGrid<R>::locations_const_iterator iter=hg.locations_begin(); 
      iter!=hg.locations_end(); ++iter) 
  {
    this->new_location(iter->first,iter->second);
  }
}
  


template<class R> inline
Geometry::HybridGridCellListSet<R>
Geometry::difference(const HybridGridCellListSet<R>& hgcl, const HybridGridMaskSet<R>& hgms) 
{
  ARIADNE_CHECK_SAME_LOCATIONS(hgcl,hgms,"HybridGridCellListSet difference(HybridGridCellListSet hgcl, HybridGridMaskSet hgms)");
  
  HybridGridCellListSet<R> result=hgcl;
  result.clear();
  for(typename HybridGridCellListSet<R>::locations_const_iterator loc_iter=result.locations_begin(); 
      loc_iter!=result.locations_end(); ++loc_iter) 
  {
    DiscreteState q=loc_iter->first;
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
  for(typename HybridGridMaskSet<R>::locations_const_iterator loc_iter=result.begin(); 
      loc_iter!=result.end(); ++loc_iter) 
  {
    DiscreteState q=loc_iter->first;
    result.new_location(q,difference(hgms1[q],hgms2[q]));
  }
  return result;
}


template<class BS> inline
Geometry::HybridGridCellListSet<typename BS::real_type>
Geometry::outer_approximation(const HybridListSet<BS>& hls, const HybridGrid<typename BS::real_type>& hg) 
{
  typedef typename BS::real_type R;
  ARIADNE_CHECK_SAME_LOCATIONS(hls,hg,"HybridGridCellListSet outer_approximation(HybridListSet hls, HybridGrid hg)");

  HybridGridCellListSet<R> result(hg);
  result.clear();
  for(typename HybridListSet<BS>::locations_const_iterator loc_iter=hls.locations_begin(); 
      loc_iter!=hls.locations_end(); ++loc_iter) 
  {
    DiscreteState q=loc_iter->first;
    result[q]=outer_approximation(hls[q],hg[q]);
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
  for(typename HybridGridCellListSet<R>::locations_const_iterator loc_iter=result.begin(); 
      loc_iter!=result.end(); ++loc_iter) 
  {
    DiscreteState q=loc_iter->first;
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
  for(typename HybridGridMaskSet<R>::locations_const_iterator loc_iter=result.begin(); 
      loc_iter!=result.end(); ++loc_iter) 
  {
    DiscreteState q=loc_iter->first;
    result.new_location(q,regular_intersection(hgms1[q],hgms2[q]));
  }
  return result;
}




}
