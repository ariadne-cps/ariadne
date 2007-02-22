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

namespace Ariadne { namespace Geometry {


template<class S> inline
HybridSet<S>::HybridSet(location_type nq) 
  : _component_sets(nq)
{
}
  
  
template<class S> inline
HybridSet<S>::HybridSet(location_type nq, const S& s) 
  : _component_sets(nq,s) 
{
}
   
 
template<class S> template<class S1> inline
HybridSet<S>::HybridSet(const HybridSet<S1>& hs)
  : _component_sets() 
{
  for(location_type q=0; q!=hs.number_of_discrete_locations(); ++q) {
    this->_component_sets.push_back(static_cast<S>(hs[q]));
  }
}


template<class S> inline
location_type 
HybridSet<S>::number_of_discrete_locations() const 
{ 
  return _component_sets.size(); 
}


template<class S> inline
S& 
HybridSet<S>::operator[](const location_type& q)
{ 
  if (q>_component_sets.size()) {
    std::ostringstream o;
    o << "Requested the component number "<<q<<" of a hybrid set having " 
      << _component_sets.size() << " components.";
    throw std::runtime_error(o.str());
  }

  return _component_sets[q]; 
}


template<class S> inline  
const S& 
HybridSet<S>::operator[](const location_type& q) const 
{ 
  if (q>_component_sets.size()) {
    std::ostringstream o;
    o << "Requested the component number "<<q<<" of a hybrid set having " 
      << _component_sets.size() << " components.";
    throw std::runtime_error(o.str());
  }

  return _component_sets[q];
}


template<class S> inline  
void 
HybridSet<S>::clear()
{ 
  for(location_type q=0; q!=this->number_of_discrete_locations(); ++q) { 
    (*this)[q].clear();
  } 
}

template<class S> inline  
tribool 
HybridSet<S>::empty() const { 
  tribool result=true; 
  for(size_type q=0; q!=this->number_of_discrete_locations(); ++q) {
    result=result && (*this)[q].empty();
    if(!result) { return result; } }
  return result;
}
  

template<class S> template<class S1> inline  
void 
HybridSet<S>::adjoin(const HybridSet<S1>& hs) {
  if(this->number_of_discrete_locations()!=hs.number_of_discrete_locations()) {
    throw std::runtime_error("Invalid number of discrete components");
  }
  for(location_type q=0; q!=this->number_of_discrete_locations(); ++q) {
    (*this)[q].adjoin(hs[q]);
  }
}



template<class R> inline  
HybridGridMaskSet<R>::HybridGridMaskSet(location_type nq, const FiniteGrid<R>& fg) 
  : HybridSet< GridMaskSet<R> >(nq,GridMaskSet<R>(fg)) 
{
}



template<class R> inline  
HybridGridCellListSet<R>::HybridGridCellListSet(location_type nq, const Grid<R>& g) 
  : HybridSet< GridCellListSet<R> >(nq,GridCellListSet<R>(g)) 
{
}
    
  /*! \brief Construct a set for \a n discrete modes, based on a list of finite grids. */
template<class R> inline  
HybridGridCellListSet<R>::HybridGridCellListSet(const HybridSet< GridMaskSet<R> >& hgms)
  : HybridSet< GridCellListSet<R> >(hgms)
{
}



template<class R> inline
HybridGridCellListSet<R>
regular_intersection(const HybridGridCellListSet<R>& hgcl, const HybridGridMaskSet<R>& hgms) 
{
  if(hgcl.number_of_discrete_locations()!=hgms.number_of_discrete_locations()) {
    throw HybridSystemError("Intersection of sets with different numbers of discrete locations");
  }
  
  HybridGridCellListSet<R> result=hgcl;
  result.clear();
  for(location_type q=0; q!=hgcl.number_of_discrete_locations(); ++q) {
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
