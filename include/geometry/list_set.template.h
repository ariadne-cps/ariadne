/***************************************************************************
 *            list_set.template.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#include <iostream>

#include "list_set.h"
#include "base/stlio.h"
#include "geometry/rectangle.h"
#include "geometry/name.h"

namespace Ariadne {

template<class BS>
void
ListSet<BS>::_instantiate()
{
}


template<class BS, class R> 
tribool 
contains(const ListSet<BS>& ls, const Point<R>& p) 
{
  tribool result=false;
  for (typename ListSet<BS>::const_iterator i=ls.begin(); i!=ls.end(); ++i) {
    //result=result || i->contains(p);
    result=result || contains(*i,p);
    if(result) { return result; }
  }
  return result;
}

template<class BS> 
Box<typename BS::real_type> 
bounding_box(const ListSet<BS>& ls)  
{
  typedef typename ListSet<BS>::real_type R;
  if(ls.size()==0) { return Box<R>(ls.dimension()); }
  //Box<R> result=(*this)[0].bounding_box();
  Box<R> result=bounding_box(ls[0]);
  for(typename ListSet<BS>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
    //Box<R> bb=iter->bounding_box();
    Box<R> bb=bounding_box(*iter);
    for(size_type i=0; i!=result.dimension(); ++i) {
      if(bb.lower_bound(i) < result.lower_bound(i)) {
        result.set_lower_bound(i,bb.lower_bound(i));
      }
      if(bb.upper_bound(i) > result.upper_bound(i)) {
        result.set_upper_bound(i,bb.upper_bound(i));
      }
    }
  }
  return result;
}


template<class BS, class R>
tribool
disjoint(const ListSet<BS>& ls,
                   const Box<R>& bx)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls,bx,"tribool disjoint(ListSet<BS> ls, Box bx)");
  tribool result=true;
  for (typename ListSet<BS>::const_iterator i=ls.begin(); i!=ls.end(); ++i) {
    result = result && disjoint(*i,bx);
    if(!result) { return result; }
  }
  return result;
}

template<class BS1, class BS2>
tribool
disjoint(const ListSet<BS1>& ls1,
                   const ListSet<BS2>& ls2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls1,ls2,"tribool disjoint(ListSet<BS> ls1, ListSet<BS> ls2)");
  tribool result=true;
  for (typename ListSet<BS1>::const_iterator i=ls1.begin(); i!=ls1.end(); ++i) {
    for (typename ListSet<BS2>::const_iterator j=ls2.begin(); j!=ls2.end(); ++j) {
      result = result && disjoint(*i,*j);
      if(!result) { return result; }
    }
  }
  return result;
}



template<class BS1, class BS2>
tribool
subset(const ListSet< BS1 >& ls,
                 const BS2& bs)
{
  tribool result=true;
  for(typename ListSet< BS1 >::const_iterator ls_iter=ls.begin();
      ls_iter!=ls.end(); ++ls_iter)
    {
      result = result && subset(*ls_iter,bs);
      if(result==false) {
        return result;
      }
    }
  return result;
}


template<class BS>
ListSet<BS>
join(const ListSet<BS>& ls1,
               const ListSet<BS>& ls2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls1,ls2,"ListSet<BS> join(ListSet<BS> ls1, ListSet<BS> ls2)");
  ListSet<BS> ds_union(ls1);
  ds_union.inplace_union(ls2);
  return ds_union;
}


template<class BS>
ListSet<BS>
open_intersection(const ListSet<BS>& ls1,
                            const ListSet<BS>& ls2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls1,ls2,"ListSet<BS> open_intersection(ListSet<BS> ls1, ListSet<BS> ls2)");
  ListSet<BS> ds_inter(ls1.dimension());
  for (size_type i=0; i<ls1.size(); i++) {
    for (size_type j=0; j<ls2.size(); j++) {
      if (!disjoint(ls1[i],ls2[j])) {
        ds_inter.push_back(open_intersection(ls1[i],ls2[j]));
      }
    }
  }
  return ds_inter;
}



template<class BS>
std::string 
ListSet<BS>::summary() const
{
  std::stringstream ss;
  ss << "ListSet<"<<name<BS>()<<">( size="<<this->size()<<", dimension="<<this->dimension()<<" )";
  return ss.str();
}

template<class BS>
std::ostream& 
ListSet<BS>::write(std::ostream& os) const
{
  const ListSet<BS>& ls=*this;
  os << "ListSet<"<<name<BS>()<<">(\n size="<<ls.size()<<",\n";
  os << "  [ ";
  if (ls.size() >0 ) {
    os << ls[0];
  }
  for (size_type i=1; i<ls.size(); i++) {
    os << ",\n    " << ls[i];
  }
  os << "\n  ]\n";
  os << " }" << std::endl;
  return os;
}


template<class BS>
std::istream& 
ListSet<BS>::read(std::istream& is)
{
  ListSet<BS>& ls=*this;
  std::vector< BS >& vec(ls._vector);
  is >> vec;
  
  if(vec.size()==0) {
    ls._dimension = 0;
  }
  else {
    ls._dimension=vec[0].dimension();
  }
  return is;
}



} // namespace Ariadne
