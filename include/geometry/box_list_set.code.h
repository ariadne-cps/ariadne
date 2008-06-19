/***************************************************************************
 *            box_list_set.code.h
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

#include "base/stlio.h"
#include "box_list_set.h"

namespace Ariadne {

template<class R>
void
BoxListSet<R>::_instantiate()
{
  BoxListSet<R>* bxls=0;
  Ariadne::open_intersection(*bxls,*bxls);
  Ariadne::subset(*bxls,*bxls);
  Ariadne::disjoint(*bxls,*bxls);
}

template<class R>
BoxListSet<R>*
BoxListSet<R>::clone() const
{
  return new BoxListSet<R>(*this);
}


template<class R> 
tribool
BoxListSet<R>::intersects(const Box<R>& r) const
{
  return !Ariadne::disjoint(*this,r);
}


template<class R> 
tribool
BoxListSet<R>::disjoint(const Box<R>& r) const
{
  return Ariadne::disjoint(*this,r);
}


template<class R>
tribool
BoxListSet<R>::superset(const Box<R>& r) const
{
  BoxListSet<R> ls;
  ls.adjoin(r);
  return Ariadne::subset(ls,*this);
}


template<class R>
tribool
BoxListSet<R>::subset(const Box<R>& r) const
{
  return Ariadne::subset(*this,r);
}


template<class R>
tribool
disjoint(const BoxListSet<R>& ls1,
                   const BoxListSet<R>& ls2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls1,ls2,"tribool disjoint(BoxListSet<R> ls1, BoxListSet<R> ls2)");
  tribool result=true;
  for (typename BoxListSet<R>::const_iterator i=ls1.begin(); i!=ls1.end(); ++i) {
    for (typename BoxListSet<R>::const_iterator j=ls2.begin(); j!=ls2.end(); ++j) {
      result = result && Ariadne::disjoint(*i,*j);
      if(!result) { return result; }
    }
  }
  return result;
}



template<class R>
tribool
subset(const BoxListSet<R>& bxls1,
                 const BoxListSet<R>& bxls2)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
tribool
subset(const BoxListSet<R>& ls,
                 const Box<R>& bx)
{
  tribool result;
  for(typename BoxListSet<R>::const_iterator ls_iter=ls.begin();
      ls_iter!=ls.end(); ++ls_iter)
    {
      result = result && subset(*ls_iter,bx);
      if(result==false) {
        return result;
      }
    }
  return result;
}



template<class R>
BoxListSet<R>
join(const BoxListSet<R>& ls1,
               const BoxListSet<R>& ls2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls1,ls2,"BoxListSet<R> join(BoxListSet<R> ls1, BoxListSet<R> ls2)");
  BoxListSet<R> result(ls1);
  result.adjoin(ls2);
  return result;
}


template<class R>
BoxListSet<R>
open_intersection(const BoxListSet<R>& ls1,
                            const BoxListSet<R>& ls2)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls1,ls2,"BoxListSet<R> open_intersection(BoxListSet<R> ls1, BoxListSet<R> ls2)");
  BoxListSet<R> ds_inter;
  for (size_type i=0; i<ls1.size(); i++) {
    for (size_type j=0; j<ls2.size(); j++) {
      if (!disjoint(ls1[i],ls2[j])) {
        ds_inter.push_back(open_intersection(ls1[i],ls2[j]));
      }
    }
  }
  return ds_inter;
}

template<class R>
BoxListSet<R>
inner_intersection(const BoxListSet<R>& ls,
                             const SetInterface<R>& ops)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls,ops,"BoxListSet<R> inner_intersection(BoxListSet<R> ls, SetInterface<R> ops)");
  BoxListSet<R> ds(ls.dimension());
  for (size_type i=0; i<ls.size(); i++) {
    if(ops.subset(ls[i].bounding_box())) {
      ds.push_back(ls[i]);
    }
  }
  return ds;
}

template<class R>
BoxListSet<R>
lower_intersection(const BoxListSet<R>& ls,
                             const SetInterface<R>& cls)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls,cls,"BoxListSet<R> inner_intersection(BoxListSet<R> ls, SetInterface<R> cls)");
  BoxListSet<R> ds(ls.dimension());
  for (size_type i=0; i<ls.size(); i++) {
    if(cls.intersects(ls[i].bounding_box())) {
      ds.push_back(ls[i]);
    }
  }
  return ds;
}

template<class R>
BoxListSet<R>
outer_intersection(const BoxListSet<R>& ls,
                             const SetInterface<R>& cps)
{
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls,cps,"BoxListSet<R> inner_intersection(BoxListSet<R> ls, SetInterface<R> cps)");
  BoxListSet<R> ds(ls.dimension());
  for (size_type i=0; i<ls.size(); i++) {
    if(cps.disjoint(ls[i].bounding_box())) {
    } else {
      ds.push_back(ls[i]);
    }
  }
  return ds;
}





template<class R>
std::string 
BoxListSet<R>::summary() const
{
  std::stringstream ss;
  ss << "BoxListSet( size="<<this->size()<<", dimension="<<this->dimension()<<" )";
  return ss.str();
}

template<class R>
std::ostream& 
BoxListSet<R>::write(std::ostream& os) const
{
  const BoxListSet<R>& ls=*this;
  os << "BoxListSet(\n size="<<ls.size()<<",\n";
  os << "  [ ";
  if (ls.size()>0 ) {
    os << ls[0];
  }
  for (size_type i=1; i<ls.size(); i++) {
    os << ",\n    " << ls[i];
  }
  os << "\n  ]\n";
  os << " }" << std::endl;
  return os;
}


template<class R>
std::istream& 
BoxListSet<R>::read(std::istream& is)
{
  BoxListSet<R>& ls=*this;
  std::vector< Box<R> >& vec(ls._vector);
  is >> vec;
  return is;
}



} // namespace Ariadne

