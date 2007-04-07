/***************************************************************************
 *            list_set.code.h
 *
 *  23 June 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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
#include "../base/stlio.h"
#include "../geometry/rectangle.h"
#include "../geometry/irregular_grid_set.h"

namespace Ariadne {
  namespace Geometry {

    template<class BS>
    ListSet<BS>*
    ListSet<BS>::clone() const
    {
      return new ListSet<BS>(*this);
    }

    
    template<class BS> 
    tribool
    ListSet<BS>::intersects(const Rectangle<R>& r) const
    {
      return !Geometry::disjoint(*this,r);
    }

    
    template<class BS> 
    tribool
    ListSet<BS>::disjoint(const Rectangle<R>& r) const
    {
      return Geometry::disjoint(*this,r);
    }

    
    template<class BS>
    tribool
    ListSet<BS>::superset(const Rectangle<R>& r) const
    {
      const ListSet< Rectangle<R> >* rls=
        dynamic_cast< const ListSet< Rectangle<R> > * >(this);
      if(rls) {
        return Geometry::subset(r,IrregularGridMaskSet<R>(*rls));
      } else {
        return indeterminate;
      }
    }

    
    template<class BS>
    tribool
    ListSet<BS>::subset(const Rectangle<R>& r) const
    {
      return Geometry::subset(*this,r);
    }

    
    template<class BS>
    tribool
    disjoint (const ListSet<BS>& ls1,
              const ListSet<BS>& ls2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(ls1,ls2,"tribool disjoint(ListSet<BS> ls1, ListSet<BS> ls2)");
      tribool result=true;
      for (typename ListSet<BS>::const_iterator i=ls1.begin(); i!=ls1.end(); ++i) {
        for (typename ListSet<BS>::const_iterator j=ls2.begin(); j!=ls2.end(); ++j) {
          result = result && disjoint(*i,*j);
          if(!result) { return result; }
        }
      }
      return result;
    }


    template<class R>
    tribool
    subset(const ListSet< Rectangle<R> >& rls1,
           const ListSet< Rectangle<R> >& rls2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(rls1,rls2,"tribool subset(ListSet<Rectangle> rls1, ListSet<Rectangle> rls2)");
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    

    template<class R, class BS>
    tribool
    disjoint(const ListSet< BS >& ls,
             const Rectangle<R>& r)
    {
      tribool result=true;
      for(typename ListSet< BS >::const_iterator ls_iter=ls.begin();
          ls_iter!=ls.end(); ++ls_iter)
      {
        result = result && Geometry::disjoint(*ls_iter,r);
        if(result==false) {
          return result;
        }
      }
      return result;
    }


    template<class R, class BS>
    tribool
    subset(const ListSet< BS >& ls,
           const Rectangle<R>& r)
    {
      tribool result=true;
      for(typename ListSet< BS >::const_iterator ls_iter=ls.begin();
          ls_iter!=ls.end(); ++ls_iter)
      {
        result = result && Geometry::subset(*ls_iter,r);
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
    void
    ListSet<BS>::_instantiate_geometry_operators()
    {
      //Rectangle<R>* r=0;
      //ListSet<BS>* ls=0;
      //disjoint(*ls,*r);
    }
    
    
    template<class BS>
    std::ostream& 
    ListSet<BS>::write(std::ostream& os) const
    {
      const ListSet<BS>& ls=*this;
      os << "ListSet<"<<BS::name()<<">(\n  size="<<ls.size()<<",\n";
      if(!ls.empty()) { os << "  front="<<ls[0]<<",\n"; }
      // Don't write list
      if(false) {
        os << "  [ ";
        if (ls.size() >0 ) {
          os << ls[0];
        }
        for (size_type i=1; i<ls.size(); i++) {
          os << ",\n    " << ls[i];
        }
        os << "\n  ]\n";
      }
      os << "}" << std::endl;
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
  

  }
}
