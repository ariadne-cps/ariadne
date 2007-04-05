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
    disjoint (const ListSet<BS>& A,
              const ListSet<BS>& B)
    {
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      tribool result=true;
      for (typename ListSet<BS>::const_iterator i=A.begin(); i!=A.end(); ++i) {
        for (typename ListSet<BS>::const_iterator j=B.begin(); j!=B.end(); ++j) {
          result = result && disjoint(*i,*j);
          if(!result) { return result; }
        }
      }
      return result;
    }


    template<class R>
    tribool
    subset(const ListSet< Rectangle<R> >& A,
           const ListSet< Rectangle<R> >& B)
    {
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
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
    join(const ListSet<BS>& A,
         const ListSet<BS>& B)
    {
      ListSet<BS> ds_union(A);
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      ds_union.inplace_union(B);
      return ds_union;
    }


    template<class BS>
    ListSet<BS>
    open_intersection(const ListSet<BS>& A,
                      const ListSet<BS>& B)
    {
      ListSet<BS> ds_inter(A.dimension());
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for (size_type i=0; i<A.size(); i++) {
        for (size_type j=0; j<B.size(); j++) {
          if (!disjoint(A[i],B[j])) {
              ds_inter.push_back(open_intersection(A[i],B[j]));
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
      const ListSet<BS>& A=*this;
      os << "ListSet<" << Numeric::name<R>() << ",BS>{\n  ";
      os << "[ ";
      if (A.size() >0 ) {
        os << A[0];
      }
      for (size_type i=1; i<A.size(); i++) {
        os << ",\n    " << A[i];
      }
      os << "\n  ]\n}" << std::endl;
      return os;
    }


    template<class BS>
    std::istream& 
    ListSet<BS>::read(std::istream& is)
    {
      ListSet<BS>& A=*this;
      std::vector< BS >& vec(A._vector);
      is >> vec;

      if(vec.size()==0) {
        A._dimension = 0;
      }
      else {
        A._dimension=vec[0].dimension();
      }
      return is;
    }
  

  }
}
