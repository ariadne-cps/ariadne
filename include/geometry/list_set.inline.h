/***************************************************************************
 *            list_set.h
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

#include "geometrical_traits.h"

namespace Ariadne {
  namespace Geometry {

    template<class BS> inline
    ListSet<BS>::ListSet()
      : _dimension(0), _vector()
    {
    }

    template<class BS> inline
    ListSet<BS>::ListSet(size_type n) 
      : _dimension(n), _vector() 
    {
    }

    template<class BS> inline
    ListSet<BS>::ListSet(const BS& A)
      : _dimension(A.dimension()), _vector()
    {
      if (A.empty()) {
        return;
      }
      _vector.push_back(A);
    }

    template<class BS> inline
    ListSet<BS>::ListSet(const ListSet<BS>& A)
      : _dimension(A.dimension()), _vector(A._vector) 
    {
    }

    template<class BS> inline
    ListSet<BS>::~ListSet() 
    {
      this->_vector.clear();
    }


    template<class BS> inline
    size_type 
    ListSet<BS>::size() const 
    {
      return this->_vector.size();
    }

    template<class BS> inline
    void 
    ListSet<BS>::push_back(const BS& bs) 
    {
      if (this->dimension()==0) { this->_dimension=bs.dimension(); }
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,bs,"void ListSet<BS>::push_back(BS bs)");
      this->_vector.push_back(bs);
    }

    template<class BS> inline
    void 
    ListSet<BS>::pop_back() 
    {
      if (this->_vector.empty()) { 
        ARIADNE_THROW(LengthError,"void ListSet<BS>::pop_back()"," empty list");
      }
      this->_vector.pop_back();
    }

    template<class BS> inline
    dimension_type 
    ListSet<BS>::dimension() const 
    {
      return this->_dimension;
    }

    template<class BS> inline
    const BS& 
    ListSet<BS>::get(size_type index) const 
    {
      ARIADNE_CHECK_ARRAY_INDEX(*this,index,"BS ListSet<BS>::get(size_type index)");
      return this->_vector[index];
    }

    template<class BS> inline
    void 
    ListSet<BS>::set(size_type index, const BS& set) 
    {
      ARIADNE_CHECK_ARRAY_INDEX(*this,index,"void ListSet<BS>::set(size_type index, BS set)");
      this->_vector[index]=set;
    }

    template<class BS> inline
    const BS& 
    ListSet<BS>::operator[](size_type index) const 
    {
      ARIADNE_CHECK_ARRAY_INDEX(*this,index,"void ListSet<BS>::operator[](size_type index)");
      return this->_vector[index];
    }


    template<class BS> inline
    const ListSet<BS>& 
    ListSet<BS>::operator=(const ListSet<BS>& A) 
    {
      if(this !=& A) {
        this->_dimension = A._dimension;
        this->_vector = A._vector;
      }
      return *this;
    }

    template<class BS> 
    template<class BSt> 
    inline
    ListSet<BS>::operator ListSet<BSt> () const 
    {
      ListSet<BSt> result(this->dimension());
      BSt bs(this->dimension());
      for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        bs=BSt(*iter);
        result.push_back(bs);
      }
      return result;
    }
      
    template<class BS> inline
    tribool 
    ListSet<BS>::contains(const Point<R>& p) const 
    {
      tribool result=false;
      for (typename ListSet<BS>::const_iterator i=this->begin(); i!=this->end(); ++i) {
        result=result || i->contains(p);
        if(result) { return result; }
      }
      return result;
    }

    template<class BS> inline
    tribool 
    ListSet<BS>::empty() const 
    {
      tribool result=true;
      for (typename ListSet<BS>::const_iterator i=this->begin(); i!=this->end(); ++i) {
        result = result && i->empty();
        if(!result) { return result; }
      }
      return result;
    }

    template<class BS> inline
    tribool 
    ListSet<BS>::bounded() const 
    {
      tribool result=true;
      for (typename ListSet<BS>::const_iterator i=this->begin(); i!=this->end(); ++i) {
        result = result && i->bounded();
        if(!result) { return result; }
      }
      return result;
    }

    template<class BS> inline
    Rectangle<typename ListSet<BS>::real_type> 
    ListSet<BS>::bounding_box() const 
    {
      if(this->empty()) { return Rectangle<R>(this->dimension()); }
      Rectangle<R> result=(*this)[0].bounding_box();
      for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        Rectangle<R> bb=iter->bounding_box();
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
        
    template<class BS> inline
    void 
    ListSet<BS>::clear() 
    { 
      this->_vector.clear();
    }
      
    template<class BS> inline
    typename ListSet<BS>::iterator 
    ListSet<BS>::begin() 
    {
      return _vector.begin();
    }

    template<class BS> inline
    typename ListSet<BS>::iterator 
    ListSet<BS>::end() 
    {
      return _vector.end();
    }

    template<class BS> inline
    typename ListSet<BS>::const_iterator 
    ListSet<BS>::begin() const 
    {
      return _vector.begin();
    }

    template<class BS> inline
    typename ListSet<BS>::const_iterator 
    ListSet<BS>::end() const 
    {
      return _vector.end();
    }

    template<class BS> inline
    void 
    ListSet<BS>::adjoin(const ListSet<BS>& ls) 
    {
      if(this->dimension()==0) { 
        this->_dimension=ls.dimension(); 
      }
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,ls,"void ListSet<BS>::adjoin(ListSet<BS> ls)");
      this->_vector.reserve(ls.size());
      for(typename ListSet<BS>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
        this->_vector.push_back(*iter);
      }
    }
      
    template<class BS> inline
    void 
    ListSet<BS>::adjoin(const BS& bs) 
    {
      if(this->dimension()==0) { 
        this->_dimension=bs.dimension(); 
      }
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,bs,"void ListSet<BS>::adjoin(BS bs)");
      if(bs.empty()) {
      } else {
        this->_vector.push_back(bs);
      }
    }


    template<class BS> template<class S> inline
    void 
    ListSet<BS>::adjoin(const S& s) 
    {
      this->adjoin(s, typename geometrical_traits<S>::set_category());
    }


    template<class BS> template<class S> inline
    void 
    ListSet<BS>::adjoin(const S& bs, basic_set_tag) 
    {
      this->adjoin(static_cast<BS>(bs));
    }


    template<class BS> template<class DS> inline
    void 
    ListSet<BS>::adjoin(const DS& ds, denotable_set_tag) 
    {
      for(typename DS::const_iterator iter=ds.begin();
          iter!=ds.end(); ++iter)
      {
        this->adjoin(static_cast<BS>(*iter));
      }
    }



      
    template<class BS> inline
    std::ostream& 
    operator<<(std::ostream& os, const ListSet<BS>& ls)
    {
      return ls.write(os);
    }

    template<class BS> inline
    std::istream& 
    operator>>(std::istream& is, ListSet<BS>& ls)
    {
      return ls.read(is);
    }

  }
}
