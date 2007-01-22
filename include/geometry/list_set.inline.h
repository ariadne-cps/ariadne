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

namespace Ariadne {
  namespace Geometry {

    template<class R, template<class> class BS> inline
    ListSet<R,BS>::ListSet()
      : _dimension(0), _vector()
    {
    }

    template<class R, template<class> class BS> inline
    ListSet<R,BS>::ListSet(size_type n) 
      : _dimension(n), _vector() 
    {
    }

    template<class R, template<class> class BS> inline
    ListSet<R,BS>::ListSet(const BS<R>& A)
      : _dimension(A.dimension()), _vector()
    {
      if (A.empty()) {
        return;
      }
      _vector.push_back(A);
    }

    template<class R, template<class> class BS> inline
    ListSet<R,BS>::ListSet(const ListSet<R,BS>& A)
      : _dimension(A.dimension()), _vector(A._vector) 
    {
    }

    template<class R, template<class> class BS> inline
    ListSet<R,BS>::~ListSet() 
    {
      this->_vector.clear();
    }


    template<class R, template<class> class BS> inline
    size_type 
    ListSet<R,BS>::size() const 
    {
      return this->_vector.size();
    }

    template<class R, template<class> class BS> inline
    void 
    ListSet<R,BS>::push_back(const BS<R>& A) 
    {
      if (this->dimension()==0) { this->_dimension=A.dimension(); }
      check_equal_dimensions(*this,A,__PRETTY_FUNCTION__);
      this->_vector.push_back(A);
    }

    template<class R, template<class> class BS> inline
    void 
    ListSet<R,BS>::pop_back() 
    {
      if (this->_vector.empty()) { 
        throw std::runtime_error("Attempting to pop from an empty ListSet");
      }
      this->_vector.pop_back();
    }

    template<class R, template<class> class BS> inline
    dimension_type 
    ListSet<R,BS>::dimension() const 
    {
      return this->_dimension;
    }

    template<class R, template<class> class BS> inline
    const BS<R>& 
    ListSet<R,BS>::get(size_type index) const 
    {
      check_array_index(*this,index,__PRETTY_FUNCTION__);
      return this->_vector[index];
    }

    template<class R, template<class> class BS> inline
    void 
    ListSet<R,BS>::set(size_type index, const BS<R>& set) 
    {
      check_array_index(*this,index,__PRETTY_FUNCTION__);
      this->_vector[index]=set;
    }

    template<class R, template<class> class BS> inline
    const BS<R>& 
    ListSet<R,BS>::operator[](size_type index) const 
    {
      check_array_index(*this,index,__PRETTY_FUNCTION__);
      return this->_vector[index];
    }


    template<class R, template<class> class BS> inline
    const ListSet<R,BS>& 
    ListSet<R,BS>::operator=(const ListSet<R,BS>& A) 
    {
      if(this !=& A) {
        this->_dimension = A._dimension;
        this->_vector = A._vector;
      }
      return *this;
    }

    template<class R, template<class> class BS> 
    template<class Rl, template<class> class BSt> 
    inline
    ListSet<R,BS>::operator ListSet<Rl,BSt> () const 
    {
      ListSet<Rl,BSt> result(this->dimension());
      BSt<Rl> bs(this->dimension());
      for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        bs=BSt<Rl>(*iter);
        result.push_back(bs);
      }
      return result;
    }
      
    template<class R, template<class> class BS> inline
    tribool 
    ListSet<R,BS>::contains(const Point<R>& p) const 
    {
      tribool result=false;
      for (typename ListSet<R,BS>::const_iterator i=this->begin(); i!=this->end(); ++i) {
        result=result || i->contains(p);
        if(result) { return result; }
      }
      return result;
    }

    template<class R, template<class> class BS> inline
    tribool 
    ListSet<R,BS>::empty() const 
    {
      tribool result=true;
      for (typename ListSet<R,BS>::const_iterator i=this->begin(); i!=this->end(); ++i) {
        result = result && i->empty();
        if(!result) { return result; }
      }
      return result;
    }

    template<class R, template<class> class BS> inline
    tribool 
    ListSet<R,BS>::bounded() const 
    {
      tribool result=true;
      for (typename ListSet<R,BS>::const_iterator i=this->begin(); i!=this->end(); ++i) {
        result = result && i->bounded();
        if(!result) { return result; }
      }
      return result;
    }

    template<class R, template<class> class BS> inline
    Rectangle<R> 
    ListSet<R,BS>::bounding_box() const 
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
        
    template<class R, template<class> class BS> inline
    void 
    ListSet<R,BS>::clear() 
    { 
      this->_vector.clear();
    }
      
    template<class R, template<class> class BS> inline
    typename ListSet<R,BS>::iterator 
    ListSet<R,BS>::begin() 
    {
      return _vector.begin();
    }

    template<class R, template<class> class BS> inline
    typename ListSet<R,BS>::iterator 
    ListSet<R,BS>::end() 
    {
      return _vector.end();
    }

    template<class R, template<class> class BS> inline
    typename ListSet<R,BS>::const_iterator 
    ListSet<R,BS>::begin() const 
    {
      return _vector.begin();
    }

    template<class R, template<class> class BS> inline
    typename ListSet<R,BS>::const_iterator 
    ListSet<R,BS>::end() const 
    {
      return _vector.end();
    }

    template<class R, template<class> class BS> inline
    void 
    ListSet<R,BS>::adjoin(const ListSet<R,BS>& A) 
    {
      if(this->dimension()==0) { 
        this->_dimension=A.dimension(); 
      }
      check_equal_dimensions(*this,A,__PRETTY_FUNCTION__);
      this->_vector.reserve(A.size());
      for(typename ListSet<R,BS>::const_iterator iter=A.begin(); iter!=A.end(); ++iter) {
        this->_vector.push_back(*iter);
      }
    }
      
    template<class R, template<class> class BS> inline
    void 
    ListSet<R,BS>::inplace_union(const ListSet<R,BS>& A) 
    {
      this->adjoin(A);
    }

    template<class R, template<class> class BS> inline
    void 
    ListSet<R,BS>::adjoin(const BS<R>& A) 
    {
      if(this->dimension()==0) { 
        this->_dimension=A.dimension(); 
      }
      check_equal_dimensions(*this,A,__PRETTY_FUNCTION__);
      if(!A.empty()) {
        this->_vector.push_back(A);
      }
    }

    template<class R, template<class> class BS> inline
    void 
    ListSet<R,BS>::inplace_union(const BS<R>& A) 
    {
      this->adjoin(A);
    }



    template<class R, template<class> class BS> inline
    std::ostream& 
    operator<<(std::ostream& os, const ListSet<R,BS>& ls)
    {
      return ls.write(os);
    }

    template<class R, template<class> class BS> inline
    std::istream& 
    operator>>(std::istream& is, ListSet<R,BS>& ls)
    {
      return ls.read(is);
    }

  }
}
