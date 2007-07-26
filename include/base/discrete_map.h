/***************************************************************************
 *            discrete_map.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or5
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
 
/*! \file discrete_map.h
 *  \brief STL style discrete maps.
 */

#ifndef ARIADNE_DISCRETE_MAP_H
#define ARIADNE_DISCRETE_MAP_H

#include <cassert>
#include <iostream>
#include <map>



namespace Ariadne {
  
  namespace Base {
    
  
  template<class M>
  class _discrete_map_reference
  {
    typedef typename M::key_type K;
    typedef typename M::data_type D;
   public:
    _discrete_map_reference(M& m, const K& k) : map(m), key(k) { }
    void operator=(const D& d) { map._set_image(key,d); }
    operator const D& () const { return static_cast<const M&>(map)[key]; }
   private:
    M& map;
    const K& key;
  };
  
  template<class M> inline 
  std::ostream& operator<<( std::ostream& os, const _discrete_map_reference<M>& dmr)
  {
    return os << static_cast<const typename M::data_type&>(dmr);
  }


  template<class M>
  class _discrete_map_value
    : public std::pair<const typename M::key_type, typename M::data_type>
  {
   public:
    const typename M::key_type& key() const { return this->first; }
    typename M::data_type& data() { return this->second; }
    const typename M::data_type& data() const { return this->second; }
  };


  template<class Iter, class Value>
  class _discrete_map_iterator
    : public Iter
  {
   public:
    template<class T> _discrete_map_iterator(const T& t) : Iter(t) { }
    Value& operator*() { return reinterpret_cast<Value&>(this->Iter::operator*()); }
    Value* operator->() { return reinterpret_cast<Value*>(this->Iter::operator->()); }
  };
  
  
    

  /*! \ingroup Storage
   *  \brief Extended functionality for maps.
   */
  template<class K, class D>
  class discrete_map
    : public std::map<K,D>
  {
    friend class _discrete_map_reference< discrete_map<K,D> >;
   public:
    typedef D data_type;
    typedef const D& const_reference;
    typedef _discrete_map_reference< discrete_map<K,D> > reference;

    typedef _discrete_map_iterator<typename std::map<K,D>::iterator, _discrete_map_value< discrete_map<K,D> > > iterator;
    typedef _discrete_map_iterator<typename std::map<K,D>::const_iterator, const _discrete_map_value< discrete_map<K,D> > > const_iterator;
    
    discrete_map() : std::map<K,D>() { }

    bool has_key(const K& k) const {
      return this->find(k)!=this->end();
    }

    void insert(const K& k, const D& d) {
      assert(!this->has_key(k));
      this->std::map<K,D>::insert(std::make_pair(k,d));
    }
    
    const_reference operator[](const K& k) const { 
      typename std::map<K,D>::const_iterator iter=this->find(k);
      if(iter==this->end()) {
        assert(this->has_key(k));
      } else {
        return iter->second;
      }
    }

    reference operator[](const K& k) { 
      return reference(*this,k);
    }

    iterator begin() { return this->std::map<K,D>::begin(); }
    iterator end() { return this->std::map<K,D>::end(); }
    const_iterator begin() const { return this->std::map<K,D>::begin(); }
    const_iterator end() const { return this->std::map<K,D>::end(); }

   private:
    void _set_image(const K& k, const D& d) {
      typename std::map<K,D>::iterator iter=this->find(k);
      if(iter==this->end()) {
        this->std::map<K,D>::insert(std::make_pair(k,d));
      } else {
        iter->second=d;
      }
    }
  };


     
  } // namespace Base
} // namespace Ariadne
  
#endif /* ARIADNE_DISCRETE_MAP_H */
