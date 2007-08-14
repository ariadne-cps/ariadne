/***************************************************************************
 *            container_utilities.h
 *
 *  Copyright  2007  Pieter Collins  pieter.collins@cwi.nl
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


#include <boost/iterator/indirect_iterator.hpp>



#ifndef ARIADNE_CONTAINER_UTILITIES_H
#define ARIADNE_CONTAINER_UTILITIES_H

namespace Ariadne {
  namespace Base {

    /*! \brief Compare two pointers by dereferencing their values and using the binary predicate object \a Compare. */
    template<class Compare> 
    class dereference_compare
    {
     public:
      typedef typename Compare::first_argument_type const* first_argument_type;
      typedef typename Compare::second_argument_type const* second_argument_type;
      typedef typename Compare::result_type result_type;
 
      result_type operator() (first_argument_type const ptr1, second_argument_type const ptr2) const {
        return Compare()(*ptr1,*ptr2); }
    };
 


    template<class Key, class Data>
    struct _map_value 
      : private std::pair<const Key,Data>
    {
     public:
      Key const& key() const { return this->std::pair<const Key,Data>::first; }
      Data & data() { return *this->std::pair<const Key,Data>::second; }
      Data const& data() const { return *this->std::pair<const Key,Data>::second; }
    };
      

    template<class Key, class Data>
    struct _reference_key_map_value 
      : private std::pair<Key const*,Data>
    {
     public:
      Key const& key() const { return *this->std::pair<Key const*,Data>::first; }
      Data & data() { return this->std::pair<Key const*,Data>::second; }
      Data const& data() const { return this->std::pair<Key const*,Data>::second; }
    };
      

    template<class Key, class Data>
    struct _reference_data_map_value 
      : private std::pair<const Key,Data*>
    {
     public:
      Key const& key() const { return this->std::pair<const Key,Data*>::first; }
      Data & data() { return *this->std::pair<const Key,Data*>::second; }
      Data const& data() const { return *this->std::pair<const Key,Data*>::second; }
    };
      

    template<class Map>
    class _map_reference
    {
      typedef typename Map::key_type Key;
      typedef typename Map::data_type Data;
     public:
      _map_reference(Map& m, Key const& k) : map(m), key(k) { }
      operator Data& () { return map.get(key); }
      Data& operator=(Data& x) { map.set(key,x); }
     private:
      Map& map;
      Key const& key;
    };


    template<class Base, class Value, class Ref>
    class _map_iterator
      : public boost::iterator_adaptor<_map_iterator<Base,Value,Ref>,
                                       Base,Value,boost::use_default,Ref>
    {
     public:
      _map_iterator() : _map_iterator::iterator_adaptor_() { }
      _map_iterator(Base iter) : _map_iterator::iterator_adaptor_(iter) { }
      Ref dereference() const { return reinterpret_cast<Ref>(*this->base_reference()); }
     private:
      friend class boost::iterator_core_access;
    };


    template<class Base, class Value, class Ref>
    class _indirect_map_iterator
      : public boost::iterator_adaptor<_map_iterator<Base,Value,Ref>,
                                       Base,Value,boost::use_default,Ref>
    {
     public:
      _indirect_map_iterator() : _indirect_map_iterator::iterator_adaptor_() { }
      _indirect_map_iterator(Base iter) : _indirect_map_iterator::iterator_adaptor_(iter) { }
      template<class Iter> _indirect_map_iterator(Iter iter) : _indirect_map_iterator::iterator_adaptor_(Base(iter.base())) { }
      Ref dereference() const { return reinterpret_cast<Ref>(*this->base_reference()); }
      template<class Iter> bool operator!=(const Iter& iter) { return this->base()!=iter.base(); }
     private:
      friend class boost::iterator_core_access;
    };


  }
}

#endif // ARIADNE_CONTAINER_UTILITIES_H
