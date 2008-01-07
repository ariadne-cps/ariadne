/***************************************************************************
 *            timed_list_set.h
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
 
/*! \file timed_list_set.h
 *  \brief A list of timed sets.
 */

#ifndef ARIADNE_TIMED_LIST_SET_H
#define ARIADNE_TIMED_LIST_SET_H

#include "timed_set.h"

namespace Ariadne {
  namespace Geometry {

    template<class BS> class ListSet;

    template<class T, class BS>
    class ListSet< TimedSet<T,BS> > 
    {
     private:
      typedef typename BS::real_type R;
      typedef TimedSet<T,BS> TBS;
      std::vector< TBS > _vector;
     public:
      typedef R real_type;
      typedef TBS basic_set_type;

      typedef typename std::vector<basic_set_type>::const_iterator const_iterator;
      typedef typename std::vector<basic_set_type>::iterator iterator;

     public:
      ListSet() : _vector() { }
      ListSet(const BS& bs) : _vector(1,TBS(0,bs)) { }
      ListSet(const TBS& tbs) : _vector(1,tbs) { }

      template<class DS> ListSet(const DS& ds) {
        _vector.reserve(ds.size()); 
        for(size_type i=0; i!=ds.size(); ++i) {
          _vector.push_back(TBS(0,ds[i]));
        } 
      }

      dimension_type dimension() const {
        if(this->_vector.empty()) { return 0; } 
        else { return this->_vector.front().dimension(); } }
      void clear() {
        this->_vector.clear(); }
      size_type size() const {
        return this->_vector.size(); }
      const TBS& operator[](size_type i) const {
        return this->_vector[i]; }
      TBS pop() {
        TBS result=_vector.back(); _vector.pop_back(); return result; }
      void adjoin(const TBS& s) {
        _vector.push_back(s); }

      void adjoin(const T& t, const ListSet<BS>& ls) {
        _vector.reserve(this->size()+ls.size()); 
        for(size_type i=0; i!=ls.size(); ++i) {
          _vector.push_back(TBS(t,ls[i]));
        } 
      }

      void adjoin(const ListSet< TimedSet<T,BS> >& tls) {
        _vector.reserve(this->size()+tls.size()); 
        for(size_type i=0; i!=tls.size(); ++i) {
          _vector.push_back(tls[i]);
        } 
      }
        
      const_iterator begin() const {
        return this->_vector.begin(); }
      const_iterator end() const {
        return this->_vector.begin(); }
      std::ostream& write(std::ostream& os) const {
        return os << "TimedListSet("<<this->_vector<<")"; }
    };
  
    template<class T, class BS> inline
    std::ostream& operator<<(std::ostream& os, const ListSet< TimedSet<T,BS> >& tls) {
      return tls.write(os); 
    }

  
  }
}


#endif /* ARIADNE_TIMED_LIST_SET_H */
