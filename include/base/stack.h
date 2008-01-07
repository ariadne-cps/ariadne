/***************************************************************************
 *            base/stack.h
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

#include "sequence_io.h"

#ifndef ARIADNE_BASE_STACK_H
#define ARIADNE_BASE_STACK_H

namespace Ariadne {
  namespace Base {

    template<class T>
    class stack
    {
     public:
      typedef typename std::vector<T>::size_type size_type;
      typedef typename std::vector<T>::const_iterator iterator;
      typedef typename std::vector<T>::const_iterator const_iterator;

      stack() { }

      template<class C> stack(const C& c) {
        _vector.reserve(c.size());
        for(typename C::const_iterator iter=c.begin();
            iter!=c.end(); ++iter) {
          _vector.push_back(T(*iter)); } }

      template<class Iter> stack(Iter curr, Iter last) {
        for( ; curr!=last; ++curr) { _vector.push_back(T(*curr)); } }

      bool empty() const { return _vector.empty(); }
      size_type size() const { return _vector.size(); }
      void push(const T& t) { _vector.push_back(t); }
      void pop(T& t) { t=_vector.back(); _vector.pop_back(); }
      T pop() { T t=_vector.back(); _vector.pop_back(); return t; }

      const_iterator begin() const { return _vector.begin(); }
      const_iterator end() const { return _vector.end(); }
     private:
      std::vector<T> _vector;
    };

    template<class T>
    std::ostream& operator<<(std::ostream& os, const stack<T>& stack) {
      return Base::write_sequence(os,stack.begin(),stack.end()); 
    }


  }
}

#endif // ARIADNE_BASE_STACK_H
