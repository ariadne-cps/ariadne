/***************************************************************************
 *            pointer.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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

/*! \file pointer.h
 *  \brief Smart pointers.
 */

#ifndef ARIADNE_POINTER_H
#define ARIADNE_POINTER_H

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

namespace Ariadne {

    using boost::shared_ptr;
    using boost::scoped_ptr;

} // namespace Ariadne

#ifdef DOXYGEN
namespace Ariadne {
  /*! \brief Reference-counted shared pointer. Initialise with a dynamically-allocated object; destruction is performed automatically. */
  template<class T> class shared_ptr { };
  /*! \brief Scope smart pointer. Initialise with a dynamically-allocated object; destruction is performed automatically when object goes out of scope. Cannot be used in containers. */
  template<class T> class scoped_ptr { };
}
#endif

namespace Ariadne {

  template<class T> 
  class copy_ptr
  {
   public:
    ~copy_ptr() { delete _ptr; }
    copy_ptr() : _ptr(new T()) { }
    copy_ptr(T* t) : _ptr(t) { }
    copy_ptr(const copy_ptr<T>& p) : _ptr(new T(*p._ptr)) { }
    copy_ptr& operator=(const copy_ptr<T>& p) { if(this->_ptr!=p._ptr) { *this->_ptr=*p._ptr; } return *this; }
    
    const T* operator->() const { return this->_ptr; }
    T* operator->() { return this->_ptr; }
    
    const T& operator*() const { return *this->_ptr; }
    T& operator*() { return *this->_ptr; }
   private:
    T* _ptr;
  };

}

#endif /* ARIADNE_POINTER_H */
