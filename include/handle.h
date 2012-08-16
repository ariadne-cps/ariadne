/***************************************************************************
 *            handle.h
 *
 *  Copyright 2011  Pieter Collins
 *
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

/*! \file handle.h
 *  \brief Generic handle base class
 */

#ifndef ARIADNE_HANDLE_H
#define ARIADNE_HANDLE_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "pointer.h"

namespace Ariadne {

template<class T> class Handle
{
  private:
    clone_on_write_ptr<T> _ptr;
  public:
    Handle() : _ptr() { }
    Handle(T* p) : _ptr(p) { }
    Handle(const T& t) : _ptr(t._clone()) { }
    Handle(std::shared_ptr<const T> p) : _ptr(p->_clone()) { }
    clone_on_write_ptr<T> managed_pointer() { return _ptr; }
    T* raw_pointer()  { return _ptr.operator->(); }
    const T* raw_pointer() const  { return _ptr.operator->(); }
    T& reference()  { return _ptr.operator*(); }
    const T& reference() const  { return _ptr.operator*(); }
    operator T& () { return _ptr->operator*(); }
    operator const T& () const { return _ptr->operator*(); }
    Handle<T>& operator=(const T& t) { _ptr(clone_on_write_ptr<T>(t._clone())); }
};

template<class T> class Handle<const T>
{
  private:
    std::shared_ptr<const T> _ptr;
  public:
    Handle() : _ptr() { }
    Handle(T* p) : _ptr(p) { }
    Handle(const T& t) : _ptr(t._clone()) { }
    Handle(std::shared_ptr<const T> p) : _ptr(p) { }
    std::shared_ptr<const T> managed_pointer() const  { return _ptr; }
    const T* raw_pointer() const  { return _ptr.operator->(); }
    const T& reference() const  { return _ptr.operator*(); }
    operator const T& () const { return _ptr->operator*(); }
    Handle<const T>& operator=(const T& t) { _ptr(std::shared_ptr<T>(t._clone())); }
};



} // namespace Ariadne

#endif /* ARIADNE_HANDLE_H */
