/***************************************************************************
 *            utility/handle.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file utility/handle.h
 *  \brief
 */



#ifndef ARIADNE_HANDLE_H
#define ARIADNE_HANDLE_H

#include "stdlib.h"

#include <memory>

#include "typedefs.h"
#include "metaprogramming.h"

#include "writable.h"

namespace Ariadne {

/************ Handle *********************************************************/

template<class T> struct HandleObject { };

template<class I> class Handle;
template<class T, class I> const T& dynamic_handle_cast(const Handle<I>& h);
template<class T, class I> T& dynamic_handle_cast(Handle<I>& h);

template<class I, class T> class Wrapper
    : public virtual I, public T
{
  protected:
    Wrapper(const T& t) : T(t) { }
    Wrapper(T&& t) : T(std::move(t)) { }
};

template<class I> class Handle {
    typedef I Interface;
  protected:
    mutable SharedPointer<I> _ptr;
  public:
    ~Handle() { }
    Handle(I* p) : _ptr(p) { }
    Handle(SharedPointer<I> p) : _ptr(p) { }
    Handle(I&& r) : _ptr(r._move()) { }
    Handle(const I& r) : _ptr(r._copy()) { }
    //Handle(const Handle<I>& h) : _ptr(h._ptr) { }
    //Handle(Handle<I>&& h) : _ptr(h._ptr) { }
    Handle(const Handle<I>& h) = default;
    Handle(Handle<I>&& h) = default;
    Handle<I>& operator=(const Handle<I>& h) = default;
    Handle<I>& operator=(Handle<I>&& h) = default;
    template<class II> Handle(const Handle<II>& h) : _ptr(h.managed_pointer()) { }
    //explicit operator I const& () const { return *_ptr; }
    operator I const& () const { return *_ptr; }
  public:
    const I& const_reference() const { return _ptr.operator*(); }
    const I& reference() const { return _ptr.operator*(); }
    I& reference() { make_unique(); return _ptr.operator*(); }
  public:
    const I* const_pointer() const { return _ptr.operator->(); }
    const I* pointer() const { return _ptr.operator->(); }
    I* pointer() { make_unique(); return _ptr.operator->(); }
  public:
    const I* raw_const_pointer() const { return _ptr.operator->(); }
    const I* raw_pointer() const { return _ptr.operator->(); }
    I* raw_pointer() { make_unique(); return _ptr.operator->(); }
  public:
    const I& ref() const { return _ptr.operator*(); }
    I& ref() { make_unique(); return _ptr.operator*(); }
  public:
    const I* ptr() const { return _ptr.operator->(); }
    I* ptr() { make_unique(); return _ptr.operator->(); }
  public:
    const SharedPointer<I> managed_const_pointer() const { return _ptr; }
    const SharedPointer<I> managed_pointer() const { return _ptr; }
    SharedPointer<I> managed_pointer() { make_unique(); return _ptr; }
  protected:
    void make_unique() { if(!_ptr.unique()) { _ptr=SharedPointer<I>(_ptr->_copy()); } }
  private:
    template<class T, class II> friend T& dynamic_handle_cast(Handle<II>& h);
    template<class T, class II> friend const T& dynamic_handle_cast(const Handle<II>& h);
};

void write_error(OutputStream& os, const WritableInterface* w, const char* i, const char* c, const char* t);

template<class T, class I> const T& dynamic_handle_cast(const Handle<I>& h) {
    const I* i=h.raw_const_pointer();
    const T* p=dynamic_cast<const Wrapper<I,T>*>(i);
    if(p) { return *p; }
    const WritableInterface* w=dynamic_cast<const WritableInterface*>(i);
    write_error(std::cerr,w,typeid(i).name(),typeid(*i).name(),typeid(T).name());
    throw std::bad_cast();
}

template<class T, class I> T& dynamic_handle_cast(Handle<I>& h) {
    const I* i=h.raw_const_pointer();
    const T* p=dynamic_cast<const Wrapper<I,T>*>(i);
    if(p) { return const_cast<T&>(*p); }
    const WritableInterface* w=dynamic_cast<const WritableInterface*>(i);
    write_error(std::cerr,w,typeid(i).name(),typeid(*i).name(),typeid(T).name());
    throw std::bad_cast();
}



} // namespace Ariadne

#endif
