/***************************************************************************
 *            utility/handle.hpp
 *
 *  Copyright  2013-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file utility/handle.hpp
 *  \brief
 */



#ifndef ARIADNE_HANDLE_HPP
#define ARIADNE_HANDLE_HPP

#include "stdlib.hpp"

#include <memory>

#include "typedefs.hpp"
#include "metaprogramming.hpp"

#include "writable.hpp"

namespace Ariadne {

/************ Handle *********************************************************/

template<class T> struct HandleObject { };

template<class I> class Handle;
template<class T, class I> const T& dynamic_handle_cast(const Handle<I>& h);
template<class T, class I> T& dynamic_handle_cast(Handle<I>& h);

template<class T, class I> class Wrapper
    : public virtual I, public T
{
  protected:
    Wrapper(const T& t) : T(t) { }
    Wrapper(T&& t) : T(std::move(t)) { }
};

template<class I> class Handle {
  public:
    typedef I Interface;
  protected:
    mutable SharedPointer<I> _ptr;
  public:
    ~Handle() { }
    explicit Handle(I* p) : _ptr(p) { }
    Handle(SharedPointer<I> p) : _ptr(p) { }
    Handle(I&& r) : _ptr(r._move()) { }
    Handle(const I& r) : _ptr(r._copy()) { }
    template<DerivedFrom<I> T> Handle(T&& t) : _ptr(new T(std::move(t))) { }
    template<DerivedFrom<I> T> requires CopyConstructible<T> Handle(const T& t)
        : _ptr(new T(t)) { }
    //Handle(const Handle<I>& h) : _ptr(h._ptr) { }
    //Handle(Handle<I>&& h) : _ptr(h._ptr) { }
    Handle(const Handle<I>& h) = default;
    Handle(Handle<I>&& h) = default;
    Handle<I>& operator=(const Handle<I>& h) = default;
    Handle<I>& operator=(Handle<I>&& h) = default;
    template<class II> explicit Handle(const Handle<II>& h) : _ptr(h.managed_pointer()) { }
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
    template<class II> static Void _make_unique(SharedPointer<II>& ptr) {
        if(ptr.use_count()>1) { ptr=SharedPointer<I>(ptr->_copy()); } }
    template<class II> static Void _make_unique(SharedPointer<const II>& ptr) { }
    void make_unique() { _make_unique(_ptr); }
  private:
    template<class T, class II> friend T& dynamic_handle_extract(Handle<II>& h);
    template<class T, class II> friend const T& dynamic_handle_extract(const Handle<II>& h);
    template<class D, class B> friend Handle<D> dynamic_handle_cast(Handle<B>& h);
    template<class D, class B> friend const Handle<D> dynamic_handle_cast(const Handle<B>& h);
};

template<class T, class... AS> Handle<T> make_handle(AS&& ... as) {
    return Handle<T>(std::make_shared<T>(as...));
}

inline void write_error(OutputStream& os, const char* f, const WritableInterface* w, const char* i, const char* c, const char* t) {
    os << "Error in " << f << ":" << std::flush;
    os << " cannot convert "; if (w) { w->_write(os); } else { os << "object"; } os << std::flush;
    os << " of static type " << c << " and dynamic type " << i << " to type " << t << std::endl;
}

template<class D, class B> D dynamic_handle_cast(B const& h) {
    typedef typename B::Interface BI;
    typedef typename D::Interface DI;
    if constexpr (Same<decltype(h.managed_pointer()),SharedPointer<BI>>) {
        SharedPointer<DI> p=std::dynamic_pointer_cast<DI>(h.managed_pointer());
        if(p) { return D(Handle<DI>(p)); }
    } else {
        SharedPointer<const DI> p=std::dynamic_pointer_cast<const DI>(h.managed_pointer());
        if(p) { return D(Handle<const DI>(p)); }
    }
    const BI* i=h.raw_pointer();
    const WritableInterface* w=dynamic_cast<const WritableInterface*>(i);
    write_error(std::cerr,"dynamic_handle_cast",w,typeid(*i).name(),typeid(i).name(),typeid(D).name());
    throw std::bad_cast();
}

template<class T, class I> const T& dynamic_handle_extract(const Handle<const I>& h) {
    const I* i=h.raw_const_pointer();
    const T* p=dynamic_cast<const T*>(i);
    if(p) { return const_cast<T&>(*p); }
    p=dynamic_cast<const Wrapper<T,I>*>(i);
    if(p) { return *p; }
    const WritableInterface* w=dynamic_cast<const WritableInterface*>(i);
    write_error(std::cerr,"dynamic_handle_extract",w,typeid(*i).name(),typeid(i).name(),typeid(T).name());
    throw std::bad_cast();
}

template<class T, class I> const T& dynamic_handle_extract(const Handle<I>& h) {
    return dynamic_handle_extract<T>(Handle<const I>(h));
}

template<class T, class I> T& dynamic_handle_extract(Handle<I>& h) {
    const I* i=h.raw_const_pointer();
    const T* p=dynamic_cast<const T*>(i);
    if(p) { return const_cast<T&>(*p); }
    p=dynamic_cast<const Wrapper<T,I>*>(i);
    if(p) { return const_cast<T&>(*p); }
    const WritableInterface* w=dynamic_cast<const WritableInterface*>(i);
    write_error(std::cerr,"dynamic_handle_extract",w,typeid(*i).name(),typeid(i).name(),typeid(T).name());
    throw std::bad_cast();
}



template<class T, class I> SharedPointer<const T> dynamic_pointer_extract(const SharedPointer<const I>& ip) {
    SharedPointer<const T> tp=dynamic_pointer_cast<const T>(ip);
    if(tp) { return tp; }
    tp=static_pointer_cast<const T>(dynamic_pointer_cast<const Wrapper<T,I>>(ip));
    if(tp) { return tp; }
    SharedPointer<const WritableInterface> wp=dynamic_pointer_cast<const WritableInterface>(ip);
    return tp;
}


} // namespace Ariadne

#endif
