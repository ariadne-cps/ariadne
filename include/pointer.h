/***************************************************************************
 *            pointer.h
 *
 *  Copyright 2008  Pieter Collins
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

/*! \file pointer.h
 *  \brief Smart pointers in Euclidean space.
 */

#ifndef ARIADNE_POINTER_H
#define ARIADNE_POINTER_H

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>


namespace Ariadne {


using boost::scoped_ptr;
using boost::shared_ptr;

template<class T>
class copy_on_write_ptr
{
  private:
    mutable int* _ref_count;
    mutable T* _ptr;
  public:
    ~copy_on_write_ptr() { this->_deallocate(); }
    template<class S> copy_on_write_ptr(S* pointer)
        : _ref_count(new int(1)), _ptr(pointer) { }
    copy_on_write_ptr(const copy_on_write_ptr<T>& other)
        : _ref_count(other._ref_count), _ptr(other._ptr) { ++(*_ref_count); }
    template<class S> copy_on_write_ptr<T>& operator=(S* pointer) {
        _deallocate(); _allocate(pointer); return *this; }
    copy_on_write_ptr<T>& operator=(const copy_on_write_ptr<T>& other) {
        if(_ptr!=other._ptr) { _deallocate(); _ptr=other._ptr;
            _ref_count=other._ref_count; ++(*_ref_count); } return *this; }
    const T& operator*() const { return *_ptr; }
    const T* operator->() const { return _ptr; }
    T& operator*() { if(*_ref_count>1) { _make_copy(); } return *_ptr; }
    T* operator->() { if(*_ref_count>1) { _make_copy(); } return _ptr; }
  public:
    int count() const { return *_ref_count; }
    const void* address() const { return _ptr; }
    const T& value() const { return *_ptr; }
  private:
    void _make_copy() {
        --(*_ref_count); _ref_count=new int(1); _ptr=new T(*_ptr); }
    void _deallocate() {
        if(*_ref_count==1) { delete _ptr; _ptr=0; delete _ref_count; }
        else { --(*_ref_count); } }
    template<class S> void _allocate(S* pointer) {
        _ref_count=new int(1); _ptr=pointer; }
};

template<class T> std::ostream&
operator<<(std::ostream& os, const copy_on_write_ptr<T>& ptr) {
    return os << "["<<ptr.count()<<"]<"<<ptr.address()<<">("<<ptr.value()<<")";
}

template<class T>
class clone_on_write_ptr
{
  private:
    mutable int* _ref_count;
    mutable T* _ptr;
  public:
    ~clone_on_write_ptr() { this->_deallocate(); }
    template<class S> clone_on_write_ptr(S* pointer)
        : _ref_count(new int(1)), _ptr(pointer) { }
    clone_on_write_ptr(const clone_on_write_ptr<T>& other)
        : _ref_count(other._ref_count), _ptr(other._ptr) { ++(*_ref_count); }
    template<class S> clone_on_write_ptr<T>& operator=(S* pointer) {
        _deallocate(); _allocate(pointer); return *this; }
    clone_on_write_ptr<T>& operator=(const clone_on_write_ptr<T>& other) {
        if(_ptr!=other._ptr) { _deallocate(); _ptr=other._ptr;
            _ref_count=other._ref_count; ++(*_ref_count); } return *this; }
    const T& operator*() const { return *_ptr; }
    const T* operator->() const { return _ptr; }
    T& operator*() { if(*_ref_count>1) { _make_clone(); } return *_ptr; }
    T* operator->() { if(*_ref_count>1) { _make_clone(); } return _ptr; }
  public:
    int count() const { return *_ref_count; }
    const void* address() const { return _ptr; }
    const T& value() const { return *_ptr; }
  private:
    void _make_clone() {
        --(*_ref_count); _ref_count=new int(1); _ptr=_ptr->clone(); }
    void _deallocate() {
        if(*_ref_count==1) { delete _ptr; _ptr=0; delete _ref_count; }
        else { --(*_ref_count); } }
    template<class S> void _allocate(S* pointer) {
        _ref_count=new int(1); _ptr=pointer; }
};

template<class T> std::ostream&
operator<<(std::ostream& os, const clone_on_write_ptr<T>& ptr) {
    return os << "["<<ptr.count()<<"]<"<<ptr.address()<<">("<<ptr.value()<<")";
}

template<class T>
class copy_on_copy_ptr
{
  private:
    mutable T* _ptr;
  public:
    ~copy_on_copy_ptr() { delete _ptr; _ptr=0; }
    template<class S> copy_on_copy_ptr(S* pointer)
        : _ptr(pointer) { }
    copy_on_copy_ptr(const copy_on_copy_ptr<T>& other)
        : _ptr(new T(other._ptr)) { }
    template<class S> copy_on_copy_ptr<T>& operator=(S* pointer) {
        delete _ptr; _ptr=pointer; return *this; }
    copy_on_copy_ptr<T>& operator=(const copy_on_copy_ptr<T>& other) {
        if(_ptr!=other._ptr) { delete _ptr; _ptr=new T(other._ptr); } return *this; }
    const T& operator*() const { return *_ptr; }
    const T* operator->() const { return _ptr; }
    T& operator*() { return *_ptr; }
    T* operator->() { return _ptr; }
  public:
    const void* address() const { return _ptr; }
    const T& value() const { return *_ptr; }
};

template<class T> std::ostream&
operator<<(std::ostream& os, const copy_on_copy_ptr<T>& ptr) {
    return os<<"<"<<ptr.address()<<">("<<ptr.value()<<")";
}

template<class T>
class clone_on_copy_ptr
{
  private:
    mutable T* _ptr;
  public:
    ~clone_on_copy_ptr() { delete _ptr; _ptr=0; }
    template<class S> clone_on_copy_ptr(S* pointer)
        : _ptr(pointer) { }
    clone_on_copy_ptr(const clone_on_copy_ptr<T>& other)
        : _ptr(other._ptr->clone()) { }
    template<class S> clone_on_copy_ptr<T>& operator=(S* pointer) {
        delete _ptr; _ptr=pointer; return *this; }
    copy_on_copy_ptr<T>& operator=(const clone_on_copy_ptr<T>& other) {
        if(_ptr!=other._ptr) { delete _ptr; _ptr=other._ptr->clone(); } return *this; }
    const T& operator*() const { return *_ptr; }
    const T* operator->() const { return _ptr; }
    T& operator*() { return *_ptr; }
    T* operator->() { return _ptr; }
  public:
    const void* address() const { return _ptr; }
    const T& value() const { return *_ptr; }
};

template<class T> std::ostream&
operator<<(std::ostream& os, const clone_on_copy_ptr<T>& ptr) {
    return os<<"<"<<ptr.address()<<">("<<ptr.value()<<")";
}


} // namespace Ariadne

#endif // ARIADNE_POINTER_H
