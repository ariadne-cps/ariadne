/***************************************************************************
 *            pybind11.hpp
 *
 *  Copyright  2014  Pieter Collins
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

/*! \file pybind11.hpp
 *  PyBind11 headers for precompiling.
 */

#ifndef ARIADNE_PYBIND11_HPP
#define ARIADNE_PYBIND11_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include <string>
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include <iostream>

#include "config.hpp"


namespace pybind11 {

template<class... BS> class bases { };
template<class T, class... BS> class class_<T,bases<BS...>> : public class_<T,BS...> {
  public:
    using class_<T,BS...>::class_;
};

class noncopyable { };
template<class T> class class_<T,noncopyable> : public class_<T> {
  public:
    using class_<T>::class_;
};

// The result of calling a method.
class method_result {
  public:
    explicit method_result() { } // : m_obj(x) { }
    explicit method_result(PyObject* x); // : m_obj(x) { }
  public:
    template <class T> operator T() { assert(false); }
        // { converter::return_from_python<T> converter; return converter(m_obj.release()); }
    template <class T> operator T&() const { assert(false); }
        // { converter::return_from_python<T&> converter; return converter(const_cast<handle<>&>(m_obj).release()); }
    template <class T> T as();
        // { converter::return_from_python<T> converter; return converter(m_obj.release()); }
    template <class T> T unchecked();
        // { return extract<T>(m_obj.get())(); }
private:
//    mutable handle<> m_obj;
};

class override_helper { // : public object {
 private:
//    override_helper(handle x) : object(x) { }
 public:
    method_result operator()() const;
        // { method_result x(PyEval_CallFunction(this->ptr(), const_cast<char*>("()"))); return x; }
};

template<class T> class wrapper {
  protected:
    //override_helper get_override(const char* str) const;
    method_result get_override(const char* str) const { return method_result(); }
};

} // namespace pybind11

#endif /* ARIADNE_PYBIND11_HPP */
