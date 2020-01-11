/***************************************************************************
 *            pybind11.hpp
 *
 *  Copyright  2014-20  Pieter Collins
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
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <string>
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include <iostream>

#include "config.hpp"


namespace pybind11 {

// Indicate bases of a pybind11 class.
template<class... BS> class bases { };
template<class T, class... BS> class class_<T,bases<BS...>> : public class_<T,BS...> {
  public:
    using class_<T,BS...>::class_;
};

// A pybind11 object convertible to C++. Useful for the result of an overloaded function.
class __attribute__ ((visibility ("default"))) result_object {
  public:
    explicit result_object() : _obj() { }
    explicit result_object(pybind11::object obj) : _obj(obj) { }
  public:
    template <class T> operator T() { return this->convert_to<T>(); }
    template <class T> operator T&() const { return const_cast<result_object*>(this)->convert_to<T&>(); }
  private:
    template<class T> T convert_to() { 
        if (pybind11::detail::cast_is_temporary_value_reference<T>::value) { 
            static detail::overload_caster_t<T> caster; 
            return detail::cast_ref<T>(std::move(_obj), caster); 
        } else {
            return detail::cast_safe<T>(std::move(_obj));
        }
    }
    pybind11::object _obj;
};

// A pybind11 function which returns an object convertible to a C++ type.
struct __attribute__ ((visibility ("default"))) override_function {
    override_function(function func) : _func(func) { }
    template<class... ARGS> result_object operator() (ARGS& ... args) { return result_object(_func(args...)); }
  private:
    function _func;
};


// Define get_override for classes wrapping interface I
template<class I> class wrapper : public I {
    const char* _pyclass_name;
  protected:
    wrapper(const char* pyclass_name) : _pyclass_name(pyclass_name) { } 
    override_function get_override(const char* pymethod) const { 
       {
            gil_scoped_acquire gil; 
            function overload = get_overload(static_cast<const I*>(this), pymethod); 
            if (overload) { 
                return override_function(overload);
            } else {
                std::string error_msg=std::string("Tried to call pure virtual function \"")+_pyclass_name+"::"+pymethod+"\"";
                pybind11_fail(error_msg.c_str()); 
            }
        }
        // FIXME: For non-abstract method, should return I::method value
    }
};

} // namespace pybind11

#endif /* ARIADNE_PYBIND11_HPP */
