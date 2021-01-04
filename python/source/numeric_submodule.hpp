/***************************************************************************
 *            numeric_submodule.hpp
 *
 *  Copyright  2021  Pieter Collins
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

#ifndef ARIADNE_PYTHON_NUMERIC_SUBMODULE_HPP
#define ARIADNE_PYTHON_NUMERIC_SUBMODULE_HPP

#include "pybind11.hpp"
#include <pybind11/operators.h>

#include "utilities.hpp"

namespace Ariadne {

template<class P> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Number<P>>& repr) {
    return os << class_name<Number<P>>()<<"("<<repr.reference()<<")"; }

template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Approximation<F>>& repr) {
    return os << class_name<F>() << "Approximation("<<repr.reference().raw()<<")"; }
template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<LowerBound<F>>& repr) {
    return os << class_name<F>() << "LowerBound("<<repr.reference().raw()<<")"; }
template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<UpperBound<F>>& repr) {
    return os << class_name<F>() << "UpperBound("<<repr.reference().raw()<<")"; }
template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Bounds<F>>& repr) {
    return os << class_name<F>() << "Bounds("<<repr.reference().lower().raw()<<","<<repr.reference().upper().raw()<<")"; }
template<class F, class FE> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Ball<F,FE>>& repr) {
    return os << class_name<F>() << "Ball("<<repr.reference().value_raw()<<","<<repr.reference().error_raw()<<")"; }
template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Value<F>>& repr) {
    return os << class_name<F>() << "Value("<<repr.reference().raw()<<")"; }
template<class FE> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Error<FE>>& repr) {
    return os << class_name<FE>() << "Error("<<repr.reference().raw()<<")"; }

}

#endif /* ARIADNE_PYTHON_NUMERIC_SUBMODULE_HPP */
