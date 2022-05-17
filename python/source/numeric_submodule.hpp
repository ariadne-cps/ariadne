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

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Dyadic>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Decimal>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Rational>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Real>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDP>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatMP>& repr);

template<class P> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Number<P>>& repr);

template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Approximation<F>>& x);
template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<LowerBound<F>>& repr);
template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<UpperBound<F>>& repr);
template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Bounds<F>>& x);
template<class F, class FE> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Ball<F,FE>>& x);
template<class F> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Value<F>>& x);
template<class FE> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Error<FE>>& repr);

OutputStream& operator<<(OutputStream& os, const PythonLiteral<FloatDP>& x);
OutputStream& operator<<(OutputStream& os, const PythonLiteral<FloatMP>& repr);
template<class F> OutputStream& operator<<(OutputStream& os, const PythonLiteral<Approximation<F>>& x);
template<class F> OutputStream& operator<<(OutputStream& os, const PythonLiteral<LowerBound<F>>& repr);
template<class F> OutputStream& operator<<(OutputStream& os, const PythonLiteral<UpperBound<F>>& repr);
template<class F> OutputStream& operator<<(OutputStream& os, const PythonLiteral<Bounds<F>>& x);
template<class F, class FE> OutputStream& operator<<(OutputStream& os, const PythonLiteral<Ball<F,FE>>& x);
template<class F> OutputStream& operator<<(OutputStream& os, const PythonLiteral<Value<F>>& x);
template<class FE> OutputStream& operator<<(OutputStream& os, const PythonLiteral<Error<FE>>& repr);


template<class X> X from_python_object_or_literal(pybind11::handle h) {
    if constexpr (Constructible<X,std::string>) {
        try { return X(pybind11::cast<std::string>(h)); } 
        catch(...) { } 
    }
    return pybind11::cast<X>(h);
}
template<> ValidatedNumber from_python_object_or_literal<ValidatedNumber>(pybind11::handle h);
template<> ValidatedUpperNumber from_python_object_or_literal<ValidatedUpperNumber>(pybind11::handle h);
template<> ValidatedLowerNumber from_python_object_or_literal<ValidatedLowerNumber>(pybind11::handle h);
template<> ApproximateNumber from_python_object_or_literal<ApproximateNumber>(pybind11::handle h);

} // namespace Ariadne

#endif /* ARIADNE_PYTHON_NUMERIC_SUBMODULE_HPP */
