/***************************************************************************
 *            symbolic/expression_patch.inl.hpp
 *
 *  Copyright  2024  Pieter Collins
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

#ifndef ARIADNE_EXPRESSION_PATCH_INL_HPP
#define ARIADNE_EXPRESSION_PATCH_INL_HPP

#include "symbolic/expression_patch.hpp"

#include "symbolic/variable.hpp"
#include "symbolic/expression.hpp"
#include "symbolic/expression_set.hpp"


namespace Ariadne {


template<class P,class T1, class T2>
Tuple<Vector<RealVariable>,MultivariateFunctionPatch<P,T1>,MultivariateFunctionPatch<P,T2>>
make_common_variables(RestrictedExpression<P,T1> ep1, RestrictedExpression<P,T2> ep2);



template<class OP> auto RestrictedExpression<ValidatedTag,RealScalar>::_apply(OP op, RestrictedExpression<P,T> const& ep1, RestrictedExpression<P,T> const& ep2) -> RestrictedExpression<P,T> {
    auto sff = make_common_variables(ep1,ep2);
    return RestrictedExpression<P,T>(std::get<0>(sff), op(std::get<1>(sff),std::get<2>(sff)));
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealScalar>::_apply(OP op, RestrictedExpression<P,T> const& ep) -> RestrictedExpression<P,T> {
    return RestrictedExpression<P,T>(ep._vars,op(ep._fp));
}



template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, RestrictedExpression<P,VT> ep1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
    auto sff = make_common_variables(ep1,ep2);
    return RestrictedExpression<P,VT>(std::get<0>(sff), op(std::get<1>(sff),std::get<2>(sff)));
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, RestrictedExpression<P,VT> ep1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,VT> {
    auto sff = make_common_variables(ep1,ep2);
    return RestrictedExpression<P,VT>(std::get<0>(sff), op(std::get<1>(sff),std::get<2>(sff)));
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, RestrictedExpression<P,ST> ep1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
    auto sff = make_common_variables(ep1,ep2);
    return RestrictedExpression<P,VT>(std::get<0>(sff), op(std::get<1>(sff),std::get<2>(sff)));
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, RestrictedExpression<P,VT> ep1, Vector<Number<P>> c2) -> RestrictedExpression<P,VT> {
    return RestrictedExpression<P,VT>(ep1._vars, op(ep1._fp,c2));
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, RestrictedExpression<P,ST> ep1, Vector<Number<P>> c2) -> RestrictedExpression<P,VT> {
    return RestrictedExpression<P,VT>(ep1._vars, op(ep1._fp,c2));
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, RestrictedExpression<P,VT> ep1, Scalar<Number<P>> c2) -> RestrictedExpression<P,VT> {
    return RestrictedExpression<P,VT>(ep1._vars, op(ep1._fp,c2));
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, Vector<Number<P>> c1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
    return RestrictedExpression<P,VT>(ep2._vars, op(c1,ep2._fp));
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, Scalar<Number<P>> c1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
    return RestrictedExpression<P,VT>(ep2._vars, op(c1,ep2._fp));
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, Vector<Number<P>> c1, RestrictedExpression<P,ST> ep2) -> RestrictedExpression<P,VT> {
    return RestrictedExpression<P,VT>(ep2._vars, op(c1,ep2._fp));
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, RestrictedExpression<P,VT> ep1, Expression<VT> e2) -> RestrictedExpression<P,VT> {
    RestrictedExpression<ValidatedTag,RealVector> ep2(ep1.domains(),Vector<Expression<Real>>(e2));
    return _apply(op,ep1,ep2);
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, RestrictedExpression<P,VT> ep1, Expression<ST> e2) -> RestrictedExpression<P,VT> {
    RestrictedExpression<ValidatedTag,RealScalar> ep2(ep1.domains(),Scalar<Expression<Real>>(e2));
    return _apply(op,ep1,ep2);
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, Expression<VT> e1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
    RestrictedExpression<ValidatedTag,RealVector> ep1(ep2.domains(),Vector<Expression<Real>>(e1));
    return _apply(op,ep1,ep2);
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, Expression<ST> e1, RestrictedExpression<P,VT> ep2) -> RestrictedExpression<P,VT> {
    RestrictedExpression<ValidatedTag,RealScalar> ep1(ep2.domains(),Scalar<Expression<Real>>(e1));
    return _apply(op,ep1,ep2);
}
template<class OP> auto RestrictedExpression<ValidatedTag,RealVector>::_apply(OP op, RestrictedExpression<P,VT> ep) -> RestrictedExpression<P,VT> {
    return RestrictedExpression<P,VT>(ep._vars, op(ep._fp));
}


} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_PATCH_INL_HPP */
