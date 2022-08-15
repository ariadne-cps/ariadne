/***************************************************************************
 *            chebyshev_model.tpl.hpp
 *
 *  Copyright 2008-18  Pieter Collins
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
 *  GNU Library General Public License for more detai1ls.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "numeric/numeric.hpp"

#include <iomanip>
#include <limits>

#include "numeric/rounding.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/expansion.hpp"
#include "algebra/series.hpp"
#include "algebra/differential.hpp"
#include "function/chebyshev_model.hpp"
#include "function/taylor_series.hpp"
#include "function/function.hpp"
#include "utility/exceptions.hpp"

#include "algebra/expansion.inl.hpp"
#include "algebra/expansion.tpl.hpp"
#include "algebra/evaluate.tpl.hpp"
#include "algebra/algebra_operations.tpl.hpp"

#include "algebra/multi_index-noaliasing.hpp"
#include "function/function_mixin.hpp"
#include "algebra/vector.hpp"

namespace Ariadne {

#warning Move to scaled_function_patch.tpl.hpp with version from taylor_model.tpl.hpp, or remove
inline Interval<FloatDP> convert_interval(Interval<FloatDP> const& ivl, DoublePrecision pr) {
    return ivl; }
inline Interval<FloatMP> convert_interval(Interval<FloatDP> const& ivl, MultiplePrecision pr) {
    return Interval<FloatMP>(FloatMP(ivl.lower_bound().raw(),pr),FloatMP(ivl.upper_bound().raw(),pr)); }

#warning Class should not be instantiated
template<> auto Error<Approximation<FloatDP>>::operator=(Positive<ValidatedUpperNumber> y) -> Error<Approximation<FloatDP>>& {
    ARIADNE_NOT_IMPLEMENTED; }
template<> auto Error<Approximation<FloatMP>>::operator=(Positive<ValidatedUpperNumber> y) -> Error<Approximation<FloatMP>>& {
    ARIADNE_NOT_IMPLEMENTED; }


template<class X, class XE> auto
MultivariateChebyshevModel<X,XE>::_factory() const -> FunctionModelFactoryInterface<P,PR,PRE>* {
    ARIADNE_NOT_IMPLEMENTED;
}


template<class X, class XE> auto
MultivariateChebyshevModel<X,XE>::_create() const -> FunctionModelInterface<P,SIG,PR,PRE>* {
    ARIADNE_NOT_IMPLEMENTED;
}


template<class X, class XE> auto
MultivariateChebyshevModel<X,XE>::result_size() const -> SizeOne {
    return SizeOne();
}

template<class X, class XE> auto
MultivariateChebyshevModel<X,XE>::domain() const -> DomainType const {
    return DomainType(this->_cp.argument_size());
}

template<class X, class XE> auto
MultivariateChebyshevModel<X,XE>::codomain() const -> CodomainType const {
    return cast_exact(IntervalRangeType(this->range()));
}

template<class X, class XE> auto
MultivariateChebyshevModel<X,XE>::clobber() -> Void {
    this->_e=0u;
}

template<class X, class XE> auto
MultivariateChebyshevModel<X,XE>::range() const -> RangeType const {
    ARIADNE_NOT_IMPLEMENTED;
}



} //namespace Ariadne


