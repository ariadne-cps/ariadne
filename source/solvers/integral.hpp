/***************************************************************************
 *            integral.hpp
 *
 *  Copyright  2019  Pieter Collins
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

#ifndef ARIADNE_INTEGRAL_HPP
#define ARIADNE_INTEGRAL_HPP

namespace Ariadne {

class IntegralInterface {
  public:
    ~IntegralInterface() = default;
    virtual ValidatedNumber integral(ValidatedScalarUnivariateFunction const& f, ValidatedNumber const& a, ValidatedNumber const& b) const = 0;
};

} // namespace Ariadne


#include "function/taylor_function.hpp"

namespace Ariadne {

template<class F> using ValidatedScalarTaylorFunctionModel = ScaledFunctionPatch<ValidatedTaylorModel<F>>;

template<class F> class TaylorIntegral : public IntegralInterface {
    Sweeper<F> _swp;
  public:
    TaylorIntegral(Sweeper<F> swp) : _swp(swp) { }
    virtual ValidatedNumber integral(ValidatedScalarUnivariateFunction const& f, ValidatedNumber const& a, ValidatedNumber const& b) const override;
    ValidatedNumber integral(ValidatedScalarMultivariateFunction const& f, ValidatedNumber const& a, ValidatedNumber const& b) const;
    Bounds<F> integral(ValidatedScalarTaylorFunctionModel<F> const& tf, Bounds<F> const& a, Bounds<F> const& b) const;
};

template<class F> ValidatedNumber TaylorIntegral<F>::integral(ValidatedScalarUnivariateFunction const& f, ValidatedNumber const& a, ValidatedNumber const& b)  const {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class F> ValidatedNumber TaylorIntegral<F>::integral(ValidatedScalarMultivariateFunction const& vf, ValidatedNumber const& a, ValidatedNumber const& b) const {
    auto pr=this->_swp.properties().precision();
    Bounds<F> ax(a,pr); Bounds<F> bx(b,pr);
    IntervalDomainType dom(cast_exact(ax.lower()),cast_exact(bx.upper()));
    BoxDomainType bxdom={dom};
    ValidatedScalarTaylorFunctionModel<F> tf(bxdom,vf,this->_swp);
    return this->integral(tf,ax,bx);
}

template<class F> Bounds<F> TaylorIntegral<F>::integral(ValidatedScalarTaylorFunctionModel<F> const& tf, Bounds<F> const& a, Bounds<F> const& b) const {
    ARIADNE_ASSERT(tf.argument_size()==1);
    auto pr=tf.properties().precision();
    Vector<Bounds<F>> va={a};
    Vector<Bounds<F>> vb={b};
    auto integral_tf = antiderivative(tf, 0);
    return integral_tf(vb)-integral_tf(va);
}

} // namespace Ariadne

#endif
