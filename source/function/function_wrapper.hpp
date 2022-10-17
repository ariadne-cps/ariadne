/***************************************************************************
 *            function_wrapper.hpp
 *
 *  Copyright  2022  Pieter Collins
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

#ifndef ARIADNE_FUNCTION_WRAPPER_HPP
#define ARIADNE_FUNCTION_WRAPPER_HPP

#include "function_concepts.hpp"
#include "function_mixin.hpp"

namespace Ariadne {

template<class F, class P, class SIG> requires AFunction<F,P,SIG> class FunctionWrapper;

template<class P, class SIG, class F> requires AFunction<F,P,SIG>
FunctionWrapper<F,P,SIG>* make_function_wrapper(F f) { return new FunctionWrapper<F,P,SIG>(f); }

template<class F, class P, class SIG> requires AFunction<F,P,SIG>
class FunctionWrapper
    : public FunctionMixin<FunctionWrapper<F,P,SIG>,P,SIG>
{
    using typename FunctionInterface<P,SIG>::ArgumentSizeType;
    using typename FunctionInterface<P,SIG>::ResultSizeType;
    using typename FunctionInterface<P,SIG>::ArgumentIndexType;
    template<class X> using Argument = typename SignatureTraits<SIG>::template Argument<X>;
    template<class X> using Result = typename SignatureTraits<SIG>::template Result<X>;
  protected:
    F _f;
  public:
    FunctionWrapper(F f) : _f(f) { }
    operator F const& () const { return this->_f; }
    virtual ArgumentSizeType argument_size() const override { return this->_f.argument_size(); }
    virtual ResultSizeType result_size() const override { return this->_f.result_size(); }
    virtual FunctionInterface<P,SIG>* _derivative(ArgumentIndexType k) const override {
        return make_function_wrapper<P,SIG>(derivative(this->_f,k)); }
    template<class T> Result<T> operator() (Argument<T> const& x) const { return this->_f(x); }
};

template<class P, class SIG> template<class F> requires (not IsFunctionClass<F,SIG>) and AFunction<F,P,SIG>
Function<P,SIG>::Function(F const& f)
    : Handle<const Interface>(make_function_wrapper<F,P>(f)) { }

} // namespace Ariadne

#endif // ARIADNE_FUNCTION_WRAPPER_HPP
