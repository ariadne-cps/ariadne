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


template<class T, class I> class Wrapper;

template<class FLT, class P, class SIG> requires AFunction<FLT,P,SIG> class FunctionWrapper;

template<class FLT, class P, class SIG> struct Wrapper<FLT,FunctionInterface<P,SIG>> : public FunctionWrapper<FLT,P,SIG> { virtual ~Wrapper() = default; using FunctionWrapper<FLT,P,SIG>::FunctionWrapper; };

template<class P, class SIG, class FLT> requires AFunction<FLT,P,SIG> FunctionWrapper<FLT,P,SIG>* make_function_wrapper(FLT f);


template<class FLT, class P, class SIG> class FunctionWrapperGetterMixin;

template<class FLT, class P, class... ARGS> class FunctionWrapperGetterMixin<FLT,P,RealScalar(ARGS...)> {
};

template<class FLT, class P, class... ARGS> class FunctionWrapperGetterMixin<FLT,P,RealVector(ARGS...)>
    : public virtual VectorOfFunctionInterface<P,ARGS...>
{
    virtual ScalarFunctionInterface<P,ARGS...>* _get(SizeType i) const override;
};



template<class FLT, class P, class SIG> requires AFunction<FLT,P,SIG>
class FunctionWrapper
    : public FunctionMixin<FunctionWrapper<FLT,P,SIG>,P,SIG>
    , public FunctionWrapperGetterMixin<FunctionWrapper<FLT,P,SIG>,P,SIG>
    , public FLT
{
    using typename FunctionInterface<P,SIG>::ArgumentSizeType;
    using typename FunctionInterface<P,SIG>::ResultSizeType;
    using typename FunctionInterface<P,SIG>::ArgumentIndexType;
    template<class X> using Argument = typename SignatureTraits<SIG>::template Argument<X>;
    template<class X> using Result = typename SignatureTraits<SIG>::template Result<X>;
  public:
    virtual ~FunctionWrapper() = default;
    FunctionWrapper(FLT f) : FLT(f) { }
    operator FLT const& () const { return *this; }
    FLT const& wrapped() const { return *this; }
    virtual ArgumentSizeType argument_size() const override { return this->wrapped().argument_size(); }
    virtual ResultSizeType result_size() const override { return this->wrapped().result_size(); }
    virtual FunctionInterface<P,SIG>* _derivative(ArgumentIndexType k) const override {
        return make_function_wrapper<P,SIG>(derivative(this->wrapped(),k)); }
    template<class T> Result<T> operator() (Argument<T> const& x) const { return this->wrapped()(x); }
    friend inline OutputStream& operator<<(OutputStream& os, FunctionWrapper<FLT,P,SIG> const& f) { return os << f.wrapped(); }
        template<class I> decltype(auto) operator[](I i) const { return this->wrapped()[i]; }
};

template<class P, class SIG> template<class FLT> requires (not IsFunctionClass<FLT,SIG>) and AFunction<FLT,P,SIG>
Function<P,SIG>::Function(FLT const& f)
    : Handle<const Interface>(make_function_wrapper<P,SIG>(f)) { }

template<class P, class SIG, class FLT> requires AFunction<FLT,P,SIG>
FunctionWrapper<FLT,P,SIG>* make_function_wrapper(FLT f) { return new Wrapper<FLT,FunctionInterface<P,SIG>>(f); }

template<class FLT, class P, class SIG> FLT const& extract(Function<P,SIG> const& f) {
    if constexpr(DerivedFrom<FLT,FunctionInterface<P,SIG>>) {
        return dynamic_cast<FLT const&>(f.reference());
    } else {
        return dynamic_cast<FunctionWrapper<FLT,P,SIG> const&>(f.reference());
    }
}

} // namespace Ariadne

#endif // ARIADNE_FUNCTION_WRAPPER_HPP
