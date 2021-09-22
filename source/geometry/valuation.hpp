/***************************************************************************
 *            valuation.hpp
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
 *  MERCHANTABILITY or FITNEVV FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file valuation.hpp
 *  \brief Valuations, which are measures taking values only on open sets.
 */

#ifndef ARIADNE_VALUATION_HPP
#define ARIADNE_VALUATION_HPP

#include "utility/handle.hpp"
#include "numeric/number.hpp"

#include "numeric/sequence.hpp"
#include "numeric/integer.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/twoexp.hpp"

#include "set_interface.hpp"
#include "interval.hpp"

namespace Ariadne {

using Void = void;

template<class P, class T> class OpenSet;

template<class P, class T> class Valuation;
template<class T> using EffectiveValuation = Valuation<EffectiveTag,T>;
template<class T> using ValidatedValuation = Valuation<ValidatedTag,T>;

template<class V, class P, class T> struct IsValuation {
    using LY=PositiveLowerNumber<P>; using OPS=typename SetTraits<T>::OpenSetType;
    template<class VV,class=decltype(declval<VV>().measure(declval<OPS>()))>
        static std::true_type test(int);
    template<class VV>
        static std::false_type test(...);
    static const bool value = decltype(test<V>(1))::value;
};

template<class V, class P, class T> concept AValuation = IsValuation<V,P,T>::value;

template<class P, class T> class ValuationInterface {
  public:
    typedef typename SetTraits<T>::OpenSetType OpenSetType;
    virtual ~ValuationInterface() = default;
    virtual PositiveLowerNumber<P> _measure(OpenSetType const&) const = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

template<class P, class T> class Valuation : public Handle<const ValuationInterface<P,T>> {
  public:
    using typename Handle<const ValuationInterface<P,T>>::Interface;
    typedef typename SetTraits<T>::OpenSetType OpenSetType;
    using Handle<Interface>::Handle;
    template<AValuation<P,T> V> Valuation(V const& s);
    PositiveLowerNumber<P> measure(OpenSetType const& ops) const { return this->reference()._measure(ops); }
    friend OutputStream& operator<<(OutputStream& os, Valuation<P,T> const& val) { return val.reference()._write(os); }
};


template<class V, class P, class T> class ValuationWrapper;

template<class V, class T> class ValuationWrapper<V,ValidatedTag,T>
    : public virtual Valuation<ValidatedTag,T>::Interface
{
    using P=ValidatedTag;
    V const& base() const { return this->_v; }
  private:
    V _v;
  public:
    typedef typename SetTraits<T>::OpenSetType OpenSetType;
    ValuationWrapper(V const& val) : _v(val) { }
    virtual PositiveLowerNumber<P> _measure(OpenSetType const& ops) const final override {
        return PositiveLowerNumber<P>(this->base().measure(ops)); }
    OutputStream& _write(OutputStream& os) const final override {
        return os << this->base(); }
};

template<class V> decltype(auto) wrap_valuation(V const& v) {
    using P=ValidatedTag; using T=Real;
    return Valuation<P,T>(std::make_shared<ValuationWrapper<V,P,T>>(v));
}

template<class P, class T> template<AValuation<P,T> V>
Valuation<P,T>::Valuation(V const& v)
    : Valuation(std::make_shared<ValuationWrapper<V,P,T>>(v)) {
}

} // namespace Ariadne

#endif
