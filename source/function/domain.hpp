/***************************************************************************
 *            function/domain.hpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
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

/*! \file function/domain.hpp
 *  \brief Interval and box domains for functions.
 */

#ifndef ARIADNE_DOMAIN_HPP
#define ARIADNE_DOMAIN_HPP

#include "geometry/interval.hpp"
#include "geometry/box.hpp"

namespace Ariadne {

//! \ingroup FunctionModule
//! \brief The type used as a bounded domain for a univariate function.
using IntervalDomainType = Interval<FloatDPValue>;
//! \ingroup FunctionModule
//! \brief The type used as a bounded domain for a multivariate function.
using BoxDomainType = Box<Interval<FloatDPValue>>;

template<class... DS> using IntersectionType = decltype(intersection(declval<DS>()...));
template<class... DS> using CartesianProductType = decltype(product(declval<DS>()...));

template<class D> using DimensionalType = decltype(declval<D>().dimension());
template<class D> using IsOneDimensional = IsSame<DimensionalType<D>,SizeOne>;
template<class D> using IsMultiDimensional = IsSame<DimensionalType<D>,SizeType>;

class ScalarDomain;
class VectorDomain;

class UnitInterval;
class UnitBox;

class RealDomain;
class EuclideanDomain;

typedef Interval<FloatDPValue> IntervalDomainType;
typedef Box<Interval<FloatDPValue>> BoxDomainType;

class ScalarDomain {
    Interval<FloatDPValue> _repr;
  public:
    typedef SizeOne DimensionType;
    ScalarDomain(RealDomain const&);
    ScalarDomain(IntervalDomainType const&);
    ScalarDomain(UnitInterval const&);

    explicit ScalarDomain(); // FIXME: Should this be allowed as a constructor for the entire space?
    explicit operator IntervalDomainType () const;

    friend Boolean operator==(ScalarDomain const& dom1, ScalarDomain const& dom2) { return dom1._repr==dom2._repr; }

    static constexpr SizeOne dimension() { return SizeOne(); }
    friend ScalarDomain intersection(ScalarDomain const& dom1, ScalarDomain const& dom2) { return ScalarDomain(intersection(dom1._repr,dom2._repr)); }
    friend Boolean subset(ScalarDomain const& dom1, ScalarDomain const& dom2) { return subset(dom1._repr,dom2._repr); }
    friend Boolean equal(ScalarDomain const& dom1, ScalarDomain const& dom2) { return equal(dom1._repr,dom2._repr); }
    friend VectorDomain product(ScalarDomain const& dom1, ScalarDomain const& dom2);
    friend OutputStream& operator<<(OutputStream& os, ScalarDomain const& dom) { return os << dom._repr; }
  private:
    friend class VectorDomain;
};

class VectorDomain {
    Box<Interval<FloatDPValue>> _repr;
  public:
    typedef SizeType DimensionType;
    VectorDomain(EuclideanDomain const&);
    VectorDomain(BoxDomainType const&);
    VectorDomain(UnitBox const&);

    explicit VectorDomain(SizeType n); // FIXME: Should this be allowed as a constructor for the entire space?
    explicit operator BoxDomainType () const;

    friend Boolean operator==(VectorDomain const& dom1, VectorDomain const& dom2) { return dom1._repr==dom2._repr; }

    SizeType dimension() const { return _repr.dimension(); }
    ScalarDomain operator[] (SizeType i) const { return ScalarDomain(_repr[i]); }
    friend Boolean subset(VectorDomain const& dom1, VectorDomain const& dom2) { return subset(dom1._repr,dom2._repr); }
    friend Boolean equal(VectorDomain const& dom1, VectorDomain const& dom2) { return equal(dom1._repr,dom2._repr); }
    friend VectorDomain intersection(VectorDomain const& dom1, VectorDomain const& dom2) { return VectorDomain(intersection(dom1._repr,dom2._repr)); }
    friend VectorDomain product(ScalarDomain const& dom1, ScalarDomain const& dom2) { return VectorDomain(_product(dom1._repr,dom2._repr)); }
    friend VectorDomain product(ScalarDomain const& dom1, VectorDomain const& dom2) { return VectorDomain(product(dom1._repr,dom2._repr)); }
    friend VectorDomain product(VectorDomain const& dom1, ScalarDomain const& dom2) { return VectorDomain(product(dom1._repr,dom2._repr)); }
    friend VectorDomain product(VectorDomain const& dom1, VectorDomain const& dom2) { return VectorDomain(product(dom1._repr,dom2._repr)); }
    friend OutputStream& operator<<(OutputStream& os, VectorDomain const& dom) { return os << dom._repr; }
  private:
    static BoxDomainType _product(IntervalDomainType const& ivl1,IntervalDomainType const& ivl2) {
        return BoxDomainType({ivl1,ivl2}); }

};


/*

class IntervalDomain : public Interval<FloatDPValue> {
    typedef FloatDPValue L; typedef FloatDPValue U;
    typedef DoublePrecision PR;
  public:
    using Interval<FloatDPValue>::Interval;
    template<class V, EnableIf<IsConstructible<U,V,PR>> =dummy> IntervalDomainType(V const& l, V const& u) : Interval<FloatDPValue>(L(l,dp),U(u,dp)) { }
    template<class V, EnableIf<IsConstructible<U,V,PR>> =dummy> IntervalDomainType(Interval<V> const& ivl) : IntervalDomainType(ivl.lower(),ivl.upper()) { }
};

class BoxDomainType : public Box<Interval<FloatDPValue>> {
  public:
    using Box<Interval<FloatDPValue>>::Box;
    IntervalDomainType operator[] (SizeType i) const { return this->Box<Interval<FloatDPValue>>::operator[](i); }
    BoxDomainType operator[] (Range rng) const { return this->Box<Interval<FloatDPValue>>::operator[](rng); }
};

class BoxDomainType {
    SharedArray<IntervalDomainType> _ivls;
    Box<Interval<FloatDPValue>> const& _box_reference() const { return reinterpret_cast<Box<Interval<FloatDPValue>>const&>(_ivls.reference()); }
  public:
    template<class... ARGS, EnableIf<IsConstructible<Box<IntervalDomainType>,ARGS...>> =dummy> BoxDomainType(ARGS... args);
    operator Box<Interval<FloatDPValue>> const& () { return this->_box_reference(); }
    //operator Box<Interval<FloatDPValue>> () { return this->_box_reference(); }
    DimensionType dimension() const { return this->_box_reference().dimension(); }
    IntervalDomainType const& operator[](SizeType i) const { return _ivls[i]; }
    friend OutputStream& operator<<(OutputStream& os, BoxDomainType const&);
};

using IntervalDomainType = IntervalDomainType;
using BoxDomainType = BoxDomainType;
*/

//! \ingroup FunctionModule \brief The domain of an entire univariate function.
class RealDomain {
  public:
    typedef SizeOne DimensionType; //!< .
    typedef IndexZero IndexType; //!< .

    constexpr RealDomain() { } //!< .
    constexpr RealDomain(SizeOne) { } //!< .
    constexpr SizeOne dimension() const { return SizeOne(); } //!< .
    operator IntervalDomainType() const { return IntervalDomainType(-inf,+inf); } //!< .
    friend EuclideanDomain product(RealDomain const& dom1, RealDomain const& dom2); //!< .
    friend RealDomain intersection(RealDomain const& dom1, RealDomain const& dom2) { return RealDomain(); } //!< .
    friend Bool operator==(RealDomain const& dom1, RealDomain const& dom2) { return true; } //!< .
    friend Bool operator==(RealDomain const& dom1, IntervalDomainType const& dom2) { return IntervalDomainType(dom1)==dom2; } //!< .
    friend Bool operator==(IntervalDomainType const& dom1, RealDomain const& dom2) { return dom1==IntervalDomainType(dom2); } //!< .
    friend OutputStream& operator<<(OutputStream& os, RealDomain const& dom) { return os << "R"; } //!< .
};

//! \ingroup FunctionModule \brief The domain of an entire multivariate function.
class EuclideanDomain {
    SizeType _dim;
  public:
<<<<<<< HEAD
    typedef SizeType DimensionType; //!< .
    typedef SizeType IndexType; //!< .

    constexpr EuclideanDomain(SizeType dim) : _dim(dim) { } //!< .
    constexpr EuclideanDomain(SizeType dim, RealDomain) : _dim(dim) { } //!< .
    constexpr SizeType dimension() const { return this->_dim; } //!< .
    constexpr RealDomain operator[](SizeType ind) { return RealDomain(); } //!< .
    operator BoxDomainType() const { return BoxDomainType(this->dimension(),IntervalDomainType(RealDomain())); } //!< .
    friend EuclideanDomain intersection(EuclideanDomain const& dom1, EuclideanDomain const& dom2) { assert(dom1==dom2); return dom1; } //!< .
    friend EuclideanDomain product(RealDomain const& dom1, RealDomain const& dom2) { return EuclideanDomain(2u); } //!< .
    friend EuclideanDomain product(RealDomain const& dom1, EuclideanDomain const& dom2) { return EuclideanDomain(1u+dom2.dimension()); } //!< .
    friend EuclideanDomain product(EuclideanDomain const& dom1, RealDomain const& dom2) { return EuclideanDomain(dom1.dimension()+1u); } //!< .
    friend EuclideanDomain product(EuclideanDomain const& dom1, EuclideanDomain const& dom2) { return EuclideanDomain(dom1.dimension()+dom2.dimension()); } //!< .
    friend Bool operator==(EuclideanDomain const& dom1, EuclideanDomain const& dom2) { return dom1.dimension() == dom2.dimension(); } //!< .
    friend Bool operator==(EuclideanDomain const& dom1, BoxDomainType const& dom2) { return BoxDomainType(dom1)==dom2; } //!< .
    friend Bool operator==(BoxDomainType const& dom1, EuclideanDomain const& dom2) { return dom1==BoxDomainType(dom2); } //!< .
    friend OutputStream& operator<<(OutputStream& os, EuclideanDomain const& dom) { return os << "R" << dom.dimension(); } //!< .
};

//! \ingroup FunctionModule \brief The signed unit interval \f$[-1:+1]\f$.
class UnitInterval {
  public:
    typedef SizeOne DimensionType; //!< .
    constexpr UnitInterval() { } //!< .
    constexpr SizeOne dimension() const { return SizeOne(); } //!< .
    operator IntervalDomainType() const { return IntervalDomainType(-1,+1); } //!< .
    friend UnitInterval intersection(UnitInterval const& dom1, UnitInterval const& dom2) { return UnitInterval(); } //!< .
    friend Bool operator==(UnitInterval const& dom1, UnitInterval const& dom2) { return true; } //!< .
    friend OutputStream& operator<<(OutputStream& os, UnitInterval const& dom) { return os << "[-1:+1]"; } //!< .
};

//! \ingroup FunctionModule \brief The signed unit interval \f$[-1:+1]^n\f$.
class UnitBox {
    SizeType _dim;
  public:
    typedef SizeType DimensionType; //!> .
    constexpr UnitBox(SizeType dim) : _dim(dim) { } //!< .
    constexpr UnitBox(SizeType dim, UnitInterval) : _dim(dim) { } //!< .
    constexpr SizeType dimension() const { return this->_dim; } //!< .
    constexpr UnitInterval operator[](SizeType ind) { return UnitInterval(); } //!< .
    operator BoxDomainType() const { return BoxDomainType(this->dimension(),IntervalDomainType(UnitInterval())); } //!< .
    friend UnitBox intersection(UnitBox const& dom1, UnitBox const& dom2) { assert(dom1==dom2); return dom1; } //!< .
    friend UnitBox product(UnitBox const& dom1, UnitBox const& dom2) { return UnitBox(dom1.dimension()+dom2.dimension()); } //!< .
    friend UnitBox product(UnitBox const& dom1, UnitInterval const& dom2) { return UnitBox(dom1.dimension()+dom2.dimension()); } //!< .
    friend Bool operator==(UnitBox const& dom1, UnitBox const& dom2) { return dom1.dimension() == dom2.dimension(); } //!< .
    friend OutputStream& operator<<(OutputStream& os, UnitBox const& dom) { return os << "[-1:+1]^" << dom.dimension(); } //!< .
};

struct ScalarTraits {
    template<class X> using Type=Scalar<X>;
    typedef Scalar<Real> Kind;
    typedef SizeOne SizeType;
    typedef IndexZero IndexType;
    using EntireDomainType = RealDomain;
    using BoundedDomainType = IntervalDomainType;
    template<class PR> using RangeType = Interval<FloatUpperBound<PR>>;
};

struct VectorTraits {
    template<class X> using Type=Vector<X>;
    typedef Vector<Real> Kind;
    typedef Ariadne::SizeType SizeType;
    typedef Ariadne::SizeType IndexType;
    using EntireDomainType = EuclideanDomain;
    using BoundedDomainType = BoxDomainType;
    template<class PR> using RangeType = Box<Interval<FloatUpperBound<PR>>>;
};

template<class R> struct DomainTraits;
template<> struct DomainTraits<RealScalar> : ScalarTraits { };
template<> struct DomainTraits<RealVector> : VectorTraits { };

template<class S> struct ElementTraits;
template<class S, class X> using ElementType = typename ElementTraits<S>::template Type<X>;

template<class UB> struct ElementTraits<Interval<UB>> : ScalarTraits { };
template<class IVL> struct ElementTraits<Box<IVL>> : VectorTraits { };
template<> struct ElementTraits<RealDomain> : ScalarTraits { };
template<> struct ElementTraits<EuclideanDomain> : VectorTraits { };
template<> struct ElementTraits<UnitInterval> : ScalarTraits { };
template<> struct ElementTraits<UnitBox> : VectorTraits { };

template<class... TS> using CartesianProductType = decltype(product(declval<TS>()...));

} // namespace Ariadne


#endif /* ARIADNE_DOMAIN_HPP */
