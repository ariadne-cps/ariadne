/***************************************************************************
 *            chebyshev_model.hpp
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
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file chebyshev_model.hpp
 *  \brief Chebyshev functions on a bounded domain with a sparse representation.
 */

#ifndef ARIADNE_CHEBYSHEV_MODEL_HPP
#define ARIADNE_CHEBYSHEV_MODEL_HPP

#include <map>

#include "utility/macros.hpp"
#include "utility/typedefs.hpp"
#include "utility/array.hpp"
#include "utility/pointer.hpp"
#include "utility/metaprogramming.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/expansion.hpp"
#include "algebra/operations.hpp"

#include "function/chebyshev_polynomial.hpp"
#include "function/function_model_mixin.hpp"

namespace Ariadne {

struct IndexLess;

template<class X, class XE=X> class MultivariateChebyshevModel;


template<class X, class XE> class MultivariateChebyshevModel
{
    using ARG=RealVector;
    using RES=RealScalar;
    using SIG=RES(ARG);

    typedef typename X::Paradigm P;
    typedef typename X::PrecisionType PR;
    typedef typename XE::PrecisionType PRE;
    using F=Float<PR>;
    using FE=Float<PRE>;

    MultivariateChebyshevPolynomial<X> _cp;
    Error<XE> _e;
  public:
    template<class XX> using ChebyshevModel = MultivariateChebyshevModel<XX>;

    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef X NumericType;
    typedef typename Expansion<MultiIndex,X>::ConstIterator ConstIterator;

    typedef BoxDomainType DomainType;
    typedef IntervalDomainType CodomainType;
    typedef Interval<UpperBound<Float<PR>>> RangeType;
    typedef PRE ErrorPrecisionType;
    typedef Float<PR> RawFloatType;
    typedef Expansion<MultiIndex,X> ExpansionType;
  public:
#warning
//    MultivariateChebyshevModel(SizeType as, PR pr) : _cp(as,X(pr)) { }

    static ChebyshevModel<X> constant(SizeType as, Number<P> y, PR pr);
    static ChebyshevModel<X> constant(SizeType as, X const& c);
    static ChebyshevModel<X> coordinate(SizeType as, SizeType i, PR pr);
    static ChebyshevModel<X> basis(SizeType as, SizeType i, SizeType k, PR pr);

    ChebyshevModel<X> create_constant(Number<P> y);

//#warning Move these to correct place
//#warning Unimplemented methods
    SizeOne result_size() const;
    SizeType argument_size() const { return this->_cp.argument_size(); }
    PR precision() const { return this->_cp.zero_coefficient().precision(); }
    X const& zero_coefficient() const { return _cp.zero_coefficient(); }
    virtual FunctionModelFactoryInterface<P,PR,PRE>* _factory() const;
//#warning Derivative should take ArgumentSizeType as second argument
    friend ChebyshevModel<X> derivative(ChebyshevModel<X> cm, SizeType j);
    friend FunctionModelCreator<FunctionModelFactory<P,PR,PRE>,ARG> factory(ChebyshevModel<X> const& cm);
    template<class PR, class PRE> MultivariateChebyshevModel(FunctionModel<P,SIG,PR,PRE> fm);
    XE error() const;
    friend auto norm(ChebyshevModel<X> const cm) { return cm.sup_norm(); }
//#warning End of new methods

    virtual DomainType const domain() const;
    virtual CodomainType const codomain() const;
    virtual FunctionModelInterface<P,SIG,PR,PRE>* _create() const;
    virtual RangeType const range() const;
    virtual Void clobber();

    MagType<X> sup_norm() const;
    friend MagType<X> sup_norm(ChebyshevModel<X> const& cm) { return cm.sup_norm(); }

    friend ChebyshevModel<X> operator*(Number<P> const& s1, ChebyshevModel<X> cm2) { return X(s1,cm2.precision()) * cm2; }
    friend OutputStream& operator<<(OutputStream& os, ChebyshevModel<X> const& cm) { return cm._write(os); }

    X operator() (Vector<X> const& x) const;
    friend X evaluate(ChebyshevModel<X> const& f, Vector<X> const& x) { return f(x); }
  private: public:
    static ChebyshevModel<X> apply(Neg, ChebyshevModel<X> cm);
    static ChebyshevModel<X> apply(Sqr, ChebyshevModel<X> cm);
    static ChebyshevModel<X> apply(Add, ChebyshevModel<X> const& cm1, ChebyshevModel<X> const& cm2);
    static ChebyshevModel<X> apply(Sub, ChebyshevModel<X> const& cm1, ChebyshevModel<X> const& cm2);
    static ChebyshevModel<X> apply(Mul, ChebyshevModel<X> const& cm1, ChebyshevModel<X> const& cm2);
    static ChebyshevModel<X> apply(Add, ChebyshevModel<X> cm1, Scalar<X> const& s2);
    static ChebyshevModel<X> apply(Mul, ChebyshevModel<X> cm1, Scalar<X> const& s2);
    OutputStream& _write(OutputStream& os) const;

    static X _eval(ConstIterator& from, ConstIterator end, Vector<X> const& s, SizeType k);
    static ChebyshevModel<X> _mul_from(ConstIterator& from1, ConstIterator end1, ConstIterator& from2, ConstIterator end2, SizeType k);
  public:
};

template<class X> struct AlgebraOperations<MultivariateChebyshevModel<X>,X> : public MultivariateChebyshevModel<X> { };

} // namespace Ariadne

#endif // ARIADNE_CHEBYSHEV_MODEL_HPP
