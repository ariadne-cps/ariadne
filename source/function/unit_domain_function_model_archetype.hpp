/***************************************************************************
 *            unit_domain_function_model_archetype.hpp
 *
 *  Copyright 2022  Pieter Collins
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

#ifndef ARIADNE_UNIT_DOMAIN_FUNCTION_MODEL_ARCHETYPE_HPP
#define ARIADNE_UNIT_DOMAIN_FUNCTION_MODEL_ARCHETYPE_HPP

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"

namespace Ariadne {

class ModelArchetype;
template<class X> X evaluate(ModelArchetype const&, Vector<X> const&);
template<class X> Vector<X> evaluate(Vector<ModelArchetype> const&, Vector<X> const&);

class ModelArchetype {

    using P = ValidatedTag;
    using RES = Real;
    using ARG = RealVector;
    using SIG = RES(ARG);

    using X = FloatDP;
    using PR = DoublePrecision;
    using PRE = DoublePrecision;

    using F = Float<PR>;
    using FE = Error<Float<PRE>>;
  public:

    typedef P Paradigm;
    typedef Bounds<X> NumericType;
    typedef Number<P> GenericNumericType;

    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;

        typedef X RawFloatType;

    typedef BoxDomainType DomainType;
    typedef IntervalDomainType CodomainType;
    typedef IntervalRangeType RangeType;

    typedef SizeType ArgumentSizeType;
    typedef SizeType ArgumentIndexType;

        typedef Expansion<MultiIndex,X> ExpansionType;
    typedef F CoefficientType;
    typedef F ValueType;
    typedef FE ErrorType;

    typedef PositiveUpperBound<F> NormType;
//#warning PropertiesType wrong
    typedef Sweeper<X> PropertiesType;

    //typedef Function<P,SIG> FunctionType;
    typedef Function<P,RealScalar(ARG)> ScalarFunctionType;
    typedef Function<P,RealVector(ARG)> VectorFunctionType;

  public:
//#warning  No default constructor!
        ModelArchetype();
    ModelArchetype(ArgumentSizeType, PropertiesType);
        ModelArchetype(ExpansionType, ErrorType, PropertiesType);
        ModelArchetype(ExpansionType, Float<PRE>, PropertiesType);
        ModelArchetype(Expansion<MultiIndex,ExactDouble>, ExactDouble, PropertiesType);

    static ModelArchetype zero(ArgumentSizeType,PropertiesType); // scaled_function_patch.tpl.hpp:249
    static ModelArchetype constant(ArgumentSizeType,GenericNumericType,PropertiesType); // scaled_function_patch.tpl.hpp:249
    static ModelArchetype scaling(ArgumentSizeType,ArgumentIndexType,IntervalDomainType,PropertiesType);
    static ModelArchetype unit_ball(ArgumentSizeType,PropertiesType);
    static Vector<ModelArchetype> constants(ArgumentSizeType,Vector<GenericNumericType>,PropertiesType);
    static Vector<ModelArchetype> scalings(DomainType,PropertiesType);

    Void simplify();
    Void simplify(PropertiesType);

    ModelArchetype& operator=(NumericType const&);
    ModelArchetype& operator=(GenericNumericType const&);

            ExpansionType& expansion();
        CoefficientType& operator[](MultiIndex const&);
    ErrorType& error();

    CodomainType codomain() const;
    RangeType range() const;
    ArgumentSizeType argument_size() const;
    DegreeType degree() const;
            SizeType number_of_nonzeros() const;
            ExpansionType const& expansion() const;
        CoefficientType const& operator[](MultiIndex const&) const;
    ValueType const& value() const;
    ErrorType const& error() const;
    PropertiesType const& properties() const;
        PrecisionType const precision() const;

    Void set_value(ValueType const& c);
    Void set_error(ErrorType const& e);
    Void set_properties(PropertiesType const& prp);
    Void clobber();

    template<class X> friend X evaluate(ModelArchetype const&, Vector<X> const&);
    template<class X> friend Vector<X> evaluate(Vector<ModelArchetype> const&, Vector<X> const&);
    static ModelArchetype apply(ScalarFunctionType const&, Vector<ModelArchetype> const&);


    friend OutputStream& operator<<(OutputStream& os, ModelArchetype const&);

    friend ModelArchetype compose(ModelArchetype const&, Vector<ModelArchetype> const&);

    friend ModelArchetype nul(ModelArchetype);
    friend ModelArchetype pos(ModelArchetype);
    friend ModelArchetype neg(ModelArchetype);
    friend ModelArchetype sqr(ModelArchetype);
    friend ModelArchetype hlf(ModelArchetype);
    friend ModelArchetype rec(ModelArchetype);
    friend ModelArchetype add(ModelArchetype const&, ModelArchetype const&);
    friend ModelArchetype sub(ModelArchetype const&, ModelArchetype const&);
    friend ModelArchetype mul(ModelArchetype const&, ModelArchetype const&);
    friend ModelArchetype div(ModelArchetype const&, ModelArchetype const&);
    friend ModelArchetype pow(ModelArchetype const&, Int n);
    friend ModelArchetype sqrt(ModelArchetype);
    friend ModelArchetype exp(ModelArchetype);
    friend ModelArchetype log(ModelArchetype);
    friend ModelArchetype sin(ModelArchetype);
    friend ModelArchetype cos(ModelArchetype);
    friend ModelArchetype tan(ModelArchetype);
    friend ModelArchetype asin(ModelArchetype);
    friend ModelArchetype acos(ModelArchetype);
    friend ModelArchetype atan(ModelArchetype);

    friend ModelArchetype max(ModelArchetype const&, ModelArchetype const&);
    friend ModelArchetype min(ModelArchetype const&, ModelArchetype const&);
    friend ModelArchetype abs(ModelArchetype const&);

    friend ModelArchetype add(ModelArchetype const&, NumericType const&);
    friend ModelArchetype sub(ModelArchetype const&, NumericType const&);
    friend ModelArchetype mul(ModelArchetype const&, NumericType const&);
    friend ModelArchetype div(ModelArchetype const&, NumericType const&);
    friend ModelArchetype add(NumericType const&, ModelArchetype const&);
    friend ModelArchetype sub(NumericType const&, ModelArchetype const&);
    friend ModelArchetype mul(NumericType const&, ModelArchetype const&);
    friend ModelArchetype div(NumericType const&, ModelArchetype const&);
    friend ModelArchetype max(ModelArchetype const&, NumericType const&);
    friend ModelArchetype min(ModelArchetype const&, NumericType const&);
    friend ModelArchetype max(NumericType const&, ModelArchetype const&);
    friend ModelArchetype min(NumericType const&, ModelArchetype const&);


    friend ModelArchetype operator+(ModelArchetype const&, ModelArchetype const&);
    friend ModelArchetype operator-(ModelArchetype const&, ModelArchetype const&);
    friend ModelArchetype operator*(ModelArchetype const&, ModelArchetype const&);
    friend ModelArchetype operator/(ModelArchetype const&, ModelArchetype const&);

    friend ModelArchetype operator+(ModelArchetype const&, NumericType const&);
    friend ModelArchetype operator-(ModelArchetype const&, NumericType const&);
    friend ModelArchetype operator*(ModelArchetype const&, NumericType const&);
    friend ModelArchetype operator/(ModelArchetype const&, NumericType const&);
    friend ModelArchetype operator+(NumericType const&, ModelArchetype const&);
    friend ModelArchetype operator-(NumericType const&, ModelArchetype const&);
    friend ModelArchetype operator*(NumericType const&, ModelArchetype const&);
    friend ModelArchetype operator/(NumericType const&, ModelArchetype const&);

    friend ModelArchetype operator+(ModelArchetype const&, GenericNumericType const&);
    friend ModelArchetype operator-(ModelArchetype const&, GenericNumericType const&);
    friend ModelArchetype operator*(ModelArchetype const&, GenericNumericType const&);
    friend ModelArchetype operator/(ModelArchetype const&, GenericNumericType const&);
    friend ModelArchetype operator+(GenericNumericType const&, ModelArchetype const&);
    friend ModelArchetype operator-(GenericNumericType const&, ModelArchetype const&);
    friend ModelArchetype operator*(GenericNumericType const&, ModelArchetype const&);
    friend ModelArchetype operator/(GenericNumericType const&, ModelArchetype const&);

    friend NumericType evaluate(ModelArchetype const&, Vector<NumericType>);
    friend Vector<NumericType> evaluate(Vector<ModelArchetype> const&, Vector<NumericType>);

    friend ModelArchetype partial_evaluate(ModelArchetype const&, SizeType j, NumericType);
    friend Vector<ModelArchetype> partial_evaluate(Vector<ModelArchetype> const&, SizeType j, NumericType);

    friend PositiveUpperBound<Float<PR>> norm(ModelArchetype const&);

    friend ModelArchetype derivative(ModelArchetype, ArgumentIndexType);
    friend ModelArchetype antiderivative(ModelArchetype, ArgumentIndexType);
    friend Covector<NumericType> gradient(ModelArchetype, Vector<NumericType>);

    friend ModelArchetype embed(SizeType,ModelArchetype,SizeType);
    friend Vector<ModelArchetype> embed(SizeType,Vector<ModelArchetype>,SizeType);
    friend Vector<ModelArchetype> combine(Vector<ModelArchetype>,Vector<ModelArchetype>);

    friend Bool inconsistent(ModelArchetype const&, ModelArchetype const&);
    friend Bool refines(ModelArchetype const&, ModelArchetype const&);
    friend ModelArchetype refinement(ModelArchetype const&, ModelArchetype const&);

};
template<Same<ModelArchetype> M, Same<Dyadic> N> M operator-(M,N);
template<Same<ModelArchetype> M, Same<Dyadic> N> M operator/(M,N);
template<Same<ModelArchetype> M, Same<Int> N> M operator*(M,N);

} // namespace Ariadne

#endif /* ARIADNE_UNIT_DOMAIN_FUNCTION_MODEL_ARCHETYPE_HPP */
