/***************************************************************************
 *            function/formula_function.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file function/formula_function.hpp
 *  \brief Formula functions
 */

#ifndef ARIADNE_FORMULA_FUNCTION_HPP
#define ARIADNE_FORMULA_FUNCTION_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "../function/function_interface.hpp"

#include "../utility/macros.hpp"
#include "../utility/pointer.hpp"
#include "../utility/container.hpp"
#include "../utility/metaprogramming.hpp"

#include "../numeric/numeric.hpp"
#include "../numeric/operators.tpl.hpp"
#include "../algebra/vector.hpp"

#include "../function/function_mixin.hpp"
#include "../function/projection.hpp"
#include "../function/formula.hpp"

namespace Ariadne {

//------------------------ Formula functions  -----------------------------------//


//! A function defined by a formula
template<class Y>
struct ScalarUnivariateFormulaFunction
    : ScalarUnivariateFunctionMixin<ScalarUnivariateFormulaFunction<Y>,InformationTag<Y>>
{
    typedef InformationTag<Y> P;
    Formula<Y> _formula;

    ScalarUnivariateFormulaFunction(const Formula<Y>& f) : _formula(f) { }
    operator Formula<Y>() const { return _formula; }

    virtual SizeOne argument_size() const final { return SizeOne(); }
    virtual SizeOne result_size() const final { return SizeOne(); }
    virtual ScalarUnivariateFormulaFunction<Y>* _derivative(SizeOne j) const final {
        return new ScalarUnivariateFormulaFunction<Y>(Ariadne::derivative(_formula,0)); }
    virtual OutputStream& _write(OutputStream& os) const final { return os << this->_formula; }
    virtual OutputStream& repr(OutputStream& os) const final { return os << "FormulaFunction("<<this->_formula<<")"; }

    template<class X> Void _compute(X& r, const X& x) const;
};



//! A function defined by a formula
template<class Y>
struct VectorUnivariateFormulaFunction
    : VectorUnivariateFunctionMixin<VectorUnivariateFormulaFunction<Y>,InformationTag<Y>>
{
    typedef InformationTag<Y> P;
    Vector<Formula<Y>> _formulae;

    VectorUnivariateFormulaFunction(const List< Formula<Y> >& f) : _formulae(f) { }
    VectorUnivariateFormulaFunction(const Vector< Formula<Y> >& f) : _formulae(f) { }

    ScalarUnivariateFormulaFunction<Y> operator[](SizeType i) const { return ScalarUnivariateFormulaFunction(_formulae[i]); }

    virtual SizeType result_size() const { return this->_formulae.size(); }
    virtual SizeOne argument_size() const { return SizeOne(); }
    virtual ScalarUnivariateFormulaFunction<Y>* _get(SizeType i) const { return new ScalarUnivariateFormulaFunction<Y>(this->_formulae[i]); }
    virtual VectorUnivariateFormulaFunction<Y>* _derivative(SizeOne k) const {
        return new VectorUnivariateFormulaFunction<Y>(Vector<Formula<Y>>(this->_formulae.size(),[&](SizeType i){return derivative(this->_formulae[i],k);})); }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_formulae; }
    virtual OutputStream& repr(OutputStream& os) const { return os << "VectorUnivariateFormulaFunction("<<this->result_size()<<","<<this->_formulae<<")"; }

    template<class X> Void _compute(Vector<X>& r, const X& x) const;
};

//! A function defined by a formula
template<class Y>
struct ScalarFormulaFunction
    : ScalarMultivariateFunctionMixin<ScalarFormulaFunction<Y>,InformationTag<Y>>
{
    typedef InformationTag<Y> P;
    SizeType _argument_size;
    Formula<Y> _formula;

    ScalarFormulaFunction(SizeType as, const Formula<Y>& f) : _argument_size(as), _formula(f) { }
    operator Formula<Y>() const { return _formula; }

    virtual SizeType argument_size() const final { return _argument_size; }
    virtual SizeOne result_size() const final { return SizeOne(); }
    virtual ScalarFormulaFunction<Y>* _derivative(SizeType j) const final { return new ScalarFormulaFunction<Y>(_argument_size,Ariadne::derivative(_formula,j)); }
    virtual OutputStream& _write(OutputStream& os) const final { return os << this->_formula; }
    virtual OutputStream& repr(OutputStream& os) const final { return os << "FormulaFunction("<<this->_argument_size<<","<<this->_formula<<")"; }

    template<class X> Void _compute(X& r, const Vector<X>& x) const;
};



//! A vector function defined by formulae
template<class Y>
struct VectorFormulaFunction
    : VectorMultivariateFunctionMixin<VectorFormulaFunction<Y>,InformationTag<Y>>
{
    SizeType _argument_size;
    Vector< Formula<Y> > _formulae;

    VectorFormulaFunction(SizeType as, const List< Formula<Y> >& f) : _argument_size(as), _formulae(f) { }
    VectorFormulaFunction(SizeType as, const Vector< Formula<Y> >& f) : _argument_size(as), _formulae(f) { }

    ScalarFormulaFunction<Y> operator[](SizeType i) const { return ScalarFormulaFunction(_argument_size,_formulae[i]); }

    virtual SizeType result_size() const { return this->_formulae.size(); }
    virtual SizeType argument_size() const { return this->_argument_size; }
    virtual ScalarFormulaFunction<Y>* _get(SizeType i) const { return new ScalarFormulaFunction<Y>(this->_argument_size,this->_formulae[i]); }
    virtual VectorFormulaFunction<Y>* _derivative(SizeType k) const {
        return new VectorFormulaFunction<Y>(this->_argument_size, Vector<Formula<Y>>(this->_formulae.size(),[&](SizeType i){return derivative(this->_formulae[i],k);})); }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_formulae; }
    virtual OutputStream& repr(OutputStream& os) const { return os << "VectorFormulaFunction("<<this->result_size()<<","<<this->argument_size()<<","<<this->_formulae<<")"; }

    template<class X> Void _compute(Vector<X>& r, const Vector<X>& x) const;
};



typedef ScalarUnivariateFormulaFunction<EffectiveNumber> EffectiveScalarUnivariateFormulaFunction;
typedef VectorUnivariateFormulaFunction<EffectiveNumber> EffectiveVectorUnivariateFormulaFunction;
typedef ScalarFormulaFunction<EffectiveNumber> EffectiveScalarFormulaFunction;
typedef VectorFormulaFunction<EffectiveNumber> EffectiveVectorFormulaFunction;

typedef Pair<Nat,EffectiveFormula> CoordinateFormulaPair;
typedef List<CoordinateFormulaPair> CoordinateFormulaPairs;

class NotFormulaFunctionException : public std::runtime_error {
  public:
    NotFormulaFunctionException(const String& str) : std::runtime_error(str) { }
};

//! \brief Returns \a true if the function \a f is syntactically constant in the indices \a is.
template<class Y> Bool is_constant_in(const ScalarFormulaFunction<Y>& f, const Set<Nat>& is) { return is_constant_in(f._formula,is); }
//! \brief Returns \a true if the function \a f is syntactically affine in the indices \a is.
template<class Y> Bool is_affine_in(const ScalarFormulaFunction<Y>& f, const Set<Nat>& is) { return is_affine_in(f._formula,is); }
//! \brief Returns \a true if the vector function \a f is syntactically affine in the indices \a is.
template<class Y> Bool is_affine_in(const VectorFormulaFunction<Y>& f, const Set<Nat>& is) { return is_affine_in(f._formulae,is); }
//! \brief Returns \a true if the vector function \a f is syntactically additive (possibly with multipliers) in the indices \a is.
template<class Y> Bool is_additive_in(const VectorFormulaFunction<Y>& f, const Set<Nat>& is) { return is_additive_in(f._formulae,is); }

template<class Y> Bool is_affine_in(const VectorMultivariateFunction<Y>& f, const Set<Nat>& is) {
    auto ff = dynamic_cast<const EffectiveVectorFormulaFunction*>(f.raw_pointer());
    if (ff == nullptr) ARIADNE_THROW(NotFormulaFunctionException,"is_affine_in(f,is)","Affinity checking currently available only for formula functions.");
    return is_affine_in(ff->_formulae,is);
}
template<class Y> Bool is_additive_in(const VectorMultivariateFunction<Y>& f, const Set<Nat>& is) {
    auto ff = dynamic_cast<const EffectiveVectorFormulaFunction*>(f.raw_pointer());
    if (ff == nullptr) ARIADNE_THROW(NotFormulaFunctionException,"is_additive_in(f,is)","Additivity checking currently available only for formula functions.");
    return is_additive_in(ff->_formulae,is);
}



} // namespace Ariadne

#endif
