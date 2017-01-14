/***************************************************************************
 *            user_function.h
 *
 *  Copyright 2008-9  Pieter Collins
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

/*! \file user_function.h
 *  \brief Wrappers for user functions
 */

#ifndef ARIADNE_USER_FUNCTION_H
#define ARIADNE_USER_FUNCTION_H

#include <cstdarg>
#include <iosfwd>

#include "function/function_interface.h"
#include "function/function_mixin.h"

#include "utility/macros.h"
#include "utility/pointer.h"
#include "utility/container.h"

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "function/taylor_model.h"
#include "algebra/differential.h"
#include "function/formula.h"

#include "function_mixin.tcc"

namespace Ariadne {

//! \related ScalarUserFunction
//! \brief Template structure containing scalar function meta-information.
template<SizeType AS, SizeType PS=0u, DegreeType SM=255u>
struct ScalarFunctionData
{
    static const SizeType argument_size() { return AS; }
    static const SizeType parameter_size() { return PS; }
    static const DegreeType smoothness() { return SM; }
};

//! \brief A wrapper for converting templated C++ functions to %Ariadne functions.
//!
//! Given a C++ class F with a (static) template method
//!   <code>template<class R, class A, class P> compute(R& r, const A& a, const P& p);</code>
//! the type <code>ScalarUserFunction<F></code> is an Ariadne function defined by \f$r=f(a)\f$.
//! The constructor for ScalarUserFunction<F> takes a Vector<ApproximateNumericType> or Vector<ValidatedNumericType> argument which is used for \a p.
//!
//! The class F must also define meta-data <c>argument_size(), parameter_size()
//! and smoothness()</c>. These are most easily defined by inheriting from the
//! <tt>ScalarFunctionData<AS,PS,SM=SMOOTH></tt> class.
//!
//! The constant \a SMOOTH is used for an arbitrarily-differentiable function.
template<class F> class ScalarUserFunction
    : public EffectiveScalarFunction
{
  private:
    class Representation
        : public ScalarFunctionMixin< Representation, EffectiveTag >
    {
      private:
        Vector<EffectiveNumericType> _p;
      public:
        Representation(const Vector<EffectiveNumericType>& p) : _p(p) { }

        template<class R, class A> inline Void _compute(R& r, const A& a) const { F::compute(r,a,_p); }

        virtual Representation* clone() const { return new Representation(*this); }

        virtual SizeType argument_size() const { return F::argument_size(); }
        virtual SizeType parameter_size() const { return F::parameter_size(); }


        virtual EffectiveScalarFunction derivative(SizeType j) const { ARIADNE_NOT_IMPLEMENTED; }

        virtual Covector<ApproximateNumericType> gradient(const Vector<ApproximateNumericType>& x) const {
            return this->evaluate(Differential<ApproximateNumericType>::variables(1u,x)).gradient(); }
        virtual Covector<ValidatedNumericType> gradient(const Vector<ValidatedNumericType>& x) const {
            return this->evaluate(Differential<ValidatedNumericType>::variables(1u,x)).gradient(); }

         virtual OutputStream& repr(OutputStream& os) const  {
            return os << "USER"; }
       virtual OutputStream& write(OutputStream& os) const  {
            return os << "ScalarUserFunction( argument_size="<<this->argument_size()<<" )"; }
    };
  public:
    ScalarUserFunction() : EffectiveScalarFunction(new Representation(Vector<EffectiveNumericType>(this->parameter_size()))) { }
    ScalarUserFunction(const Vector<ExactNumericType>& p) : EffectiveScalarFunction(new Representation(Vector<EffectiveNumericType>(p))) { }
    ScalarUserFunction(const Vector<EffectiveNumericType>& p) : EffectiveScalarFunction(new Representation(p)) { }

    SizeType parameter_size() const { return F().parameter_size(); }
    //const Vector<ValidatedNumericType> parameters() const { return _p; }
};




//! \brief A convenience class defining the meta-data of an %Ariadne function.
template<SizeType RS, SizeType AS, SizeType PS=0u, DegreeType SM=255u>
class VectorFunctionData
{
  public:
    //!
    static const SizeType result_size() { return RS; }
    //!
    static const SizeType argument_size() { return AS; }
    //!
    static const SizeType parameter_size() { return PS; }
    //!
    static const DegreeType smoothness() { return SM; }
};




//! \brief A wrapper for converting templated C++ functions to %Ariadne functions.
//!
//! Given a C++ class F with a (static) template method
//!   <code>template<class R, class A, class P> compute(R& r, const A& a, const P& p);</code>
//! the type <code>VectorUserFunction<F></code> is an Ariadne function defined by \f$r=f(a)\f$.
//! The constructor for VectorUserFunction<F> takes a Vector<EffectiveNumericType> argument which is used for \a p.
//!
//! The class F must also define meta-data <c>result_size(), argument_size(), parameter_size()
//! and smoothness()</c> as static functions. These are most easily defined by inheriting from the
//! <tt>VectorFunctionData<RS,AS,PS,SM=SMOOTH></tt> class.
//!
//! The constant \a SMOOTH is used for an arbitrarily-differentiable function.
template<class F> class VectorUserFunction
    : public EffectiveVectorFunction
{
  private:
    class Representation
        : public VectorFunctionMixin< Representation, EffectiveTag >
    {
      public:
        Representation(const Vector<EffectiveNumericType>& p) : _p(p) { }

        virtual Representation* clone() const { return new Representation(*this); }

        virtual SizeType result_size() const { return F::result_size(); }
        virtual SizeType argument_size() const { return F::argument_size(); }
        virtual SizeType parameter_size() const { return F::parameter_size(); }

        template<class R, class A> inline Void _compute(R& r, const A& a) const { F::compute(r,a,_p); }

        virtual Matrix<ApproximateNumericType> jacobian(const Vector<ApproximateNumericType>& x) const {
            return Ariadne::jacobian(this->evaluate(Differential<ApproximateNumericType>::variables(1u,x))); }
        virtual Matrix<ValidatedNumericType> jacobian(const Vector<ValidatedNumericType>& x) const {
            return Ariadne::jacobian(this->evaluate(Differential<ValidatedNumericType>::variables(1u,x))); }

        virtual EffectiveScalarFunctionInterface* _get(SizeType i) const { ARIADNE_NOT_IMPLEMENTED; }
        virtual EffectiveScalarFunction operator[](SizeType i) const { ARIADNE_NOT_IMPLEMENTED; }

        // TODO: Find a better way for writing functions which can handle transformations which may not have a
        // write() method or operator<<.
        virtual OutputStream& write(OutputStream& os) const  {
            return os << "VectorUserFunction( result_size="<<this->result_size()<<", argument_size="<<this->argument_size()<<" )"; }

        Vector<EffectiveNumericType> _p;
    };
  public:
    //VectorUserFunction() : VectorFunction(new Representation(Vector<EffectiveNumericType>(this->parameter_size()))) { }
    //VectorUserFunction(const Vector<ApproximateNumericType>& p) : VectorFunction(new Representation(Vector<EffectiveNumericType>(p))) { }
    //VectorUserFunction(const Vector<ValidatedNumericType>& p) : VectorFunction(new Representation(Vector<EffectiveNumericType>(p))) { }
    VectorUserFunction(const Vector<EffectiveNumericType>& p) : EffectiveVectorFunction(new Representation(p)) { }
    const Vector<EffectiveNumericType>& parameters() const { return dynamic_cast<const Representation*>(this->raw_pointer())->_p; }
};


} // namespace Ariadne

#endif /* ARIADNE_USER_FUNCTION_H */
