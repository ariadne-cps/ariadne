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

#include "function_interface.h"
#include "function_mixin.h"

#include "macros.h"
#include "pointer.h"
#include "container.h"

#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "taylor_model.h"
#include "differential.h"
#include "formula.h"

#include "../src/function_mixin.tcc"

namespace Ariadne {

//! \related ScalarUserFunction
//! \brief Template structure containing scalar function meta-information.
template<uint AS, uint PS=0u, uint SM=255u>
struct ScalarFunctionData
{
    static const uint argument_size() { return AS; }
    static const uint parameter_size() { return PS; }
    static const uint smoothness() { return SM; }
};

//! \brief A wrapper for converting templated C++ functions to %Ariadne functions.
//!
//! Given a C++ class T with a (static) template method
//!   <code>template<class R, class A, class P> compute(R& r, const A& a, const P& p);</code>
//! the type <code>ScalarUserFunction<T></code> is an Ariadne function defined by \f$r=f(a)\f$.
//! The constructor for ScalarUserFunction<T> takes a Vector<ApproximateNumber> or Vector<ValidatedNumber> argument which is used for \a p.
//!
//! The class T must also define meta-data <c>argument_size(), parameter_size()
//! and smoothness()</c>. These are most easily defined by inheriting from the
//! <tt>ScalarFunctionData<AS,PS,SM=SMOOTH></tt> class.
//!
//! The constant \a SMOOTH is used for an arbitrarily-differentiable function.
template<class T> class ScalarUserFunction
    : public EffectiveScalarFunction
{
  private:
    class Representation
        : public ScalarFunctionMixin< Representation, EffectiveNumber >
    {
      private:
        Vector<EffectiveNumber> _p;
      public:
        typedef uint SizeType;
        Representation(const Vector<EffectiveNumber>& p) : _p(p) { }

        template<class R, class A> inline void _compute(R& r, const A& a) const { T::compute(r,a,_p); }

        virtual Representation* clone() const { return new Representation(*this); }

        virtual SizeType argument_size() const { return T::argument_size(); }
        virtual SizeType parameter_size() const { return T::parameter_size(); }


        virtual EffectiveScalarFunction derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }

        virtual Vector<ApproximateNumber> gradient(const Vector<ApproximateNumber>& x) const {
            return this->evaluate(Differential<ApproximateNumber>::variables(1u,x)).gradient(); }
        virtual Vector<ValidatedNumber> gradient(const Vector<ValidatedNumber>& x) const {
            return this->evaluate(Differential<ValidatedNumber>::variables(1u,x)).gradient(); }

         virtual std::ostream& repr(std::ostream& os) const  {
            return os << "USER"; }
       virtual std::ostream& write(std::ostream& os) const  {
            return os << "ScalarUserFunction( argument_size="<<this->argument_size()<<" )"; }
    };
  public:
    ScalarUserFunction() : EffectiveScalarFunction(new Representation(Vector<EffectiveNumber>(this->parameter_size()))) { }
    ScalarUserFunction(const Vector<ExactNumber>& p) : EffectiveScalarFunction(new Representation(Vector<EffectiveNumber>(p))) { }
    ScalarUserFunction(const Vector<EffectiveNumber>& p) : EffectiveScalarFunction(new Representation(p)) { }

    uint parameter_size() const { return T().parameter_size(); }
    //const Vector<ValidatedNumber> parameters() const { return _p; }
};




//! \brief A convenience class defining the meta-data of an %Ariadne function.
template<uint RS, uint AS, uint PS=0u, uint SM=255u>
class VectorFunctionData
{
  public:
    //!
    static const uint result_size() { return RS; }
    //!
    static const uint argument_size() { return AS; }
    //!
    static const uint parameter_size() { return PS; }
    //!
    static const uint smoothness() { return SM; }
};




//! \brief A wrapper for converting templated C++ functions to %Ariadne functions.
//!
//! Given a C++ class T with a (static) template method
//!   <code>template<class R, class A, class P> compute(R& r, const A& a, const P& p);</code>
//! the type <code>VectorUserFunction<T></code> is an Ariadne function defined by \f$r=f(a)\f$.
//! The constructor for VectorUserFunction<T> takes a Vector<EffectiveNumber> argument which is used for \a p.
//!
//! The class T must also define meta-data <c>result_size(), argument_size(), parameter_size()
//! and smoothness()</c> as static functions. These are most easily defined by inheriting from the
//! <tt>VectorFunctionData<RS,AS,PS,SM=SMOOTH></tt> class.
//!
//! The constant \a SMOOTH is used for an arbitrarily-differentiable function.
template<class T> class VectorUserFunction
    : public EffectiveVectorFunction
{
  private:
    class Representation
        : public VectorFunctionMixin< Representation, EffectiveNumber >
    {
      public:
        typedef uint SizeType;
        Representation(const Vector<EffectiveNumber>& p) : _p(p) { }

        virtual Representation* clone() const { return new Representation(*this); }

        virtual SizeType result_size() const { return T::result_size(); }
        virtual SizeType argument_size() const { return T::argument_size(); }
        virtual SizeType parameter_size() const { return T::parameter_size(); }

        template<class R, class A> inline void _compute(R& r, const A& a) const { T::compute(r,a,_p); }

        virtual Matrix<ApproximateNumber> jacobian(const Vector<ApproximateNumber>& x) const {
            return Ariadne::jacobian(this->evaluate(Differential<ApproximateNumber>::variables(1u,x))); }
        virtual Matrix<ValidatedNumber> jacobian(const Vector<ValidatedNumber>& x) const {
            return Ariadne::jacobian(this->evaluate(Differential<ValidatedNumber>::variables(1u,x))); }

        virtual EffectiveScalarFunctionInterface* _get(uint i) const { ARIADNE_NOT_IMPLEMENTED; }
        virtual EffectiveScalarFunction operator[](uint i) const { ARIADNE_NOT_IMPLEMENTED; }

        // TODO: Find a better way for writing functions which can handle transformations which may not have a
        // write() method or operator<<.
        virtual std::ostream& write(std::ostream& os) const  {
            return os << "VectorUserFunction( result_size="<<this->result_size()<<", argument_size="<<this->argument_size()<<" )"; }

        Vector<EffectiveNumber> _p;
    };
  public:
    //VectorUserFunction() : VectorFunction(new Representation(Vector<EffectiveNumber>(this->parameter_size()))) { }
    //VectorUserFunction(const Vector<ApproximateNumber>& p) : VectorFunction(new Representation(Vector<EffectiveNumber>(p))) { }
    //VectorUserFunction(const Vector<ValidatedNumber>& p) : VectorFunction(new Representation(Vector<EffectiveNumber>(p))) { }
    VectorUserFunction(const Vector<EffectiveNumber>& p) : EffectiveVectorFunction(new Representation(p)) { }
    const Vector<EffectiveNumber>& parameters() const { return dynamic_cast<const Representation*>(this->raw_pointer())->_p; }
};


} // namespace Ariadne

#endif /* ARIADNE_USER_FUNCTION_H */
