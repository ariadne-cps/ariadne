/***************************************************************************
 *            function.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "numeric/numeric.hpp"

#include "numeric/operators.hpp"
#include "algebra/differential.hpp"
#include "algebra/algebra.hpp"
#include "function/formula.hpp"
#include "function/taylor_model.hpp"

#include "function/function.hpp"

#include "function/function_mixin.hpp"
#include "function/function_mixin.tpl.hpp"

#include "function/function_model.hpp"

#include "function/symbolic_function.hpp"

#include "expression/expression.hpp"


namespace Ariadne {

template<class T> inline String class_name() { return "Unknown"; }

template<class T> inline StringType str(const T& t) {
    StringStream ss; ss << t; return ss.str(); }

// Templated conversions dynamically checked at runtime
template<class R, class A, EnableIf<IsSame<R,A>> =dummy> R const& checked_same(A const& a) { return a; }
template<class R, class A, DisableIf<IsSame<R,A>> =dummy> R const& checked_same(A const& a) {
    ARIADNE_THROW(std::runtime_error,"checked_same<R,A> with R="<<class_name<R>()<<", A="<<class_name<A>(),"argument "<<a<<" does not have the same type as result."); }
template<class R, class A, EnableIf<IsConvertible<A,R>> =dummy> R checked_convert(A&& a) { return a; }
template<class R, class A, DisableIf<IsConvertible<A,R>> =dummy> R checked_convert(A&& a) {
    ARIADNE_THROW(std::runtime_error,"checked_convert<R,A> with R="<<class_name<R>()<<", A="<<class_name<A>(),"argument "<<a<<" is not convertible to result."); }
template<class R, class A, EnableIf<IsConstructible<R,A>> =dummy> R checked_construct(A const& a) { return R(a); }
template<class R, class A, DisableIf<IsConstructible<R,A>> =dummy> R checked_construct(A const& a) {
    ARIADNE_THROW(std::runtime_error,"checked_construct<R,A> with R="<<class_name<R>()<<", A="<<class_name<A>(),"argument "<<a<<" is not explicitly convertible to result."); }



//------------------------ Formula Function ----------------------------------//

template<class P, class Y> ScalarFunction<P,IntervalDomainType> make_formula_function(IntervalDomainType dom, Scalar<Formula<Y>> const& e) {
    assert(false);
}

template<class P, class Y> VectorFunction<P,IntervalDomainType> make_formula_function(IntervalDomainType dom, Vector<Formula<Y>> const& e) {
    assert(false);
}

template<class P, class Y> ScalarFunction<P,BoxDomainType> make_formula_function(BoxDomainType dom, Scalar<Formula<Y>> const& e) {
    return ScalarFunction<P,BoxDomainType>(new ScalarFormulaFunction<Y>(dom.dimension(),e));
}

template<class P, class Y> VectorFunction<P,BoxDomainType> make_formula_function(BoxDomainType dom, Vector<Formula<Y>> const& e) {
    return VectorFunction<P,BoxDomainType>(new VectorFormulaFunction<Y>(dom.dimension(),e));
}

//------------------------ Function ----------------------------------//

namespace {
OutputStream& operator<<(OutputStream& os, SizeOne so) { return os << "1u"; }
OutputStream& operator<<(OutputStream& os, RealDomain const& dom) { return os << "R"; }
OutputStream& operator<<(OutputStream& os, EuclideanDomain const& dom) { return os << "R" << dom.dimension(); }

template<class D, class DD> D make_domain(DD dom);
template<> IntervalDomainType make_domain<IntervalDomainType,BoxDomainType>(BoxDomainType dom) { throw std::runtime_error(""); }
template<> BoxDomainType make_domain<BoxDomainType,IntervalDomainType>(IntervalDomainType dom) { throw std::runtime_error(""); }
template<> IntervalDomainType make_domain<IntervalDomainType,IntervalDomainType>(IntervalDomainType dom) { return dom; }
template<> BoxDomainType make_domain<BoxDomainType,BoxDomainType>(BoxDomainType dom) { return dom; }

template<class P, class D, class DD> ScalarFunction<P,D> make_zero_function(SizeOne rs, DD dom) {
    return FunctionConstructors<P>::zero(make_domain<D>(dom)); }
template<class P, class D, class DD> VectorFunction<P,D> make_zero_function(SizeType rs, DD dom) {
    return  FunctionConstructors<P>::zeros(rs,make_domain<D>(dom)); }
}

template<class P, class D, class C> Function<P,D,C>::Function() : _ptr() {
}

template<class P, class D, class C> Function<P,D,C>::Function(EuclideanDomain dom) {
    ResultSizeType rs=ResultSizeType(); BoxDomainType bx_dom=dom;
    (*this) = make_zero_function<P,D>(rs,bx_dom);
}

template<class P, class D, class C> Function<P,D,C>::Function(ResultSizeType rs, EuclideanDomain dom) {
    BoxDomainType const& bx_dom=dom;
    (*this) = make_zero_function<P,D>(rs,bx_dom);
}

template<class P, class D, class C> Function<P,D,C>::Function(DomainType dom) {
    ResultSizeType rs=ResultSizeType(); (*this) = make_zero_function<P,D>(rs,dom);
}

template<class P, class D, class C> Function<P,D,C>::Function(ResultSizeType rs, DomainType dom) {
    (*this) = make_zero_function<P,D>(rs,dom);
}

template<class P, class D, class C> Function<P,D,C>::Function(ResultSizeType rs, ScalarFunction<P,D> sf)
    : Function(Vector<ScalarFunction<P,D>>(SizeType(rs),sf)) {
}

template<class P, class D, class C> Function<P,D,C>::Function(InitializerList<ScalarFunction<P,D>> const& lsf)
    : Function(Vector<ScalarFunction<P,D>>(lsf)) {
}

template<class P, class D, class C> Function<P,D,C>::Function(List<ScalarFunction<P,D>> const& lsf)
    : Function(Vector<ScalarFunction<P,D>>(lsf)) {
}

template<class P, class D, class C> Function<P,D,C>::Function(DomainType dom, Result<Formula<Y>>const& e) {
    *this = make_formula_function<P>(dom,e);
}

template<class P, class D, class C> struct MakeVectorFunction;
template<class P, class D> struct MakeVectorFunction<P,D,IntervalDomainType> {
    Function<P,D,IntervalDomainType> create(Vector<ScalarFunction<P,D>> const& lsf) {
        ARIADNE_FAIL_MSG("Cannot construct scalar function from list."); }
};
template<class P, class D> struct MakeVectorFunction<P,D,BoxDomainType> {
    Function<P,D,BoxDomainType> create(Vector<ScalarFunction<P,D>> const& lsf) {
        return Function<P,D,BoxDomainType>(std::make_shared<VectorOfScalarFunction<P,D>>(lsf)); }
};

template<class P, class D, class C> Function<P,D,C> make_vector_function(Vector<ScalarFunction<P,D>> const& lsf) {
    return MakeVectorFunction<P,D,C>().create(lsf);
}

template<class P, class D, class C> Function<P,D,C>::Function(Vector<ScalarFunction<P,D>> const& vsf)
    : Function<P,D,C>(make_vector_function<P,D,C>(vsf)) {
}


//------------------------ Function Constructors -----------------------------------//

template<class P> ScalarFunction<P> FunctionConstructors<P>::zero(BoxDomainType dom) {
    return ConstantFunction<Y>(dom, Y(0));
}


template<class P> ScalarFunction<P> FunctionConstructors<P>::constant(BoxDomainType dom, NumericType c) {
    return ConstantFunction<Y>(dom, c);
}

template<class P> ScalarFunction<P> FunctionConstructors<P>::coordinate(BoxDomainType dom, SizeType j) {
    return CoordinateFunction<P>(dom, j);
}

template<class P> List<ScalarFunction<P>> FunctionConstructors<P>::coordinates(BoxDomainType dom) {
    List<ScalarFunction<P>> r; r.reserve(dom.dimension());
    for(SizeType j=0; j!=dom.dimension(); ++j) { r.append(coordinate(dom,j)); }
    return std::move(r);
}

template<class P> VectorFunction<P> FunctionConstructors<P>::zeros(SizeType rs, BoxDomainType dom) {
    return VectorFunction<P>(new VectorOfScalarFunction<P>(rs,zero(dom)));
}

template<class P> VectorFunction<P> FunctionConstructors<P>::identity(BoxDomainType dom) {
    SizeType n=dom.dimension();
    ScalarFunction<P> z=ScalarFunction<P,BoxDomainType>::zero(dom);
    VectorOfScalarFunction<P>* res = new VectorOfScalarFunction<P>(n,z);
    for(SizeType i=0; i!=n; ++i) {
        res->_vec[i]=ScalarFunction<P>::coordinate(dom,i);
    }
    return VectorFunction<P>(res);
}


template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::zero(IntervalDomainType dom) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::constant(IntervalDomainType dom, NumericType c) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::coordinate(IntervalDomainType dom, SizeType j) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class P> VectorUnivariateFunction<P> FunctionConstructors<P>::zeros(SizeType rs, IntervalDomainType dom) {
    ARIADNE_NOT_IMPLEMENTED;
}


template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::identity(IntervalDomainType dom) {
    ARIADNE_NOT_IMPLEMENTED;
}


template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::zero() {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::constant(NumericType c) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::coordinate() {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class P> VectorUnivariateFunction<P> FunctionConstructors<P>::zeros(SizeType rs) {
    ARIADNE_NOT_IMPLEMENTED;
}


template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::identity() {
    ARIADNE_NOT_IMPLEMENTED;
}


template<class P> ScalarFunction<P> FunctionConstructors<P>::zero(SizeType as) {
    ScalarFunction<P> sf(new ScalarFormulaFunction<Y>(as,Formula<Y>::zero()));
    return sf;
}

template<class P> ScalarFunction<P> FunctionConstructors<P>::constant(SizeType as, NumericType c) {
    return ScalarFunction<P>(new ScalarFormulaFunction<Y>(as,Formula<Y>::constant(c)));
}

template<class P> ScalarFunction<P> FunctionConstructors<P>::coordinate(SizeType as, SizeType j) {
    return ScalarFunction<P>(new ScalarFormulaFunction<Y>(as,Formula<Y>::coordinate(j)));
}

template<class P> VectorFunction<P> FunctionConstructors<P>::zeros(SizeType rs, SizeType as) {
    VectorOfScalarFunction<P>* res = new VectorOfScalarFunction<P>(rs,as);
    for(SizeType i=0; i!=rs; ++i) {
        res->_vec[i]=ScalarFunction<P>::zero(as);
    }
    return VectorFunction<P>(res);
}

template<class P> List<ScalarFunction<P>> FunctionConstructors<P>::coordinates(SizeType as) {
    List<ScalarFunction<P>> r; r.reserve(as);
    for(SizeType j=0; j!=as; ++j) { r.append(coordinate(as,j)); }
    return std::move(r);
}

template<class P> VectorFunction<P> FunctionConstructors<P>::identity(SizeType n) {
    ScalarFunction<P> z=ScalarFunction<P,BoxDomainType>::zero(n);
    VectorOfScalarFunction<P>* res = new VectorOfScalarFunction<P>(n,n);
    for(SizeType i=0; i!=n; ++i) {
        res->_vec[i]=ScalarFunction<P>::coordinate(n,i);
    }
    return VectorFunction<P>(res);
}

template<class P> VectorFunction<P,BoxDomainType> FunctionConstructors<P>::constant(BoxDomainType dom, Vector<NumericType> c) {
    SizeType n=c.size();
    ScalarFunction<P> z=ScalarFunction<P,BoxDomainType>::zero(dom);
    VectorOfScalarFunction<P>* res = new VectorOfScalarFunction<P>(n,z);
    for(SizeType i=0; i!=n; ++i) {
        res->_vec[i]=ScalarFunction<P>::constant(dom,c[i]);
    }
    return VectorFunction<P>(res);
}

template<class P> VectorFunction<P,IntervalDomainType> FunctionConstructors<P>::constant(IntervalDomainType dom, Vector<NumericType> c) {
    SizeType n=c.size();
    ScalarFunction<P,IntervalDomainType> z=ScalarFunction<P,IntervalDomainType>::zero(dom);
    VectorOfScalarFunction<P,IntervalDomainType>* res = new VectorOfScalarFunction<P,IntervalDomainType>(n,z);
    for(SizeType i=0; i!=n; ++i) {
        res->_vec[i]=ScalarFunction<P,IntervalDomainType>::constant(dom,c[i]);
    }
    return VectorFunction<P,IntervalDomainType>(res);
}


template class FunctionConstructors<ApproximateTag>;
template class FunctionConstructors<ValidatedTag>;
template class FunctionConstructors<EffectiveTag>;

//------------------------ Converting to and from Expression classes to Function classes via Formula -----------------------------------//

template<class X> class Expression;
template<class X> class Variable;
template<class X> class Space;

SizeType dimension(const Space<Real>& spc);

Formula<EffectiveNumber> make_formula(const Expression<Real>& expr, const Variable<Real>& var);
Formula<EffectiveNumber> make_formula(const Expression<Real>& expr, const Space<Real>& spc);
Vector<Formula<EffectiveNumber>> make_formula(const Vector<Expression<Real>>& e, const Space<Real>& spc);


EffectiveScalarUnivariateFunction make_formula_function(SizeOne as, const Formula<EffectiveNumber>& fm) {
    return EffectiveScalarUnivariateFunction(RealDomain(),fm); }
EffectiveScalarFunction make_formula_function(SizeType as, const Formula<EffectiveNumber>& fm) {
    return EffectiveScalarFunction(EuclideanDomain(as),fm); }
EffectiveVectorFunction make_formula_function(SizeType as, const Vector<Formula<EffectiveNumber>>& fm) {
    return EffectiveVectorFunction(EuclideanDomain(as),fm); }

EffectiveScalarUnivariateFunction make_function(const Variable<Real>& var, const Expression<Real>& expr) {
    return make_formula_function(SizeOne(),make_formula(expr,var)); }
EffectiveScalarFunction make_function(const Space<Real>& spc, const Expression<Real>& expr) {
    return make_formula_function(dimension(spc),make_formula(expr,spc)); }
EffectiveVectorFunction make_function(const Space<Real>& spc, const Vector<Expression<Real>>& expr) {
    return make_formula_function(dimension(spc),make_formula(expr,spc)); }

[[deprecated]]
EffectiveScalarFunction make_function(const Expression<Real>& expr, const Space<Real>& spc) {
    return make_function(spc,expr); }

Formula<Real> make_formula(const EffectiveScalarFunction& f) {
    Vector<Formula<Real>> x(f.argument_size());
    for(SizeType i=0; i!=x.size(); ++i) {
        x[i]=Formula<Real>::coordinate(i);
    }
    return f(x);
}

Vector<Formula<Real>> make_formula(const EffectiveVectorFunction& f) {
    const VectorFunctionInterface<EffectiveTag>& fi=f;
    const VectorFormulaFunction<Real>* ff;
    const VectorOfScalarFunction<EffectiveTag>* vf;
    if( (vf=dynamic_cast<const VectorOfScalarFunction<EffectiveTag>*>(&fi)) ) {
        Vector< Formula<Real> > r(vf->result_size());
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=make_formula((*vf)[i]); }
        return r;
    } else if( (ff=dynamic_cast<const VectorFormulaFunction<Real>*>(&fi)) ) {
        return ff->_formulae;
    } else {
        ARIADNE_FAIL_MSG("Cannot compute formula for function "<<f<<"\n");
    }
}

//------------------------ Scalar Function ----------------------------------//

//------------------------ Vector Function ----------------------------------//

template<class P, class D, class C> ScalarFunction<P,D> Function<P,D,C>::get(SizeType i) const {
    ARIADNE_ASSERT((IsSame<ResultSizeType,SizeType>::value));
    const VectorOfFunctionInterface<P,D>* vfp = dynamic_cast<const VectorOfFunctionInterface<P,D>*>(this->raw_pointer());
    if(!vfp) { std::cerr<<"\nCannot get element of "<<*this<<"\n  of type "<<typeid(this->raw_pointer()).name()<<":"<<typeid(this->reference()).name()<<"\n\n"; }
    return ScalarFunction<P,D>(SharedPointer<ScalarFunctionInterface<P,D>>(vfp->_get(i)));
}

template<class P, class D, class C> Void Function<P,D,C>::set(SizeType i, ScalarFunction<P,D> sf) {
    ARIADNE_ASSERT((IsSame<ResultSizeType,SizeType>::value));
    const VectorOfScalarFunction<P,D>& cvf = dynamic_cast<const VectorOfScalarFunction<P,D>&>(this->_ptr.operator*());
    VectorOfScalarFunction<P,D>& vf = const_cast<VectorOfScalarFunction<P,D>&>(cvf);
    vf[i]=sf;
}

/*
template<class P, class D> VectorFunction<P,D>::VectorFunction()
    : Function<P,D,C>(new VectorOfScalarFunction<P,D>(0u,ScalarFunction<P,D>()))
{
}

template<class P, class D> VectorFunction<P,D>::VectorFunction(SizeType rs, SizeType as)
    : Function<P,D,C>(new VectorOfScalarFunction<P,D>(rs,ScalarFunction<P,D>::zero(as)))
{
}

template<class P, class D> VectorFunction<P,D>::VectorFunction(const InitializerList<ScalarFunction<P,D>>& lsf)
    : VectorFunction<P,D>(List<ScalarFunction<P,D>>(lsf)) { }

template<class P, class D> VectorFunction<P,D>::VectorFunction(const List<ScalarFunction<P,D>>& lsf) {
    ARIADNE_ASSERT(lsf.size()>0);
    SizeType as=lsf[0].argument_size();
    VectorOfScalarFunction<P,D>* new_ptr=new VectorOfScalarFunction<P,D>(lsf.size(),as);
    for(SizeType i=0; i!=lsf.size(); ++i) {
        new_ptr->set(i,lsf[i]);
    }
    this->_ptr=std::shared_ptr< const VectorFunctionInterface<P,D> >(new_ptr);
}

template<class P, class D> VectorFunction<P,D>::VectorFunction(const Vector<ScalarFunction<P,D>>& vsf) {
    VectorOfScalarFunction<P,D>* new_ptr=new VectorOfScalarFunction<P,D>(vsf);
    this->_ptr=std::shared_ptr< const VectorFunctionInterface<P,D> >(new_ptr);
}

template<class P, class D> VectorFunction<P,D>::VectorFunction(SizeType as, const List<Formula<X>>& le) {
    ARIADNE_ASSERT(le.size()>0);
    VectorOfScalarFunction<P,D>* new_ptr=new VectorOfScalarFunction<P,D>(le.size(),as);
    for(SizeType i=0; i!=le.size(); ++i) {
        new_ptr->set(i,ScalarFormulaFunction<X>(as,le[i]));
    }
    this->_ptr=std::shared_ptr< const VectorFunctionInterface<P,D> >(new_ptr);
}

template<class P, class D> VectorFunction<P,D>::VectorFunction(SizeType rs, ScalarFunction<P,D> const& sf)
    : Function<P,D,C>(new VectorOfScalarFunction<P,D>(rs,sf))
{
}


template<class P, class D> ScalarFunction<P,D> VectorFunction<P,D>::get(SizeType i) const {
    VectorOfFunctionInterface<P,D> const* vfi=dynamic_cast<VectorOfFunctionInterface<P,D>const*>(this->raw_pointer());
    ARIADNE_PRECONDITION(vfi);
    return vfi->_get(i);
}

template<class P, class D> Void VectorFunction<P,D>::set(SizeType i, ScalarFunction<P,D> const& sf) {
    const VectorOfFunctionInterface<P,D>* cvfi=dynamic_cast<const VectorOfFunctionInterface<P,D>*>(this->_ptr.operator->());
    VectorOfFunctionInterface<P,D>* vfi=const_cast<VectorOfFunctionInterface<P,D>*>(cvfi);
    ARIADNE_PRECONDITION(vfi);
    vfi->_set(i,sf);
}

template class VectorFunction<ApproximateTag>;
template class VectorFunction<ValidatedTag>;
template class VectorFunction<EffectiveTag>;

*/

//------------------------ Instantiate functions -----------------------------------//

template class Function<ApproximateTag,IntervalDomainType,IntervalDomainType>;
template class Function<ApproximateTag,IntervalDomainType,BoxDomainType>;
template class Function<ApproximateTag,BoxDomainType,IntervalDomainType>;
template class Function<ApproximateTag,BoxDomainType,BoxDomainType>;

template class Function<ValidatedTag,IntervalDomainType,IntervalDomainType>;
template class Function<ValidatedTag,IntervalDomainType,BoxDomainType>;
template class Function<ValidatedTag,BoxDomainType,IntervalDomainType>;
template class Function<ValidatedTag,BoxDomainType,BoxDomainType>;

template class Function<EffectiveTag,IntervalDomainType,IntervalDomainType>;
template class Function<EffectiveTag,IntervalDomainType,BoxDomainType>;
template class Function<EffectiveTag,BoxDomainType,IntervalDomainType>;
template class Function<EffectiveTag,BoxDomainType,BoxDomainType>;




//------------------------ Scalar arithmetic operators -----------------------------------//





template<class OP> inline
EffectiveScalarFunction make_unary_function(OP op, const EffectiveScalarFunction& f) {
    return EffectiveScalarFunction(new EffectiveUnaryFunction(op.code(),f)); }

template<class OP> inline
EffectiveScalarFunction make_binary_function(OP op, const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2) {
    return EffectiveScalarFunction(new EffectiveBinaryFunction(op.code(),f1,f2)); }

template<class OP> inline
EffectiveScalarFunction make_binary_function(OP op, const EffectiveScalarFunction& f1, const Int& n2) {
    return EffectiveScalarFunction(new EffectiveGradedFunction(op.code(),f1,n2)); }


EffectiveScalarFunction operator+(const EffectiveScalarFunction& f)
{
    return f;
}

EffectiveScalarFunction operator-(const EffectiveScalarFunction& f)
{
    const EffectiveScalarFormulaFunction* ff=dynamic_cast<const EffectiveScalarFormulaFunction*>(f.raw_pointer());
    if(ff) { return make_formula_function(ff->_argument_size,-ff->_formula); }
    return make_unary_function(Neg(),f);
}

EffectiveScalarFunction operator+(const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2)
{
    const EffectiveScalarFormulaFunction* e1=dynamic_cast<const EffectiveScalarFormulaFunction*>(f1.raw_pointer());
    const EffectiveScalarFormulaFunction* e2=dynamic_cast<const EffectiveScalarFormulaFunction*>(f2.raw_pointer());
    if(e1 && e2 && e1->_argument_size==e2->_argument_size) {
        return make_formula_function(e1->_argument_size,e1->_formula+e2->_formula);
    }
    return make_binary_function(Add(),f1,f2);
}

EffectiveScalarFunction operator-(const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2)
{
    const EffectiveScalarFormulaFunction* e1=dynamic_cast<const EffectiveScalarFormulaFunction*>(f1.raw_pointer());
    const EffectiveScalarFormulaFunction* e2=dynamic_cast<const EffectiveScalarFormulaFunction*>(f2.raw_pointer());
    if(e1 && e2 && e1->_argument_size==e2->_argument_size) {
        return make_formula_function(e1->_argument_size,e1->_formula-e2->_formula);
    }
    return make_binary_function(Sub(),f1,f2);
}

EffectiveScalarFunction operator*(const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2)
{
    const EffectiveScalarFormulaFunction* e1=dynamic_cast<const EffectiveScalarFormulaFunction*>(f1.raw_pointer());
    const EffectiveScalarFormulaFunction* e2=dynamic_cast<const EffectiveScalarFormulaFunction*>(f2.raw_pointer());
    if(e1 && e2 && e1->_argument_size==e2->_argument_size) {
        return make_formula_function(e1->_argument_size,e1->_formula*e2->_formula);
    }
    return make_binary_function(Mul(),f1,f2);
}

EffectiveScalarFunction operator/(const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2)
{
    const EffectiveScalarFormulaFunction* e1=dynamic_cast<const EffectiveScalarFormulaFunction*>(f1.raw_pointer());
    const EffectiveScalarFormulaFunction* e2=dynamic_cast<const EffectiveScalarFormulaFunction*>(f2.raw_pointer());
    if(e1 && e2 && e1->_argument_size==e2->_argument_size) {
        return make_formula_function(e1->_argument_size,e1->_formula/e2->_formula);
    }
    return make_binary_function(Div(),f1,f2);
}


EffectiveScalarFunction operator+(const EffectiveScalarFunction& f1, const EffectiveNumber& s2)
{
    const EffectiveScalarFormulaFunction* e1=dynamic_cast<const EffectiveScalarFormulaFunction*>(f1.raw_pointer());
    if(e1) { return make_formula_function(e1->_argument_size,e1->_formula+s2); }
    return f1+EffectiveScalarFunction::constant(f1.argument_size(),s2);
}

EffectiveScalarFunction operator-(const EffectiveScalarFunction& f1, const EffectiveNumber& s2)
{
    const EffectiveScalarFormulaFunction* e1=dynamic_cast<const EffectiveScalarFormulaFunction*>(f1.raw_pointer());
    if(e1) { return make_formula_function(e1->_argument_size,e1->_formula-s2); }
    return f1-EffectiveScalarFunction::constant(f1.argument_size(),s2);
}

EffectiveScalarFunction operator*(const EffectiveScalarFunction& f1, const EffectiveNumber& s2)
{
    const EffectiveScalarFormulaFunction* e1=dynamic_cast<const EffectiveScalarFormulaFunction*>(f1.raw_pointer());
    if(e1) { return make_formula_function(e1->_argument_size,e1->_formula*s2); }
    return f1*EffectiveScalarFunction::constant(f1.argument_size(),s2);
}

EffectiveScalarFunction operator/(const EffectiveScalarFunction& f1, const EffectiveNumber& s2)
{
    const EffectiveScalarFormulaFunction* e1=dynamic_cast<const EffectiveScalarFormulaFunction*>(f1.raw_pointer());
    if(e1) { return make_formula_function(e1->_argument_size,e1->_formula/s2); }
    return f1/EffectiveScalarFunction::constant(f1.argument_size(),s2);
}

EffectiveScalarFunction operator/(const EffectiveScalarFunction& f1, Int s2)
{
    return f1/EffectiveNumber(s2);
}

EffectiveScalarFunction operator+(const EffectiveNumber& s1, const EffectiveScalarFunction& f2)
{
    return f2+s1;
}

EffectiveScalarFunction operator-(const EffectiveNumber& s1, const EffectiveScalarFunction& f2)
{
    const EffectiveScalarFormulaFunction* e2=dynamic_cast<const EffectiveScalarFormulaFunction*>(f2.raw_pointer());
    if(e2) { return make_formula_function(e2->_argument_size,s1-e2->_formula); }
    return EffectiveScalarFunction::constant(f2.argument_size(),s1)-f2;
}

EffectiveScalarFunction operator*(const EffectiveNumber& s1, const EffectiveScalarFunction& f2)
{
    return f2*s1;
}


EffectiveScalarFunction operator/(const EffectiveNumber& s1, const EffectiveScalarFunction& f2)
{
    const EffectiveScalarFormulaFunction* e2=dynamic_cast<const EffectiveScalarFormulaFunction*>(f2.raw_pointer());
    if(e2) { return make_formula_function(e2->_argument_size,s1/e2->_formula); }
    return EffectiveScalarFunction::constant(f2.argument_size(),s1)/f2;
}

EffectiveScalarFunction pow(const EffectiveScalarFunction& f, SizeType m)
{
    const EffectiveScalarFormulaFunction* e=dynamic_cast<const EffectiveScalarFormulaFunction*>(f.raw_pointer());
    if(e) { return make_formula_function(e->_argument_size,pow(e->_formula,m)); }
    return make_binary_function(Pow(),f,m);
}


EffectiveScalarFunction pow(const EffectiveScalarFunction& f, Int n)
{
    const EffectiveScalarFormulaFunction* e=dynamic_cast<const EffectiveScalarFormulaFunction*>(f.raw_pointer());
    if(e) { return make_formula_function(e->_argument_size, pow(e->_formula,n)); }
    return make_binary_function(Pow(),f,n);
}

EffectiveScalarFunction neg(const EffectiveScalarFunction& f) {
    return make_unary_function(Neg(),f); }

EffectiveScalarFunction rec(const EffectiveScalarFunction& f) {
    return make_unary_function(Rec(),f); }

EffectiveScalarFunction sqr(const EffectiveScalarFunction& f) {
    return make_unary_function(Sqr(),f); }

EffectiveScalarFunction sqrt(const EffectiveScalarFunction& f) {
    return make_unary_function(Sqrt(),f); }

EffectiveScalarFunction exp(const EffectiveScalarFunction& f) {
    return make_unary_function(Exp(),f); }

EffectiveScalarFunction log(const EffectiveScalarFunction& f) {
    return make_unary_function(Log(),f); }

EffectiveScalarFunction sin(const EffectiveScalarFunction& f) {
    return make_unary_function(Sin(),f); }

EffectiveScalarFunction cos(const EffectiveScalarFunction& f) {
    return make_unary_function(Cos(),f); }

EffectiveScalarFunction tan(const EffectiveScalarFunction& f) {
    return make_unary_function(Tan(),f); }

EffectiveScalarFunction atan(const EffectiveScalarFunction& f) {
    return make_unary_function(Atan(),f); }

// Deprecated
EffectiveScalarFunction operator+(const Int& s1, const EffectiveScalarFunction& f2) { return EffectiveNumericType(s1)+f2; }
EffectiveScalarFunction operator+(const EffectiveScalarFunction& f1, const Int& s2) { return f1+EffectiveNumericType(s2); }
EffectiveScalarFunction operator-(const Int& s1, const EffectiveScalarFunction& f2) { return EffectiveNumericType(s1)-f2; }
EffectiveScalarFunction operator-(const EffectiveScalarFunction& f1, const Int& s2) { return f1-EffectiveNumericType(s2); }
EffectiveScalarFunction operator*(const Int& s1, const EffectiveScalarFunction& f2) { return EffectiveNumericType(s1)*f2; }
EffectiveScalarFunction operator*(const EffectiveScalarFunction& f1, const Int& s2) { return f1*EffectiveNumericType(s2); }
EffectiveScalarFunction operator/(const Int& s1, const EffectiveScalarFunction& f2) { return EffectiveNumericType(s1)/f2; }
EffectiveScalarFunction operator/(const EffectiveScalarFunction& f1, const Int& s2) { return f1/EffectiveNumericType(s2); }

ValidatedScalarFunction operator+(const Int& s1, const ValidatedScalarFunction& f2) { return ValidatedNumber(s1)+f2; }
ValidatedScalarFunction operator+(const ValidatedScalarFunction& f1, const Int& s2) { return f1+ValidatedNumber(s2); }
ValidatedScalarFunction operator-(const Int& s1, const ValidatedScalarFunction& f2) { return ValidatedNumber(s1)-f2; }
ValidatedScalarFunction operator-(const ValidatedScalarFunction& f1, const Int& s2) { return f1-ValidatedNumber(s2); }
ValidatedScalarFunction operator*(const Int& s1, const ValidatedScalarFunction& f2) { return ValidatedNumber(s1)*f2; }
ValidatedScalarFunction operator*(const ValidatedScalarFunction& f1, const Int& s2) { return f1*ValidatedNumber(s2); }
ValidatedScalarFunction operator/(const Int& s1, const ValidatedScalarFunction& f2) { return ValidatedNumber(s1)/f2; }
ValidatedScalarFunction operator/(const ValidatedScalarFunction& f1, const Int& s2) { return f1/ValidatedNumber(s2); }







//------------------------ Vector function operators -------------------------------//

EffectiveVectorFunction operator*(const EffectiveScalarFunction& f, const Vector<EffectiveNumber>& e) {
    for(SizeType i=0; i!=e.size(); ++i) { ARIADNE_ASSERT(decide(e[i]==EffectiveNumber(0)) || decide(e[i]==EffectiveNumber(1))); }
    VectorFunction<EffectiveTag> r(e.size(),f.domain());
    for(SizeType i=0; i!=e.size(); ++i) {
        if(decide(e[i]==EffectiveNumber(1))) { r.set(i,f); }
    }
    return r;
}

EffectiveVectorFunction operator+(const EffectiveVectorFunction& f1, const EffectiveVectorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    EffectiveVectorFunction r(f1.result_size(),f1.domain());
    for(SizeType i=0; i!=r.result_size(); ++i) {
        r.set(i,f1[i]+f2[i]);
    }
    return r;
}

EffectiveVectorFunction operator-(const EffectiveVectorFunction& f1, const EffectiveVectorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    EffectiveVectorFunction r(f1.result_size(),f1.domain());
    for(SizeType i=0; i!=r.result_size(); ++i) {
        r.set(i,f1[i]-f2[i]);
    }
    return r;
}

EffectiveVectorFunction operator*(const EffectiveVectorFunction& vf, const EffectiveScalarFunction& sf) {
    ARIADNE_ASSERT(vf.argument_size()==sf.argument_size());
    EffectiveVectorFunction r(vf.result_size(),vf.domain());
    for(SizeType i=0; i!=r.result_size(); ++i) {
        r.set(i,vf[i]*sf);
    }
    return r;
}

EffectiveVectorFunction operator*(const EffectiveScalarFunction& sf, const EffectiveVectorFunction& vf) {
    ARIADNE_ASSERT(sf.argument_size()==vf.argument_size());
    EffectiveVectorFunction r(vf.result_size(),vf.domain());
    for(SizeType i=0; i!=r.result_size(); ++i) {
        r.set(i,sf*vf[i]);
    }
    return r;
}

EffectiveVectorFunction operator*(const EffectiveNumber& c, const EffectiveVectorFunction& vf) {
    EffectiveVectorFunction r(vf.result_size(),vf.domain());
    for(SizeType i=0; i!=r.result_size(); ++i) {
        r.set(i,c*vf[i]);
    }
    return r;
}



EffectiveVectorFunction join(const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2) {
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    EffectiveVectorFunction r(2,f1.domain());
    r.set(0,f1);
    r.set(1,f2);
    return r;
}

EffectiveVectorFunction join(const EffectiveVectorFunction& f1, const EffectiveScalarFunction& f2) {
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    EffectiveVectorFunction r(f1.result_size()+1u,f1.domain());
    for(SizeType i=0u; i!=f1.result_size(); ++i) {
        r.set(i,f1.get(i));
    }
    r.set(f1.result_size(),f2);
    return r;
}

EffectiveVectorFunction join(const EffectiveScalarFunction& f1, const EffectiveVectorFunction& f2) {
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    EffectiveVectorFunction r(f2.result_size()+1u,f1.domain());
    r.set(0u,f1);
    for(SizeType i=0u; i!=f2.result_size(); ++i) {
        r.set(i+1u,f2.get(i));
    }
    return r;
}

EffectiveVectorFunction join(const EffectiveVectorFunction& f1, const EffectiveVectorFunction& f2) {
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    EffectiveVectorFunction r(f1.result_size()+f2.result_size(),f1.domain());
    for(SizeType i=0u; i!=f1.result_size(); ++i) {
        r.set(i,f1.get(i));
    }
    for(SizeType i=0u; i!=f2.result_size(); ++i) {
        r.set(i+f1.result_size(),f2.get(i));
    }
    return r;
}

ValidatedVectorFunction operator-(const ValidatedVectorFunction& f1, const ValidatedVectorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    ValidatedVectorFunction r(f1.result_size(),f1.domain());
    for(SizeType i=0; i!=r.result_size(); ++i) {
        r.set(i,f1[i]-f2[i]);
    }
    return r;
}


EffectiveScalarFunction embed(SizeType as1, const EffectiveScalarFunction& f2, SizeType as3) {
    return EffectiveScalarFunction(new ScalarEmbeddedFunction<EffectiveTag>(as1,f2,as3));
}

EffectiveVectorFunction embed(SizeType as1, const EffectiveVectorFunction& f2, SizeType as3) {
    return EffectiveVectorFunction(new VectorEmbeddedFunction<EffectiveTag>(as1,f2,as3));
}

EffectiveScalarFunction compose(const EffectiveScalarFunction& f, const EffectiveVectorFunction& g) {
    ARIADNE_ASSERT(f.argument_size()==g.result_size());
    return EffectiveScalarFunction(new ScalarComposedFunction<EffectiveTag>(f,g));
}

EffectiveVectorFunction compose(const EffectiveVectorFunction& f, const EffectiveVectorFunction& g) {
    ARIADNE_ASSERT(f.argument_size()==g.result_size());
    return EffectiveVectorFunction(new VectorComposedFunction<EffectiveTag>(f,g));
}

EffectiveScalarFunction lie_derivative(const EffectiveScalarFunction& g, const EffectiveVectorFunction& f) {
    ARIADNE_ASSERT_MSG(g.argument_size()==f.result_size(),"f="<<f<<", g="<<g<<"\n");
    ARIADNE_ASSERT_MSG(f.result_size()==f.argument_size(),"f="<<f<<", g="<<g<<"\n");
    ARIADNE_ASSERT_MSG(f.result_size()>0,"f="<<f<<", g="<<g<<"\n");

    try {
        EffectiveScalarFunction r=g.derivative(0)*f[0];
        for(SizeType i=1; i!=g.argument_size(); ++i) {
            r=r+g.derivative(i)*f[i];
        }
        return r;
    }
    catch(...) {
        ARIADNE_FAIL_MSG("Failed to compute Lie derivative of "<<g<<" under vector field "<<f<<"\n");
    }
}



//------------------------ Validated function operators -------------------------------//


namespace {

template<class OP> ValidatedScalarFunction apply(OP op, ValidatedScalarFunction const& f) {
    std::shared_ptr<ValidatedScalarFunctionModel64Interface const>
        fp=std::dynamic_pointer_cast<ValidatedScalarFunctionModel64Interface const>(f.managed_pointer());
    if(fp) {
        ValidatedScalarFunctionModel64 fm(fp); return op(fm);
    }
    return ValidatedScalarFunction(new UnaryFunction<ValidatedTag>(op.code(),f));
}

template<class OP> ValidatedScalarFunction apply(OP op, ValidatedScalarFunction const& f1, ValidatedScalarFunction const& f2) {
    std::shared_ptr<ValidatedScalarFunctionModel64Interface const>
        f1p=std::dynamic_pointer_cast<ValidatedScalarFunctionModel64Interface const>(f1.managed_pointer());
    std::shared_ptr<ValidatedScalarFunctionModel64Interface const>
        f2p=std::dynamic_pointer_cast<ValidatedScalarFunctionModel64Interface const>(f2.managed_pointer());
    if(f1p && f2p) {
        ValidatedScalarFunctionModel64 f1m(f1p); ValidatedScalarFunctionModel64 f2m(f2p); return op(f1m,f2m);
    } else if(f1p) {
        ValidatedScalarFunctionModel64 f1m(f1p); return op(f1m,factory(f1m).create(f2));
    } else if(f2p) {
        ValidatedScalarFunctionModel64 f2m(f2p); return op(factory(f2m).create(f1),f2m);
    } else {
        return ValidatedScalarFunction(new BinaryFunction<ValidatedTag>(op.code(),f1,f2));
    }
}

template<class OP> ValidatedScalarFunction apply(OP op, ValidatedScalarFunction const& f1, ValidatedNumber const& c2) {
    std::shared_ptr<ValidatedScalarFunctionModel64Interface const>
        f1p=std::dynamic_pointer_cast<ValidatedScalarFunctionModel64Interface const>(f1.managed_pointer());
    if(f1p) {
        ValidatedScalarFunctionModel64 f1m=f1p; return op(f1,c2);
    } else {
        return ValidatedScalarFunction(new BinaryFunction<ValidatedTag>(op.code(),f1,f1.create_constant(c2)));
    }
}

template<class OP> ValidatedScalarFunction apply(OP op, ValidatedNumber const& c1, ValidatedScalarFunction const& f2) {
    std::shared_ptr<ValidatedScalarFunctionModel64Interface const>
        f2p=std::dynamic_pointer_cast<ValidatedScalarFunctionModel64Interface const>(f2.managed_pointer());
    if(f2p) {
        ValidatedScalarFunctionModel64 f2m=f2p; return op(c1,f2m);
    } else {
        return ValidatedScalarFunction(new BinaryFunction<ValidatedTag>(op.code(),f2.create_constant(c1),f2));
    }
}

ValidatedScalarFunction apply(Pow op, ValidatedScalarFunction const& f, Int n) {
    std::shared_ptr<ValidatedScalarFunctionModel64Interface const>
        fp=std::dynamic_pointer_cast<ValidatedScalarFunctionModel64Interface const>(f.managed_pointer());
    if(fp) {
        ValidatedScalarFunctionModel64 fm=fp; return op(fm,n);
    } else {
        return ValidatedScalarFunction(new GradedFunction<ValidatedTag>(op.code(),f,n));
    }
}

} // namespace



ValidatedScalarFunction operator+(ValidatedScalarFunction const& f) {
    return apply(Pos(),f);
}

ValidatedScalarFunction operator-(ValidatedScalarFunction const& f) {
    return apply(Neg(),f);
}

ValidatedScalarFunction operator+(ValidatedScalarFunction const& f1, ValidatedScalarFunction const& f2) {
    return apply(Add(),f1,f2);
}

ValidatedScalarFunction operator-(ValidatedScalarFunction const& f1, ValidatedScalarFunction const& f2) {
    return apply(Sub(),f1,f2);
}

ValidatedScalarFunction operator*(ValidatedScalarFunction const& f1, ValidatedScalarFunction const& f2) {
    return apply(Mul(),f1,f2);
}

ValidatedScalarFunction operator/(ValidatedScalarFunction const& f1, ValidatedScalarFunction const& f2) {
    return apply(Div(),f1,f2);
}

ValidatedScalarFunction operator+(ValidatedScalarFunction const& f1, ValidatedNumber const& c2) {
    return apply(Add(),f1,c2);
}

ValidatedScalarFunction operator-(ValidatedScalarFunction const& f1, ValidatedNumber const& c2) {
    return apply(Sub(),f1,c2);
}

ValidatedScalarFunction operator*(ValidatedScalarFunction const& f1, ValidatedNumber const& c2) {
    return apply(Mul(),f1,c2);
}

ValidatedScalarFunction operator/(ValidatedScalarFunction const& f1, ValidatedNumber const& c2) {
    return apply(Div(),f1,c2);
}

ValidatedScalarFunction operator+(ValidatedNumber const& c1, ValidatedScalarFunction const& f2) {
    return apply(Add(),c1,f2);
}

ValidatedScalarFunction operator-(ValidatedNumber const& c1, ValidatedScalarFunction const& f2) {
    return apply(Sub(),c1,f2);
}

ValidatedScalarFunction operator*(ValidatedNumber const& c1, ValidatedScalarFunction const& f2) {
    return apply(Mul(),c1,f2);
}

ValidatedScalarFunction operator/(ValidatedNumber const& c1, ValidatedScalarFunction const& f2) {
    return apply(Div(),c1,f2);
}


ValidatedScalarFunction& operator+=(ValidatedScalarFunction& f1, const ValidatedScalarFunction& f2) {
    return f1=f1+f2;
}

ValidatedScalarFunction& operator-=(ValidatedScalarFunction& f1, const ValidatedScalarFunction& f2) {
    return f1=f1-f2;
}

ValidatedScalarFunction& operator+=(ValidatedScalarFunction& f1, const ValidatedNumber& c2) {
    return f1=f1+c2;
}

ValidatedScalarFunction& operator-=(ValidatedScalarFunction& f1, const ValidatedNumber& c2) {
    return f1=f1-c2;
}

ValidatedScalarFunction& operator*=(ValidatedScalarFunction& f1, const ValidatedNumber& c2) {
    return f1=f1*c2;
}

ValidatedScalarFunction& operator/=(ValidatedScalarFunction& f1, const ValidatedNumber& c2) {
    return f1=f1/c2;
}


ValidatedScalarFunction add(ValidatedScalarFunction const& f1, const ValidatedNumber& c2) {
    return apply(Add(),f1,c2);
}

ValidatedScalarFunction sub(ValidatedScalarFunction const& f1, const ValidatedNumber& c2) {
    return apply(Sub(),f1,c2);
}

ValidatedScalarFunction mul(ValidatedScalarFunction const& f1, const ValidatedNumber& c2) {
    return apply(Mul(),f1,c2);
}

ValidatedScalarFunction div(ValidatedScalarFunction const& f1, const ValidatedNumber& c2) {
    return apply(Div(),f1,c2);
}

ValidatedScalarFunction pos(ValidatedScalarFunction const& f) {
    return apply(Pos(),f);
}

ValidatedScalarFunction neg(ValidatedScalarFunction const& f) {
    return apply(Neg(),f);
}

ValidatedScalarFunction pow(ValidatedScalarFunction const& f, Int n) {
    return apply(Pow(),f,n);
}

ValidatedScalarFunction sqr(ValidatedScalarFunction const& f) {
    return apply(Sqr(),f);
}

ValidatedScalarFunction rec(ValidatedScalarFunction const& f) {
    return apply(Rec(),f);
}

ValidatedScalarFunction sqrt(ValidatedScalarFunction const& f) {
    return apply(Sqrt(),f);
}

ValidatedScalarFunction exp(ValidatedScalarFunction const& f) {
    return apply(Exp(),f);
}

ValidatedScalarFunction log(ValidatedScalarFunction const& f) {
    return apply(Log(),f);
}

ValidatedScalarFunction sin(ValidatedScalarFunction const& f) {
    return apply(Sin(),f);
}

ValidatedScalarFunction cos(ValidatedScalarFunction const& f) {
    return apply(Cos(),f);
}

ValidatedScalarFunction tan(ValidatedScalarFunction const& f) {
    return apply(Tan(),f);
}

ValidatedScalarFunction atan(ValidatedScalarFunction const& f) {
    return apply(Atan(),f);
}


/*
ValidatedVectorFunction operator-(ValidatedVectorFunction const& f1, ValidatedVectorFunction const& f2) {
    std::shared_ptr<ValidatedVectorFunctionModel64Interface const> f1p=std::dynamic_pointer_cast<ValidatedVectorFunctionModel64Interface const>(f1.managed_pointer());
    std::shared_ptr<ValidatedVectorFunctionModel64Interface const> f2p=std::dynamic_pointer_cast<ValidatedVectorFunctionModel64Interface const>(f2.managed_pointer());
    if(f1p && f2p) {
        return ValidatedVectorFunctionModel64(*f1p) - ValidatedVectorFunctionModel64(*f2p);
    } else if(f1p) {
        return ValidatedVectorFunctionModel64(*f1p) - f2.reference();
    } else if(f2p) {
        return f1.reference() - ValidatedVectorFunctionModel64(*f2p);
    } else {
        VectorOfScalarFunction<ValidatedTag> r(f1.result_size(),ValidatedScalarFunction(f1.argument_size()));
        for(SizeType i=0; i!=r.result_size(); ++i) {
            r[i]=f1[i]-f2[i];
        }
        return r;
    }
}
*/

ValidatedScalarFunction compose(const ValidatedScalarFunction& f, const ValidatedVectorFunction& g) {
    ARIADNE_ASSERT(f.argument_size()==g.result_size());
    std::shared_ptr<ValidatedVectorFunctionModel64Interface const>
        gp=std::dynamic_pointer_cast<ValidatedVectorFunctionModel64Interface const>(g.managed_pointer());
    if(gp) {
        return compose(f,ValidatedVectorFunctionModel64(gp->_clone()));
    } else {
        return ValidatedScalarFunction(new ScalarComposedFunction<ValidatedTag>(f,g));
    }
}

ValidatedVectorFunction compose(const ValidatedVectorFunction& f, const ValidatedVectorFunction& g) {
    ARIADNE_ASSERT(f.argument_size()==g.result_size());
    std::shared_ptr<ValidatedVectorFunctionModel64Interface const>
        gp=std::dynamic_pointer_cast<ValidatedVectorFunctionModel64Interface const>(g.managed_pointer());
    if(gp) {
        return compose(f,ValidatedVectorFunctionModel64(gp->_clone()));
    } else {
        return ValidatedVectorFunction(new VectorComposedFunction<ValidatedTag>(f,g));
    }
}

ValidatedVectorFunction join(ValidatedVectorFunction const& f1, const ValidatedVectorFunction& f2) {
    std::shared_ptr<ValidatedVectorFunctionModel64Interface const>
        f1p=std::dynamic_pointer_cast<ValidatedVectorFunctionModel64Interface const>(f1.managed_pointer());
    std::shared_ptr<ValidatedVectorFunctionModel64Interface const>
        f2p=std::dynamic_pointer_cast<ValidatedVectorFunctionModel64Interface const>(f2.managed_pointer());
    if(f1p && f2p) {
        ValidatedVectorFunctionModel64 f1m(f1p); ValidatedVectorFunctionModel64 f2m(f2p); return join(f1m,f2m);
    } else if(f1p) {
        ValidatedVectorFunctionModel64 f1m(f1p); ValidatedVectorFunctionModel64 f2m=factory(f1m).create(f2); return join(f1m,f2m);
    } else if(f2p) {
        ValidatedVectorFunctionModel64 f2m(f2p); ValidatedVectorFunctionModel64 f1m=factory(f2m).create(f1); return join(f1m,f2m);
    } else {
        return ValidatedVectorFunction(new JoinedFunction<ValidatedTag>(f1,f2));
    }
}

ValidatedVectorFunction join(ValidatedVectorFunction const& f1, const ValidatedScalarFunction& f2) {
    std::shared_ptr<ValidatedVectorFunctionModel64Interface const>
        f1p=std::dynamic_pointer_cast<ValidatedVectorFunctionModel64Interface const>(f1.managed_pointer());
    std::shared_ptr<ValidatedScalarFunctionModel64Interface const>
        f2p=std::dynamic_pointer_cast<ValidatedScalarFunctionModel64Interface const>(f2.managed_pointer());
    if(f1p && f2p) {
        ValidatedVectorFunctionModel64 f1m(f1p); ValidatedScalarFunctionModel64 f2m(f2p); return join(f1m,f2m);
    } else if(f1p) {
        ValidatedVectorFunctionModel64 f1m(f1p); ValidatedScalarFunctionModel64 f2m=factory(f1m).create(f2); return join(f1m,f2m);
    } else if(f2p) {
        ValidatedScalarFunctionModel64 f2m(f2p); ValidatedVectorFunctionModel64 f1m=factory(f2m).create(f1); return join(f1m,f2m);
    } else {
        VectorOfScalarFunction<ValidatedTag> r(f1.result_size()+1u,f1.domain());
        for(SizeType i=0; i!=f1.result_size(); ++i) { r[i]=f1[i]; }
        r[f1.result_size()]=f2;
        return r;
    }
}

ValidatedVectorFunction join(ValidatedScalarFunction const& f1, const ValidatedVectorFunction& f2) {
    std::shared_ptr<ValidatedScalarFunctionModel64Interface const>
        f1p=std::dynamic_pointer_cast<ValidatedScalarFunctionModel64Interface const>(f1.managed_pointer());
    std::shared_ptr<ValidatedVectorFunctionModel64Interface const>
        f2p=std::dynamic_pointer_cast<ValidatedVectorFunctionModel64Interface const>(f2.managed_pointer());
    if(f1p && f2p) {
        ValidatedScalarFunctionModel64 f1m(f1p); ValidatedVectorFunctionModel64 f2m(f2p); return join(f1m,f2m);
    } else if(f1p) {
        ValidatedScalarFunctionModel64 f1m(f1p); ValidatedVectorFunctionModel64 f2m=factory(f1m).create(f2); return join(f1m,f2m);
    } else if(f2p) {
        ValidatedVectorFunctionModel64 f2m(f2p); ValidatedScalarFunctionModel64 f1m=factory(f2m).create(f1); return join(f1m,f2m);
    } else {
        VectorOfScalarFunction<ValidatedTag> r(f1.result_size()+1u,f1.domain());
        r[0u]=f1;
        for(SizeType i=0; i!=f1.result_size(); ++i) { r[i+1]=f2[i]; }
        return r;
    }
}

ValidatedVectorFunction join(ValidatedScalarFunction const& f1, const ValidatedScalarFunction& f2) {
    std::shared_ptr<ValidatedScalarFunctionModel64Interface const>
        f1p=std::dynamic_pointer_cast<ValidatedScalarFunctionModel64Interface const>(f1.managed_pointer());
    std::shared_ptr<ValidatedScalarFunctionModel64Interface const>
        f2p=std::dynamic_pointer_cast<ValidatedScalarFunctionModel64Interface const>(f2.managed_pointer());
    if(f1p && f2p) {
        ValidatedScalarFunctionModel64 f1m(f1p); ValidatedScalarFunctionModel64 f2m(f2p); return join(f1m,f2m);
    } else if(f1p) {
        ValidatedScalarFunctionModel64 f1m(f1p); ValidatedScalarFunctionModel64 f2m=factory(f1m).create(f2); return join(f1m,f2m);
    } else if(f2p) {
        ValidatedScalarFunctionModel64 f2m(f2p); ValidatedScalarFunctionModel64 f1m=factory(f2m).create(f1); return join(f1m,f2m);
    } else {
        VectorOfScalarFunction<ValidatedTag> r(2u,f1.domain());
        r[0u]=f1;
        r[1u]=f2;
        return r;
    }
}


UpperIntervalType evaluate_range(ScalarFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<UpperIntervalType>(f(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
Vector<UpperIntervalType> evaluate_range(VectorFunction<ValidatedTag>const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Vector<UpperIntervalType>>(f(reinterpret_cast<Vector<ValidatedNumericType>const&>(x))); }
Vector<Differential<UpperIntervalType>> derivative_range(VectorFunction<ValidatedTag>const& f, const Vector<Differential<UpperIntervalType>>& x) {
    return static_cast<Vector<Differential<UpperIntervalType>>>(f(reinterpret_cast<Vector<Differential<ValidatedNumericType>>const&>(x))); }
Covector<UpperIntervalType> gradient_range(ValidatedScalarFunction const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Covector<UpperIntervalType>>(static_cast<Covector<ValidatedNumericType>>(gradient(f,reinterpret_cast<Vector<ValidatedNumericType>const&>(x)))); }
Matrix<UpperIntervalType> jacobian_range(ValidatedVectorFunction const& f, const Vector<UpperIntervalType>& x) {
    return static_cast<Matrix<UpperIntervalType>>(static_cast<Matrix<ValidatedNumericType>>(jacobian(f,reinterpret_cast<Vector<ValidatedNumericType>const&>(x)))); }


} // namespace Ariadne
