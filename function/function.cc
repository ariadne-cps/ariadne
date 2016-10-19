/***************************************************************************
 *            function.cc
 *
 *  Copyright 2008-12  Pieter Collins
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

#include "numeric/numeric.h"

#include "algebra/differential.h"
#include "numeric/operators.h"
#include "expression/variables.h"
#include "expression/expression.h"
#include "expression/space.h"
#include "expression/formula.h"
#include "algebra/algebra.h"
#include "function/taylor_model.h"

#include "function/function.h"

#include "function/function_mixin.h"
#include "function_mixin.tcc"

#include "function/function_model.h"

#include "function/symbolic_function.h"
#include "function/taylor_function.h"


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

// Functions for converting from Expression classes to Function classes via Formula
template<class X> class Expression;
template<class X> class Variable;
template<class X> class Space;
SizeType dimension(const Space<Real>& spc);
SizeType len(const List< Variable<Real> >& vars);
Formula<Real> formula(const Expression<Real>& expr, const Space<Real>& spc);
Formula<Real> formula(const Expression<Real>& expr, const List< Variable<Real> >& vars);

EffectiveScalarFunction make_function(const Expression<Real>& expr, const Space<Real>& spc) {
    return EffectiveScalarFunction(EuclideanDomain(dimension(spc)),formula(expr,spc)); }
EffectiveScalarFunction make_function(const Expression<Real>& expr, const List< Variable<Real> >& vars) {
    return EffectiveScalarFunction(EuclideanDomain(len(vars)),formula(expr,vars)); }

OutputStream& operator<<(OutputStream& os, Expression<Real> const&);

EffectiveScalarUnivariateFunction make_function(const Variable<Real>& var, const Expression<Real>& expr) {
    Space<Real> spc={var}; return EffectiveScalarUnivariateFunction(RealDomain(),formula(expr,spc)); }
EffectiveScalarFunction make_function(const Space<Real>& spc, const Expression<Real>& expr) {
    return EffectiveScalarFunction(EuclideanDomain(dimension(spc)),formula(expr,spc)); }
EffectiveVectorFunction make_function(const Space<Real>& spc, const Vector<Expression<Real>>& expr) {
    List<EffectiveScalarFunction> lsf;
    for(SizeType i=0; i!=expr.size(); ++i) {
        lsf.append(make_function(spc,expr[i]));
    }
    return EffectiveVectorFunction(lsf);
}


//------------------------ Vector of Scalar functions  -----------------------------------//


template<class P, class D=BoxDomain> class NonResizableScalarFunction : public ScalarFunction<P,D> {
  public:
    NonResizableScalarFunction<P,D>& operator=(const ScalarFunction<P,D>& f) {
        ARIADNE_ASSERT(this->argument_size()==f.argument_size());
        this->ScalarFunction<P,D>::operator=(f);
        return *this;
    }
};

template<class P, class D=BoxDomain>
struct VectorOfScalarFunction
    : VectorFunctionMixin<VectorOfScalarFunction<P,D>,P,D>
    , public virtual VectorOfFunctionInterface<P,D>
{
    typedef D DomainType;
    VectorOfScalarFunction(SizeType rs, SizeType as)
        : VectorOfScalarFunction(rs, BoxDomain(as,IntervalDomain(-inf,+inf))) { }
    VectorOfScalarFunction(SizeType rs, DomainType dom)
        : _dom(dom), _vec(rs,ScalarFunction<P,D>(dom)) { }
    VectorOfScalarFunction(SizeType rs, const ScalarFunction<P,D>& f)
        : _dom(f.domain()), _vec(rs,f) { }
    VectorOfScalarFunction(const Vector<ScalarFunction<P,D>>& vsf)
        : _dom(vsf.zero_element().domain()), _vec(vsf) { }

    Void set(SizeType i, ScalarFunction<P,D> f) {
        if(this->argument_size()==0u) { this->_dom=f.domain(); }
        ARIADNE_ASSERT(f.argument_size()==this->argument_size());
        this->_vec[i]=f; }
    ScalarFunction<P,D> get(SizeType i) const {
        return this->_vec[i]; }

    virtual SizeType result_size() const final {
        return _vec.size(); }
    virtual SizeType argument_size() const final {
        return dimension(_dom); }
    virtual DomainType const domain() const final {
        return _dom; }

    virtual ScalarFunctionInterface<P,D>* _get(SizeType i) const final {
        return this->_vec[i].raw_pointer()->_clone(); }
    virtual Void _set(SizeType i, const ScalarFunctionInterface<P,D>* sf) final {
        this->_vec[i]=ScalarFunction<P,D>(sf->_clone()); }
    virtual VectorFunctionInterface<P,D>* _derivative(SizeType i) const {
        ARIADNE_NOT_IMPLEMENTED; }

    const ScalarFunction<P,D> operator[](SizeType i) const {
        return this->_vec[i]; }

    NonResizableScalarFunction<P,D>& operator[](SizeType i) {
        return static_cast<NonResizableScalarFunction<P,D>&>(this->_vec[i]); }

    virtual OutputStream& write(OutputStream& os) const {
        os << "[";
        for(SizeType i=0; i!=this->_vec.size(); ++i) {
            if(i!=0) { os << ","; }
            this->_vec[i].raw_pointer()->write(os); }
        return os << "]"; }

    virtual OutputStream& repr(OutputStream& os) const {
        //os << "VoSF[R" << this->argument_size() << "->R" << this->result_size() << "]";
        os << "[";
        for(SizeType i=0; i!=this->_vec.size(); ++i) {
            if(i!=0) { os << ","; }
            this->_vec[i].raw_pointer()->repr(os); }
        return os << "]"; }

    template<class X> inline Void _compute(Vector<X>& r, const ElementType<D,X>& x) const {
        r=Vector<X>(this->_vec.size(),zero_element(x));
        for(SizeType i=0; i!=r.size(); ++i) {
            r[i]=_vec[i].evaluate(x); } }

    DomainType _dom;
    Vector< ScalarFunction<P,D> > _vec;

};


template<class X>
struct FunctionElement
    : ScalarFunctionMixin<FunctionElement<X>,X>
{
    FunctionElement(const VectorFunction<X>& f, SizeType i)
        : _f(f), _i(i) { ARIADNE_ASSERT(i<f.result_size()); }

    virtual SizeType argument_size() const { return _f.argument_size(); }
    virtual OutputStream& write(OutputStream& os) const { return os<<_f<<"["<<_i<<"]"; }
    virtual ScalarFunctionInterface<X>* _derivative(SizeType j) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class XX> inline Void _compute(XX& r, const Vector<XX>& x) const {
        r=this->_f.evaluate(x)[_i]; }

    VectorFunction<X> _f;
    SizeType _i;
};



//------------------------ Conversion to a formula -----------------------------------//

Formula<Real> formula(const EffectiveScalarFunction& f) {
    const ScalarFunctionInterface<EffectiveTag>* fptr=f.raw_pointer();
    const ScalarFormulaFunction<Real>* ff=dynamic_cast<const ScalarFormulaFunction<Real>*>(fptr);
    if(ff) { return ff->_formula; }
    const RealConstantFunction* cf=dynamic_cast<const RealConstantFunction*>(fptr);
    const EffectiveCoordinateFunction* indf=dynamic_cast<const EffectiveCoordinateFunction*>(fptr);
    const EffectiveUnaryFunction* uf=dynamic_cast<const EffectiveUnaryFunction*>(fptr);
    const EffectiveBinaryFunction* bf=dynamic_cast<const EffectiveBinaryFunction*>(fptr);
    const EffectivePowerFunction* pf=dynamic_cast<const EffectivePowerFunction*>(fptr);
    if(cf) { return Formula<Real>::constant(cf->_value); }
    if(indf) { return Formula<Real>::coordinate(indf->_index); }
    if(uf) { return make_formula(uf->_op,formula(uf->_arg)); }
    if(bf) { return make_formula(bf->_op,formula(bf->_arg1),formula(bf->_arg2)); }
    if(pf) { return make_formula(pf->_op,formula(pf->_arg1),pf->_arg2); }
    ARIADNE_FAIL_MSG("Cannot compute formula for function "<<f<<"\n");
}

Vector< Formula<Real> > formula(const EffectiveVectorFunction& f) {
    const VectorFunctionInterface<EffectiveTag>& fi=f;
    const VectorFormulaFunction<Real>* ff;
    const VectorOfScalarFunction<EffectiveTag>* vf;
    if( (vf=dynamic_cast<const VectorOfScalarFunction<EffectiveTag>*>(&fi)) ) {
        Vector< Formula<Real> > r(vf->result_size());
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=formula((*vf)[i]); }
        return r;
    } else if( (ff=dynamic_cast<const VectorFormulaFunction<Real>*>(&fi)) ) {
        return ff->_formulae;
    } else {
        ARIADNE_FAIL_MSG("Cannot compute formula for function "<<f<<"\n");
    }
}



//------------------------ Function Constructors -----------------------------------//

template<class P> ScalarFunction<P> FunctionConstructors<P>::zero(BoxDomain dom) {
    return ConstantFunction<Y>(dom, Y(0));
}


template<class P> ScalarFunction<P> FunctionConstructors<P>::constant(BoxDomain dom, NumericType c) {
    return ConstantFunction<Y>(dom, c);
}

template<class P> ScalarFunction<P> FunctionConstructors<P>::coordinate(BoxDomain dom, SizeType j) {
    return CoordinateFunction<P>(dom, j);
}

template<class P> List<ScalarFunction<P>> FunctionConstructors<P>::coordinates(BoxDomain dom) {
    List<ScalarFunction<P>> r; r.reserve(dom.dimension());
    for(SizeType j=0; j!=dom.dimension(); ++j) { r.append(coordinate(dom,j)); }
    return std::move(r);
}

template<class P> VectorFunction<P> FunctionConstructors<P>::zeros(SizeType rs, BoxDomain dom) {
    return VectorFunction<P>(new VectorOfScalarFunction<P>(rs,zero(dom)));
}

template<class P> VectorFunction<P> FunctionConstructors<P>::identity(BoxDomain dom) {
    SizeType n=dom.dimension();
    ScalarFunction<P> z=ScalarFunction<P,BoxDomain>::zero(dom);
    VectorOfScalarFunction<P>* res = new VectorOfScalarFunction<P>(n,z);
    for(SizeType i=0; i!=n; ++i) {
        res->_vec[i]=ScalarFunction<P>::coordinate(dom,i);
    }
    return VectorFunction<P>(res);
}


template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::zero(IntervalDomain dom) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::constant(IntervalDomain dom, CanonicalNumericType<P> c) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::coordinate(IntervalDomain dom, SizeType j) {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class P> VectorUnivariateFunction<P> FunctionConstructors<P>::zeros(SizeType rs, IntervalDomain dom) {
    ARIADNE_NOT_IMPLEMENTED;
}


template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::identity(IntervalDomain dom) {
    ARIADNE_NOT_IMPLEMENTED;
}


template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::zero() {
    ARIADNE_NOT_IMPLEMENTED;
}

template<class P> ScalarUnivariateFunction<P> FunctionConstructors<P>::constant(CanonicalNumericType<P> c) {
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

template<class P> ScalarFunction<P> FunctionConstructors<P>::constant(SizeType as, CanonicalNumericType<P> c) {
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
    ScalarFunction<P> z=ScalarFunction<P,BoxDomain>::zero(n);
    VectorOfScalarFunction<P>* res = new VectorOfScalarFunction<P>(n,n);
    for(SizeType i=0; i!=n; ++i) {
        res->_vec[i]=ScalarFunction<P>::coordinate(n,i);
    }
    return VectorFunction<P>(res);
}

template<class P> VectorFunction<P,BoxDomain> FunctionConstructors<P>::constant(BoxDomain dom, Vector<NumericType> c) {
    SizeType n=c.size();
    ScalarFunction<P> z=ScalarFunction<P,BoxDomain>::zero(dom);
    VectorOfScalarFunction<P>* res = new VectorOfScalarFunction<P>(n,z);
    for(SizeType i=0; i!=n; ++i) {
        res->_vec[i]=ScalarFunction<P>::constant(dom,c[i]);
    }
    return VectorFunction<P>(res);
}

template<class P> VectorFunction<P,IntervalDomain> FunctionConstructors<P>::constant(IntervalDomain dom, Vector<NumericType> c) {
    SizeType n=c.size();
    ScalarFunction<P,IntervalDomain> z=ScalarFunction<P,IntervalDomain>::zero(dom);
    VectorOfScalarFunction<P,IntervalDomain>* res = new VectorOfScalarFunction<P,IntervalDomain>(n,z);
    for(SizeType i=0; i!=n; ++i) {
        res->_vec[i]=ScalarFunction<P,IntervalDomain>::constant(dom,c[i]);
    }
    return VectorFunction<P,IntervalDomain>(res);
}


template class FunctionConstructors<ApproximateTag>;
template class FunctionConstructors<ValidatedTag>;
template class FunctionConstructors<EffectiveTag>;

//------------------------ Function ----------------------------------//

namespace {
OutputStream& operator<<(OutputStream& os, SizeOne so) { return os << "1u"; }
OutputStream& operator<<(OutputStream& os, RealDomain const& dom) { return os << "R"; }
OutputStream& operator<<(OutputStream& os, EuclideanDomain const& dom) { return os << "R" << dom.dimension(); }

template<class D, class DD> D make_domain(DD dom);
template<> IntervalDomain make_domain<IntervalDomain,BoxDomain>(BoxDomain dom) { throw std::runtime_error(""); }
template<> BoxDomain make_domain<BoxDomain,IntervalDomain>(IntervalDomain dom) { throw std::runtime_error(""); }
template<> IntervalDomain make_domain<IntervalDomain,IntervalDomain>(IntervalDomain dom) { return dom; }
template<> BoxDomain make_domain<BoxDomain,BoxDomain>(BoxDomain dom) { return dom; }

template<class P, class D, class DD> ScalarFunction<P,D> make_zero_function(SizeOne rs, DD dom) {
    return FunctionConstructors<P>::zero(make_domain<D>(dom)); }
template<class P, class D, class DD> VectorFunction<P,D> make_zero_function(SizeType rs, DD dom) {
    return  FunctionConstructors<P>::zeros(rs,make_domain<D>(dom)); }
}



template<class P, class D, class C> Function<P,D,C>::Function() : _ptr() {
}

template<class P, class D, class C> Function<P,D,C>::Function(EuclideanDomain dom) {
    ResultSizeType rs; BoxDomain bx_dom=dom;
    (*this) = make_zero_function<P,D>(rs,bx_dom);
}

template<class P, class D, class C> Function<P,D,C>::Function(ResultSizeType rs, EuclideanDomain dom) {
    BoxDomain const& bx_dom=dom;
    (*this) = make_zero_function<P,D>(rs,bx_dom);
}

template<class P, class D, class C> Function<P,D,C>::Function(EuclideanDomain dom, Result<Formula<Y>>const& e) {
    BoxDomain const& bx_dom=dom;
    (*this) = Function<P,D,C>(make_domain<D>(bx_dom),e);
}

template<class P, class D, class C> Function<P,D,C>::Function(DomainType dom) {
    ResultSizeType rs; (*this) = make_zero_function<P,D>(rs,dom);
}

template<class P, class D, class C> Function<P,D,C>::Function(ResultSizeType rs, DomainType dom) {
    (*this) = make_zero_function<P,D>(rs,dom);
}

template<class P, class Y> ScalarFunction<P,IntervalDomain> make_formula_function(IntervalDomain dom, Scalar<Formula<Y>> const& e) {
    assert(false);
}

template<class P, class Y> VectorFunction<P,IntervalDomain> make_formula_function(IntervalDomain dom, Vector<Formula<Y>> const& e) {
    assert(false);
}

template<class P, class Y> ScalarFunction<P,BoxDomain> make_formula_function(BoxDomain dom, Scalar<Formula<Y>> const& e) {
    return ScalarFunction<P,BoxDomain>(new ScalarFormulaFunction<Y>(dom.dimension(),e));
}

template<class P, class Y> VectorFunction<P,BoxDomain> make_formula_function(BoxDomain dom, Vector<Formula<Y>> const& e) {
    return VectorFunction<P,BoxDomain>(new VectorFormulaFunction<Y>(dom.dimension(),e));
}

template<class P, class D, class C> Function<P,D,C>::Function(DomainType dom, Result<Formula<Y>>const& e) {
    *this = make_formula_function<P>(dom,e);
}

template<class P, class D, class C> Function<P,D,C>::Function(EuclideanDomain dom, List<Formula<Y>> const& e)
    : Function(dom, checked_construct<Result<Formula<Y>>>(e)) {
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

template<class P, class D, class C> struct MakeVectorFunction;
template<class P, class D> struct MakeVectorFunction<P,D,IntervalDomain> {
    Function<P,D,IntervalDomain> create(Vector<ScalarFunction<P,D>> const& lsf) {
        ARIADNE_FAIL_MSG("Cannot construct scalar function from list."); }
};
template<class P, class D> struct MakeVectorFunction<P,D,BoxDomain> {
    Function<P,D,BoxDomain> create(Vector<ScalarFunction<P,D>> const& lsf) {
        return Function<P,D,BoxDomain>(std::make_shared<VectorOfScalarFunction<P,D>>(lsf)); }
};

template<class P, class D, class C> Function<P,D,C> make_vector_function(Vector<ScalarFunction<P,D>> const& lsf) {
    return MakeVectorFunction<P,D,C>().create(lsf);
}

template<class P, class D, class C> Function<P,D,C>::Function(Vector<ScalarFunction<P,D>> const& vsf)
    : Function<P,D,C>(make_vector_function<P,D,C>(vsf)) {
}

//------------------------ Scalar Function ----------------------------------//

/*
template<class P, class D> ScalarFunction<P,D>::ScalarFunction(DomainType dom)
    : ScalarFunction<P,D>(dom,Formula<X>::zero())
{
}

template<class P, class D> ScalarFunction<P,D>::ScalarFunction(DomainType dom, Formula<X> const& e)
    : Function<P,D,C>(new ScalarFormulaFunction<X>(dom.dimension(),e))
{
}
*/

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

template class Function<ApproximateTag,IntervalDomain,IntervalDomain>;
template class Function<ApproximateTag,IntervalDomain,BoxDomain>;
template class Function<ApproximateTag,BoxDomain,IntervalDomain>;
template class Function<ApproximateTag,BoxDomain,BoxDomain>;

template class Function<ValidatedTag,IntervalDomain,IntervalDomain>;
template class Function<ValidatedTag,IntervalDomain,BoxDomain>;
template class Function<ValidatedTag,BoxDomain,IntervalDomain>;
template class Function<ValidatedTag,BoxDomain,BoxDomain>;

template class Function<EffectiveTag,IntervalDomain,IntervalDomain>;
template class Function<EffectiveTag,IntervalDomain,BoxDomain>;
template class Function<EffectiveTag,BoxDomain,IntervalDomain>;
template class Function<EffectiveTag,BoxDomain,BoxDomain>;




//------------------------ Scalar arithmetic operators -----------------------------------//





template<class OP> inline
EffectiveScalarFunction make_unary_function(OP op, const EffectiveScalarFunction& f) {
    return EffectiveScalarFunction(new EffectiveUnaryFunction(op.code(),f)); }

template<class OP> inline
EffectiveScalarFunction make_binary_function(OP op, const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2) {
    return EffectiveScalarFunction(new EffectiveBinaryFunction(op.code(),f1,f2)); }

template<class OP> inline
EffectiveScalarFunction make_binary_function(OP op, const EffectiveScalarFunction& f1, const Int& n2) {
    return EffectiveScalarFunction(new EffectivePowerFunction(op.code(),f1,n2)); }


EffectiveScalarFunction operator+(const EffectiveScalarFunction& f)
{
    return f;
}

EffectiveScalarFunction operator-(const EffectiveScalarFunction& f)
{
    const RealScalarFormulaFunction* ff=dynamic_cast<const RealScalarFormulaFunction*>(f.raw_pointer());
    if(ff) { return function(ff->_argument_size,-ff->_formula); }
    return make_unary_function(Neg(),f);
}

EffectiveScalarFunction operator+(const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    const RealScalarFormulaFunction* e2=dynamic_cast<const RealScalarFormulaFunction*>(f2.raw_pointer());
    if(e1 && e2 && e1->_argument_size==e2->_argument_size) {
        return function(e1->_argument_size,e1->_formula+e2->_formula);
    }
    return make_binary_function(Add(),f1,f2);
}

EffectiveScalarFunction operator-(const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    const RealScalarFormulaFunction* e2=dynamic_cast<const RealScalarFormulaFunction*>(f2.raw_pointer());
    if(e1 && e2 && e1->_argument_size==e2->_argument_size) {
        return function(e1->_argument_size,e1->_formula-e2->_formula);
    }
    return make_binary_function(Sub(),f1,f2);
}

EffectiveScalarFunction operator*(const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    const RealScalarFormulaFunction* e2=dynamic_cast<const RealScalarFormulaFunction*>(f2.raw_pointer());
    if(e1 && e2 && e1->_argument_size==e2->_argument_size) {
        return function(e1->_argument_size,e1->_formula*e2->_formula);
    }
    return make_binary_function(Mul(),f1,f2);
}

EffectiveScalarFunction operator/(const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    const RealScalarFormulaFunction* e2=dynamic_cast<const RealScalarFormulaFunction*>(f2.raw_pointer());
    if(e1 && e2 && e1->_argument_size==e2->_argument_size) {
        return function(e1->_argument_size,e1->_formula/e2->_formula);
    }
    return make_binary_function(Div(),f1,f2);
}


EffectiveScalarFunction operator+(const EffectiveScalarFunction& f1, const Real& s2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    if(e1) { return function(e1->_argument_size,e1->_formula+s2); }
    return f1+EffectiveScalarFunction::constant(f1.argument_size(),s2);
}

EffectiveScalarFunction operator-(const EffectiveScalarFunction& f1, const Real& s2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    if(e1) { return function(e1->_argument_size,e1->_formula-s2); }
    return f1-EffectiveScalarFunction::constant(f1.argument_size(),s2);
}

EffectiveScalarFunction operator*(const EffectiveScalarFunction& f1, const Real& s2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    if(e1) { return function(e1->_argument_size,e1->_formula*s2); }
    return f1*EffectiveScalarFunction::constant(f1.argument_size(),s2);
}

EffectiveScalarFunction operator/(const EffectiveScalarFunction& f1, const Real& s2)
{
    const RealScalarFormulaFunction* e1=dynamic_cast<const RealScalarFormulaFunction*>(f1.raw_pointer());
    if(e1) { return function(e1->_argument_size,e1->_formula/s2); }
    return f1/EffectiveScalarFunction::constant(f1.argument_size(),s2);
}

EffectiveScalarFunction operator/(const EffectiveScalarFunction& f1, Int s2)
{
    return f1/Real(s2);
}

EffectiveScalarFunction operator+(const Real& s1, const EffectiveScalarFunction& f2)
{
    return f2+s1;
}

EffectiveScalarFunction operator-(const Real& s1, const EffectiveScalarFunction& f2)
{
    const RealScalarFormulaFunction* e2=dynamic_cast<const RealScalarFormulaFunction*>(f2.raw_pointer());
    if(e2) { return function(e2->_argument_size,s1-e2->_formula); }
    return EffectiveScalarFunction::constant(f2.argument_size(),s1)-f2;
}

EffectiveScalarFunction operator*(const Real& s1, const EffectiveScalarFunction& f2)
{
    return f2*s1;
}


EffectiveScalarFunction operator/(const Real& s1, const EffectiveScalarFunction& f2)
{
    const RealScalarFormulaFunction* e2=dynamic_cast<const RealScalarFormulaFunction*>(f2.raw_pointer());
    if(e2) { return function(e2->_argument_size,s1/e2->_formula); }
    return EffectiveScalarFunction::constant(f2.argument_size(),s1)/f2;
}

EffectiveScalarFunction pow(const EffectiveScalarFunction& f, SizeType m)
{
    const RealScalarFormulaFunction* e=dynamic_cast<const RealScalarFormulaFunction*>(f.raw_pointer());
    if(e) { return function(e->_argument_size,pow(e->_formula,m)); }
    return make_binary_function(Pow(),f,m);
}


EffectiveScalarFunction pow(const EffectiveScalarFunction& f, Int n)
{
    const RealScalarFormulaFunction* e=dynamic_cast<const RealScalarFormulaFunction*>(f.raw_pointer());
    if(e) { return function(e->_argument_size, pow(e->_formula,n)); }
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

ValidatedScalarFunction operator+(const Int& s1, const ValidatedScalarFunction& f2) { return ValidatedNumericType(s1)+f2; }
ValidatedScalarFunction operator+(const ValidatedScalarFunction& f1, const Int& s2) { return f1+ValidatedNumericType(s2); }
ValidatedScalarFunction operator-(const Int& s1, const ValidatedScalarFunction& f2) { return ValidatedNumericType(s1)-f2; }
ValidatedScalarFunction operator-(const ValidatedScalarFunction& f1, const Int& s2) { return f1-ValidatedNumericType(s2); }
ValidatedScalarFunction operator*(const Int& s1, const ValidatedScalarFunction& f2) { return ValidatedNumericType(s1)*f2; }
ValidatedScalarFunction operator*(const ValidatedScalarFunction& f1, const Int& s2) { return f1*ValidatedNumericType(s2); }
ValidatedScalarFunction operator/(const Int& s1, const ValidatedScalarFunction& f2) { return ValidatedNumericType(s1)/f2; }
ValidatedScalarFunction operator/(const ValidatedScalarFunction& f1, const Int& s2) { return f1/ValidatedNumericType(s2); }







//------------------------ Vector function operators -------------------------------//

EffectiveVectorFunction operator*(const EffectiveScalarFunction& f, const Vector<Real>& e) {
    for(SizeType i=0; i!=e.size(); ++i) { ARIADNE_ASSERT(decide(e[i]==Real(0)) || decide(e[i]==Real(1))); }
    VectorFunction<EffectiveTag> r(e.size(),f.domain());
    for(SizeType i=0; i!=e.size(); ++i) {
        if(decide(e[i]==Real(1))) { r.set(i,f); }
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

EffectiveVectorFunction operator*(const Real& c, const EffectiveVectorFunction& vf) {
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



//------------------------ ValidatedNumber function operators -------------------------------//


namespace {

template<class OP> ValidatedScalarFunction apply(OP op, ValidatedScalarFunction const& f) {
    std::shared_ptr<ValidatedScalarFunctionModelInterface const>
        fp=std::dynamic_pointer_cast<ValidatedScalarFunctionModelInterface const>(f.managed_pointer());
    if(fp) {
        ValidatedScalarFunctionModel fm(fp); return op(fm);
    }
    return ValidatedScalarFunction(new UnaryFunction<ValidatedTag>(op.code(),f));
}

template<class OP> ValidatedScalarFunction apply(OP op, ValidatedScalarFunction const& f1, ValidatedScalarFunction const& f2) {
    std::shared_ptr<ValidatedScalarFunctionModelInterface const>
        f1p=std::dynamic_pointer_cast<ValidatedScalarFunctionModelInterface const>(f1.managed_pointer());
    std::shared_ptr<ValidatedScalarFunctionModelInterface const>
        f2p=std::dynamic_pointer_cast<ValidatedScalarFunctionModelInterface const>(f2.managed_pointer());
    if(f1p && f2p) {
        ValidatedScalarFunctionModel f1m(f1p); ValidatedScalarFunctionModel f2m(f2p); return op(f1m,f2m);
    } else if(f1p) {
        ValidatedScalarFunctionModel f1m(f1p); return op(f1m,f1m.create(f2));
    } else if(f2p) {
        ValidatedScalarFunctionModel f2m(f2p); return op(f2m.create(f1),f2m);
    } else {
        return ValidatedScalarFunction(new BinaryFunction<ValidatedTag>(op.code(),f1,f2));
    }
}

template<class OP> ValidatedScalarFunction apply(OP op, ValidatedScalarFunction const& f1, ValidatedNumericType const& c2) {
    std::shared_ptr<ValidatedScalarFunctionModelInterface const>
        f1p=std::dynamic_pointer_cast<ValidatedScalarFunctionModelInterface const>(f1.managed_pointer());
    if(f1p) {
        ValidatedScalarFunctionModel f1m=f1p; return op(f1m,c2);
    } else {
        return ValidatedScalarFunction(new BinaryFunction<ValidatedTag>(op.code(),f1,f1.create_constant(c2)));
    }
}

template<class OP> ValidatedScalarFunction apply(OP op, ValidatedNumericType const& c1, ValidatedScalarFunction const& f2) {
    std::shared_ptr<ValidatedScalarFunctionModelInterface const>
        f2p=std::dynamic_pointer_cast<ValidatedScalarFunctionModelInterface const>(f2.managed_pointer());
    if(f2p) {
        ValidatedScalarFunctionModel f2m=f2p; return op(c1,f2m);
    } else {
        return ValidatedScalarFunction(new BinaryFunction<ValidatedTag>(op.code(),f2.create_constant(c1),f2));
    }
}

ValidatedScalarFunction apply(Pow op, ValidatedScalarFunction const& f, Int n) {
    std::shared_ptr<ValidatedScalarFunctionModelInterface const>
        fp=std::dynamic_pointer_cast<ValidatedScalarFunctionModelInterface const>(f.managed_pointer());
    if(fp) {
        ValidatedScalarFunctionModel fm=fp; return op(fm,n);
    } else {
        return ValidatedScalarFunction(new PowerFunction<ValidatedTag>(op.code(),f,n));
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

ValidatedScalarFunction operator+(ValidatedScalarFunction const& f1, ValidatedNumericType const& c2) {
    return apply(Add(),f1,c2);
}

ValidatedScalarFunction operator-(ValidatedScalarFunction const& f1, ValidatedNumericType const& c2) {
    return apply(Sub(),f1,c2);
}

ValidatedScalarFunction operator*(ValidatedScalarFunction const& f1, ValidatedNumericType const& c2) {
    return apply(Mul(),f1,c2);
}

ValidatedScalarFunction operator/(ValidatedScalarFunction const& f1, ValidatedNumericType const& c2) {
    return apply(Div(),f1,c2);
}

ValidatedScalarFunction operator+(ValidatedNumericType const& c1, ValidatedScalarFunction const& f2) {
    return apply(Add(),c1,f2);
}

ValidatedScalarFunction operator-(ValidatedNumericType const& c1, ValidatedScalarFunction const& f2) {
    return apply(Sub(),c1,f2);
}

ValidatedScalarFunction operator*(ValidatedNumericType const& c1, ValidatedScalarFunction const& f2) {
    return apply(Mul(),c1,f2);
}

ValidatedScalarFunction operator/(ValidatedNumericType const& c1, ValidatedScalarFunction const& f2) {
    return apply(Div(),c1,f2);
}


ValidatedScalarFunction& operator+=(ValidatedScalarFunction& f1, const ValidatedScalarFunction& f2) {
    return f1=f1+f2;
}

ValidatedScalarFunction& operator-=(ValidatedScalarFunction& f1, const ValidatedScalarFunction& f2) {
    return f1=f1-f2;
}

ValidatedScalarFunction& operator+=(ValidatedScalarFunction& f1, const ValidatedNumericType& c2) {
    return f1=f1+c2;
}

ValidatedScalarFunction& operator-=(ValidatedScalarFunction& f1, const ValidatedNumericType& c2) {
    return f1=f1-c2;
}

ValidatedScalarFunction& operator*=(ValidatedScalarFunction& f1, const ValidatedNumericType& c2) {
    return f1=f1*c2;
}

ValidatedScalarFunction& operator/=(ValidatedScalarFunction& f1, const ValidatedNumericType& c2) {
    return f1=f1/c2;
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
    std::shared_ptr<ValidatedVectorFunctionModelInterface const> f1p=std::dynamic_pointer_cast<ValidatedVectorFunctionModelInterface const>(f1.managed_pointer());
    std::shared_ptr<ValidatedVectorFunctionModelInterface const> f2p=std::dynamic_pointer_cast<ValidatedVectorFunctionModelInterface const>(f2.managed_pointer());
    if(f1p && f2p) {
        return ValidatedVectorFunctionModel(*f1p) - ValidatedVectorFunctionModel(*f2p);
    } else if(f1p) {
        return ValidatedVectorFunctionModel(*f1p) - f2.reference();
    } else if(f2p) {
        return f1.reference() - ValidatedVectorFunctionModel(*f2p);
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
    std::shared_ptr<ValidatedVectorFunctionModelInterface const>
        gp=std::dynamic_pointer_cast<ValidatedVectorFunctionModelInterface const>(g.managed_pointer());
    if(gp) {
        return compose(f,ValidatedVectorFunctionModel(gp->_clone()));
    } else {
        return ValidatedScalarFunction(new ScalarComposedFunction<ValidatedTag>(f,g));
    }
}

ValidatedVectorFunction compose(const ValidatedVectorFunction& f, const ValidatedVectorFunction& g) {
    ARIADNE_ASSERT(f.argument_size()==g.result_size());
    std::shared_ptr<ValidatedVectorFunctionModelInterface const>
        gp=std::dynamic_pointer_cast<ValidatedVectorFunctionModelInterface const>(g.managed_pointer());
    if(gp) {
        return compose(f,ValidatedVectorFunctionModel(gp->_clone()));
    } else {
        return ValidatedVectorFunction(new VectorComposedFunction<ValidatedTag>(f,g));
    }
}

ValidatedVectorFunction join(ValidatedVectorFunction const& f1, const ValidatedVectorFunction& f2) {
    VectorTaylorFunction const*
        f1p=dynamic_cast<VectorTaylorFunction const*>(f1.raw_pointer());
    VectorTaylorFunction const*
        f2p=dynamic_cast<VectorTaylorFunction const*>(f2.raw_pointer());
    if(f1p && f2p) {
        return join(*f1p,*f2p);
    }
    VectorOfScalarFunction<ValidatedTag> r(f1.result_size()+f2.result_size(),f1.domain());
    for(SizeType i=0; i!=f1.result_size(); ++i) { r[i]=f1[i]; }
    for(SizeType i=0; i!=f2.result_size(); ++i) { r[i+f1.result_size()]=f2[i]; }
    return r;
}

ValidatedVectorFunction join(ValidatedVectorFunction const& f1, const ValidatedScalarFunction& f2) {
    VectorTaylorFunction const*
        f1p=dynamic_cast<VectorTaylorFunction const*>(f1.raw_pointer());
    ScalarTaylorFunction const*
        f2p=dynamic_cast<ScalarTaylorFunction const*>(f2.raw_pointer());
    if(f1p && f2p) {
        return join(*f1p,*f2p);
    }
    VectorOfScalarFunction<ValidatedTag> r(f1.result_size()+1u,f1.domain());
    for(SizeType i=0; i!=f1.result_size(); ++i) { r[i]=f1[i]; }
    r[f1.result_size()]=f2;
    return r;
}




RealExpression evaluate(EffectiveScalarFunction const& f, Vector<RealVariable> const& vars) {
    typedef Algebra<EffectiveNumericType> EffectiveAlgebra;
    EffectiveAlgebra az(RealExpression::constant(0));
    Vector<EffectiveAlgebra> va(vars.size(),az);
    for(SizeType i=0; i!=va.size(); ++i) { va[i]=EffectiveAlgebra(RealExpression(vars[i])); }
    //Vector<EffectiveAlgebra> va=Vector<EffectiveAlgebra>(Vector<RealExpression>(vars));
    EffectiveAlgebra fa=f(va);
    return fa.template extract<RealExpression>();
}





} // namespace Ariadne
