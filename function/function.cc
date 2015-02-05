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
#include "expression/operators.h"
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

template<class T> inline StringType str(const T& t) {
    StringStream ss; ss << t; return ss.str(); }


// Functions for converting from Expression classes to Function classes via Formula
template<class X> class Expression;
template<class X> class Variable;
template<class X> class Space;
SizeType dimension(const Space<Real>& spc);
SizeType len(const List< Variable<Real> >& vars);
Formula<Real> formula(const Expression<Real>& expr, const Space<Real>& spc);
Formula<Real> formula(const Expression<Real>& expr, const List< Variable<Real> >& vars);
EffectiveScalarFunction make_function(const Expression<Real>& expr, const Space<Real>& spc) {
    return EffectiveScalarFunction(dimension(spc),formula(expr,spc)); }
EffectiveScalarFunction make_function(const Expression<Real>& expr, const List< Variable<Real> >& vars) {
    return EffectiveScalarFunction(len(vars),formula(expr,vars)); }



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
    : VectorFunctionMixin<VectorOfScalarFunction<P,D>,P>
{
    typedef D DomainType;
    VectorOfScalarFunction(SizeType rs, SizeType as)
        : VectorOfScalarFunction(rs, BoxDomain(as,IntervalDomain(-inf,+inf))) { }
    VectorOfScalarFunction(SizeType rs, DomainType dom)
        : _as(dom.dimension()), _vec(rs,ScalarFunction<P,D>(dom)) { }
    VectorOfScalarFunction(SizeType rs, const ScalarFunction<P,D>& f)
        : _as(f.argument_size()), _vec(rs,f) { }
    VectorOfScalarFunction(const Vector<ScalarFunction<P>>& vsf)
        : _as(0u), _vec(vsf) { if(vsf.size()!=0) { _as=vsf[0].argument_size(); } }

    Void set(SizeType i, const ScalarFunction<P,D>& f) {
        if(this->argument_size()==0u) { this->_as=f.argument_size(); }
        ARIADNE_ASSERT(f.argument_size()==this->argument_size());
        this->_vec[i]=f; }
    ScalarFunction<P,D> get(SizeType i) const {
        return this->_vec[i]; }

    virtual SizeType result_size() const {
        return _vec.size(); }
    virtual SizeType argument_size() const {
        return _as; }

    virtual ScalarFunctionInterface<P,D>* _get(SizeType i) const {
        return static_cast<const ScalarFunctionInterface<P,D>&>(this->_vec[i])._clone(); }
    virtual Void _set(SizeType i, const ScalarFunctionInterface<P,D>* fptr) {
        this->_vec[i]=ScalarFunction<P,D>(fptr->_clone()); }
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

    template<class X> inline Void _compute(Vector<X>& r, const Vector<X>& x) const {
        r=Vector<X>(this->_vec.size(),x.zero_element());
        for(SizeType i=0; i!=r.size(); ++i) {
            r[i]=_vec[i].evaluate(x); } }

    SizeType _as;
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



//------------------------ Function -----------------------------------//

template<class P, class D, class C> ScalarFunction<P> Function<P,D,C>::zero(SizeType n) {
    return ScalarFunction<P>(new ScalarFormulaFunction<Y>(n,Formula<Y>::zero()));
}

template<class P, class D, class C> ScalarFunction<P> Function<P,D,C>::constant(SizeType n, CanonicalNumberType<P> c) {
    return ScalarFunction<P>(new ScalarFormulaFunction<Y>(n,Formula<Y>::constant(c)));
}

template<class P, class D, class C> ScalarFunction<P> Function<P,D,C>::coordinate(SizeType n, SizeType j) {
    return ScalarFunction<P>(new ScalarFormulaFunction<Y>(n,Formula<Y>::coordinate(j)));
}

template<class P, class D, class C> VectorFunction<P> Function<P,D,C>::zeros(SizeType rs, SizeType as) {
    VectorOfScalarFunction<P>* res = new VectorOfScalarFunction<P>(rs,as);
    for(SizeType i=0; i!=rs; ++i) {
        res->_vec[i]=ScalarFunction<P>::zero(as);
    }
    return VectorFunction<P>(res);
}

template<class P, class D, class C> VectorFunction<P> Function<P,D,C>::identity(SizeType n) {
    VectorOfScalarFunction<P>* res = new VectorOfScalarFunction<P>(n,n);
    for(SizeType i=0; i!=n; ++i) {
        res->_vec[i]=ScalarFunction<P>::coordinate(n,i);
    }
    return VectorFunction<P>(res);
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

ValidatedScalarFunction operator-(ValidatedScalarFunction const& f1, ValidatedScalarFunction const& f2) {
    std::shared_ptr<ValidatedScalarFunctionModelInterface const> f1p=std::dynamic_pointer_cast<ValidatedScalarFunctionModelInterface const>(f1.managed_pointer());
    std::shared_ptr<ValidatedScalarFunctionModelInterface const> f2p=std::dynamic_pointer_cast<ValidatedScalarFunctionModelInterface const>(f2.managed_pointer());
    if(f1p && f2p) {
        return ValidatedScalarFunctionModel(*f1p) - ValidatedScalarFunctionModel(*f2p);
    }
    return ValidatedScalarFunction(new BinaryFunction<ValidatedTag>(OperatorCode::SUB,f1,f2));
}

ValidatedScalarFunction operator-(ValidatedScalarFunction const& f, ValidatedNumber const& c) {
    std::shared_ptr<ValidatedScalarFunctionModelInterface const> fp=std::dynamic_pointer_cast<ValidatedScalarFunctionModelInterface const>(f.managed_pointer());
    if(fp) { return ValidatedScalarFunctionModel(*fp) - c; }
    std::shared_ptr<EffectiveScalarFunctionInterface const> rfp=std::dynamic_pointer_cast<EffectiveScalarFunctionInterface const>(f.managed_pointer());
    if(rfp && c.lower().raw()==c.upper().raw()) { return EffectiveScalarFunction(*rfp) - ExactFloat(c.lower().raw()); }
    return ValidatedScalarFunction(new BinaryFunction<ValidatedTag>(OperatorCode::SUB,f,ValidatedScalarFunction::constant(f.argument_size(),c)));
}

ValidatedScalarFunction operator-(ValidatedNumber const& c, ValidatedScalarFunction const& f) {
    std::shared_ptr<ValidatedScalarFunctionModelInterface const> fp=std::dynamic_pointer_cast<ValidatedScalarFunctionModelInterface const>(f.managed_pointer());
    if(fp) { return c - ValidatedScalarFunctionModel(*fp); }
    return ValidatedScalarFunction(new BinaryFunction<ValidatedTag>(OperatorCode::SUB,ValidatedScalarFunction::constant(f.argument_size(),c),f));
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
    return compose(f,ValidatedVectorFunctionModel(dynamic_cast<ValidatedVectorFunctionModelInterface const&>(*g.raw_pointer())));
}

ValidatedVectorFunction compose(const ValidatedVectorFunction& f, const ValidatedVectorFunction& g) {
    ARIADNE_ASSERT(f.argument_size()==g.result_size());
    return compose(f,ValidatedVectorFunctionModel(dynamic_cast<ValidatedVectorFunctionModelInterface const&>(*g.raw_pointer())));
}

ValidatedVectorFunction join(ValidatedVectorFunction const& f1, const ValidatedVectorFunction& f2) {
    if(dynamic_cast<VectorTaylorFunction const*>(f1.raw_pointer()) && dynamic_cast<VectorTaylorFunction const*>(f2.raw_pointer())) {
        return join(dynamic_cast<VectorTaylorFunction const&>(*f1.raw_pointer()),dynamic_cast<VectorTaylorFunction const&>(*f2.raw_pointer()));
    }
    VectorOfScalarFunction<ValidatedTag> r(f1.result_size()+f2.result_size(),f1.domain());
    for(SizeType i=0; i!=f1.result_size(); ++i) { r[i]=f1[i]; }
    for(SizeType i=0; i!=f2.result_size(); ++i) { r[i+f1.result_size()]=f2[i]; }
    return r;
}









} // namespace Ariadne
