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
#include "function/function.h"
#include "expression/operators.h"
#include "expression/formula.h"
#include "algebra/algebra.h"

#include "function/function_mixin.h"
#include "function_mixin.tcc"

#include "function/function_model.h"

#include "function/symbolic_function.h"
#include "function/taylor_function.h"


namespace Ariadne {

template<class T> inline std::string str(const T& t) {
    std::stringstream ss; ss << t; return ss.str(); }

typedef uint Nat;
typedef int Int;

// Functions for converting from Expression classes to Function classes via Formula
template<class X> class Expression;
template<class X> class Variable;
template<class X> class Space;
Nat dimension(const Space<Real>& spc);
Nat len(const List< Variable<Real> >& vars);
Formula<Real> formula(const Expression<Real>& expr, const Space<Real>& spc);
Formula<Real> formula(const Expression<Real>& expr, const List< Variable<Real> >& vars);
EffectiveScalarFunction make_function(const Expression<Real>& expr, const Space<Real>& spc) {
    return EffectiveScalarFunction(dimension(spc),formula(expr,spc)); }
EffectiveScalarFunction make_function(const Expression<Real>& expr, const List< Variable<Real> >& vars) {
    return EffectiveScalarFunction(len(vars),formula(expr,vars)); }



//------------------------ FunctionInterface forwarded functions  -----------------------------------//

Vector<ApproximateNumber> ScalarFunctionInterface<ApproximateTag>::gradient(const Vector<ApproximateNumber>& x) const {
    return this->evaluate(Differential<ApproximateNumber>::variables(1u,x)).gradient(); }
Vector<ValidatedNumber> ScalarFunctionInterface<ValidatedTag>::gradient(const Vector<ValidatedNumber>& x) const {
    return this->evaluate(Differential<ValidatedNumber>::variables(1u,x)).gradient(); }
Differential<ApproximateNumber> ScalarFunctionInterface<ApproximateTag>::differential(const Vector<ApproximateNumber>& x, Nat d) const {
    return this->evaluate(Differential<ApproximateNumber>::variables(d,x)); }
Differential<ValidatedNumber> ScalarFunctionInterface<ValidatedTag>::differential(const Vector<ValidatedNumber>& x, Nat d) const {
    return this->evaluate(Differential<ValidatedNumber>::variables(d,x)); }
Matrix<ApproximateNumber> VectorFunctionInterface<ApproximateTag>::jacobian(const Vector<ApproximateNumber>& x) const {
    return this->evaluate(Differential<ApproximateNumber>::variables(1u,x)).jacobian(); }
Matrix<ValidatedNumber> VectorFunctionInterface<ValidatedTag>::jacobian(const Vector<ValidatedNumber>& x) const {
    return this->evaluate(Differential<ValidatedNumber>::variables(1u,x)).jacobian(); }
Vector< Differential<ApproximateNumber> > VectorFunctionInterface<ApproximateTag>::differentials(const Vector<ApproximateNumber>& x, Nat d) const {
    return this->evaluate(Differential<ApproximateNumber>::variables(d,x)); }
Vector< Differential<ValidatedNumber> > VectorFunctionInterface<ValidatedTag>::differentials(const Vector<ValidatedNumber>& x, Nat d) const {
    return this->evaluate(Differential<ValidatedNumber>::variables(d,x)); }

//------------------------ Vector of Scalar functions  -----------------------------------//


template<class X> class NonResizableScalarFunction : public ScalarFunction<X> {
  public:
    NonResizableScalarFunction<X>& operator=(const ScalarFunction<X>& f) {
        ARIADNE_ASSERT(this->argument_size()==f.argument_size());
        this->ScalarFunction<X>::operator=(f);
        return *this;
    }
};

template<class X>
struct VectorOfScalarFunction
    : VectorFunctionMixin<VectorOfScalarFunction<X>,X>
{
    VectorOfScalarFunction(uint rs, uint as)
        : _as(as), _vec(rs,ScalarFunction<X>(as)) { }
    VectorOfScalarFunction(uint rs, const ScalarFunction<X>& f)
        : _as(f.argument_size()), _vec(rs,f) { }

    void set(uint i, const ScalarFunction<X>& f) {
        if(this->argument_size()==0u) { this->_as=f.argument_size(); }
        ARIADNE_ASSERT(f.argument_size()==this->argument_size());
        this->_vec[i]=f; }
    ScalarFunction<X> get(uint i) const {
        return this->_vec[i]; }

    virtual uint result_size() const {
        return _vec.size(); }
    virtual uint argument_size() const {
        return _as; }

    virtual ScalarFunctionInterface<X>* _get(uint i) const { return static_cast<const ScalarFunctionInterface<X>&>(this->_vec[i])._clone(); }

    const ScalarFunction<X> operator[](uint i) const {
        return this->_vec[i]; }

    NonResizableScalarFunction<X>& operator[](uint i) {
        return static_cast<NonResizableScalarFunction<X>&>(this->_vec[i]); }

    virtual std::ostream& write(std::ostream& os) const {
        os << "[";
        for(uint i=0; i!=this->_vec.size(); ++i) {
            if(i!=0) { os << ","; }
            this->_vec[i].raw_pointer()->write(os); }
        return os << "]"; }

    virtual std::ostream& repr(std::ostream& os) const {
        //os << "VoSF[R" << this->argument_size() << "->R" << this->result_size() << "]";
        os << "[";
        for(uint i=0; i!=this->_vec.size(); ++i) {
            if(i!=0) { os << ","; }
            this->_vec[i].raw_pointer()->repr(os); }
        return os << "]"; }

    template<class XX> inline void _compute(Vector<XX>& r, const Vector<XX>& x) const {
        r=Vector<XX>(this->_vec.size(),x.zero_element());
        for(uint i=0; i!=r.size(); ++i) {
            r[i]=_vec[i].evaluate(x); } }

    uint _as;
    Vector< ScalarFunction<X> > _vec;

};


template<class X>
struct FunctionElement
    : ScalarFunctionMixin<FunctionElement<X>,X>
{
    FunctionElement(const VectorFunction<X>& f, Nat i)
        : _f(f), _i(i) { ARIADNE_ASSERT(i<f.result_size()); }

    virtual Nat argument_size() const { return _f.argument_size(); }
    virtual std::ostream& write(std::ostream& os) const { return os<<_f<<"["<<_i<<"]"; }
    virtual ScalarFunctionInterface<X>* _derivative(uint j) const { ARIADNE_NOT_IMPLEMENTED; }

    template<class XX> inline void _compute(XX& r, const Vector<XX>& x) const {
        r=this->_f.evaluate(x)[_i]; }

    VectorFunction<X> _f;
    Nat _i;
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
        for(uint i=0; i!=r.size(); ++i) { r[i]=formula((*vf)[i]); }
        return r;
    } else if( (ff=dynamic_cast<const VectorFormulaFunction<Real>*>(&fi)) ) {
        return ff->_formulae;
    } else {
        ARIADNE_FAIL_MSG("Cannot compute formula for function "<<f<<"\n");
    }
}



//------------------------ ScalarFunction -----------------------------------//

template<class P> ScalarFunction<P>::ScalarFunction()
    : _ptr(new ScalarFormulaFunction<X>(0u,Formula<X>::zero()))
{
}

template<class P> ScalarFunction<P>::ScalarFunction(Nat n)
    : _ptr(new ScalarFormulaFunction<X>(n,Formula<X>::zero()))
{
}

template<class P> ScalarFunction<P>::ScalarFunction(Nat n, Formula<X> e)
    : _ptr(new ScalarFormulaFunction<X>(n,e))
{
}

template<class P> ScalarFunction<P> ScalarFunction<P>::zero(Nat n)
{
    return ScalarFunction<P>(new ScalarFormulaFunction<X>(n,Formula<X>::zero()));
}

template<class P> ScalarFunction<P> ScalarFunction<P>::constant(Nat n, CanonicalNumberType<P> c)
{
    return ScalarFunction<P>(new ScalarFormulaFunction<X>(n,Formula<X>::constant(c)));
}

template<class P> ScalarFunction<P> ScalarFunction<P>::coordinate(Nat n, Nat j)
{
    return ScalarFunction<P>(new ScalarFormulaFunction<X>(n,Formula<X>::coordinate(j)));
}

template<class P> List< ScalarFunction<P> > ScalarFunction<P>::coordinates(Nat n)
{
    List< ScalarFunction<P> > lsf; lsf.reserve(n);
    for(Nat i=0; i!=n; ++i) { lsf.append(ScalarFunction<P>::coordinate(n,i)); }
    return lsf;
}

template class ScalarFunction<ApproximateTag>;
template class ScalarFunction<ValidatedTag>;
template class ScalarFunction<EffectiveTag>;



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

EffectiveScalarFunction operator/(const EffectiveScalarFunction& f1, int s2)
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

EffectiveScalarFunction pow(const EffectiveScalarFunction& f, Nat m)
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






//------------------------ Vector Function ----------------------------------//


template<class P> VectorFunction<P>::VectorFunction()
    : _ptr(new VectorOfScalarFunction<P>(0u,ScalarFunction<P>()))
{
}

template<class P> VectorFunction<P>::VectorFunction(Nat rs, Nat as)
    : _ptr(new VectorOfScalarFunction<P>(rs,ScalarFunction<P>::zero(as)))
{
}

template<class P> VectorFunction<P>::VectorFunction(std::initializer_list< ScalarFunction<P> > lsf)
{
    ARIADNE_ASSERT(lsf.size()>0);
    Nat as=lsf.begin()->argument_size();
    VectorOfScalarFunction<P>* new_ptr=new VectorOfScalarFunction<P>(lsf.size(),as);
    for(uint i=0; i!=lsf.size(); ++i) {
        new_ptr->set(i,*(lsf.begin()+i));
    }
    this->_ptr=std::shared_ptr< const VectorFunctionInterface<P> >(new_ptr);
}

template<class P> VectorFunction<P>::VectorFunction(const List< ScalarFunction<P> >& lsf)
{
    ARIADNE_ASSERT(lsf.size()>0);
    Nat as=lsf[0].argument_size();
    VectorOfScalarFunction<P>* new_ptr=new VectorOfScalarFunction<P>(lsf.size(),as);
    for(uint i=0; i!=lsf.size(); ++i) {
        new_ptr->set(i,lsf[i]);
    }
    this->_ptr=std::shared_ptr< const VectorFunctionInterface<P> >(new_ptr);
}

template<class P> VectorFunction<P>::VectorFunction(Nat as, const List< Formula<X> >& le)
{
    ARIADNE_ASSERT(le.size()>0);
    VectorOfScalarFunction<P>* new_ptr=new VectorOfScalarFunction<P>(le.size(),as);
    for(uint i=0; i!=le.size(); ++i) {
        new_ptr->set(i,ScalarFormulaFunction<X>(as,le[i]));
    }
    this->_ptr=std::shared_ptr< const VectorFunctionInterface<P> >(new_ptr);
}

template<class P> VectorFunction<P>::VectorFunction(Nat rs, ScalarFunction<P> sf)
    : _ptr(new VectorOfScalarFunction<P>(rs,sf))
{
}

template<class P> VectorFunction<P> VectorFunction<P>::zeros(Nat rs, Nat as)
{
    VectorOfScalarFunction<P>* res = new VectorOfScalarFunction<P>(rs,as);
    for(uint i=0; i!=rs; ++i) {
        res->_vec[i]=ScalarFunction<P>::zero(as);
    }
    return VectorFunction<P>(res);
}

template<class P> VectorFunction<P> VectorFunction<P>::identity(Nat n)
{
    VectorOfScalarFunction<P>* res = new VectorOfScalarFunction<P>(n,n);
    for(uint i=0; i!=n; ++i) {
        res->_vec[i]=ScalarFunction<P>::coordinate(n,i);
    }
    return VectorFunction<P>(res);
}

template<class P> Void VectorFunction<P>::set(Nat i, ScalarFunction<P> sf)
{
    if(!this->_ptr.unique()) {
        ARIADNE_WARN("VectorFunction<P>::set(Nat, ScalarFunction<P>): Cloning shared pointer.\n");
        this->_ptr=std::shared_ptr< const VectorFunctionInterface<P> >(this->_ptr->_clone());
    }
    VectorFunctionInterface<P>* vf=const_cast<VectorFunctionInterface<P>*>(this->raw_pointer());
    VectorOfScalarFunction<P>* vsf=dynamic_cast<VectorOfScalarFunction<P>*>(vf);
    if(vsf) {
        vsf->set(i,sf);
    } else {
        ARIADNE_THROW(std::runtime_error,"Void VectorFunction<P>::set(Nat i, ScalarFunction<P> sf)","Cannot set element of "<<*vf<<"\n");
    }
}

template class VectorFunction<ApproximateTag>;
template class VectorFunction<ValidatedTag>;
template class VectorFunction<EffectiveTag>;



//------------------------ Vector function operators -------------------------------//

EffectiveVectorFunction operator*(const EffectiveScalarFunction& f, const Vector<Real>& e) {
    for(uint i=0; i!=e.size(); ++i) { ARIADNE_ASSERT(e[i]==Real(0) || e[i]==Real(1)); }
    EffectiveVectorFunction r(e.size(),f.argument_size());
    for(uint i=0; i!=e.size(); ++i) {
        if(e[i]==Real(1)) { r.set(i,f); }
    }
    return r;
}

EffectiveVectorFunction operator+(const EffectiveVectorFunction& f1, const EffectiveVectorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    EffectiveVectorFunction r(f1.result_size(),f1.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,f1[i]+f2[i]);
    }
    return r;
}

EffectiveVectorFunction operator-(const EffectiveVectorFunction& f1, const EffectiveVectorFunction& f2) {
    ARIADNE_ASSERT(f1.result_size()==f2.result_size());
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    EffectiveVectorFunction r(f1.result_size(),f1.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,f1[i]-f2[i]);
    }
    return r;
}

EffectiveVectorFunction operator*(const EffectiveVectorFunction& vf, const EffectiveScalarFunction& sf) {
    ARIADNE_ASSERT(vf.argument_size()==sf.argument_size());
    EffectiveVectorFunction r(vf.result_size(),vf.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,vf[i]*sf);
    }
    return r;
}

EffectiveVectorFunction operator*(const EffectiveScalarFunction& sf, const EffectiveVectorFunction& vf) {
    ARIADNE_ASSERT(sf.argument_size()==vf.argument_size());
    EffectiveVectorFunction r(vf.result_size(),vf.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,sf*vf[i]);
    }
    return r;
}

EffectiveVectorFunction operator*(const Real& c, const EffectiveVectorFunction& vf) {
    EffectiveVectorFunction r(vf.result_size(),vf.argument_size());
    for(uint i=0; i!=r.result_size(); ++i) {
        r.set(i,c*vf[i]);
    }
    return r;
}



EffectiveVectorFunction join(const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2) {
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    EffectiveVectorFunction r(2,f1.argument_size());
    r.set(0,f1);
    r.set(1,f2);
    return r;
}

EffectiveVectorFunction join(const EffectiveVectorFunction& f1, const EffectiveScalarFunction& f2) {
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    EffectiveVectorFunction r(f1.result_size()+1u,f1.argument_size());
    for(uint i=0u; i!=f1.result_size(); ++i) {
        r.set(i,f1.get(i));
    }
    r.set(f1.result_size(),f2);
    return r;
}

EffectiveVectorFunction join(const EffectiveScalarFunction& f1, const EffectiveVectorFunction& f2) {
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    EffectiveVectorFunction r(f2.result_size()+1u,f1.argument_size());
    r.set(0u,f1);
    for(uint i=0u; i!=f2.result_size(); ++i) {
        r.set(i+1u,f2.get(i));
    }
    return r;
}

EffectiveVectorFunction join(const EffectiveVectorFunction& f1, const EffectiveVectorFunction& f2) {
    ARIADNE_ASSERT(f1.argument_size()==f2.argument_size());
    EffectiveVectorFunction r(f1.result_size()+f2.result_size(),f1.argument_size());
    for(uint i=0u; i!=f1.result_size(); ++i) {
        r.set(i,f1.get(i));
    }
    for(uint i=0u; i!=f2.result_size(); ++i) {
        r.set(i+f1.result_size(),f2.get(i));
    }
    return r;
}

EffectiveScalarFunction embed(Nat as1, const EffectiveScalarFunction& f2, Nat as3) {
    return EffectiveScalarFunction(new ScalarEmbeddedFunction<EffectiveTag>(as1,f2,as3));
}

EffectiveVectorFunction embed(Nat as1, const EffectiveVectorFunction& f2, Nat as3) {
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
        for(uint i=1; i!=g.argument_size(); ++i) {
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
    return ValidatedScalarFunction(new BinaryFunction<ValidatedTag>(SUB,f1,f2));
}

ValidatedScalarFunction operator-(ValidatedScalarFunction const& f, ValidatedNumber const& c) {
    std::shared_ptr<ValidatedScalarFunctionModelInterface const> fp=std::dynamic_pointer_cast<ValidatedScalarFunctionModelInterface const>(f.managed_pointer());
    if(fp) { return ValidatedScalarFunctionModel(*fp) - c; }
    std::shared_ptr<EffectiveScalarFunctionInterface const> rfp=std::dynamic_pointer_cast<EffectiveScalarFunctionInterface const>(f.managed_pointer());
    if(rfp && c.lower().value()==c.upper().value()) { return EffectiveScalarFunction(*rfp) - ExactFloat(c.lower().value()); }
    return ValidatedScalarFunction(new BinaryFunction<ValidatedTag>(SUB,f,ValidatedScalarFunction::constant(f.argument_size(),c)));
}

ValidatedScalarFunction operator-(ValidatedNumber const& c, ValidatedScalarFunction const& f) {
    std::shared_ptr<ValidatedScalarFunctionModelInterface const> fp=std::dynamic_pointer_cast<ValidatedScalarFunctionModelInterface const>(f.managed_pointer());
    if(fp) { return c - ValidatedScalarFunctionModel(*fp); }
    return ValidatedScalarFunction(new BinaryFunction<ValidatedTag>(SUB,ValidatedScalarFunction::constant(f.argument_size(),c),f));
}

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
        for(uint i=0; i!=r.result_size(); ++i) {
            r[i]=f1[i]-f2[i];
        }
        return r;
    }
}

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
    VectorOfScalarFunction<ValidatedTag> r(f1.result_size()+f2.result_size(),f1.argument_size());
    for(uint i=0; i!=f1.result_size(); ++i) { r[i]=f1[i]; }
    for(uint i=0; i!=f2.result_size(); ++i) { r[i+f1.result_size()]=f2[i]; }
    return r;
}









} // namespace Ariadne
