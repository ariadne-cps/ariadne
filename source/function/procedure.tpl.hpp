/***************************************************************************
 *            procedure.tcc
 *
 *  Copyright 2010-17  Pieter Collins
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

#include <iostream>

#include "../utility/container.hpp"
#include "../algebra/vector.hpp"

#include "../numeric/operators.hpp"
#include "../function/formula.hpp"
#include "../algebra/expansion.hpp"
#include "../geometry/box.hpp"

#include "../algebra/evaluate.tpl.hpp"

// FIXME: Added to prevent compilation error in Clang++-5.0. Should not be necessary.
#include "../function/taylor_model.hpp"
#include "../algebra/algebra.hpp"

namespace Ariadne {

extern template class Procedure<ApproximateNumber>;
extern template class Procedure<ValidatedNumber>;

template<class P> Procedure<Number<P>> make_procedure(ScalarMultivariateFunction<P> const& f) {
    typedef Number<P> Y;
    Formula<Y> e=f(Formula<Y>::identity(f.argument_size()));
    return Procedure<Y>(e);
}

template<class Y> SizeType _convert(List<ProcedureInstruction>& p, List<Y>& c, const Formula<Y>& f, Map<const Void*,SizeType>& cache) {
    if(cache.has_key(f.node_ptr())) { return cache[f.node_ptr()]; }
    switch(f.kind()) { // Can't use simple evaluate (above) as we need to pass the cache to subformulae
        case OperatorKind::COORDINATE: p.append(ProcedureInstruction(OperatorCode::IND,f.ind())); break;
        case OperatorKind::NULLARY: p.append(ProcedureInstruction(OperatorCode::CNST,c.size())); c.append(f.val()); break;
        case OperatorKind::BINARY:
            p.append(ProcedureInstruction(f.op(),_convert(p,c,f.arg1(),cache),_convert(p,c,f.arg2(),cache))); break;
        case OperatorKind::UNARY:
            p.append(ProcedureInstruction(f.op(),_convert(p,c,f.arg(),cache))); break;
        case OperatorKind::SCALAR:
            c.append(f.cnst()); p.append(ProcedureInstruction(f.op(),c.size()-1u,_convert(p,c,f.arg(),cache))); break;
        case OperatorKind::GRADED:
            p.append(ProcedureInstruction(f.op(),_convert(p,c,f.arg(),cache),f.num())); break;
        default: ARIADNE_FAIL_MSG("Unrecognised operator "<<f.op()<<" of kind "<<f.kind()<<" in formula "<<f);
    }
    const SizeType num=p.size()-1;
    cache.insert(f.node_ptr(),num);
    return num;
}

template<class Y>
Void _write(OutputStream& os, const List<ProcedureInstruction>& p, const List<Y>& c)
{
    for(SizeType i=0; i!=p.size(); ++i) {
        const ProcedureInstruction& instruction=p[i];
        os << "v[" << i << "]=";
        switch(instruction.op) {
            case OperatorCode::CNST:
                os << c[instruction.arg()]; break;
            case OperatorCode::IND:
                os << "x[" << instruction.arg() << "]"; break;
            case OperatorCode::ADD: case OperatorCode::SUB: case OperatorCode::MUL: case OperatorCode::DIV:
            case OperatorCode::MAX: case OperatorCode::MIN:
                 os << "v[" << instruction.arg1() << "]" << symbol(instruction.op) << "v[" << instruction.arg2() << "]"; break;
            case OperatorCode::SADD: case OperatorCode::SSUB: case OperatorCode::SMUL: case OperatorCode::SDIV:
                os << c[instruction.arg1()] << symbol(instruction.op) << "v[" << instruction.arg2() << "]"; break;
            case OperatorCode::POW:
                os<<"pow(v[" << instruction.arg() << "],"<<instruction.np()<<")"; break;
            case OperatorCode::POS: case OperatorCode::NEG:
                os << symbol(instruction.op) << "v[" << instruction.arg() << "]"; break;
            case OperatorCode::ABS: case OperatorCode::REC: case OperatorCode::SQR: case OperatorCode::SQRT:
            case OperatorCode::EXP: case OperatorCode::LOG: case OperatorCode::SIN: case OperatorCode::COS: case OperatorCode::TAN:
            case OperatorCode::ASIN: case OperatorCode::ACOS: case OperatorCode::ATAN:
                os << instruction.op << "(v[" << instruction.arg() << "])"; break;
            default:   ARIADNE_FAIL_MSG("Unrecognised operator "<<instruction.op);
        }
        os << "; ";
    }
}


//template<class Y, class X> Formula<Y> to_formula(const Expansion<MultiIndex,X>& e) {
//    return horner_evaluate(e,Formula<Y>::identity(e.argument_size()));
//};


template<class Y>
Procedure<Y>::Procedure()
{
}

template<class Y>
Procedure<Y>::Procedure(const Formula<Y>& f)
{
    Map<const Void*, SizeType> ind;
    _convert(this->_instructions,this->_constants,f, ind);
}


template<class Y> template<class X, EnableIf<IsConvertible<X,Y>>>
Procedure<Y>::Procedure(const Expansion<MultiIndex,X>& e)
    : Procedure(horner_evaluate(e,Formula<Y>::identity(e.argument_size())))
{
}


template<class Y> OutputStream& Procedure<Y>::_write(OutputStream& os) const {
    os<<"Procedure( ";
    Ariadne::_write(os,this->_instructions,this->_constants);
    os << "r=v[" << this->_instructions.size()-1u << "] )";
    return os;
}

template<class Y> Procedure<Y>& operator+=(Procedure<Y>& f, const Y& c) {
    f._constants.append(c);
    f.new_unary_instruction(OperatorCode::CNST,f._constants.size()-1);
    f.new_binary_instruction(OperatorCode::ADD,f._instructions.size()-1,f._instructions.size()-2);
    return f;
}


template<class Y>
Vector<Procedure<Y>>::Vector(const Vector<Formula<Y>>& f)
    : _results(f.size(),0u)
{
    Map<const Void*, SizeType> ind;
    for(SizeType i=0; i!=f.size(); ++i) {
        _convert(this->_instructions,this->_constants,f[i], ind);
    }
    for(Nat i=0; i!=f.size(); ++i) {
        this->_results[i]=ind[f[i].node_ptr()];
    }
}

template<class Y>
Vector<Procedure<Y>>::Vector(const Procedure<Y>& f)
    : _constants(f._constants)
    , _instructions(f._instructions)
    , _results(1u,f._instructions.size()-1u)
{
}

template<class Y>
OutputStream& Vector<Procedure<Y>>::_write(OutputStream& os) const {
    os<<"Procedure( ";
    Ariadne::_write(os,this->_instructions,this->_constants);
    os << "r=v" << this->_results << " )";
    return os;
}



// \related Procedure \brief Evaluate a function \a f defined by an algorithmic procedure
template<class X, class Y> Void _execute(List<X>& v, const List<ProcedureInstruction>& p, const List<Y>& c, const Vector<X>& x)
{
    ARIADNE_ASSERT(v.size()==p.size());
    for(SizeType i=0; i!=p.size(); ++i) {
        const ProcedureInstruction& instruction=p[i];
        SizeType a=instruction.arg(); SizeType a1=instruction.arg1(); SizeType a2=instruction.arg2(); Int n=instruction.np();
        switch(instruction.op) {
//            case OperatorCode::CNST: { X r=x.zero_element(); r=c[instruction.arg]; v[i]=r; } break;
            case OperatorCode::CNST: v[i]=(x.zero_element()+c[a]); break;
            case OperatorCode::IND:  v[i]=(x[a]); break;
            case OperatorCode::ADD:  v[i]=(v[a1]+v[a2]); break;
            case OperatorCode::SUB:  v[i]=(v[a1]-v[a2]); break;
            case OperatorCode::MUL:  v[i]=(v[a1]*v[a2]); break;
            case OperatorCode::DIV:  v[i]=(v[a1]/v[a2]); break;
            case OperatorCode::SADD:  v[i]=(c[a1]+v[a2]); break;
            case OperatorCode::SSUB:  v[i]=(c[a1]-v[a2]); break;
            case OperatorCode::SMUL:  v[i]=(c[a1]*v[a2]); break;
            case OperatorCode::SDIV:  v[i]=(c[a1]/v[a2]); break;
            case OperatorCode::POW:  v[i]=(pow(v[a],n)); break;
            case OperatorCode::ABS:  v[i]=(abs(v[a])); break;
            case OperatorCode::POS:  v[i]=(pos(v[a])); break;
            case OperatorCode::NEG:  v[i]=(neg(v[a])); break;
            case OperatorCode::REC:  v[i]=(rec(v[a])); break;
            case OperatorCode::SQR:  v[i]=(sqr(v[a])); break;
            case OperatorCode::SQRT: v[i]=(sqrt(v[a])); break;
            case OperatorCode::EXP:  v[i]=(exp(v[a])); break;
            case OperatorCode::LOG:  v[i]=(log(v[a])); break;
            case OperatorCode::SIN:  v[i]=(sin(v[a])); break;
            case OperatorCode::COS:  v[i]=(cos(v[a])); break;
            case OperatorCode::TAN:  v[i]=(tan(v[a])); break;
            case OperatorCode::ATAN:  v[i]=(atan(v[a])); break;
            default:   ARIADNE_FAIL_MSG("Unrecognised operator "<<instruction.op);
        }
    }
}


template<class X, class Y> Void _backpropagate(Vector<X>& x, List<X>& v, const List<ProcedureInstruction>& p, const List<Y>& c)
{
    FloatDPValue inf_(inf);

    ARIADNE_ASSERT(v.size()==p.size());
    SizeType r=p.size();
    while(r!=0u) {
        --r;
        SizeType a=p[r].arg(); SizeType a1=p[r].arg1(); SizeType a2=p[r].arg2(); Int n=p[r].np();
        switch(p[r].op) {
            case OperatorCode::CNST: break;
            case OperatorCode::IND:  restrict(x[a],v[r]); break;
            case OperatorCode::ADD:  restrict(v[a1],v[r]-v[a2]); restrict(v[a2],v[r]-v[a1]); break;
            case OperatorCode::SUB:  restrict(v[a1],v[a2]+v[r]); restrict(v[a2],v[a1]-v[r]); break;
            case OperatorCode::MUL:  restrict(v[a1],v[r]/v[a2]); restrict(v[a2],v[r]/v[a1]); break;
            case OperatorCode::DIV:  restrict(v[a1],v[a2]*v[r]); restrict(v[a2],v[a1]/v[r]); break;
            case OperatorCode::MAX:  restrict(v[a1],max(v[r],v[a2])); restrict(v[a2],max(v[r],v[a1])); break;
            case OperatorCode::MIN:  restrict(v[a1],min(v[r],v[a2])); restrict(v[a2],min(v[r],v[a1])); break;
            case OperatorCode::SADD: restrict(v[a2],v[r]-c[a1]); break;
            case OperatorCode::SSUB: restrict(v[a2],c[a1]-v[r]); break;
            case OperatorCode::SMUL: restrict(v[a2],v[r]/c[a1]); break;
            case OperatorCode::SDIV: restrict(v[a2],c[a1]/v[r]); break;
            case OperatorCode::POS:  restrict(v[a],v[r]); break;
            case OperatorCode::NEG:  restrict(v[a],neg(v[r])); break;
            case OperatorCode::REC:  restrict(v[a],rec(v[r])); break;
            case OperatorCode::SQR:  restrict(v[a],sqrt(v[r])); break;
            case OperatorCode::POW:  restrict(v[a],exp(log(v[r])/n)); break;
            case OperatorCode::SQRT: restrict(v[a],sqr(v[r])); break;
            case OperatorCode::EXP:  restrict(v[a],log(v[r])); break;
            case OperatorCode::LOG:  restrict(v[a],exp(v[r])); break;
            // FIXME: restricting asin/acos of interval should not be dependent on any branch of asin/acos
            case OperatorCode::SIN:  restrict(v[a],asin(v[r])); break;
            case OperatorCode::COS:  restrict(v[a],acos(v[r])); break;
            case OperatorCode::TAN:  restrict(v[a],atan(v[r])); break;
            case OperatorCode::ATAN: restrict(v[a],tan(v[r])); break;
            case OperatorCode::EQ:   restrict(v[a1],v[r]); restrict(v[a2],v[r]); break;
            case OperatorCode::LEQ:  restrict(v[a1],X(-inf_,v[a2].upper())); restrict(v[a1],X(v[a2].lower(),+inf_)); break;
            default: ARIADNE_THROW(std::runtime_error,"_propagate(Vector<X>,List<X>,List<ProcedureInstruction>)","Unhandled operator "<<p[r].op<<" at instruction "<<r<<"\n");
        }
    }
    // POSTCONDITION: No nan's get propagated to x
}


template<class X, class Y> Covector<X> gradient(Procedure<Y> const& f, Vector<X> const& x) {
    X z=x.zero_element();
    auto& p=f._instructions;
    auto& c=f._constants;

    List<X> v(p.size(),z);
    List<X> dv(v.size(),z);
    Covector<X> dfx(x.size(),z);

    execute(v,f,x);

    SizeType r=p.size();
    dv[r-1]=1;
    while(r!=0u) {
        --r;
        SizeType a=p[r].arg(); SizeType a1=p[r].arg1(); SizeType a2=p[r].arg2(); Int n=p[r].np();
        switch(p[r].op) {
            case OperatorCode::CNST: break;
            case OperatorCode::IND:  dfx[a]+=dv[r]; break;
            case OperatorCode::ADD:  dv[a1]+=dv[r]; dv[a2]+=dv[r]; break;
            case OperatorCode::SUB:  dv[a1]+=dv[r]; dv[a2]-=dv[r]; break;
            case OperatorCode::MUL:  dv[a1]+=dv[r]*v[a2]; dv[a2]+=dv[r]*v[a1]; break;
            case OperatorCode::DIV:  dv[a1]+=dv[r]*rec(v[a2]); dv[a2]-=dv[r]*v[a1]*sqr(rec(v[a2])); break;
            case OperatorCode::SADD: dv[a2]+=dv[r]; break;
            case OperatorCode::SSUB: dv[a2]-=dv[r]; break;
            case OperatorCode::SMUL: dv[a2]+=dv[r]*c[a1]; break;
            case OperatorCode::SDIV: dv[a2]-=dv[r]*c[a1]*sqr(rec(v[a2])); break;
            case OperatorCode::POS:  dv[a]+=Pos().derivative(v[a],dv[r]); break;
            case OperatorCode::NEG:  dv[a]+=Neg().derivative(v[a],dv[r]); break;
            case OperatorCode::REC:  dv[a]+=Rec().derivative(v[a],dv[r]); break;
            case OperatorCode::SQR:  dv[a]+=Sqr().derivative(v[a],dv[r]); break;
            case OperatorCode::POW:  dv[a]+=Pow().derivative(v[a],dv[r],n); break;
            case OperatorCode::SQRT: dv[a]+=Sqrt().derivative(v[a],dv[r]); break;
            case OperatorCode::EXP:  dv[a]+=Exp().derivative(v[a],dv[r]); break;
            case OperatorCode::LOG:  dv[a]+=Log().derivative(v[a],dv[r]); break;
            case OperatorCode::SIN:  dv[a]+=Sin().derivative(v[a],dv[r]); break;
            case OperatorCode::COS:  dv[a]+=Cos().derivative(v[a],dv[r]); break;
            case OperatorCode::TAN:  dv[a]+=Tan().derivative(v[a],dv[r]); break;
            case OperatorCode::ATAN: dv[a]+=Atan().derivative(v[a],dv[r]); break;
            default: ARIADNE_THROW(std::runtime_error,"gradient(...)","Unhandled operator "<<p[r].op<<" at instruction "<<r<<"\n");
        }
    }
    return dfx;
}

} // namespace Ariadne

#include "../algebra/fixed_differential.hpp"

namespace Ariadne {

template<class X, class Y> X hessian(Procedure<Y> const& f, Vector<X> const& x, Vector<X> const& s) {
    auto da=UnivariateSecondDifferential<X>::variable(s.zero_element());
    Vector<UnivariateSecondDifferential<X>> dxpas=x+da*s;
    UnivariateSecondDifferential<X> dfxpas=evaluate(f,dxpas);
    return dfxpas.hessian();
}

}

