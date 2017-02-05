/***************************************************************************
 *            procedure.tcc
 *
 *  Copyright 2010-17  Pieter Collins
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


#include <iostream>

#include "utility/container.h"
#include "algebra/vector.h"

#include "numeric/operators.h"
#include "function/formula.h"
#include "algebra/evaluate.h"
#include "algebra/expansion.h"
#include "geometry/box.h"

namespace Ariadne {

extern template class Procedure<ApproximateNumber>;
extern template class Procedure<ValidatedNumber>;

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
                os << c[instruction.arg]; break;
            case OperatorCode::IND:
                os << "x[" << instruction.arg << "]"; break;
            case OperatorCode::ADD: case OperatorCode::SUB: case OperatorCode::MUL: case OperatorCode::DIV:
            case OperatorCode::MAX: case OperatorCode::MIN:
                 os << "v[" << instruction.arg1 << "]" << symbol(instruction.op) << "v[" << instruction.arg2 << "]"; break;
            case OperatorCode::SADD: case OperatorCode::SSUB: case OperatorCode::SMUL: case OperatorCode::SDIV:
                os << c[instruction.arg1] << symbol(instruction.op) << "v[" << instruction.arg2 << "]"; break;
            case OperatorCode::POW:
                os<<"pow(v[" << instruction.arg << "],"<<instruction.np<<")"; break;
            case OperatorCode::POS: case OperatorCode::NEG:
                os << symbol(instruction.op) << "v[" << instruction.arg << "]"; break;
            case OperatorCode::ABS: case OperatorCode::REC: case OperatorCode::SQR: case OperatorCode::SQRT:
            case OperatorCode::EXP: case OperatorCode::LOG: case OperatorCode::SIN: case OperatorCode::COS: case OperatorCode::TAN:
            case OperatorCode::ASIN: case OperatorCode::ACOS: case OperatorCode::ATAN:
                os << instruction.op << "(v[" << instruction.arg << "])"; break;
            default:   ARIADNE_FAIL_MSG("Unrecognised operator "<<instruction.op);
        }
        os << "; ";
    }
}


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
        switch(instruction.op) {
//            case OperatorCode::CNST: { X r=x.zero_element(); r=c[instruction.arg]; v[i]=r; } break;
            case OperatorCode::CNST: v[i]=(x.zero_element()+c[instruction.arg]); break;
            case OperatorCode::IND:  v[i]=(x[instruction.arg]); break;
            case OperatorCode::ADD:  v[i]=(v[instruction.arg1]+v[instruction.arg2]); break;
            case OperatorCode::SUB:  v[i]=(v[instruction.arg1]-v[instruction.arg2]); break;
            case OperatorCode::MUL:  v[i]=(v[instruction.arg1]*v[instruction.arg2]); break;
            case OperatorCode::DIV:  v[i]=(v[instruction.arg1]/v[instruction.arg2]); break;
            case OperatorCode::SADD:  v[i]=(c[instruction.arg1]+v[instruction.arg2]); break;
            case OperatorCode::SSUB:  v[i]=(c[instruction.arg1]-v[instruction.arg2]); break;
            case OperatorCode::SMUL:  v[i]=(c[instruction.arg1]*v[instruction.arg2]); break;
            case OperatorCode::SDIV:  v[i]=(c[instruction.arg1]/v[instruction.arg2]); break;
            case OperatorCode::POW:  v[i]=(pow(v[instruction.arg],instruction.np)); break;
            case OperatorCode::ABS:  v[i]=(abs(v[instruction.arg])); break;
            case OperatorCode::POS:  v[i]=(pos(v[instruction.arg])); break;
            case OperatorCode::NEG:  v[i]=(neg(v[instruction.arg])); break;
            case OperatorCode::REC:  v[i]=(rec(v[instruction.arg])); break;
            case OperatorCode::SQR:  v[i]=(sqr(v[instruction.arg])); break;
            case OperatorCode::SQRT: v[i]=(sqrt(v[instruction.arg])); break;
            case OperatorCode::EXP:  v[i]=(exp(v[instruction.arg])); break;
            case OperatorCode::LOG:  v[i]=(log(v[instruction.arg])); break;
            case OperatorCode::SIN:  v[i]=(sin(v[instruction.arg])); break;
            case OperatorCode::COS:  v[i]=(cos(v[instruction.arg])); break;
            case OperatorCode::TAN:  v[i]=(tan(v[instruction.arg])); break;
            case OperatorCode::ATAN:  v[i]=(atan(v[instruction.arg])); break;
            default:   ARIADNE_FAIL_MSG("Unrecognised operator "<<instruction.op);
        }
    }
}


template<class X, class Y> Void _backpropagate(Vector<X>& x, List<X>& v, const List<ProcedureInstruction>& p, const List<Y>& c)
{
    Float64Value infty(inf);

    ARIADNE_ASSERT(v.size()==p.size());
    SizeType r=p.size();
    while(r!=0u) {
        --r;
        SizeType a=p[r].arg; SizeType a1=p[r].arg1; SizeType a2=p[r].arg2;
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
            case OperatorCode::POW:  restrict(v[a],exp(log(v[r])/p[r].np)); break;
            case OperatorCode::SQRT: restrict(v[a],sqr(v[r])); break;
            case OperatorCode::EXP:  restrict(v[a],log(v[r])); break;
            case OperatorCode::LOG:  restrict(v[a],exp(v[r])); break;
            case OperatorCode::SIN:  restrict(v[a],asin(v[r])); break;
            case OperatorCode::COS:  restrict(v[a],acos(v[r])); break;
            case OperatorCode::TAN:  restrict(v[a],atan(v[r])); break;
            case OperatorCode::ATAN: restrict(v[a],tan(v[r])); break;
            case OperatorCode::EQ:   restrict(v[a1],v[r]); restrict(v[a2],v[r]); break;
            case OperatorCode::LEQ:  restrict(v[a1],X(-infty,v[a2].upper())); restrict(v[a1],X(v[a2].lower(),+infty)); break;
            default: ARIADNE_THROW(std::runtime_error,"_propagate(Vector<X>,List<X>,List<ProcedureInstruction>)","Unhandled operator "<<p[r].op<<" at instruction "<<r<<"\n");
        }
    }
    // POSTCONDITION: No nan's get propagated to x
}

Void simple_hull_reduce(UpperBoxType& dom, const ValidatedProcedure& f, IntervalDomainType codom);
Void simple_hull_reduce(UpperBoxType& dom, const Vector<ValidatedProcedure>& f, BoxDomainType codom);

} // namespace Ariadne
