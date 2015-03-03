/***************************************************************************
 *            procedure.h
 *
 *  Copyright 2010  Pieter Collins
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


/*! \file procedure.h
 *  \brief Procedure to compute a real function
 */

#ifndef ARIADNE_PROCEDURE_H
#define ARIADNE_PROCEDURE_H

#include <iostream>

#include "utility/container.h"
#include "algebra/vector.h"

#include "geometry/box.h"

#include "expression/operators.h"
#include "expression/formula.h"
#include "algebra/expansion.h"

namespace Ariadne {

typedef Nat Nat;
template<class X> class Formula;
template<class X> class Graded;

template<class X> class Procedure;
typedef Procedure<ApproximateNumericType> ApproximateProcedure;
typedef Procedure<ValidatedNumericType> ValidatedProcedure;

struct ProcedureInstruction
{
    explicit ProcedureInstruction(OperatorCode o, SizeType a) : op(o), arg(a) { }
    explicit ProcedureInstruction(OperatorCode o, SizeType a1, SizeType a2) : op(o), arg1(a1), arg2(a2) { }
    explicit ProcedureInstruction(OperatorCode o, SizeType a, Int n) : op(o), arg(a), np(n) { }
    OperatorCode op;
    union {
        struct { SizeType arg; Int np; };
        struct { SizeType arg1; SizeType arg2; };
    };
};

//! \brief An algorithmic procedure for computing a function.
//!
//! A Procedure is more efficient to compute than a Formula, since common
//! subexpressions have already been removed. However, it is also more
//! difficult to manipulate, so it should usually only be used when no further
//! manipulations of the function are possible.
template<class X>
class Procedure {
  public:
    explicit Procedure();
    explicit Procedure(const Formula<X>& f);
    explicit Procedure(const Expansion<X>& f);
  public:
    List<X> _constants;
    List<ProcedureInstruction> _instructions;
  public:
    Void new_unary_instruction(OperatorCode o, SizeType a) { _instructions.append(ProcedureInstruction(o,a)); }
    Void new_binary_instruction(OperatorCode o, SizeType a1, SizeType a2) { _instructions.append(ProcedureInstruction(o,a1,a2)); }
    Void new_power_instruction(OperatorCode o, SizeType a, Int n) { _instructions.append(ProcedureInstruction(o,a,n)); }
  public:
};

//! \brief An algorithmic procedure for computing a function.
//!
//! A Procedure is more efficient to compute than a Formula, since common
//! subexpressions have already been removed. However, it is also more
//! difficult to manipulate, so it should usually only be used when no further
//! manipulations of the function are possible.
template<class X>
class Vector< Procedure<X> > {
  public:
    explicit Vector(const Vector< Formula<X> >& f);
    explicit Vector(const Procedure<X>& p);
  public:
    List<X> _constants;
    List<ProcedureInstruction> _instructions;
    Vector<Nat> _results;
  public:
    Nat result_size() const { return _results.size(); }
    Nat temporaries_size() const { return _instructions.size(); }
    Void new_instruction(OperatorCode o, SizeType a) { _instructions.append(ProcedureInstruction(o,a)); }
    Void new_instruction(OperatorCode o, SizeType a, Int n) { _instructions.append(ProcedureInstruction(o,a,n)); }
    Void new_instruction(OperatorCode o, SizeType a1, SizeType a2) { _instructions.append(ProcedureInstruction(o,a1,a2)); }
    Void set_return(SizeType i, SizeType a) { _results[i]=a; }
};

// \related Procedure \brief Evaluate a function \a f defined by an algorithmic procedure.
template<class X, class T> Void _execute(List<T>& v, const List<ProcedureInstruction>& p, const List<X>& c, const Vector<T>& x)
{
    T z=x.zero_element();
    for(SizeType i=0; i!=p.size(); ++i) {
        const ProcedureInstruction& instruction=p[i];
        switch(instruction.op) {
            case OperatorCode::CNST: v.append(z+static_cast<T>(c[instruction.arg])); break;
            case OperatorCode::IND:  v.append(x[instruction.arg]); break;
            case OperatorCode::ADD:  v.append(v[instruction.arg1]+v[instruction.arg2]); break;
            case OperatorCode::SUB:  v.append(v[instruction.arg1]-v[instruction.arg2]); break;
            case OperatorCode::MUL:  v.append(v[instruction.arg1]*v[instruction.arg2]); break;
            case OperatorCode::DIV:  v.append(v[instruction.arg1]/v[instruction.arg2]); break;
            case OperatorCode::POW:  v.append(pow(v[instruction.arg],instruction.np)); break;
            case OperatorCode::ABS:  v.append(abs(v[instruction.arg])); break;
            case OperatorCode::POS:  v.append(pos(v[instruction.arg])); break;
            case OperatorCode::NEG:  v.append(neg(v[instruction.arg])); break;
            case OperatorCode::REC:  v.append(rec(v[instruction.arg])); break;
            case OperatorCode::SQR:  v.append(sqr(v[instruction.arg])); break;
            case OperatorCode::SQRT: v.append(sqrt(v[instruction.arg])); break;
            case OperatorCode::EXP:  v.append(exp(v[instruction.arg])); break;
            case OperatorCode::LOG:  v.append(log(v[instruction.arg])); break;
            case OperatorCode::SIN:  v.append(sin(v[instruction.arg])); break;
            case OperatorCode::COS:  v.append(cos(v[instruction.arg])); break;
            case OperatorCode::TAN:  v.append(tan(v[instruction.arg])); break;
            case OperatorCode::ATAN:  v.append(atan(v[instruction.arg])); break;
            default:   ARIADNE_FAIL_MSG("Unrecognised operator "<<instruction.op);
        }
    }
}

// \related Procedure \brief Evaluate a function \a f defined by an algorithmic procedure.
template<class X, class T> Void _compute(List<T>& v, const List<ProcedureInstruction>& p, const List<X>& c, const Vector<T>& x)
{
    T z=x.zero_element();
    ARIADNE_ASSERT(v.size()==p.size());
    for(SizeType i=0; i!=p.size(); ++i) {
        const ProcedureInstruction& instruction=p[i];
        switch(instruction.op) {
            case OperatorCode::CNST: v[i]=(z+c[instruction.arg]); break;
            case OperatorCode::IND:  v[i]=(x[instruction.arg]); break;
            case OperatorCode::ADD:  v[i]=(v[instruction.arg1]+v[instruction.arg2]); break;
            case OperatorCode::SUB:  v[i]=(v[instruction.arg1]-v[instruction.arg2]); break;
            case OperatorCode::MUL:  v[i]=(v[instruction.arg1]*v[instruction.arg2]); break;
            case OperatorCode::DIV:  v[i]=(v[instruction.arg1]/v[instruction.arg2]); break;
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

template<class T> Void _propagate(Vector<T>& x, List<T>& v, const List<ProcedureInstruction>& p)
{
    ExactFloat64 infty(inf);

    ARIADNE_ASSERT(v.size()==p.size());
    SizeType r=p.size();
    while(r!=0u) {
        --r;
        SizeType a=p[r].arg; SizeType a1=p[r].arg1; SizeType a2=p[r].arg2;
        switch(p[r].op) {
            case OperatorCode::CNST: break;
            case OperatorCode::IND: restrict(x[a],v[r]); break;
            case OperatorCode::ADD: restrict(v[a1],v[r]-v[a2]); restrict(v[a2],v[r]-v[a1]); break;
            case OperatorCode::SUB: restrict(v[a1],v[a2]+v[r]); restrict(v[a2],v[a1]-v[r]); break;
            case OperatorCode::MUL: restrict(v[a1],v[r]/v[a2]); restrict(v[a2],v[r]/v[a1]); break;
            case OperatorCode::DIV: restrict(v[a1],v[a2]*v[r]); restrict(v[a2],v[a1]/v[r]); break;
            case OperatorCode::MAX: restrict(v[a1],max(v[r],v[a2])); restrict(v[a2],max(v[r],v[a1])); break;
            case OperatorCode::POS: restrict(v[a],v[r]); break;
            case OperatorCode::NEG: restrict(v[a],neg(v[r])); break;
            case OperatorCode::REC: restrict(v[a],rec(v[r])); break;
            case OperatorCode::SQR: restrict(v[a],sqrt(v[r])); break;
            case OperatorCode::POW: restrict(v[a],exp(log(v[r])/p[r].np)); break;
            case OperatorCode::SQRT: restrict(v[a],sqr(v[r])); break;
            case OperatorCode::EXP: restrict(v[a],log(v[r])); break;
            case OperatorCode::LOG: restrict(v[a],exp(v[r])); break;
            case OperatorCode::SIN: restrict(v[a],asin(v[r])); break;
            case OperatorCode::COS: restrict(v[a],acos(v[r])); break;
            case OperatorCode::TAN: restrict(v[a],atan(v[r])); break;
            case OperatorCode::ATAN: restrict(v[a],tan(v[r])); break;
            case OperatorCode::EQ: restrict(v[a1],v[r]); restrict(v[a2],v[r]); break;
            case OperatorCode::LEQ: restrict(v[a1],T(-infty,v[a2].upper())); restrict(v[a1],T(v[a2].lower(),+infty)); break;
            default: ARIADNE_THROW(std::runtime_error,"_propagate(Vector<T>,List<T>,List<ProcedureInstruction>)","Unhandled operator "<<p[r].op<<" at instruction "<<r<<"\n");
        }
    }
    // POSTCONDITION: No nan's get propagated to x
}

template<class X> SizeType _convert(List<ProcedureInstruction>& p, List<X>& c, const Formula<X>& f, Map<const Void*,SizeType>& cache) {
    if(cache.has_key(f.node_ptr())) { return cache[f.node_ptr()]; }
    switch(f.kind()) { // Can't use simple evaluate (above) as we need to pass the cache to subformulae
        case OperatorKind::COORDINATE: p.append(ProcedureInstruction(OperatorCode::IND,f.ind())); break;
        case OperatorKind::NULLARY: p.append(ProcedureInstruction(OperatorCode::CNST,c.size())); c.append(f.val()); break;
        case OperatorKind::BINARY:
            p.append(ProcedureInstruction(f.op(),_convert(p,c,f.arg1(),cache),_convert(p,c,f.arg2(),cache))); break;
        case OperatorKind::UNARY:
            p.append(ProcedureInstruction(f.op(),_convert(p,c,f.arg(),cache))); break;
        case OperatorKind::SCALAR:
            p.append(ProcedureInstruction(f.op(),_convert(p,c,f.arg(),cache),f.num())); break;
        default: ARIADNE_FAIL_MSG("Unrecognised operator "<<f.op()<<" of kind "<<f.kind()<<" in formula "<<f);
    }
    const SizeType num=p.size()-1;
    cache.insert(f.node_ptr(),num);
    return num;
}

template<class X>
Void _write(OutputStream& os, const List<ProcedureInstruction>& p, const List<X>& c)
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
                os << "v[" << instruction.arg1 << "]" << symbol(instruction.op) << "v[" << instruction.arg2 << "]"; break;
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


template<class X>
Procedure<X>::Procedure()
{
}

template<class X>
Procedure<X>::Procedure(const Formula<X>& f)
{
    Map<const Void*, SizeType> ind;
    _convert(this->_instructions,this->_constants,f, ind);
}

template<class X>
Procedure<X>::Procedure(const Expansion<X>& e)
{
    Vector< Formula<X> > id=Formula<X>::identity(e.argument_size());
    Formula<X> f=horner_evaluate(e,id);
    Map<const Void*, SizeType> ind;
    _convert(this->_instructions,this->_constants,f, ind);
}

// \related Procedure \brief Evaluate a function \a f defined by an algorithmic procedure.
template<class X, class T> inline T evaluate(const Procedure<X>& f, const Vector<T>& x)
{
    List<T> e;
    _execute(e,f._instructions,f._constants,x);
    return e.back();
}

template<class X> Procedure<X>& operator+=(Procedure<X>& f, const X& c) {
    f._constants.append(c);
    f.new_unary_instruction(OperatorCode::CNST,f._constants.size()-1);
    f.new_binary_instruction(OperatorCode::ADD,f._instructions.size()-1,f._instructions.size()-2);
    return f;
}

template<class X>
OutputStream& operator<<(OutputStream& os, const Procedure<X> f) {
    os<<"Procedure( ";
    _write(os,f._instructions,f._constants);
    os << "r=v[" << f._instructions.size()-1u << "] )";
    return os;
}


template<class X>
Vector< Procedure<X> >::Vector(const Vector< Formula<X> >& f)
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

template<class X>
Vector< Procedure<X> >::Vector(const Procedure<X>& f)
    : _instructions(f._instructions)
    , _constants(f._constants)
    , _results(1u,f._instructions.size()-1u)
{
}

// \related Procedure \brief Evaluate a function \a f defined by an algorithmic procedure.
template<class X, class T> inline Vector<T> evaluate(const Vector< Procedure<X> >& f, const Vector<T>& x)
{
    List<T> v;
    _execute(v,f._instructions,f._constants,x);
    Vector<T> r(f.result_size());
    for(SizeType i=0; i!=r.size(); ++i) {
        r[i]=v[f._results[i]];
    }
    return r;
}

template<class X>
OutputStream& operator<<(OutputStream& os, const Vector< Procedure<X> > f) {
    os<<"Procedure( ";
    _write(os,f._instructions,f._constants);
    os << "r=v" << f._results << " )";
    return os;
}


} // namespace Ariadne


#include "numeric/numeric.h"

namespace Ariadne {

// NOTE: Ordering of r and x is important, since nan elements of r are preserved,
// but nan elements of x do not affect r
inline
Void restrict(UpperIntervalType& r, const UpperIntervalType& x) {
    r.set_lower(max(r.lower(),x.lower()));
    r.set_upper(min(r.upper(),x.upper()));
};

template<class X>
Void simple_hull_reduce(UpperBoxType& dom, const Procedure<X>& f, ExactIntervalType codom)
{
    const List<ProcedureInstruction>& p=f._instructions;
    const List<X>& c=f._constants;
    List<UpperIntervalType> v; v.reserve(p.size());

    _execute(v,p,c,dom);
    restrict(v.back(),codom);
    _propagate(dom,v,p);
}

template<class X>
Void simple_hull_reduce(UpperBoxType& dom, const Vector< Procedure<X> >& f, ExactBoxType codom)
{
    const List<ProcedureInstruction>& p=f._instructions;
    const List<X>& c=f._constants;
    List<UpperIntervalType> v; v.reserve(p.size());

    ARIADNE_ASSERT(codom.size()==f._results.size());

    _execute(v,p,c,dom);
    for(SizeType i=0; i!=codom.size(); ++i) {
        restrict(v[f._results[i]],codom[i]);
    }
    _propagate(dom,v,p);
}

} // namespace Ariadne


#endif // ARIADNE_PROCEDURE_H
