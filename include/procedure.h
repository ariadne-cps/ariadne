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

#include "container.h"
#include "vector.h"

#include "operators.h"

namespace Ariadne {

template<class X> class Formula;

struct ProcedureInstruction
{
    explicit ProcedureInstruction(Operator o, uint a) : op(o), arg(a) { }
    explicit ProcedureInstruction(Operator o, uint a1, uint a2) : op(o), arg1(a1), arg2(a2) { }
    explicit ProcedureInstruction(Operator o, uint a, int n) : op(o), arg(a), np(n) { }
    Operator op;
    union {
        struct { uint arg; int np; };
        struct { uint arg1; uint arg2; };
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
    explicit Procedure(const Formula<X>& f);
  public:
    List<X> _constants;
    List<ProcedureInstruction> _instructions;
  public:
    void new_instruction(Operator o, uint a) { _instructions.append(ProcedureInstruction(o,a)); }
    void new_instruction(Operator o, uint a, int n) { _instructions.append(ProcedureInstruction(o,a,n)); }
    void new_instruction(Operator o, uint a1, uint a2) { _instructions.append(ProcedureInstruction(o,a1,a2)); }
  public:
};

//! \related Procedure \brief Evaluate a function \a f defined by an algorithmic procedure.
template<class X, class T> T evaluate(const Procedure<X>& f, const Vector<T>& x)
{
    T z=x[0]*0;
    List<T> v; v.reserve(f._instructions.size());
    for(uint i=0; i!=f._instructions.size(); ++i) {
        const ProcedureInstruction& instruction=f._instructions[i];
        switch(instruction.op) {
            case CNST: v.append(z+f._constants[instruction.arg]); break;
            case IND:  v.append(x[instruction.arg]); break;
            case ADD:  v.append(v[instruction.arg1]+v[instruction.arg2]); break;
            case SUB:  v.append(v[instruction.arg1]-v[instruction.arg2]); break;
            case MUL:  v.append(v[instruction.arg1]*v[instruction.arg2]); break;
            case DIV:  v.append(v[instruction.arg1]/v[instruction.arg2]); break;
            case POW:  v.append(pow(v[instruction.arg],instruction.np)); break;
            case ABS:  v.append(abs(v[instruction.arg])); break;
            case POS:  v.append(pos(v[instruction.arg])); break;
            case NEG:  v.append(neg(v[instruction.arg])); break;
            case REC:  v.append(rec(v[instruction.arg])); break;
            case SQR:  v.append(sqr(v[instruction.arg])); break;
            case SQRT: v.append(sqrt(v[instruction.arg])); break;
            case EXP:  v.append(exp(v[instruction.arg])); break;
            case LOG:  v.append(log(v[instruction.arg])); break;
            case SIN:  v.append(sin(v[instruction.arg])); break;
            case COS:  v.append(cos(v[instruction.arg])); break;
            case TAN:  v.append(tan(v[instruction.arg])); break;
            default:   ARIADNE_FAIL_MSG("Unrecognised operator "<<instruction.op);
        }
    }
    return v.back();
}

template<class X>
std::ostream& operator<<(std::ostream& os, const Procedure<X> f) {
    os<<"Procedure( ";
    for(uint i=0; i!=f._instructions.size(); ++i) {
        const ProcedureInstruction& instruction=f._instructions[i];
        os << "v[" << i << "]=";
        switch(instruction.op) {
            case CNST:
                os << f._constants[instruction.arg]; break;
            case IND:
                os << "x[" << instruction.arg << "]"; break;
            case ADD: case SUB: case MUL: case DIV:
                os << "v[" << instruction.arg1 << "]" << symbol(instruction.op) << "v[" << instruction.arg2 << "]"; break;
            case POW:
                os<<"pow(v[" << instruction.arg << "],"<<instruction.np<<")"; break;
            case POS: case NEG:
                os << symbol(instruction.op) << "v[" << instruction.arg << "]"; break;
            case ABS: case REC: case SQR: case SQRT:
            case EXP: case LOG: case SIN: case COS: case TAN:
                os << instruction.op << "(v[" << instruction.arg << "])"; break;
            default:   ARIADNE_FAIL_MSG("Unrecognised operator "<<instruction.op);
        }
        os << "; ";
    }
    os << "r=v["<<f._instructions.size()-1<<"] )";
    return os;
}

template<class X> uint _convert(Procedure<X>& p, const FormulaNode<X>* f, Map<const FormulaNode<X>*,uint>& ind) {
    if(ind.has_key(f)) { return ind[f]; }
    switch(f->op) { // Can't use simple evaluate (above) as we need to pass the cache to subformulae
        case CNST: p.new_instruction(CNST,p._constants.size()); p._constants.append(*f->val); break;
        case IND: p.new_instruction(IND,f->ind); break;
        case ADD: case SUB: case MUL: case DIV:
            p.new_instruction(f->op,_convert(p,f->arg1,ind),_convert(p,f->arg2,ind)); break;
        case POW:
            p.new_instruction(f->op,_convert(p,f->arg,ind),f->np); break;
        case NEG: case REC: case SQR: case SQRT:
        case EXP: case LOG: case SIN: case COS: case TAN:
            p.new_instruction(f->op,_convert(p,f->arg,ind)); break;
        default: assert(false);
    }
    const uint num=p._instructions.size()-1;
    ind.insert(f,num);
    return num;
}

template<class X>
Procedure<X>::Procedure(const Formula<X>& f)
{
    Map<const FormulaNode<X>*, uint> ind;
    _convert(*this,f._root.operator->(), ind);
}




} // namespace Ariadne


#endif // ARIADNE_PROCEDURE_H
