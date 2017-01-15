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

#include "numeric/operators.h"
#include "function/formula.h"
#include "algebra/expansion.h"

namespace Ariadne {

template<class Y> class Formula;
template<class X> class Graded;

template<class Y> class Procedure;
typedef Procedure<ApproximateNumber> ApproximateProcedure;
typedef Procedure<ValidatedNumber> ValidatedProcedure;

template<class Y, class X> Formula<Y> to_formula(const Expansion<X>& e) {
    return horner_evaluate(e,Formula<Y>::identity(e.argument_size()));
};


Void simple_hull_reduce(UpperBoxType& dom, const ValidatedProcedure& f, ExactIntervalType codom);
Void simple_hull_reduce(UpperBoxType& dom, const Vector<ValidatedProcedure>& f, ExactBoxType codom);

struct ProcedureInstruction {
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
template<class Y>
class Procedure {
  public:
    explicit Procedure<Y>();
    explicit Procedure<Y>(const Formula<Y>& f);
    template<class X, EnableIf<IsConvertible<X,Y>> =dummy> explicit Procedure<Y>(const Expansion<X>& e)
        : Procedure(to_formula<Y>(e)) { }
    friend OutputStream& operator<<(OutputStream& os, Procedure<Y> const& p) { return p._write(os); }
  public:
   template<class X, class YY> friend X evaluate(const Procedure<YY>& p, const Vector<X>& x);
  private: public:
    List<Y> _constants;
    List<ProcedureInstruction> _instructions;
  public:
    Void new_constant(Y const& c) { _constants.append(c); }
    Void new_unary_instruction(OperatorCode o, SizeType a) { _instructions.append(ProcedureInstruction(o,a)); }
    Void new_binary_instruction(OperatorCode o, SizeType a1, SizeType a2) { _instructions.append(ProcedureInstruction(o,a1,a2)); }
    Void new_scalar_instruction(OperatorCode o, SizeType c1, SizeType a2) { _instructions.append(ProcedureInstruction(o,c1,a2)); }
    Void new_graded_instruction(OperatorCode o, SizeType a, Int n) { _instructions.append(ProcedureInstruction(o,a,n)); }
  private:
    OutputStream& _write(OutputStream& os) const;
};

//! \brief An algorithmic procedure for computing a function.
//!
//! A Procedure is more efficient to compute than a Formula, since common
//! subexpressions have already been removed. However, it is also more
//! difficult to manipulate, so it should usually only be used when no further
//! manipulations of the function are possible.
template<class Y>
class Vector<Procedure<Y>> {
  public:
    explicit Vector<Procedure<Y>>(const Vector<Formula<Y>>& f);
    explicit Vector<Procedure<Y>>(const Procedure<Y>& p);
    friend OutputStream& operator<<(OutputStream& os, Vector<Procedure<Y>> const& p) { return p._write(os); }
  public:
    SizeType result_size() const { return _results.size(); }
    SizeType temporaries_size() const { return _instructions.size(); }
    Void new_instruction(OperatorCode o, SizeType a) { _instructions.append(ProcedureInstruction(o,a)); }
    Void new_instruction(OperatorCode o, SizeType a, Int n) { _instructions.append(ProcedureInstruction(o,a,n)); }
    Void new_instruction(OperatorCode o, SizeType a1, SizeType a2) { _instructions.append(ProcedureInstruction(o,a1,a2)); }
    Void set_return(SizeType i, SizeType a) { _results[i]=a; }
  public:
    List<Y> _constants;
    List<ProcedureInstruction> _instructions;
    Vector<SizeType> _results;
  private:
    OutputStream& _write(OutputStream& os) const;
};

template<class X, class Y> Void _execute(List<X>& v, const List<ProcedureInstruction>& p, const List<Y>& c, const Vector<X>& x);
template<class X, class Y> Void _backpropagate(Vector<X>& x, List<X>& v, const List<ProcedureInstruction>& p, const List<Y>& c);

// \related Procedure \brief Evaluate a function \a p defined by an algorithmic procedure.
template<class X, class Y> inline X evaluate(const Procedure<Y>& p, const Vector<X>& x) {
    List<X> t(p._instructions.size(),x.zero_element());
    _execute(t,p._instructions,p._constants,x);
    return std::move(t.back());
}

// \related Procedure \brief Evaluate a function \a p defined by an algorithmic procedure.
template<class X, class Y> inline Vector<X> evaluate(const Vector<Procedure<Y>>& p, const Vector<X>& x) {
    List<X> t(p._instructions.size(),x.zero_element());
    _execute(t,p._instructions,p._constants,x);
    Vector<X> r(p.result_size());
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=std::move(t[p._results[i]]); }
    return r;
}

// \related Procedure \brief Execute the procedure \a p storing the results in \a t.
template<class X, class Y> Void execute(List<X>& t, const Procedure<Y>& p, const Vector<X>& x) {
    _execute(t,p._instructions,p._constants,x);
}

template<class X, class Y> Void execute(List<X>& t, const Vector<Procedure<Y>>& p, const Vector<X>& x) {
    _execute(t,p._instructions,p._constants,x);
}

// \related Backpropagate the results of a validated procedure to the inputs \a x.
template<class X, class Y> Void backpropagate(List<X>& t, const Procedure<Y>& p, Vector<X>& x) {
    _backpropagate(t,p._instructions,p._constants,x);
}

template<class X, class Y> Void backpropagate(List<X>& t, const Vector<Procedure<Y>>& p, Vector<X>& x) {
    _backpropagate(t,p._instructions,p._constants,x);
}


} // namespace Ariadne

#endif // ARIADNE_PROCEDURE_H
