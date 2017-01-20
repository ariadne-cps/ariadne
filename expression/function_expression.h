/***************************************************************************
 *            function_expression.h
 *
 *  Copyright 2008-15  Pieter Collins
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

/*! \file function_expression.h
 *  \brief Expressions created as wrappers around functions
 */

#ifndef ARIADNE_FUNCTION_EXPRESSION_H
#define ARIADNE_FUNCTION_EXPRESSION_H

#include "expression/expression.h"
#include "expression/expression_set.h"
#include "expression/space.h"

#include "function/function.h"
#include "function/formula.h"
#include "function/function_model.h"
#include "function/taylor_function.h"

namespace Ariadne {

template<class P, class D, class C> class FunctionExpression;
template<class P, class D=BoxDomain> using ScalarFunctionExpression = FunctionExpression<P,D,IntervalDomain>;
template<class P, class D=BoxDomain> using VectorFunctionExpression = FunctionExpression<P,D,BoxDomain>;

using ValidatedScalarFunctionExpression = FunctionExpression<ValidatedTag,BoxDomain,IntervalDomain>;
using ValidatedVectorFunctionExpression = FunctionExpression<ValidatedTag,BoxDomain,BoxDomain>;

using VariableIntervalDomain = ExactFloat64VariableInterval;
using VariablesBoxDomain = ExactFloat64VariablesBox;

inline BoxDomain preimage(Projection const& prj, BoxDomain const& bx) {
    BoxDomain dbx(prj.argument_size(),IntervalDomain(-infty,+infty));
    for(SizeType i=0; i!=prj.result_size(); ++i) {
        IntervalDomain& dbxj=dbx[prj.index(i)]; dbxj=intersection(dbxj,bx[i]); }
    return dbx;
}

VariableIntervalDomain _make_domain(RealVariable const& arg, IntervalDomain const& dom) {
    return VariableIntervalDomain(arg,dom);
}

VariablesBoxDomain _make_domain(Vector<RealVariable> const& args, BoxDomain const& dom) {
    Map<RealVariable,Float64ExactInterval> domain;
    for(SizeType i=0; i!=args.size(); ++i) { domain[args[i]]=dom[i]; }
    return domain;
}

template<class P,class C> class FunctionExpression<P,BoxDomain,C> {
    typedef BoxDomain D;
  public:
    template<class Y> using Argument = typename ElementTraits<D>::template Type<Y>;
    template<class Y> using Result = typename ElementTraits<C>::template Type<Y>;
  private:
    Argument<RealVariable> _vars;
    Function<P,D,C> _f;
  public:
    FunctionExpression(Function<P,D,C> const& f, Argument<RealVariable> const& vars) : _vars(vars), _f(f) {
        ARIADNE_PRECONDITION(_f.argument_size()==_vars.size()); }
    Function<P,D,C> const& function() const { return this->_f; }
    Argument<RealVariable> const& variables() const { return this->_vars; }
    Set<RealVariable> arguments() const { Set<RealVariable> args; for(SizeType i=0; i!=_vars.size(); ++i) { args.insert(_vars[i]); } return args; }
    auto domain() const -> decltype(_make_domain(_vars,_f.domain())) { return _make_domain(_vars,_f.domain()); }
    template<class PP> friend OutputStream& operator<<(OutputStream& os, ScalarFunctionExpression<PP> const& sfe);
};

template<class P> OutputStream& operator<<(OutputStream& os, ScalarFunctionExpression<P> const& e) {
    os << e._f;
    for(SizeType i=0; i!=e._vars.size(); ++i) {
        os << (i==0?'[':',') << "x[" << i << "]=" << e._vars[i];
    } return os << "]";
}

template<class P> ScalarFunctionExpression<P> FunctionFacade<P,BoxDomain,IntervalDomain>::operator() (Vector<RealVariable> const& vars) const {
    return ScalarFunctionExpression<P>(static_cast<Function<P,BoxDomain,IntervalDomain>const&>(*this),vars);
}

template<class P> ScalarFunctionExpression<P> evaluate(ScalarFunction<P> const& f, Vector<RealVariable> const& vars) {
    return ScalarFunctionExpression<P>(f,vars);
}

template<class F> TaylorModel<ValidatedTag,F> compose(const TaylorModel<ValidatedTag,F>& x, Projection const& prj) {
    TaylorModel<ValidatedTag,F> r(prj.argument_size(),x.sweeper());
    r.expansion().reserve(x.number_of_nonzeros());
    MultiIndex ra(r.argument_size());
    for(auto xiter=x.begin(); xiter!=x.end(); ++xiter) {
        MultiIndex const& xa=xiter->key();
        auto const& xc=xiter->data();
        for(SizeType i=0; i!=r.argument_size(); ++i) {
            ra[i]=xa[prj.index(i)];
        }
        r.expansion().append(ra,xc);
    }
    r.cleanup();
    return std::move(r);
}

template<class M> FunctionPatch<M> compose(FunctionPatch<M> const& f, Projection const& prj) {
    SizeType as=prj.argument_size();
    auto f_dom=f.domain();
    BoxDomain dom=preimage(prj,f_dom);

    bool has_strict_subdomain = false;
    for(SizeType i=0; i!=prj.result_size(); ++i) {
        if(not same(f_dom[i],dom[prj.index(i)])) { has_strict_subdomain=true; }
    }
    if (has_strict_subdomain) {
        Vector<FunctionPatch<M>> id=f.create_coordinates(dom);
        return compose(f,prj(id));
    } else {
        return FunctionPatch<M>(dom,compose(f.model(),prj));
    }
}

template<class P> ScalarFunctionModel<P> compose(ScalarFunctionModel<P> const& f, Projection const& prj) {
    auto fpp = std::dynamic_pointer_cast<const ScalarTaylorFunction>(f.managed_pointer());
    if(fpp) {
        return compose(ScalarTaylorFunction(fpp),prj);
    }
    SizeType as=prj.argument_size();
    auto f_dom=f.domain();
    BoxDomain dom=preimage(prj,f_dom);
    Vector<ScalarFunctionModel<P>> id=f.create_coordinates(dom);
    Vector<ScalarFunctionModel<P>> pid=prj(id);
    VectorFunctionModel<P> vpid(pid.array());
    return compose(f,vpid);
}

template<class P> ScalarFunction<P> compose(ScalarFunction<P> const& f, Projection const& prj) {
    auto fmp = std::dynamic_pointer_cast<const ScalarFunctionModelInterface<P>>(f.managed_pointer());
    if(fmp) {
        return compose(ScalarFunctionModel<P>(fmp),prj);
    }
    SizeType as=prj.argument_size();
    auto f_dom=f.domain();
    BoxDomain dom=preimage(prj,f_dom);
    Vector<ScalarFunction<P>> id(ScalarFunction<P>::coordinates(dom));
    VectorFunction<P> pid=prj(id);
    return compose(f,pid);
}

RealScalarFunction compose(RealScalarFunction const& f, Projection const& prj) {
    SizeType as=prj.argument_size();
    auto f_dom=f.domain();
    BoxDomain dom=preimage(prj,f_dom);
    Vector<RealScalarFunction> id(RealScalarFunction::coordinates(dom));
    RealVectorFunction pid=prj(id);
    return compose(f,pid);
}


template<class P> ScalarFunction<P> make_function(RealSpace args, ScalarFunctionExpression<P> const& expr) {
    ARIADNE_PRECONDITION(subset(expr.arguments(),args.variables()));
    auto expr_vars=expr.variables();
    auto expr_func=expr.function();

    Array<SizeType> indices(expr_vars.size());
    for(SizeType i=0; i!=expr_vars.size(); ++i) {
        indices[i]=args.index(expr_vars[i]);
    }
    Projection projection(args.size(),indices);

    return compose(expr_func,projection);
}

} // namespace Ariadne

#endif /* ARIADNE_FUNCTION_EXPRESSION_H */
