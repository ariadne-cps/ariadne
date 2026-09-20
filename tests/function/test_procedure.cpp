/***************************************************************************
 *            test_procedure.cpp
 *
 *  Copyright  2010-20  Pieter Collins
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

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <type_traits>

#include "config.hpp"

#include "function/procedure.hpp"
#include "numeric/numeric.hpp"
#include "geometry/interval.hpp"

namespace Ariadne {

// Test-side overloads used to instantiate the generic backward contractors.
// The UpperIntervalType overload mirrors the production implementation in
// procedure.cpp; the FloatDPBounds overload is used only by the isolated Leq
// test, since Leq is not a ProcedureInstruction operator.
inline Void restrict(UpperIntervalType& r, UpperIntervalType const& x) {
    r.set_lower_bound(max(r.lower_bound(),x.lower_bound()));
    r.set_upper_bound(min(r.upper_bound(),x.upper_bound()));
}
inline Void restrict(FloatDPBounds& r, FloatDPBounds const& x) {
    r=refinement(r,x);
}

} // namespace Ariadne

#include "function/procedure.tpl.hpp"
#include "algebra/vector.hpp"
#include "algebra/expansion.hpp"
#include "algebra/evaluate.hpp"
#include "function/formula.hpp"

#include "../test.hpp"

#include "algebra/covector.hpp"
#include "function/function.hpp"
#include "algebra/expansion.inl.hpp"

using namespace std;
using namespace Ariadne;

template<class X> decltype(auto) mag(Covector<X> const& u) { return norm(transpose(u)); }

class TestProcedure
{
    DoublePrecision pr;
  public:
    TestProcedure();
    Void test();
  private:
    Void test_construct_from_formula();
    Void test_construct_from_expansion();
    Void test_evaluate();
    Void test_propagate();
    Void test_backward_contractor_soundness();
    Void test_backward_contractor_witness_preservation();
    Void test_backward_contractor_generated_witnesses();
    Void test_leq_backpropagate_soundness();
    Void test_derivative();
};

TestProcedure::TestProcedure()
{
}

Void TestProcedure::test()
{
    ARIADNE_TEST_CALL(test_construct_from_formula());
    ARIADNE_TEST_CALL(test_construct_from_expansion());
    ARIADNE_TEST_CALL(test_evaluate());
    ARIADNE_TEST_CALL(test_propagate());
    ARIADNE_TEST_CALL(test_backward_contractor_soundness());
    ARIADNE_TEST_CALL(test_backward_contractor_witness_preservation());
    ARIADNE_TEST_CALL(test_backward_contractor_generated_witnesses());
    ARIADNE_TEST_CALL(test_leq_backpropagate_soundness());
    ARIADNE_TEST_CALL(test_derivative());
}


Void TestProcedure::test_construct_from_formula()
{
    ApproximateNumber c(2.0_x);
    ApproximateFormula o(ApproximateFormula::constant(1.0_x));
    ApproximateFormula x(ApproximateFormula::coordinate(0));
    ApproximateFormula y(ApproximateFormula::coordinate(1));

    ApproximateFormula xs=pow(x,2);
    ApproximateFormula ys=pow(y,2);

    ARIADNE_TEST_PRINT(y*c);
    ARIADNE_TEST_PRINT(c*y);
    Vector<ApproximateFormula> f={sqrt(xs+c*ys), atan(y/x), xs*c-ys};
    ARIADNE_TEST_PRINT(f);

    ApproximateProcedure p0(argument_size=2u,f[0]);
    ARIADNE_TEST_EQUALS(p0.argument_size(),2u);
    ARIADNE_TEST_PRINT(p0);
    Vector<ApproximateProcedure> p(argument_size=2u,f);
    ARIADNE_TEST_PRINT(p);
    ARIADNE_TEST_EQUALS(p.result_size(),3u);
    ARIADNE_TEST_EQUALS(p.argument_size(),2u);
}


Void TestProcedure::test_construct_from_expansion()
{
    {
        Expansion<MultiIndex,FloatDPApproximation> e({ {{0,0},1.0_x}, {{1,0},2.0_x}, {{0,2},3.0_x}, {{1,4},4.0_x} },pr);
        ARIADNE_TEST_PRINT(e);
        e.reverse_lexicographic_sort();
        ARIADNE_TEST_PRINT(e);
        Procedure<ApproximateNumber> p(e);
        ARIADNE_TEST_PRINT(p);
        Vector<FloatDPApproximation> x({2.0_x,3.0_x},pr);
        ARIADNE_TEST_PRINT(simple_evaluate(e,x));
        ARIADNE_TEST_PRINT(evaluate(p,x));
        ARIADNE_TEST_EQUAL(evaluate(p,x),simple_evaluate(e,x));
    }

    {
        Expansion<MultiIndex,FloatDPApproximation> e({ {{0,0},1.0_x}, {{1,0},2.0_x}, {{0,1},3.0_x}, {{2,0},4.0_x}, {{1,1},5.0_x}, {{0,2},6.0_x} },pr);
        e.reverse_lexicographic_sort();
        Procedure<ApproximateNumber> p(e);
        ARIADNE_TEST_PRINT(p);
        Vector<FloatDPApproximation> x({2.0_x,3.0_x},pr);
        ARIADNE_TEST_EQUAL(evaluate(p,x),simple_evaluate(e,x));
    }
}


Void TestProcedure::test_evaluate()
{
    ApproximateProcedure p(2);
    p.new_instruction(Var(),0ul);
    p.new_instruction(Var(),1ul);
    p.new_instruction(Pow(),0ul,2);
    p.new_instruction(Sqr(),1ul);
    p.new_instruction(Add(),2ul,3ul);
    p.new_constant(9.0_x);
    p.new_instruction_scalar(Mul(),0ul,4ul);
    p.new_instruction(Sqrt(),5ul);
    ARIADNE_TEST_PRINT(p);

    Vector<FloatDPApproximation> x({3.0_x,4.0_x},pr);
    ARIADNE_TEST_PRINT(x);

    ARIADNE_TEST_EQUALS(evaluate(p,x),15.0_x);
}

Void TestProcedure::test_propagate()
{
    {
        ValidatedProcedure p(2);
        p.new_instruction(Var(),0u);
        p.new_instruction(Var(),1u);
        p.new_instruction(Sqr(),0u);
        p.new_instruction(Sqr(),1u);
        p.new_instruction(Add(),2u,3u);
        p.new_constant(1.125_decimal);
        p.new_instruction_scalar(Mul(),0ul,4ul);
        p.new_instruction(Sqrt(),5u);
        ARIADNE_TEST_PRINT(p);

        UpperBoxType x=ExactBoxType{ {0.25_x,2.0_x}, {0.5_x,3.0_x} };
        ARIADNE_TEST_PRINT(x);

        simple_hull_reduce(x,p,ExactIntervalType(1,1));
        ARIADNE_TEST_PRINT(x);
        simple_hull_reduce(x,p,ExactIntervalType(1,1));
        ARIADNE_TEST_PRINT(x);
    }

    ValidatedFormula x(ValidatedFormula::coordinate(0));
    ValidatedFormula y(ValidatedFormula::coordinate(1));

    Vector<ValidatedFormula> ff={sqrt(sqr(x)+sqr(y)),2*x-y};
    ARIADNE_TEST_PRINT(ff);
    Vector<ValidatedProcedure> pp(argument_size=2,ff);
    ARIADNE_TEST_PRINT(pp);
    UpperBoxType xx=ExactBoxType{ {0.125_x,2.0_x}, {0.25_x,3.0_x} };
    ExactBoxType cc=ExactBoxType{ {1.0_x,1.0_x}, {1.0_x,1.0_x} };
    ARIADNE_TEST_PRINT(xx);
    ARIADNE_TEST_PRINT(evaluate(pp,cast_vector(xx)));
    simple_hull_reduce(xx,pp,cc);
    ARIADNE_TEST_PRINT(xx);
}

Void TestProcedure::test_backward_contractor_soundness()
{
    auto check_periodic_witness = [](auto op, ExactIntervalType domain) {
        ValidatedProcedure p(1);
        p.new_instruction(Var(),0u);
        p.new_instruction(op,0u);
        UpperBoxType x=ExactBoxType({domain});
        simple_hull_reduce(x,p,ExactIntervalType(0,0));
        ARIADNE_TEST_ASSERT(!x[0].is_empty());
    };

    // Each interval contains a non-principal zero of the corresponding
    // periodic function. Backward propagation must not discard that branch.
    check_periodic_witness(Sin(),ExactIntervalType(3.125_x,3.25_x));
    check_periodic_witness(Cos(),ExactIntervalType(4.5_x,5.0_x));
    check_periodic_witness(Tan(),ExactIntervalType(3.125_x,3.25_x));

    {
        ValidatedProcedure p(1);
        p.new_instruction(Var(),0u);
        p.new_instruction(Sqr(),0u);
        UpperBoxType x=ExactBoxType({ExactIntervalType(-1,-1)});
        simple_hull_reduce(x,p,ExactIntervalType(1,1));
        ARIADNE_TEST_ASSERT(!x[0].is_empty());
    }

    {
        ValidatedProcedure p(1);
        p.new_instruction(Var(),0u);
        p.new_instruction(Pow(),0u,2);
        UpperBoxType x=ExactBoxType({ExactIntervalType(-1,-1)});
        simple_hull_reduce(x,p,ExactIntervalType(1,1));
        ARIADNE_TEST_ASSERT(!x[0].is_empty());
    }

    {
        ValidatedProcedure p(2);
        p.new_instruction(Var(),0u);
        p.new_instruction(Var(),1u);
        p.new_instruction(Max(),0u,1u);
        UpperBoxType x=ExactBoxType({ExactIntervalType(0,0),ExactIntervalType(1,1)});
        simple_hull_reduce(x,p,ExactIntervalType(1,1));
        ARIADNE_TEST_ASSERT(!x[0].is_empty());
        ARIADNE_TEST_ASSERT(!x[1].is_empty());
    }

    {
        ValidatedProcedure p(2);
        p.new_instruction(Var(),0u);
        p.new_instruction(Var(),1u);
        p.new_instruction(Min(),0u,1u);
        UpperBoxType x=ExactBoxType({ExactIntervalType(1,1),ExactIntervalType(0,0)});
        simple_hull_reduce(x,p,ExactIntervalType(0,0));
        ARIADNE_TEST_ASSERT(!x[0].is_empty());
        ARIADNE_TEST_ASSERT(!x[1].is_empty());
    }

    {
        ValidatedProcedure p(1);
        p.new_instruction(Var(),0u);
        p.new_constant(1.0_x);
        p.new_instruction_scalar(Max(),0u,0u);
        UpperBoxType x=ExactBoxType({ExactIntervalType(0,0)});
        simple_hull_reduce(x,p,ExactIntervalType(1,1));
        ARIADNE_TEST_ASSERT(!x[0].is_empty());
    }

    {
        ValidatedProcedure p(1);
        p.new_instruction(Var(),0u);
        p.new_constant(0.0_x);
        p.new_instruction_scalar(Min(),0u,0u);
        UpperBoxType x=ExactBoxType({ExactIntervalType(1,1)});
        simple_hull_reduce(x,p,ExactIntervalType(0,0));
        ARIADNE_TEST_ASSERT(!x[0].is_empty());
    }

    // Conservative periodic backward propagation must still propagate an
    // empty forward image so that an impossible codomain empties the domain.
    {
        ValidatedProcedure p(1);
        p.new_instruction(Var(),0u);
        p.new_instruction(Sin(),0u);
        UpperBoxType x=ExactBoxType({ExactIntervalType(0,0)});
        simple_hull_reduce(x,p,ExactIntervalType(2,2));
        ARIADNE_TEST_ASSERT(x[0].is_empty());
    }
}


Void TestProcedure::test_backward_contractor_witness_preservation()
{
    auto U = [](ExactIntervalType const& x) { return UpperIntervalType(x); };
    auto E = [](auto l, auto u) { return ExactIntervalType(l,u); };
    auto S = [](auto x) { return ExactIntervalType(x,x); };

    auto preserves = [](UpperIntervalType const& contracted, ExactIntervalType const& witness) {
        return intersect(contracted,witness);
    };

    auto check_unary = [&](auto op, ExactIntervalType domain, ExactIntervalType output, ExactIntervalType witness) {
        UpperIntervalType contracted=U(domain);
        backpropagate(U(output),op,contracted);
        ARIADNE_TEST_ASSERT(preserves(contracted,witness));
    };

    auto check_binary = [&](auto op,
                            ExactIntervalType lhs_domain, ExactIntervalType rhs_domain,
                            ExactIntervalType output,
                            ExactIntervalType lhs_witness, ExactIntervalType rhs_witness) {
        UpperIntervalType lhs=U(lhs_domain);
        UpperIntervalType rhs=U(rhs_domain);
        backpropagate(U(output),op,lhs,rhs);
        ARIADNE_TEST_ASSERT(preserves(lhs,lhs_witness));
        ARIADNE_TEST_ASSERT(preserves(rhs,rhs_witness));
    };

    auto check_power = [&](ExactIntervalType domain, Int exponent,
                           ExactIntervalType output, ExactIntervalType witness) {
        UpperIntervalType contracted=U(domain);
        backpropagate(U(output),Pow(),contracted,exponent);
        ARIADNE_TEST_ASSERT(preserves(contracted,witness));
    };

    // Unary inverse contractors.
    check_unary(Pos(), S(2), S(2), S(2));
    check_unary(Neg(), S(2), S(-2), S(2));
    check_unary(Rec(), S(2), S(0.5_x), S(2));

    check_unary(Sqr(), S(-2), S(4), S(-2));
    check_unary(Sqr(), S( 2), S(4), S(2));
    check_unary(Sqr(), S( 0), S(0), S(0));

    check_unary(Sqrt(), S(4), S(2), S(4));
    check_unary(Sqrt(), S(0), S(0), S(0));
    check_unary(Exp(), S(0), S(1), S(0));
    check_unary(Log(), S(1), S(0), S(1));

    // Principal and non-principal periodic witnesses.
    check_unary(Sin(), S(0), S(0), S(0));
    check_unary(Cos(), S(0), S(1), S(0));
    check_unary(Tan(), S(0), S(0), S(0));
    check_unary(Sin(), E(3.125_x,3.25_x), S(0), E(3.125_x,3.25_x));
    check_unary(Cos(), E(4.5_x,5.0_x), S(0), E(4.5_x,5.0_x));
    check_unary(Tan(), E(3.125_x,3.25_x), S(0), E(3.125_x,3.25_x));

    // Inverse trigonometric functions on their natural domains.
    check_unary(Asin(), S(0), S(0), S(0));
    check_unary(Acos(), S(1), S(0), S(1));
    check_unary(Atan(), S(0), S(0), S(0));

    // Binary contractors, including signs and zero products.
    check_binary(Add(), S(2), S(3), S(5), S(2), S(3));
    check_binary(Sub(), S(2), S(3), S(-1), S(2), S(3));
    check_binary(Mul(), S(2), S(3), S(6), S(2), S(3));
    check_binary(Mul(), S(-2), S(3), S(-6), S(-2), S(3));
    check_binary(Mul(), S(0), S(3), S(0), S(0), S(3));
    check_binary(Div(), S(6), S(3), S(2), S(6), S(3));
    check_binary(Div(), S(-6), S(3), S(-2), S(-6), S(3));

    check_binary(Max(), S(0), S(1), S(1), S(0), S(1));
    check_binary(Max(), S(1), S(0), S(1), S(1), S(0));
    check_binary(Max(), S(1), S(1), S(1), S(1), S(1));
    check_binary(Min(), S(0), S(1), S(0), S(0), S(1));
    check_binary(Min(), S(1), S(0), S(0), S(1), S(0));
    check_binary(Min(), S(0), S(0), S(0), S(0), S(0));

    // Integer powers: both branches for even powers, signed odd powers,
    // zero exponent, and negative exponents.
    check_power(S(-2),  2, S(4), S(-2));
    check_power(S( 2),  2, S(4), S(2));
    check_power(S(-2),  3, S(-8), S(-2));
    check_power(S( 2),  3, S(8), S(2));
    check_power(S( 2),  0, S(1), S(2));
    check_power(S( 2), -1, S(0.5_x), S(2));
    check_power(S(-2), -2, S(0.25_x), S(-2));

    // Non-degenerate domains exercise endpoint propagation.
    check_unary(Pos(), E(1,3), E(1,3), S(2));
    check_unary(Neg(), E(1,3), E(-3,-1), S(2));
    check_unary(Rec(), E(1,4), E(0.25_x,1.0_x), S(2));
    check_unary(Sqr(), E(-3,3), E(1,4), S(-2));
    check_unary(Sqr(), E(-3,3), E(1,4), S(2));
    check_unary(Sqrt(), E(0,9), E(1,2), S(4));
    check_unary(Exp(), E(-1,2), E(0.5_x,2.0_x), S(0));
    check_unary(Log(), E(0.5_x,3.0_x), E(-0.5_x,1.0_x), S(1));
    check_unary(Sin(), E(3,4), E(-0.5_x,0.5_x), E(3.125_x,3.25_x));
    check_unary(Cos(), E(4,5), E(-0.5_x,0.5_x), E(4.5_x,5.0_x));
    check_unary(Tan(), E(3,4), E(-0.5_x,0.5_x), E(3.125_x,3.25_x));
    check_unary(Asin(), E(-1,1), E(-0.5_x,0.5_x), S(0));
    check_unary(Acos(), E(-1,1), E(0,2), S(1));
    check_unary(Atan(), E(-2,2), E(-1,1), S(0));

    check_binary(Add(), E(0,4), E(1,5), E(4,6), S(2), S(3));
    check_binary(Sub(), E(0,4), E(1,5), E(-2,0), S(2), S(3));
    check_binary(Mul(), E(-3,-1), E(2,4), E(-8,-4), S(-2), S(3));
    check_binary(Div(), E(4,8), E(2,4), E(1,3), S(6), S(3));
    check_binary(Max(), E(-2,2), E(0,4), E(1,3), S(1), S(2));
    check_binary(Min(), E(-2,2), E(0,4), E(-1,1), S(0), S(2));

    check_power(E(-3,3), 2, E(1,4), S(-2));
    check_power(E(-3,3), 2, E(1,4), S(2));
    check_power(E(-3,3), 3, E(-9,-1), S(-2));
    check_power(E(1,4), -1, E(0.25_x,1.0_x), S(2));

    // Scalar-left overloads are distinct implementations.
    {
        UpperIntervalType a=U(S(3));
        backpropagate(U(S(5)),Add(),ValidatedNumber(2),a);
        ARIADNE_TEST_ASSERT(preserves(a,S(3)));
    }
    {
        UpperIntervalType a=U(S(3));
        backpropagate(U(S(-1)),Sub(),ValidatedNumber(2),a);
        ARIADNE_TEST_ASSERT(preserves(a,S(3)));
    }
    {
        UpperIntervalType a=U(S(3));
        backpropagate(U(S(6)),Mul(),ValidatedNumber(2),a);
        ARIADNE_TEST_ASSERT(preserves(a,S(3)));
    }
    {
        UpperIntervalType a=U(S(2));
        backpropagate(U(S(3)),Div(),ValidatedNumber(6),a);
        ARIADNE_TEST_ASSERT(preserves(a,S(2)));
    }
    {
        UpperIntervalType a=U(S(0));
        backpropagate(U(S(1)),Max(),ValidatedNumber(1),a);
        ARIADNE_TEST_ASSERT(preserves(a,S(0)));
    }
    {
        UpperIntervalType a=U(S(1));
        backpropagate(U(S(0)),Min(),ValidatedNumber(0),a);
        ARIADNE_TEST_ASSERT(preserves(a,S(1)));
    }

    // Scalar-right overloads are not used by ScalarProcedureInstruction today,
    // but remain part of the generic backpropagation API.
    {
        UpperIntervalType a=U(S(3));
        backpropagate(U(S(5)),Add(),a,ValidatedNumber(2));
        ARIADNE_TEST_ASSERT(preserves(a,S(3)));
    }
    {
        UpperIntervalType a=U(S(3));
        backpropagate(U(S(1)),Sub(),a,ValidatedNumber(2));
        ARIADNE_TEST_ASSERT(preserves(a,S(3)));
    }
    {
        UpperIntervalType a=U(S(3));
        backpropagate(U(S(6)),Mul(),a,ValidatedNumber(2));
        ARIADNE_TEST_ASSERT(preserves(a,S(3)));
    }
    {
        UpperIntervalType a=U(S(6));
        backpropagate(U(S(3)),Div(),a,ValidatedNumber(2));
        ARIADNE_TEST_ASSERT(preserves(a,S(6)));
    }
    {
        UpperIntervalType a=U(S(0));
        backpropagate(U(S(1)),Max(),a,ValidatedNumber(1));
        ARIADNE_TEST_ASSERT(preserves(a,S(0)));
    }
    {
        UpperIntervalType a=U(S(1));
        backpropagate(U(S(0)),Min(),a,ValidatedNumber(0));
        ARIADNE_TEST_ASSERT(preserves(a,S(1)));
    }
}

Void TestProcedure::test_backward_contractor_generated_witnesses()
{
    auto U = [](ExactIntervalType const& x) { return UpperIntervalType(x); };
    auto E = [](auto l, auto u) { return ExactIntervalType(l,u); };
    auto S = [](auto x) { return ExactIntervalType(x,x); };

    auto preserves = [](UpperIntervalType const& contracted, ExactIntervalType const& witness) {
        return intersect(contracted,witness);
    };

    auto check_unary = [&](auto op, ExactIntervalType domain, ExactIntervalType witness) {
        UpperIntervalType contracted=U(domain);
        UpperIntervalType witness_interval=U(witness);
        UpperIntervalType output=op(witness_interval);
        backpropagate(output,op,contracted);
        ARIADNE_TEST_ASSERT(preserves(contracted,witness));
    };

    auto check_binary = [&](auto op,
                            ExactIntervalType lhs_domain, ExactIntervalType rhs_domain,
                            ExactIntervalType lhs_witness, ExactIntervalType rhs_witness) {
        UpperIntervalType lhs=U(lhs_domain);
        UpperIntervalType rhs=U(rhs_domain);
        UpperIntervalType output=op(U(lhs_witness),U(rhs_witness));
        backpropagate(output,op,lhs,rhs);
        ARIADNE_TEST_ASSERT(preserves(lhs,lhs_witness));
        ARIADNE_TEST_ASSERT(preserves(rhs,rhs_witness));
    };

    auto check_power = [&](ExactIntervalType domain, ExactIntervalType witness, Int exponent) {
        UpperIntervalType contracted=U(domain);
        UpperIntervalType output=pow(U(witness),exponent);
        backpropagate(output,Pow(),contracted,exponent);
        ARIADNE_TEST_ASSERT(preserves(contracted,witness));
    };

    // Generic unary grid over exactly representable witnesses.
    for (auto witness : {S(-2),S(-1),S(-0.5_x),S(0),S(0.5_x),S(1),S(2)}) {
        check_unary(Pos(),E(-3,3),witness);
        check_unary(Neg(),E(-3,3),witness);
        check_unary(Sqr(),E(-3,3),witness);
        check_unary(Exp(),E(-3,3),witness);
        check_unary(Sin(),E(-3,3),witness);
        check_unary(Cos(),E(-3,3),witness);
        check_unary(Tan(),E(-3,3),witness);
        check_unary(Atan(),E(-3,3),witness);
    }

    for (auto witness : {S(-2),S(-1),S(-0.5_x),S(0.5_x),S(1),S(2)}) {
        check_unary(Rec(),E(-3,3),witness);
        for (Int exponent=-3; exponent<=4; ++exponent) {
            check_power(E(-3,3),witness,exponent);
        }
    }

    for (auto witness : {S(0),S(0.5_x),S(1),S(2),S(4)}) {
        check_unary(Sqrt(),E(0,4),witness);
    }

    for (auto witness : {S(0.5_x),S(1),S(2)}) {
        check_unary(Log(),E(0.5_x,2),witness);
    }

    for (auto witness : {S(-1),S(-0.5_x),S(0),S(0.5_x),S(1)}) {
        check_unary(Asin(),E(-1,1),witness);
        check_unary(Acos(),E(-1,1),witness);
    }

    // Exhaustive small binary grid. Every output is evaluated from the
    // singleton witnesses, so each tested tuple is valid by construction.
    for (auto lhs_witness : {S(-2),S(-1),S(0),S(1),S(2)}) {
        for (auto rhs_witness : {S(-2),S(-1),S(0),S(1),S(2)}) {
            check_binary(Add(),E(-3,3),E(-3,3),lhs_witness,rhs_witness);
            check_binary(Sub(),E(-3,3),E(-3,3),lhs_witness,rhs_witness);
            check_binary(Mul(),E(-3,3),E(-3,3),lhs_witness,rhs_witness);
            check_binary(Max(),E(-3,3),E(-3,3),lhs_witness,rhs_witness);
            check_binary(Min(),E(-3,3),E(-3,3),lhs_witness,rhs_witness);
        }
    }

    for (auto lhs_witness : {S(-2),S(-1),S(0),S(1),S(2)}) {
        for (auto rhs_witness : {S(-2),S(-1),S(1),S(2)}) {
            check_binary(Div(),E(-3,3),E(-3,3),lhs_witness,rhs_witness);
        }
    }
}


Void TestProcedure::test_leq_backpropagate_soundness()
{
    // Leq is currently not representable by ProcedureInstruction, so this
    // overload is not on the simple_hull_reduce execution path.
    ARIADNE_TEST_CONCEPT((not std::is_constructible_v<BinaryElementaryOperator,Leq>));
    ARIADNE_TEST_CONCEPT((not std::is_constructible_v<ProcedureInstruction,Leq,SizeType,SizeType>));

    // Independently exercise the overload itself. The valid witness
    // a1=0, a2=1 satisfies a1<=a2 and must not be removed by propagation.
    FloatDPBounds r(0,dp);
    FloatDPBounds a1(0,0,dp);
    FloatDPBounds a2(1,1,dp);

    backpropagate(r,Leq(),a1,a2);

    ARIADNE_TEST_ASSERT(not inconsistent(a1,FloatDPBounds(0,dp)));
    ARIADNE_TEST_ASSERT(not inconsistent(a2,FloatDPBounds(1,dp)));
}


Void TestProcedure::test_derivative()
{
    typedef ApproximateTag P;
    typedef Number<P> Y;
    typedef FloatDPApproximation X;

    Y c(2);
    Formula<Y> o(Formula<Y>::constant(1));
    Formula<Y> x(Formula<Y>::coordinate(0));
    Formula<Y> y(Formula<Y>::coordinate(1));

    auto xs=sqr(x); auto ys=sqr(y);
    Formula<Y> e=sqrt(xs+c*ys)+sin(y/x)*xs*c-ys;
    ARIADNE_TEST_PRINT(e);
    ScalarMultivariateFunction<P> f(EuclideanDomain(2),e);
    ARIADNE_TEST_PRINT(f);
    Procedure<Y> p(f);
    ARIADNE_TEST_PRINT(p);

    X zero(0,dp);

    Vector<X> q({2,1},dp);
    Vector<X> s({-1,3},dp);
    Vector<Differential<X>> ds=Differential<X>::variable(1,2,zero,0)*s+q;
    ARIADNE_TEST_WITHIN(gradient(p,q),f.gradient(q),8e-16);
    ARIADNE_TEST_EQUALS(hessian(p,q,s),f(ds).hessian().get(0,0));
}

Int main() {
    TestProcedure().test();
    return ARIADNE_TEST_FAILURES;
}

