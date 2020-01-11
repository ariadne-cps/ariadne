/***************************************************************************
 *            check_number.hpp
 *
 *  Copyright  2009-20  Pieter Collins
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

#include "../test.hpp"
#include "../utility.hpp"
#include "numeric/operators.hpp"

using namespace Ariadne;

namespace Ariadne {


template<class T1, class T2> struct IsEqualityComparible : HasOperator<Equal,T1,T2> { };
template<class T1, class T2> struct IsLessThanCompartible : HasOperator<Less,T1,T2> { };

template<class X, class N, class R=Return<DontCare>> struct HasPow : HasOperator<Pow,X,N,R> { };


template<class N> class CheckNumericTypeConcepts {
  public:
    void check_signed_concept() {
        ARIADNE_TEST_STATIC_ASSERT(HasOperatorReturning<N,OperatorPlus,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperatorReturning<N,OperatorMinus,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Pos,N,Return<N>>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Neg,N,Return<N>>);
    }

    void check_semiring_concept() {
        ARIADNE_TEST_STATIC_ASSERT(IsDefaultConstructible<N>);
        ARIADNE_TEST_STATIC_ASSERT(IsConstructible<N,uint>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorPlus,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorTimes,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Pow,N,uint>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Sqr,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Add,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Mul,N,N>);
    }

    void check_ring_concept() {
        check_semiring_concept();
        check_signed_concept();
        ARIADNE_TEST_STATIC_ASSERT(IsConstructible<N,int>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorMinus,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorMinus,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Neg,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Sub,N,N>);
    }

    void check_field_concept() {
        check_ring_concept();
        if (IsSame<Paradigm<N>,ApproximateTag>::value) {
            ARIADNE_TEST_STATIC_ASSERT(IsConstructible<N,double>);
        }
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<OperatorDivides,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Rec,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Div,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Pow,N,int>);
    };

    void check_transcendental_concept() {
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Sqrt,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Exp,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Log,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Sin,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Cos,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Tan,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Atan,N>);
    }

    void check_absolute_concept() {
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Abs,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Max,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Min,N,N>);
    }

    void check_equality_comparible_concept() {
        ARIADNE_TEST_STATIC_ASSERT(IsEqualityComparible<N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Equal,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Unequal,N,N>);
    }

    void check_order_concept() {
        ARIADNE_TEST_STATIC_ASSERT(IsLessThanCompartible<N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Less,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Leq,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Gtr,N,N>);
        ARIADNE_TEST_STATIC_ASSERT(HasOperator<Geq,N,N>);
    }

    void check_comparible_concept() {
        check_equality_comparible_concept();
        check_order_concept();
    }
};

} // namespace Ariadne

