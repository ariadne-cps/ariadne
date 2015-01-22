/***************************************************************************
 *            check_numeric.cc
 *
 *  Copyright 2006-14  Alberto Casagrande, Pieter Collins
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


#include "utility/module.h"

#include "config.h"
#include "numeric/logical.h"
#include "numeric/number.h"
#include "numeric/float.h"
#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/real.h"
#include "numeric/float64.h"
#include "numeric/floatmp.h"

#include "test.h"
#include "utility.h"
#include "check_number.h"

using namespace std;
using namespace Ariadne;

namespace Ariadne {
template<class T> String class_name();

namespace Detail {

template<class OP, class X1, class X2, class = Fallback> struct SafeTypedef {
    typedef Fallback Type; };
template<class OP, class X1, class X2> struct SafeTypedef<OP,X1,X2, EnableIf<IsConvertible<decltype(declval<OP>()(declval<X1>(),declval<X2>())),DontCare>,Fallback>> {
    typedef decltype(declval<OP>()(declval<X1>(),declval<X2>()) ) Type; };

template<class X, class = Fallback> struct SafeNegationTypedef { typedef Fallback Type; };
template<class X> struct SafeNegationTypedef<X, EnableIf<IsConvertible<decltype(-declval<X>()),DontCare>,Fallback>> { typedef decltype(-declval<X>()) Type; };
template<class X1, class X2, class = Fallback> struct SafeSumTypedef { typedef Fallback Type; };
template<class X1, class X2> struct SafeSumTypedef<X1,X2, EnableIf<IsConvertible<decltype(declval<X1>()+declval<X2>()),DontCare>,Fallback>> { typedef decltype(declval<X1>()+declval<X2>()) Type; };
template<class X1, class X2, class = Fallback> struct SafeDifferenceTypedef { typedef Fallback Type; };
template<class X1, class X2> struct SafeDifferenceTypedef<X1,X2, EnableIf<IsConvertible<decltype(declval<X1>()-declval<X2>()),DontCare>,Fallback>> { typedef decltype(declval<X1>()-declval<X2>()) Type; };
template<class X1, class X2, class = Fallback> struct SafeProductTypedef { typedef Fallback Type; };
template<class X1, class X2> struct SafeProductTypedef<X1,X2, EnableIf<IsConvertible<decltype(declval<X1>()*declval<X2>()),DontCare>,Fallback>> { typedef decltype(declval<X1>()*declval<X2>()) Type; };
template<class X1, class X2, class = Fallback> struct SafeQuotientTypedef { typedef Fallback Type; };
template<class X1, class X2> struct SafeQuotientTypedef<X1,X2, EnableIf<IsConvertible<decltype(declval<X1>()/declval<X2>()),DontCare>,Fallback>> { typedef decltype(declval<X1>()/declval<X2>()) Type; };

template<class X1, class X2, class = Fallback> struct SafeEqualsTypedef { typedef Fallback Type; };
template<class X1, class X2> struct SafeEqualsTypedef<X1,X2, EnableIf<IsConvertible<decltype(declval<X1>()==declval<X2>()),DontCare>,Fallback>> { typedef decltype(declval<X1>()==declval<X2>()) Type; };
template<class X1, class X2, class = Fallback> struct SafeLessTypedef { typedef Fallback Type; };
template<class X1, class X2> struct SafeLessTypedef<X1,X2, EnableIf<IsConvertible<decltype(declval<X1>()< declval<X2>()),DontCare>,Fallback>> { typedef decltype(declval<X1>()< declval<X2>()) Type;
};

}

template<class OP, class X1, class X2=X1> using SafeType = typename Detail::SafeTypedef<OP,X1,X2>::Type;

template<class X> using SafeNegationType = typename Detail::SafeNegationTypedef<X>::Type;
template<class X1, class X2=X1> using SafeSumType = typename Detail::SafeSumTypedef<X1,X2>::Type;
template<class X1, class X2=X1> using SafeDifferenceType = typename Detail::SafeDifferenceTypedef<X1,X2>::Type;
template<class X1, class X2=X1> using SafeProductType = typename Detail::SafeProductTypedef<X1,X2>::Type;
template<class X1, class X2=X1> using SafeQuotientType = typename Detail::SafeQuotientTypedef<X1,X2>::Type;

template<class X1, class X2=X1> using SafeArithmeticType = SafeSumType<SafeProductType<X1,X2>>;
template<class X1, class X2=X1> using SafeEqualsType = typename Detail::SafeEqualsTypedef<X1,X2>::Type;
template<class X1, class X2=X1> using SafeLessType = typename Detail::SafeLessTypedef<X1,X2>::Type;

template<class X1, class X2> struct IsEquivalent : And<IsConvertible<X1,X2>,IsConvertible<X2,X1>> { };
template<class T, class U> struct IsNotConvertible : Not<IsConvertible<T,U>> { };
template<class T, class U> struct IsExplicitlyConvertible : IsConstructible<U,T> { };
template<class T, class U> struct IsNotExplicitlyConvertible : Not<IsConstructible<U,T>> { };
template<class T, class U> struct IsStrictlyConvertible : And<IsConstructible<U,T>,IsNotConvertible<T,U>> { };

struct Plus { template<class A1, class A2> auto operator() (A1 const& a1, A2 const& a2) -> decltype(operator+(a1,a2)) { return operator+(a1,a2); } };
struct Minus { template<class A1, class A2> auto operator() (A1 const& a1, A2 const& a2) -> decltype(operator-(a1,a2)) { return operator-(a1,a2); } };
struct Times { template<class A1, class A2> auto operator() (A1 const& a1, A2 const& a2) -> decltype(operator*(a1,a2)) { return operator*(a1,a2); } };
struct Divides { template<class A1, class A2> auto operator() (A1 const& a1, A2 const& a2) -> decltype(operator/(a1,a2)) { return operator/(a1,a2); } };

template<class OP> String op_name();
template<> String op_name<Add>() { return "add"; }
template<> String op_name<Sub>() { return "sub"; }
template<> String op_name<Mul>() { return "mul"; }
template<> String op_name<Div>() { return "div"; }

template<> String op_name<Plus>() { return "operator+"; }
template<> String op_name<Minus>() { return "operator-"; }
template<> String op_name<Times>() { return "operator*"; }
template<> String op_name<Divides>() { return "operator/"; }

void output_check_result(Bool p, String e, String r, String op, String a1, String a2) {
    std::cout << "Checking " << op << "(" << a1 << "," << a2 << ") -> " << e << ":";
    if(p) {
        std::cout << " OK\n";
    } else {
        std::cout << " Failed! Actual type is " << r << "\n";
        std::cerr << "ERROR: " << op << "(" << a1 << "," << a2 << ") -> " << r << "; expected " << e << "\n";
    }
}

// Check if the result of OP(A1,A2) has expected type E
template<class E, class OP, class A1, class A2> void chk() {
    typedef SafeType<OP,A1,A2> R;
    Bool p=IsSame<R,E>();
    output_check_result(p,class_name<E>(),class_name<R>(),op_name<OP>(),class_name<A1>(),class_name<A2>());
}

template<class T, class F> void check_explicitly_constructible() {
    std::cout << "Checking " << class_name<T>() << " is constructible from " << class_name<F>() << ":";
    if(IsConstructible<T,F>::value and not IsConvertible<F,T>::value) {
        std::cout << " OK\n";
    } else {
        std::cout << " Failed!\n";
        if(IsConvertible<F,T>::value) {
            std::cerr << "ERROR: " << class_name<T>() << " is convertible from  " << class_name<F>() << ", but only explict constructor allowed.\n";
        } else {
            std::cerr << "ERROR: " << class_name<T>() << " is not constructible from " << class_name<F>() << ".\n";
        }
    }
}

template<class T, class F> void check_convertible() {
    std::cout << "Checking " << class_name<T>() << " is convertible from " << class_name<F>() << ":";
    if(IsConvertible<F,T>::value) {
        std::cout << " OK\n";
    } else {
        std::cout << " Failed!\n";
        std::cerr << "ERROR: " << class_name<T>() << " is not convertible from " << class_name<F>() << ".\n";
    }
}

template<class T, class F> void check_not_constructible() {
    std::cout << "Checking " << class_name<T>() << " is not constructible from " << class_name<F>() << ":";
    if(not IsConstructible<T,F>::value) {
        std::cout << " OK\n";
    } else {
        std::cout << " Failed!\n";
        std::cerr << "ERROR: " << class_name<T>() << " is constructible from " << class_name<F>() << "\n";
    }
}

struct Cid { } ; struct Eid { }; struct Nid { };
template<class R, class OP, class X> void chk() { return chk<R,X>(OP()); }
template<class R, class X> void chk(Cid) { check_convertible<R,X>(); }
template<class R, class X> void chk(Eid) { check_explicitly_constructible<R,X>(); }
template<class R, class X> void chk(Nid) { check_not_constructible<R,X>(); }


#define ARIADNE_CLASS_NAME(Class) \
    template<> String class_name<Class>() { return #Class; }

ARIADNE_CLASS_NAME(Fallback);

ARIADNE_CLASS_NAME(Approximate);
ARIADNE_CLASS_NAME(Lower);
ARIADNE_CLASS_NAME(Upper);
ARIADNE_CLASS_NAME(Bounded);
ARIADNE_CLASS_NAME(Validated);
ARIADNE_CLASS_NAME(Exact);

ARIADNE_CLASS_NAME(Logical<Approximate>);
ARIADNE_CLASS_NAME(Logical<Lower>);
ARIADNE_CLASS_NAME(Logical<Upper>);
ARIADNE_CLASS_NAME(Logical<Validated>);
ARIADNE_CLASS_NAME(Logical<Effective>);
ARIADNE_CLASS_NAME(Logical<Exact>);

ARIADNE_CLASS_NAME(bool);
ARIADNE_CLASS_NAME(uint);
ARIADNE_CLASS_NAME(Fuzzy);
ARIADNE_CLASS_NAME(NegSierpinski);
ARIADNE_CLASS_NAME(Sierpinski);
ARIADNE_CLASS_NAME(Tribool);
ARIADNE_CLASS_NAME(Boolean);

ARIADNE_CLASS_NAME(Integer);
ARIADNE_CLASS_NAME(Rational);
ARIADNE_CLASS_NAME(Real);

//ARIADNE_CLASS_NAME(ApproximateNumber)
//ARIADNE_CLASS_NAME(LowerNumber)
//ARIADNE_CLASS_NAME(UpperNumber)
//ARIADNE_CLASS_NAME(ValidatedNumber)
//ARIADNE_CLASS_NAME(EffectiveNumber)
//ARIADNE_CLASS_NAME(ExactNumber)

ARIADNE_CLASS_NAME(ApproximateFloat64);
ARIADNE_CLASS_NAME(LowerFloat64);
ARIADNE_CLASS_NAME(UpperFloat64);
ARIADNE_CLASS_NAME(ValidatedFloat64);
ARIADNE_CLASS_NAME(ErrorFloat64);
ARIADNE_CLASS_NAME(ExactFloat64);

ARIADNE_CLASS_NAME(ApproximateFloat);
ARIADNE_CLASS_NAME(LowerFloat);
ARIADNE_CLASS_NAME(UpperFloat);
ARIADNE_CLASS_NAME(ValidatedFloat);
//ARIADNE_CLASS_NAME(ErrorFloat);
ARIADNE_CLASS_NAME(ExactFloat);

#undef ARIADNE_CLASS_NAME

typedef int I; typedef double D;
typedef Integer Z; typedef Rational Q; typedef Real R;
typedef ExactNumber ExN; typedef EffectiveNumber EfN; typedef ValidatedNumber VaN;
typedef UpperNumber UpN; typedef LowerNumber LoN; typedef ApproximateNumber ApN;
typedef ExactFloat ExF; typedef MetricFloat MeF; typedef BoundedFloat BoF;
typedef UpperFloat UpF; typedef LowerFloat LoF; typedef ApproximateFloat ApF;

typedef decltype(declval<ExF>() + declval<ExF>()) VaF;
typedef decltype(declval<UpF>() * declval<UpF>()) PrF;
//typedef ValidatedFloat64Type VaF;

} // namespace Ariadne

class CheckNumeric
{
  public:
    void check();
  private:
    typedef SafeSumType<ExactFloat,ExactFloat> ExactFloatArithmeticType;

    void notifications();
    void check_conversions();
    void check_addition();
    void check_subtraction();
    void check_multiplication();
    void check_division();
    void check_comparison();
    void check_equality();
};

void CheckNumeric::check()
{
    ARIADNE_TEST_CALL(notifications());
    ARIADNE_TEST_CALL(check_conversions());
    ARIADNE_TEST_CALL(check_addition());
    ARIADNE_TEST_WARN("Skipping check of sub, mul, div and cmp functions.");
    return;
    ARIADNE_TEST_CALL(check_subtraction());
    ARIADNE_TEST_CALL(check_multiplication());
    ARIADNE_TEST_CALL(check_division());
    ARIADNE_TEST_CALL(check_comparison());
    ARIADNE_TEST_CALL(check_equality());
}

String to_str(bool b) { return b?"true":"false"; }

void CheckNumeric::notifications()
{

    // Operations on ExactFloat: display what is being used.
    ARIADNE_TEST_NOTIFY(String("Conversion double -> ApproximateNumber: ")+to_str(IsConvertible<double,ApproximateNumber>::value));
    ARIADNE_TEST_NOTIFY(String("Conversion double -> ApproximateFloat64: ")+to_str(IsConvertible<double,ApproximateFloat64>::value));
    ARIADNE_TEST_NOTIFY(String("Conversion double -> ExactFloat64: ")+to_str(IsConvertible<double,ExactFloat64>::value));
    ARIADNE_TEST_NOTIFY(String("Construction double -> ExactNumber: ")+to_str(IsConstructible<ExactNumber,double>::value));
    ARIADNE_TEST_NOTIFY(String("Construction double -> ExactFloat64: ")+to_str(IsConstructible<ExactFloat64,double>::value));
    ARIADNE_TEST_NOTIFY(String("Conversion int -> ExactFloat64: ")+to_str(IsConvertible<int,ExactFloat64>::value));
    ARIADNE_TEST_NOTIFY(String("Construction Integer -> ExactFloat64: ")+to_str(IsConstructible<ExactFloat64,Integer>::value));
    ARIADNE_TEST_NOTIFY((String("ExactFloat + ExactFloat -> ")+class_name<SafeSumType<ExactFloat,ExactFloat>>()));
    ARIADNE_TEST_NOTIFY((String("ExactFloat + double -> ")+class_name<SafeSumType<ExactFloat,double>>()));
    ARIADNE_TEST_NOTIFY((String("ExactNumber + double -> ")+class_name<SafeSumType<ExactNumber,double>>()));
    ARIADNE_TEST_NOTIFY((String("ApproximateNumber + double -> ")+class_name<SafeSumType<ApproximateNumber,double>>()));
    ARIADNE_TEST_NOTIFY((String("MetricFloat + BoundedFloat -> ")+class_name<SafeSumType<MetricFloat,BoundedFloat>>()));
    ARIADNE_TEST_NOTIFY((String("UpperNumber * UpperNumber -> ")+class_name<SafeProductType<UpperNumber,UpperNumber>>()));
    ARIADNE_TEST_NOTIFY((String("UpperFloat * UpperFloat -> ")+class_name<SafeProductType<UpperFloat,UpperFloat>>()));
    ARIADNE_TEST_NOTIFY((String("Rational == double -> ")+class_name<SafeEqualsType<Rational,double>>()));
    ARIADNE_TEST_NOTIFY((String("Rational < double -> ")+class_name<SafeLessType<Rational,double>>()));

    ARIADNE_TEST_STATIC_ASSERT(IsStronger<Paradigm<ExactFloatArithmeticType>,Validated>)
}

void CheckNumeric::check_conversions() {
    // Cid is conversion, Eid is explicit construction, Nid is no construction

    // Conversions involving int
    chk<ErrorFloat,Cid,uint>(); chk<ErrorFloat,Eid,int>(); chk<ErrorFloat,Eid,double>();

    // Concrete numbers
    chk<D,Cid,D>(); chk<D,Cid,I>(); chk<D,Nid,Z>(); chk<D,Nid,Q>(); chk<D,Nid,R>();
    chk<I,Cid,D>(); chk<I,Cid,I>(); chk<I,Nid,Z>(); chk<I,Nid,Q>(); chk<I,Nid,R>();
    chk<Z,Nid,D>(); chk<Z,Cid,I>(); chk<Z,Cid,Z>(); chk<Z,Nid,Q>(); chk<Z,Nid,R>();
    chk<Q,Nid,D>(); chk<Q,Cid,I>(); chk<Q,Cid,Z>(); chk<Q,Cid,Q>(); chk<Q,Nid,R>();
    chk<R,Eid,D>(); chk<R,Cid,I>(); chk<R,Cid,Z>(); chk<R,Cid,Q>(); chk<R,Cid,R>();

    // Generic numbers
    chk<ExN,Cid,ExN>(); chk<ExN,Nid,EfN>(); chk<ExN,Nid,VaN>(); chk<ExN,Nid,UpN>(); chk<ExN,Nid,LoN>(); chk<ExN,Nid,ApN>();
    chk<EfN,Cid,ExN>(); chk<EfN,Cid,EfN>(); chk<EfN,Nid,VaN>(); chk<EfN,Nid,UpN>(); chk<EfN,Nid,LoN>(); chk<EfN,Nid,ApN>();
    chk<VaN,Cid,ExN>(); chk<VaN,Cid,EfN>(); chk<VaN,Cid,VaN>(); chk<VaN,Nid,UpN>(); chk<VaN,Nid,LoN>(); chk<VaN,Nid,ApN>();
    chk<UpN,Cid,ExN>(); chk<UpN,Cid,EfN>(); chk<UpN,Cid,VaN>(); chk<UpN,Cid,UpN>(); chk<UpN,Nid,LoN>(); chk<UpN,Nid,ApN>();
    chk<LoN,Cid,ExN>(); chk<LoN,Cid,EfN>(); chk<LoN,Cid,VaN>(); chk<LoN,Nid,UpN>(); chk<LoN,Cid,LoN>(); chk<LoN,Nid,ApN>();
    chk<ApN,Cid,ExN>(); chk<ApN,Cid,EfN>(); chk<ApN,Cid,VaN>(); chk<ApN,Cid,UpN>(); chk<ApN,Cid,LoN>(); chk<ApN,Cid,ApN>();

    // Mixed Generic # Concrete numbers
    chk<ExN,Nid,D>(); chk<ExN,Cid,I>(); chk<ExN,Cid,Z>(); chk<ExN,Cid,Q>(); chk<ExN,Nid,R>();
    chk<EfN,Nid,D>(); chk<EfN,Cid,I>(); chk<EfN,Cid,Z>(); chk<EfN,Cid,Q>(); chk<EfN,Cid,R>();
    chk<VaN,Nid,D>(); chk<VaN,Cid,I>(); chk<VaN,Cid,Z>(); chk<VaN,Cid,Q>(); chk<VaN,Cid,R>();
    chk<UpN,Nid,D>(); chk<UpN,Cid,I>(); chk<UpN,Cid,Z>(); chk<UpN,Cid,Q>(); chk<UpN,Cid,R>();
    chk<LoN,Nid,D>(); chk<LoN,Cid,I>(); chk<LoN,Cid,Z>(); chk<LoN,Cid,Q>(); chk<LoN,Cid,R>();
    chk<ApN,Cid,D>(); chk<ApN,Cid,I>(); chk<ApN,Cid,Z>(); chk<ApN,Cid,Q>(); chk<ApN,Cid,R>();

    // Mixed Concrete # Generic numbers
    chk<D,Nid,ExN>(); chk<D,Nid,EfN>(); chk<D,Nid,VaN>(); chk<D,Nid,UpN>(); chk<D,Nid,LoN>(); chk<D,Nid,ApN>();
    chk<I,Nid,ExN>(); chk<I,Nid,EfN>(); chk<I,Nid,VaN>(); chk<I,Nid,UpN>(); chk<I,Nid,LoN>(); chk<I,Nid,ApN>();
    chk<Z,Nid,ExN>(); chk<Z,Nid,EfN>(); chk<Z,Nid,VaN>(); chk<Z,Nid,UpN>(); chk<Z,Nid,LoN>(); chk<Z,Nid,ApN>();
    chk<Q,Nid,ExN>(); chk<Q,Nid,EfN>(); chk<Q,Nid,VaN>(); chk<Q,Nid,UpN>(); chk<Q,Nid,LoN>(); chk<Q,Nid,ApN>();
    chk<R,Nid,ExN>(); chk<R,Nid,EfN>(); chk<R,Nid,VaN>(); chk<R,Nid,UpN>(); chk<R,Nid,LoN>(); chk<R,Nid,ApN>();

    // Float numbers
    chk<ExF,Cid,ExF>(); chk<ExF,Nid,MeF>(); chk<ExF,Nid,BoF>(); chk<ExF,Nid,UpF>(); chk<ExF,Nid,LoF>(); chk<ExF,Nid,ApF>();
    chk<MeF,Cid,ExF>(); chk<MeF,Cid,MeF>(); chk<MeF,Cid,BoF>(); chk<MeF,Nid,UpF>(); chk<MeF,Nid,LoF>(); chk<MeF,Nid,ApF>();
    chk<BoF,Cid,ExF>(); chk<BoF,Cid,MeF>(); chk<BoF,Cid,BoF>(); chk<BoF,Nid,UpF>(); chk<BoF,Nid,LoF>(); chk<BoF,Nid,ApF>();
    chk<UpF,Cid,ExF>(); chk<UpF,Cid,MeF>(); chk<UpF,Cid,BoF>(); chk<UpF,Cid,UpF>(); chk<UpF,Nid,LoF>(); chk<UpF,Nid,ApF>();
    chk<LoF,Cid,ExF>(); chk<LoF,Cid,MeF>(); chk<LoF,Cid,BoF>(); chk<LoF,Nid,UpF>(); chk<LoF,Cid,LoF>(); chk<LoF,Nid,ApF>();
    chk<ApF,Cid,ExF>(); chk<ApF,Cid,MeF>(); chk<ApF,Cid,BoF>(); chk<ApF,Cid,UpF>(); chk<ApF,Cid,LoF>(); chk<ApF,Cid,ApF>();

    // Mixed Float # Generic and Generic # Float
    chk<ExF,Nid,ExN>(); chk<ExF,Nid,EfN>(); chk<ExF,Nid,VaN>(); chk<ExF,Nid,UpN>(); chk<ExF,Nid,LoN>(); chk<ExF,Nid,ApN>();
    chk<MeF,Eid,ExN>(); chk<MeF,Eid,EfN>(); chk<MeF,Eid,VaN>(); chk<MeF,Nid,UpN>(); chk<MeF,Nid,LoN>(); chk<MeF,Nid,ApN>();
    chk<BoF,Eid,ExN>(); chk<BoF,Eid,EfN>(); chk<BoF,Eid,VaN>(); chk<BoF,Nid,UpN>(); chk<BoF,Nid,LoN>(); chk<BoF,Nid,ApN>();
    chk<UpF,Eid,ExN>(); chk<UpF,Eid,EfN>(); chk<UpF,Eid,VaN>(); chk<UpF,Eid,UpN>(); chk<UpF,Nid,LoN>(); chk<UpF,Nid,ApN>();
    chk<LoF,Eid,ExN>(); chk<LoF,Eid,EfN>(); chk<LoF,Eid,VaN>(); chk<LoF,Nid,UpN>(); chk<LoF,Eid,LoN>(); chk<LoF,Nid,ApN>();
    chk<ApF,Eid,ExN>(); chk<ApF,Eid,EfN>(); chk<ApF,Eid,VaN>(); chk<ApF,Eid,UpN>(); chk<ApF,Eid,LoN>(); chk<ApF,Eid,ApN>();

    chk<ExN,Cid,ExF>(); chk<ExN,Nid,MeF>(); chk<ExN,Nid,BoF>(); chk<ExN,Nid,UpF>(); chk<ExN,Nid,LoF>(); chk<ExN,Nid,ApF>();
    chk<EfN,Cid,ExF>(); chk<EfN,Nid,MeF>(); chk<EfN,Nid,BoF>(); chk<EfN,Nid,UpF>(); chk<EfN,Nid,LoF>(); chk<EfN,Nid,ApF>();
    chk<VaN,Cid,ExF>(); chk<VaN,Cid,MeF>(); chk<VaN,Cid,BoF>(); chk<VaN,Nid,UpF>(); chk<VaN,Nid,LoF>(); chk<VaN,Nid,ApF>();
    chk<UpN,Cid,ExF>(); chk<UpN,Cid,MeF>(); chk<UpN,Cid,BoF>(); chk<UpN,Cid,UpF>(); chk<UpN,Nid,LoF>(); chk<UpN,Nid,ApF>();
    chk<LoN,Cid,ExF>(); chk<LoN,Cid,MeF>(); chk<LoN,Cid,BoF>(); chk<LoN,Nid,UpF>(); chk<LoN,Cid,LoF>(); chk<LoN,Nid,ApF>();
    chk<ApN,Cid,ExF>(); chk<ApN,Cid,MeF>(); chk<ApN,Cid,BoF>(); chk<ApN,Cid,UpF>(); chk<ApN,Cid,LoF>(); chk<ApN,Cid,ApF>();

    // Mixed Float # Concrete and Concrete # Float
    chk<ExF,Eid,D>(); chk<ExF,Cid,I>(); chk<ExF,Eid,Z>(); chk<ExF,Nid,Q>(); chk<ExF,Nid,R>();
    chk<MeF,Eid,D>(); chk<MeF,Cid,I>(); chk<MeF,Eid,Z>(); chk<MeF,Eid,Q>(); chk<MeF,Eid,R>();
    chk<BoF,Eid,D>(); chk<BoF,Cid,I>(); chk<BoF,Eid,Z>(); chk<BoF,Eid,Q>(); chk<BoF,Eid,R>();
    chk<UpF,Eid,D>(); chk<UpF,Cid,I>(); chk<UpF,Eid,Z>(); chk<UpF,Eid,Q>(); chk<UpF,Eid,R>();
    chk<LoF,Eid,D>(); chk<LoF,Cid,I>(); chk<LoF,Eid,Z>(); chk<LoF,Eid,Q>(); chk<LoF,Eid,R>();
    chk<ApF,Cid,D>(); chk<ApF,Cid,I>(); chk<ApF,Eid,Z>(); chk<ApF,Eid,Q>(); chk<ApF,Eid,R>();

    chk<D,Nid,ExF>(); chk<D,Nid,MeF>(); chk<D,Nid,BoF>(); chk<D,Nid,UpF>(); chk<D,Nid,LoF>(); chk<D,Nid,ApF>();
    chk<I,Nid,ExF>(); chk<I,Nid,MeF>(); chk<I,Nid,BoF>(); chk<I,Nid,UpF>(); chk<I,Nid,LoF>(); chk<I,Nid,ApF>();
    chk<Z,Nid,ExF>(); chk<Z,Nid,MeF>(); chk<Z,Nid,BoF>(); chk<Z,Nid,UpF>(); chk<Z,Nid,LoF>(); chk<Z,Nid,ApF>();
    chk<Q,Eid,ExF>(); chk<Q,Nid,MeF>(); chk<Q,Nid,BoF>(); chk<Q,Nid,UpF>(); chk<Q,Nid,LoF>(); chk<Q,Nid,ApF>();
    chk<R,Eid,ExF>(); chk<R,Nid,MeF>(); chk<R,Nid,BoF>(); chk<R,Nid,UpF>(); chk<R,Nid,LoF>(); chk<R,Nid,ApF>();
}

void CheckNumeric::check_addition() {
    // Concrete numbers
    chk< D ,Add,D,D>(); chk<D,Add,D,I>(); chk<ApN,Add,D,Z>(); chk<ApN,Add,D,Q>(); chk<ApN,Add,D,R>();
    chk< D ,Add,I,D>(); chk<I,Add,I,I>(); chk<Z,Add,I,Z>(); chk<Q,Add,I,Q>(); chk<R,Add,I,R>();
    chk<ApN,Add,Z,D>(); chk<Z,Add,Z,I>(); chk<Z,Add,Z,Z>(); chk<Q,Add,Z,Q>(); chk<R,Add,Z,R>();
    chk<ApN,Add,Q,D>(); chk<Q,Add,Q,I>(); chk<Q,Add,Q,Z>(); chk<Q,Add,Q,Q>(); chk<R,Add,Q,R>();
    chk<ApN,Add,R,D>(); chk<R,Add,R,I>(); chk<R,Add,R,Z>(); chk<R,Add,R,Q>(); chk<R,Add,R,R>();

    // Generic numbers
    chk<ExN,Add,ExN,ExN>(); chk<EfN,Add,ExN,EfN>(); chk<VaN,Add,ExN,VaN>(); chk<UpN,Add,ExN,UpN>(); chk<LoN,Add,ExN,LoN>(); chk<ApN,Add,ExN,ApN>();
    chk<EfN,Add,EfN,ExN>(); chk<EfN,Add,EfN,EfN>(); chk<VaN,Add,EfN,VaN>(); chk<UpN,Add,EfN,UpN>(); chk<LoN,Add,EfN,LoN>(); chk<ApN,Add,EfN,ApN>();
    chk<VaN,Add,VaN,ExN>(); chk<VaN,Add,VaN,EfN>(); chk<VaN,Add,VaN,VaN>(); chk<UpN,Add,VaN,UpN>(); chk<LoN,Add,VaN,LoN>(); chk<ApN,Add,VaN,ApN>();
    chk<UpN,Add,UpN,ExN>(); chk<UpN,Add,UpN,EfN>(); chk<UpN,Add,UpN,VaN>(); chk<UpN,Add,UpN,UpN>(); chk<ApN,Add,UpN,LoN>(); chk<ApN,Add,UpN,ApN>();
    chk<LoN,Add,LoN,ExN>(); chk<LoN,Add,LoN,EfN>(); chk<LoN,Add,LoN,VaN>(); chk<ApN,Add,LoN,UpN>(); chk<LoN,Add,LoN,LoN>(); chk<ApN,Add,LoN,ApN>();
    chk<ApN,Add,ApN,ExN>(); chk<ApN,Add,ApN,EfN>(); chk<ApN,Add,ApN,VaN>(); chk<ApN,Add,ApN,UpN>(); chk<ApN,Add,ApN,LoN>(); chk<ApN,Add,ApN,ApN>();

    // Mixed Generic # Concrete numbers
    chk<ApN,Add,ExN,D>(); chk<ExN,Add,ExN,I>(); chk<ExN,Add,ExN,Z>(); chk<ExN,Add,ExN,Q>(); chk<EfN,Add,ExN,R>();
    chk<ApN,Add,EfN,D>(); chk<EfN,Add,EfN,I>(); chk<EfN,Add,EfN,Z>(); chk<EfN,Add,EfN,Q>(); chk<EfN,Add,EfN,R>();
    chk<ApN,Add,VaN,D>(); chk<VaN,Add,VaN,I>(); chk<VaN,Add,VaN,Z>(); chk<VaN,Add,VaN,Q>(); chk<VaN,Add,VaN,R>();
    chk<ApN,Add,UpN,D>(); chk<UpN,Add,UpN,I>(); chk<UpN,Add,UpN,Z>(); chk<UpN,Add,UpN,Q>(); chk<UpN,Add,UpN,R>();
    chk<ApN,Add,LoN,D>(); chk<LoN,Add,LoN,I>(); chk<LoN,Add,LoN,Z>(); chk<LoN,Add,LoN,Q>(); chk<LoN,Add,LoN,R>();
    chk<ApN,Add,ApN,D>(); chk<ApN,Add,ApN,I>(); chk<ApN,Add,ApN,Z>(); chk<ApN,Add,ApN,Q>(); chk<ApN,Add,ApN,R>();

    // Mixed Concrete # Generic numbers
    chk<ApN,Add,D,ExN>(); chk<ApN,Add,D,EfN>(); chk<ApN,Add,D,VaN>(); chk<ApN,Add,D,UpN>(); chk<ApN,Add,D,LoN>(); chk<ApN,Add,D,ApN>();
    chk<ExN,Add,I,ExN>(); chk<EfN,Add,I,EfN>(); chk<VaN,Add,I,VaN>(); chk<UpN,Add,I,UpN>(); chk<LoN,Add,I,LoN>(); chk<ApN,Add,I,ApN>();
    chk<ExN,Add,Z,ExN>(); chk<EfN,Add,Z,EfN>(); chk<VaN,Add,Z,VaN>(); chk<UpN,Add,Z,UpN>(); chk<LoN,Add,Z,LoN>(); chk<ApN,Add,Z,ApN>();
    chk<ExN,Add,Q,ExN>(); chk<EfN,Add,Q,EfN>(); chk<VaN,Add,Q,VaN>(); chk<UpN,Add,Q,UpN>(); chk<LoN,Add,Q,LoN>(); chk<ApN,Add,Q,ApN>();
    chk<EfN,Add,R,ExN>(); chk<EfN,Add,R,EfN>(); chk<VaN,Add,R,VaN>(); chk<UpN,Add,R,UpN>(); chk<LoN,Add,R,LoN>(); chk<ApN,Add,R,ApN>();

    // Float numbers
    chk<VaF,Add,ExF,ExF>(); chk<MeF,Add,ExF,MeF>(); chk<BoF,Add,ExF,BoF>(); chk<UpF,Add,ExF,UpF>(); chk<LoF,Add,ExF,LoF>(); chk<ApF,Add,ExF,ApF>();
    chk<MeF,Add,MeF,ExF>(); chk<MeF,Add,MeF,MeF>(); chk<VaF,Add,MeF,BoF>(); chk<UpF,Add,MeF,UpF>(); chk<LoF,Add,MeF,LoF>(); chk<ApF,Add,MeF,ApF>();
    chk<BoF,Add,BoF,ExF>(); chk<VaF,Add,BoF,MeF>(); chk<BoF,Add,BoF,BoF>(); chk<UpF,Add,BoF,UpF>(); chk<LoF,Add,BoF,LoF>(); chk<ApF,Add,BoF,ApF>();
    chk<UpF,Add,UpF,ExF>(); chk<UpF,Add,UpF,MeF>(); chk<UpF,Add,UpF,BoF>(); chk<UpF,Add,UpF,UpF>(); chk<ApF,Add,UpF,LoF>(); chk<ApF,Add,UpF,ApF>();
    chk<LoF,Add,LoF,ExF>(); chk<LoF,Add,LoF,MeF>(); chk<LoF,Add,LoF,BoF>(); chk<ApF,Add,LoF,UpF>(); chk<LoF,Add,LoF,LoF>(); chk<ApF,Add,LoF,ApF>();
    chk<ApF,Add,ApF,ExF>(); chk<ApF,Add,ApF,MeF>(); chk<ApF,Add,ApF,BoF>(); chk<ApF,Add,ApF,UpF>(); chk<ApF,Add,ApF,LoF>(); chk<ApF,Add,ApF,ApF>();

    // Mixed Float # Generic and Generic # Float
    chk<VaF,Add,ExF,ExN>(); chk<VaF,Add,ExF,EfN>(); chk<VaF,Add,ExF,VaN>(); chk<UpF,Add,ExF,UpN>(); chk<LoF,Add,ExF,LoN>(); chk<ApF,Add,ExF,ApN>();
    chk<MeF,Add,MeF,ExN>(); chk<MeF,Add,MeF,EfN>(); chk<MeF,Add,MeF,VaN>(); chk<UpF,Add,MeF,UpN>(); chk<LoF,Add,MeF,LoN>(); chk<ApF,Add,MeF,ApN>();
    chk<BoF,Add,BoF,ExN>(); chk<BoF,Add,BoF,EfN>(); chk<BoF,Add,BoF,VaN>(); chk<UpF,Add,BoF,UpN>(); chk<LoF,Add,BoF,LoN>(); chk<ApF,Add,BoF,ApN>();
    chk<UpF,Add,UpF,ExN>(); chk<UpF,Add,UpF,EfN>(); chk<UpF,Add,UpF,VaN>(); chk<UpF,Add,UpF,UpN>(); chk<ApF,Add,UpF,LoN>(); chk<ApF,Add,UpF,ApN>();
    chk<LoF,Add,LoF,ExN>(); chk<LoF,Add,LoF,EfN>(); chk<LoF,Add,LoF,VaN>(); chk<ApF,Add,LoF,UpN>(); chk<LoF,Add,LoF,LoN>(); chk<ApF,Add,LoF,ApN>();
    chk<ApF,Add,ApF,ExN>(); chk<ApF,Add,ApF,EfN>(); chk<ApF,Add,ApF,VaN>(); chk<ApF,Add,ApF,UpN>(); chk<ApF,Add,ApF,LoN>(); chk<ApF,Add,ApF,ApN>();

    chk<VaF,Add,ExN,ExF>(); chk<MeF,Add,ExN,MeF>(); chk<BoF,Add,ExN,BoF>(); chk<UpF,Add,ExN,UpF>(); chk<LoF,Add,ExN,LoF>(); chk<ApF,Add,ExN,ApF>();
    chk<VaF,Add,EfN,ExF>(); chk<MeF,Add,EfN,MeF>(); chk<BoF,Add,EfN,BoF>(); chk<UpF,Add,EfN,UpF>(); chk<LoF,Add,EfN,LoF>(); chk<ApF,Add,EfN,ApF>();
    chk<VaF,Add,VaN,ExF>(); chk<MeF,Add,VaN,MeF>(); chk<BoF,Add,VaN,BoF>(); chk<UpF,Add,VaN,UpF>(); chk<LoF,Add,VaN,LoF>(); chk<ApF,Add,VaN,ApF>();
    chk<UpF,Add,UpN,ExF>(); chk<UpF,Add,UpN,MeF>(); chk<UpF,Add,UpN,BoF>(); chk<UpF,Add,UpN,UpF>(); chk<ApF,Add,UpN,LoF>(); chk<ApF,Add,UpN,ApF>();
    chk<LoF,Add,LoN,ExF>(); chk<LoF,Add,LoN,MeF>(); chk<LoF,Add,LoN,BoF>(); chk<ApF,Add,LoN,UpF>(); chk<LoF,Add,LoN,LoF>(); chk<ApF,Add,LoN,ApF>();
    chk<ApF,Add,ApN,ExF>(); chk<ApF,Add,ApN,MeF>(); chk<ApF,Add,ApN,BoF>(); chk<ApF,Add,ApN,UpF>(); chk<ApF,Add,ApN,LoF>(); chk<ApF,Add,ApN,ApF>();

    // Mixed Float # Concrete and Concrete # Float
    chk<ApF,Add,ExF,D>(); chk<VaF,Add,ExF,I>(); chk<VaF,Add,ExF,Z>(); chk<VaF,Add,ExF,Q>(); chk<VaF,Add,ExF,R>();
    chk<ApF,Add,MeF,D>(); chk<MeF,Add,MeF,I>(); chk<MeF,Add,MeF,Z>(); chk<MeF,Add,MeF,Q>(); chk<MeF,Add,MeF,R>();
    chk<ApF,Add,BoF,D>(); chk<BoF,Add,BoF,I>(); chk<BoF,Add,BoF,Z>(); chk<BoF,Add,BoF,Q>(); chk<BoF,Add,BoF,R>();
    chk<ApF,Add,UpF,D>(); chk<UpF,Add,UpF,I>(); chk<UpF,Add,UpF,Z>(); chk<UpF,Add,UpF,Q>(); chk<UpF,Add,UpF,R>();
    chk<ApF,Add,LoF,D>(); chk<LoF,Add,LoF,I>(); chk<LoF,Add,LoF,Z>(); chk<LoF,Add,LoF,Q>(); chk<LoF,Add,LoF,R>();
    chk<ApF,Add,ApF,D>(); chk<ApF,Add,ApF,I>(); chk<ApF,Add,ApF,Z>(); chk<ApF,Add,ApF,Q>(); chk<ApF,Add,ApF,R>();

    chk<ApF,Add,D,ExF>(); chk<ApF,Add,D,MeF>(); chk<ApF,Add,D,BoF>(); chk<ApF,Add,D,UpF>(); chk<ApF,Add,D,LoF>(); chk<ApF,Add,D,ApF>();
    chk<VaF,Add,I,ExF>(); chk<MeF,Add,I,MeF>(); chk<BoF,Add,I,BoF>(); chk<UpF,Add,I,UpF>(); chk<LoF,Add,I,LoF>(); chk<ApF,Add,I,ApF>();
    chk<VaF,Add,Z,ExF>(); chk<MeF,Add,Z,MeF>(); chk<BoF,Add,Z,BoF>(); chk<UpF,Add,Z,UpF>(); chk<LoF,Add,Z,LoF>(); chk<ApF,Add,Z,ApF>();
    chk<VaF,Add,Q,ExF>(); chk<MeF,Add,Q,MeF>(); chk<BoF,Add,Q,BoF>(); chk<UpF,Add,Q,UpF>(); chk<LoF,Add,Q,LoF>(); chk<ApF,Add,Q,ApF>();
    chk<VaF,Add,R,ExF>(); chk<MeF,Add,R,MeF>(); chk<BoF,Add,R,BoF>(); chk<UpF,Add,R,UpF>(); chk<LoF,Add,R,LoF>(); chk<ApF,Add,R,ApF>();
}

void CheckNumeric::check_subtraction() {

    // Concrete numbers
    chk< D ,Sub,D,D>(); chk<D,Sub,D,I>(); chk<ApN,Sub,D,Z>(); chk<ApN,Sub,D,Q>(); chk<ApN,Sub,D,R>();
    chk< D ,Sub,I,D>(); chk<I,Sub,I,I>(); chk<Z,Sub,I,Z>(); chk<Q,Sub,I,Q>(); chk<R,Sub,I,R>();
    chk<ApN,Sub,Z,D>(); chk<Z,Sub,Z,I>(); chk<Z,Sub,Z,Z>(); chk<Q,Sub,Z,Q>(); chk<R,Sub,Z,R>();
    chk<ApN,Sub,Q,D>(); chk<Q,Sub,Q,I>(); chk<Q,Sub,Q,Z>(); chk<Q,Sub,Q,Q>(); chk<R,Sub,Q,R>();
    chk<ApN,Sub,R,D>(); chk<R,Sub,R,I>(); chk<R,Sub,R,Z>(); chk<R,Sub,R,Q>(); chk<R,Sub,R,R>();

    // Generic numbers
    chk<ExN,Sub,ExN,ExN>(); chk<EfN,Sub,ExN,EfN>(); chk<VaN,Sub,ExN,VaN>(); chk<LoN,Sub,ExN,UpN>(); chk<UpN,Sub,ExN,LoN>(); chk<ApN,Sub,ExN,ApN>();
    chk<EfN,Sub,EfN,ExN>(); chk<EfN,Sub,EfN,EfN>(); chk<VaN,Sub,EfN,VaN>(); chk<LoN,Sub,EfN,UpN>(); chk<UpN,Sub,EfN,LoN>(); chk<ApN,Sub,EfN,ApN>();
    chk<VaN,Sub,VaN,ExN>(); chk<VaN,Sub,VaN,EfN>(); chk<VaN,Sub,VaN,VaN>(); chk<LoN,Sub,VaN,UpN>(); chk<UpN,Sub,VaN,LoN>(); chk<ApN,Sub,VaN,ApN>();
    chk<UpN,Sub,UpN,ExN>(); chk<UpN,Sub,UpN,EfN>(); chk<UpN,Sub,UpN,VaN>(); chk<ApN,Sub,UpN,UpN>(); chk<UpN,Sub,UpN,LoN>(); chk<ApN,Sub,UpN,ApN>();
    chk<LoN,Sub,LoN,ExN>(); chk<LoN,Sub,LoN,EfN>(); chk<LoN,Sub,LoN,VaN>(); chk<LoN,Sub,LoN,UpN>(); chk<ApN,Sub,LoN,LoN>(); chk<ApN,Sub,LoN,ApN>();
    chk<ApN,Sub,ApN,ExN>(); chk<ApN,Sub,ApN,EfN>(); chk<ApN,Sub,ApN,VaN>(); chk<ApN,Sub,ApN,UpN>(); chk<ApN,Sub,ApN,LoN>(); chk<ApN,Sub,ApN,ApN>();

    // Mixed Generic # Concrete numbers
    chk<ApN,Sub,ExN,D>(); chk<ExN,Sub,ExN,I>(); chk<ExN,Sub,ExN,Z>(); chk<ExN,Sub,ExN,Q>(); chk<EfN,Sub,ExN,R>();
    chk<ApN,Sub,EfN,D>(); chk<EfN,Sub,EfN,I>(); chk<EfN,Sub,EfN,Z>(); chk<EfN,Sub,EfN,Q>(); chk<EfN,Sub,EfN,R>();
    chk<ApN,Sub,VaN,D>(); chk<VaN,Sub,VaN,I>(); chk<VaN,Sub,VaN,Z>(); chk<VaN,Sub,VaN,Q>(); chk<VaN,Sub,VaN,R>();
    chk<ApN,Sub,UpN,D>(); chk<UpN,Sub,UpN,I>(); chk<UpN,Sub,UpN,Z>(); chk<UpN,Sub,UpN,Q>(); chk<UpN,Sub,UpN,R>();
    chk<ApN,Sub,LoN,D>(); chk<LoN,Sub,LoN,I>(); chk<LoN,Sub,LoN,Z>(); chk<LoN,Sub,LoN,Q>(); chk<LoN,Sub,LoN,R>();
    chk<ApN,Sub,ApN,D>(); chk<ApN,Sub,ApN,I>(); chk<ApN,Sub,ApN,Z>(); chk<ApN,Sub,ApN,Q>(); chk<ApN,Sub,ApN,R>();

    // Mixed Concrete # Generic numbers
    chk<ApN,Sub,D,ExN>(); chk<ApN,Sub,D,EfN>(); chk<ApN,Sub,D,VaN>(); chk<ApN,Sub,D,UpN>(); chk<ApN,Sub,D,LoN>(); chk<ApN,Sub,D,ApN>();
    chk<ExN,Sub,I,ExN>(); chk<EfN,Sub,I,EfN>(); chk<VaN,Sub,I,VaN>(); chk<LoN,Sub,I,UpN>(); chk<UpN,Sub,I,LoN>(); chk<ApN,Sub,I,ApN>();
    chk<ExN,Sub,Z,ExN>(); chk<EfN,Sub,Z,EfN>(); chk<VaN,Sub,Z,VaN>(); chk<LoN,Sub,Z,UpN>(); chk<UpN,Sub,Z,LoN>(); chk<ApN,Sub,Z,ApN>();
    chk<ExN,Sub,Q,ExN>(); chk<EfN,Sub,Q,EfN>(); chk<VaN,Sub,Q,VaN>(); chk<LoN,Sub,Q,UpN>(); chk<UpN,Sub,Q,LoN>(); chk<ApN,Sub,Q,ApN>();
    chk<EfN,Sub,R,ExN>(); chk<EfN,Sub,R,EfN>(); chk<VaN,Sub,R,VaN>(); chk<LoN,Sub,R,UpN>(); chk<UpN,Sub,R,LoN>(); chk<ApN,Sub,R,ApN>();

    // Float numbers
    chk<VaF,Sub,ExF,ExF>(); chk<MeF,Sub,ExF,MeF>(); chk<BoF,Sub,ExF,BoF>(); chk<LoF,Sub,ExF,UpF>(); chk<UpF,Sub,ExF,LoF>(); chk<ApF,Sub,ExF,ApF>();
    chk<MeF,Sub,MeF,ExF>(); chk<MeF,Sub,MeF,MeF>(); chk<VaF,Sub,MeF,BoF>(); chk<LoF,Sub,MeF,UpF>(); chk<UpF,Sub,MeF,LoF>(); chk<ApF,Sub,MeF,ApF>();
    chk<BoF,Sub,BoF,ExF>(); chk<VaF,Sub,BoF,MeF>(); chk<BoF,Sub,BoF,BoF>(); chk<LoF,Sub,BoF,UpF>(); chk<UpF,Sub,BoF,LoF>(); chk<ApF,Sub,BoF,ApF>();
    chk<UpF,Sub,UpF,ExF>(); chk<UpF,Sub,UpF,MeF>(); chk<UpF,Sub,UpF,BoF>(); chk<ApF,Sub,UpF,UpF>(); chk<UpF,Sub,UpF,LoF>(); chk<ApF,Sub,UpF,ApF>();
    chk<LoF,Sub,LoF,ExF>(); chk<LoF,Sub,LoF,MeF>(); chk<LoF,Sub,LoF,BoF>(); chk<LoF,Sub,LoF,UpF>(); chk<ApF,Sub,LoF,LoF>(); chk<ApF,Sub,LoF,ApF>();
    chk<ApF,Sub,ApF,ExF>(); chk<ApF,Sub,ApF,MeF>(); chk<ApF,Sub,ApF,BoF>(); chk<ApF,Sub,ApF,UpF>(); chk<ApF,Sub,ApF,LoF>(); chk<ApF,Sub,ApF,ApF>();

    // Mixed Float # Generic and Generic # Float
    chk<VaF,Sub,ExF,ExN>(); chk<VaF,Sub,ExF,EfN>(); chk<VaF,Sub,ExF,VaN>(); chk<LoF,Sub,ExF,UpN>(); chk<UpF,Sub,ExF,LoN>(); chk<ApF,Sub,ExF,ApN>();
    chk<MeF,Sub,MeF,ExN>(); chk<MeF,Sub,MeF,EfN>(); chk<MeF,Sub,MeF,VaN>(); chk<LoF,Sub,MeF,UpN>(); chk<UpF,Sub,MeF,LoN>(); chk<ApF,Sub,MeF,ApN>();
    chk<BoF,Sub,BoF,ExN>(); chk<BoF,Sub,BoF,EfN>(); chk<BoF,Sub,BoF,VaN>(); chk<LoF,Sub,BoF,UpN>(); chk<UpF,Sub,BoF,LoN>(); chk<ApF,Sub,BoF,ApN>();
    chk<UpF,Sub,UpF,ExN>(); chk<UpF,Sub,UpF,EfN>(); chk<UpF,Sub,UpF,VaN>(); chk<ApF,Sub,UpF,UpN>(); chk<UpF,Sub,UpF,LoN>(); chk<ApF,Sub,UpF,ApN>();
    chk<LoF,Sub,LoF,ExN>(); chk<LoF,Sub,LoF,EfN>(); chk<LoF,Sub,LoF,VaN>(); chk<LoF,Sub,LoF,UpN>(); chk<ApF,Sub,LoF,LoN>(); chk<ApF,Sub,LoF,ApN>();
    chk<ApF,Sub,ApF,ExN>(); chk<ApF,Sub,ApF,EfN>(); chk<ApF,Sub,ApF,VaN>(); chk<ApF,Sub,ApF,UpN>(); chk<ApF,Sub,ApF,LoN>(); chk<ApF,Sub,ApF,ApN>();

    chk<VaF,Sub,ExN,ExF>(); chk<MeF,Sub,ExN,MeF>(); chk<BoF,Sub,ExN,BoF>(); chk<LoF,Sub,ExN,UpF>(); chk<UpF,Sub,ExN,LoF>(); chk<ApF,Sub,ExN,ApF>();
    chk<VaF,Sub,EfN,ExF>(); chk<MeF,Sub,EfN,MeF>(); chk<BoF,Sub,EfN,BoF>(); chk<LoF,Sub,EfN,UpF>(); chk<UpF,Sub,EfN,LoF>(); chk<ApF,Sub,EfN,ApF>();
    chk<VaF,Sub,VaN,ExF>(); chk<MeF,Sub,VaN,MeF>(); chk<BoF,Sub,VaN,BoF>(); chk<LoF,Sub,VaN,UpF>(); chk<UpF,Sub,VaN,LoF>(); chk<ApF,Sub,VaN,ApF>();
    chk<UpF,Sub,UpN,ExF>(); chk<UpF,Sub,UpN,MeF>(); chk<UpF,Sub,UpN,BoF>(); chk<ApF,Sub,UpN,UpF>(); chk<UpF,Sub,UpN,LoF>(); chk<ApF,Sub,UpN,ApF>();
    chk<LoF,Sub,LoN,ExF>(); chk<LoF,Sub,LoN,MeF>(); chk<LoF,Sub,LoN,BoF>(); chk<LoF,Sub,LoN,UpF>(); chk<ApF,Sub,LoN,LoF>(); chk<ApF,Sub,LoN,ApF>();
    chk<ApF,Sub,ApN,ExF>(); chk<ApF,Sub,ApN,MeF>(); chk<ApF,Sub,ApN,BoF>(); chk<ApF,Sub,ApN,UpF>(); chk<ApF,Sub,ApN,LoF>(); chk<ApF,Sub,ApN,ApF>();

    // Mixed Float # Concrete and Concrete # Float
    chk<ApF,Sub,ExF,D>(); chk<VaF,Sub,ExF,I>(); chk<VaF,Sub,ExF,Z>(); chk<VaF,Sub,ExF,Q>(); chk<VaF,Sub,ExF,R>();
    chk<ApF,Sub,MeF,D>(); chk<MeF,Sub,MeF,I>(); chk<MeF,Sub,MeF,Z>(); chk<MeF,Sub,MeF,Q>(); chk<MeF,Sub,MeF,R>();
    chk<ApF,Sub,BoF,D>(); chk<BoF,Sub,BoF,I>(); chk<BoF,Sub,BoF,Z>(); chk<BoF,Sub,BoF,Q>(); chk<BoF,Sub,BoF,R>();
    chk<ApF,Sub,UpF,D>(); chk<UpF,Sub,UpF,I>(); chk<UpF,Sub,UpF,Z>(); chk<UpF,Sub,UpF,Q>(); chk<UpF,Sub,UpF,R>();
    chk<ApF,Sub,LoF,D>(); chk<LoF,Sub,LoF,I>(); chk<LoF,Sub,LoF,Z>(); chk<LoF,Sub,LoF,Q>(); chk<LoF,Sub,LoF,R>();
    chk<ApF,Sub,ApF,D>(); chk<ApF,Sub,ApF,I>(); chk<ApF,Sub,ApF,Z>(); chk<ApF,Sub,ApF,Q>(); chk<ApF,Sub,ApF,R>();

    chk<ApF,Sub,D,ExF>(); chk<ApF,Sub,D,MeF>(); chk<ApF,Sub,D,BoF>(); chk<ApF,Sub,D,UpF>(); chk<ApF,Sub,D,LoF>(); chk<ApF,Sub,D,ApF>();
    chk<VaF,Sub,I,ExF>(); chk<MeF,Sub,I,MeF>(); chk<BoF,Sub,I,BoF>(); chk<LoF,Sub,I,UpF>(); chk<UpF,Sub,I,LoF>(); chk<ApF,Sub,I,ApF>();
    chk<VaF,Sub,Z,ExF>(); chk<MeF,Sub,Z,MeF>(); chk<BoF,Sub,Z,BoF>(); chk<LoF,Sub,Z,UpF>(); chk<UpF,Sub,Z,LoF>(); chk<ApF,Sub,Z,ApF>();
    chk<VaF,Sub,Q,ExF>(); chk<MeF,Sub,Q,MeF>(); chk<BoF,Sub,Q,BoF>(); chk<LoF,Sub,Q,UpF>(); chk<UpF,Sub,Q,LoF>(); chk<ApF,Sub,Q,ApF>();
    chk<VaF,Sub,R,ExF>(); chk<MeF,Sub,R,MeF>(); chk<BoF,Sub,R,BoF>(); chk<LoF,Sub,R,UpF>(); chk<UpF,Sub,R,LoF>(); chk<ApF,Sub,R,ApF>();
}

void CheckNumeric::check_multiplication() {

    // Concrete numbers
    chk< D ,Mul,D,D>(); chk<D,Mul,D,I>(); chk<ApN,Mul,D,Z>(); chk<ApN,Mul,D,Q>(); chk<ApN,Mul,D,R>();
    chk< D ,Mul,I,D>(); chk<I,Mul,I,I>(); chk<Z,Mul,I,Z>(); chk<Q,Mul,I,Q>(); chk<R,Mul,I,R>();
    chk<ApN,Mul,Z,D>(); chk<Z,Mul,Z,I>(); chk<Z,Mul,Z,Z>(); chk<Q,Mul,Z,Q>(); chk<R,Mul,Z,R>();
    chk<ApN,Mul,Q,D>(); chk<Q,Mul,Q,I>(); chk<Q,Mul,Q,Z>(); chk<Q,Mul,Q,Q>(); chk<R,Mul,Q,R>();
    chk<ApN,Mul,R,D>(); chk<R,Mul,R,I>(); chk<R,Mul,R,Z>(); chk<R,Mul,R,Q>(); chk<R,Mul,R,R>();

    // Generic numbers
    chk<ExN,Mul,ExN,ExN>(); chk<EfN,Mul,ExN,EfN>(); chk<VaN,Mul,ExN,VaN>(); chk<UpN,Mul,ExN,UpN>(); chk<LoN,Mul,ExN,LoN>(); chk<ApN,Mul,ExN,ApN>();
    chk<EfN,Mul,EfN,ExN>(); chk<EfN,Mul,EfN,EfN>(); chk<VaN,Mul,EfN,VaN>(); chk<UpN,Mul,EfN,UpN>(); chk<LoN,Mul,EfN,LoN>(); chk<ApN,Mul,EfN,ApN>();
    chk<VaN,Mul,VaN,ExN>(); chk<VaN,Mul,VaN,EfN>(); chk<VaN,Mul,VaN,VaN>(); chk<UpN,Mul,VaN,UpN>(); chk<LoN,Mul,VaN,LoN>(); chk<ApN,Mul,VaN,ApN>();
    chk<UpN,Mul,UpN,ExN>(); chk<UpN,Mul,UpN,EfN>(); chk<UpN,Mul,UpN,VaN>(); chk<UpN,Mul,UpN,UpN>(); chk<ApN,Mul,UpN,LoN>(); chk<ApN,Mul,UpN,ApN>();
    chk<LoN,Mul,LoN,ExN>(); chk<LoN,Mul,LoN,EfN>(); chk<LoN,Mul,LoN,VaN>(); chk<ApN,Mul,LoN,UpN>(); chk<LoN,Mul,LoN,LoN>(); chk<ApN,Mul,LoN,ApN>();
    chk<ApN,Mul,ApN,ExN>(); chk<ApN,Mul,ApN,EfN>(); chk<ApN,Mul,ApN,VaN>(); chk<ApN,Mul,ApN,UpN>(); chk<ApN,Mul,ApN,LoN>(); chk<ApN,Mul,ApN,ApN>();

    // Mixed Generic # Concrete numbers
    chk<ApN,Mul,ExN,D>(); chk<ExN,Mul,ExN,I>(); chk<ExN,Mul,ExN,Z>(); chk<ExN,Mul,ExN,Q>(); chk<EfN,Mul,ExN,R>();
    chk<ApN,Mul,EfN,D>(); chk<EfN,Mul,EfN,I>(); chk<EfN,Mul,EfN,Z>(); chk<EfN,Mul,EfN,Q>(); chk<EfN,Mul,EfN,R>();
    chk<ApN,Mul,VaN,D>(); chk<VaN,Mul,VaN,I>(); chk<VaN,Mul,VaN,Z>(); chk<VaN,Mul,VaN,Q>(); chk<VaN,Mul,VaN,R>();
    chk<ApN,Mul,UpN,D>(); chk<UpN,Mul,UpN,I>(); chk<UpN,Mul,UpN,Z>(); chk<UpN,Mul,UpN,Q>(); chk<UpN,Mul,UpN,R>();
    chk<ApN,Mul,LoN,D>(); chk<LoN,Mul,LoN,I>(); chk<LoN,Mul,LoN,Z>(); chk<LoN,Mul,LoN,Q>(); chk<LoN,Mul,LoN,R>();
    chk<ApN,Mul,ApN,D>(); chk<ApN,Mul,ApN,I>(); chk<ApN,Mul,ApN,Z>(); chk<ApN,Mul,ApN,Q>(); chk<ApN,Mul,ApN,R>();

    // Mixed Concrete # Generic numbers
    chk<ApN,Mul,D,ExN>(); chk<ApN,Mul,D,EfN>(); chk<ApN,Mul,D,VaN>(); chk<ApN,Mul,D,UpN>(); chk<ApN,Mul,D,LoN>(); chk<ApN,Mul,D,ApN>();
    chk<ExN,Mul,I,ExN>(); chk<EfN,Mul,I,EfN>(); chk<VaN,Mul,I,VaN>(); chk<UpN,Mul,I,UpN>(); chk<LoN,Mul,I,LoN>(); chk<ApN,Mul,I,ApN>();
    chk<ExN,Mul,Z,ExN>(); chk<EfN,Mul,Z,EfN>(); chk<VaN,Mul,Z,VaN>(); chk<UpN,Mul,Z,UpN>(); chk<LoN,Mul,Z,LoN>(); chk<ApN,Mul,Z,ApN>();
    chk<ExN,Mul,Q,ExN>(); chk<EfN,Mul,Q,EfN>(); chk<VaN,Mul,Q,VaN>(); chk<UpN,Mul,Q,UpN>(); chk<LoN,Mul,Q,LoN>(); chk<ApN,Mul,Q,ApN>();
    chk<EfN,Mul,R,ExN>(); chk<EfN,Mul,R,EfN>(); chk<VaN,Mul,R,VaN>(); chk<UpN,Mul,R,UpN>(); chk<LoN,Mul,R,LoN>(); chk<ApN,Mul,R,ApN>();

    // Float numbers
    chk<VaF,Mul,ExF,ExF>(); chk<MeF,Mul,ExF,MeF>(); chk<BoF,Mul,ExF,BoF>(); chk<UpF,Mul,ExF,UpF>(); chk<LoF,Mul,ExF,LoF>(); chk<ApF,Mul,ExF,ApF>();
    chk<MeF,Mul,MeF,ExF>(); chk<MeF,Mul,MeF,MeF>(); chk<VaF,Mul,MeF,BoF>(); chk<UpF,Mul,MeF,UpF>(); chk<LoF,Mul,MeF,LoF>(); chk<ApF,Mul,MeF,ApF>();
    chk<BoF,Mul,BoF,ExF>(); chk<VaF,Mul,BoF,MeF>(); chk<BoF,Mul,BoF,BoF>(); chk<UpF,Mul,BoF,UpF>(); chk<LoF,Mul,BoF,LoF>(); chk<ApF,Mul,BoF,ApF>();
    chk<UpF,Mul,UpF,ExF>(); chk<UpF,Mul,UpF,MeF>(); chk<UpF,Mul,UpF,BoF>(); chk<UpF,Mul,UpF,UpF>(); chk<ApF,Mul,UpF,LoF>(); chk<ApF,Mul,UpF,ApF>();
    chk<LoF,Mul,LoF,ExF>(); chk<LoF,Mul,LoF,MeF>(); chk<LoF,Mul,LoF,BoF>(); chk<ApF,Mul,LoF,UpF>(); chk<LoF,Mul,LoF,LoF>(); chk<ApF,Mul,LoF,ApF>();
    chk<ApF,Mul,ApF,ExF>(); chk<ApF,Mul,ApF,MeF>(); chk<ApF,Mul,ApF,BoF>(); chk<ApF,Mul,ApF,UpF>(); chk<ApF,Mul,ApF,LoF>(); chk<ApF,Mul,ApF,ApF>();

    // Mixed Float # Generic and Generic # Float
    chk<VaF,Mul,ExF,ExN>(); chk<VaF,Mul,ExF,EfN>(); chk<VaF,Mul,ExF,VaN>(); chk<UpF,Mul,ExF,UpN>(); chk<LoF,Mul,ExF,LoN>(); chk<ApF,Mul,ExF,ApN>();
    chk<MeF,Mul,MeF,ExN>(); chk<MeF,Mul,MeF,EfN>(); chk<MeF,Mul,MeF,VaN>(); chk<UpF,Mul,MeF,UpN>(); chk<LoF,Mul,MeF,LoN>(); chk<ApF,Mul,MeF,ApN>();
    chk<BoF,Mul,BoF,ExN>(); chk<BoF,Mul,BoF,EfN>(); chk<BoF,Mul,BoF,VaN>(); chk<UpF,Mul,BoF,UpN>(); chk<LoF,Mul,BoF,LoN>(); chk<ApF,Mul,BoF,ApN>();
    chk<UpF,Mul,UpF,ExN>(); chk<UpF,Mul,UpF,EfN>(); chk<UpF,Mul,UpF,VaN>(); chk<UpF,Mul,UpF,UpN>(); chk<ApF,Mul,UpF,LoN>(); chk<ApF,Mul,UpF,ApN>();
    chk<LoF,Mul,LoF,ExN>(); chk<LoF,Mul,LoF,EfN>(); chk<LoF,Mul,LoF,VaN>(); chk<ApF,Mul,LoF,UpN>(); chk<LoF,Mul,LoF,LoN>(); chk<ApF,Mul,LoF,ApN>();
    chk<ApF,Mul,ApF,ExN>(); chk<ApF,Mul,ApF,EfN>(); chk<ApF,Mul,ApF,VaN>(); chk<ApF,Mul,ApF,UpN>(); chk<ApF,Mul,ApF,LoN>(); chk<ApF,Mul,ApF,ApN>();

    chk<VaF,Mul,ExN,ExF>(); chk<MeF,Mul,ExN,MeF>(); chk<BoF,Mul,ExN,BoF>(); chk<UpF,Mul,ExN,UpF>(); chk<LoF,Mul,ExN,LoF>(); chk<ApF,Mul,ExN,ApF>();
    chk<VaF,Mul,EfN,ExF>(); chk<MeF,Mul,EfN,MeF>(); chk<BoF,Mul,EfN,BoF>(); chk<UpF,Mul,EfN,UpF>(); chk<LoF,Mul,EfN,LoF>(); chk<ApF,Mul,EfN,ApF>();
    chk<VaF,Mul,VaN,ExF>(); chk<MeF,Mul,VaN,MeF>(); chk<BoF,Mul,VaN,BoF>(); chk<UpF,Mul,VaN,UpF>(); chk<LoF,Mul,VaN,LoF>(); chk<ApF,Mul,VaN,ApF>();
    chk<UpF,Mul,UpN,ExF>(); chk<UpF,Mul,UpN,MeF>(); chk<UpF,Mul,UpN,BoF>(); chk<UpF,Mul,UpN,UpF>(); chk<ApF,Mul,UpN,LoF>(); chk<ApF,Mul,UpN,ApF>();
    chk<LoF,Mul,LoN,ExF>(); chk<LoF,Mul,LoN,MeF>(); chk<LoF,Mul,LoN,BoF>(); chk<ApF,Mul,LoN,UpF>(); chk<LoF,Mul,LoN,LoF>(); chk<ApF,Mul,LoN,ApF>();
    chk<ApF,Mul,ApN,ExF>(); chk<ApF,Mul,ApN,MeF>(); chk<ApF,Mul,ApN,BoF>(); chk<ApF,Mul,ApN,UpF>(); chk<ApF,Mul,ApN,LoF>(); chk<ApF,Mul,ApN,ApF>();

    // Mixed Float # Concrete and Concrete # Float
    chk<ApF,Mul,ExF,D>(); chk<VaF,Mul,ExF,I>(); chk<VaF,Mul,ExF,Z>(); chk<VaF,Mul,ExF,Q>(); chk<VaF,Mul,ExF,R>();
    chk<ApF,Mul,MeF,D>(); chk<MeF,Mul,MeF,I>(); chk<MeF,Mul,MeF,Z>(); chk<MeF,Mul,MeF,Q>(); chk<MeF,Mul,MeF,R>();
    chk<ApF,Mul,BoF,D>(); chk<BoF,Mul,BoF,I>(); chk<BoF,Mul,BoF,Z>(); chk<BoF,Mul,BoF,Q>(); chk<BoF,Mul,BoF,R>();
    chk<ApF,Mul,UpF,D>(); chk<UpF,Mul,UpF,I>(); chk<UpF,Mul,UpF,Z>(); chk<UpF,Mul,UpF,Q>(); chk<UpF,Mul,UpF,R>();
    chk<ApF,Mul,LoF,D>(); chk<LoF,Mul,LoF,I>(); chk<LoF,Mul,LoF,Z>(); chk<LoF,Mul,LoF,Q>(); chk<LoF,Mul,LoF,R>();
    chk<ApF,Mul,ApF,D>(); chk<ApF,Mul,ApF,I>(); chk<ApF,Mul,ApF,Z>(); chk<ApF,Mul,ApF,Q>(); chk<ApF,Mul,ApF,R>();

    chk<ApF,Mul,D,ExF>(); chk<ApF,Mul,D,MeF>(); chk<ApF,Mul,D,BoF>(); chk<ApF,Mul,D,UpF>(); chk<ApF,Mul,D,LoF>(); chk<ApF,Mul,D,ApF>();
    chk<VaF,Mul,I,ExF>(); chk<MeF,Mul,I,MeF>(); chk<BoF,Mul,I,BoF>(); chk<UpF,Mul,I,UpF>(); chk<LoF,Mul,I,LoF>(); chk<ApF,Mul,I,ApF>();
    chk<VaF,Mul,Z,ExF>(); chk<MeF,Mul,Z,MeF>(); chk<BoF,Mul,Z,BoF>(); chk<UpF,Mul,Z,UpF>(); chk<LoF,Mul,Z,LoF>(); chk<ApF,Mul,Z,ApF>();
    chk<VaF,Mul,Q,ExF>(); chk<MeF,Mul,Q,MeF>(); chk<BoF,Mul,Q,BoF>(); chk<UpF,Mul,Q,UpF>(); chk<LoF,Mul,Q,LoF>(); chk<ApF,Mul,Q,ApF>();
    chk<VaF,Mul,R,ExF>(); chk<MeF,Mul,R,MeF>(); chk<BoF,Mul,R,BoF>(); chk<UpF,Mul,R,UpF>(); chk<LoF,Mul,R,LoF>(); chk<ApF,Mul,R,ApF>();

}

void CheckNumeric::check_division() {

    // Concrete numbers
    chk< D ,Div,D,D>(); chk<D,Div,D,I>(); chk<ApN,Div,D,Z>(); chk<ApN,Div,D,Q>(); chk<ApN,Div,D,R>();
    chk< D ,Div,I,D>(); chk<I,Div,I,I>(); chk<Q,Div,I,Z>(); chk<Q,Div,I,Q>(); chk<R,Div,I,R>();
    chk<ApN,Div,Z,D>(); chk<Q,Div,Z,I>(); chk<Q,Div,Z,Z>(); chk<Q,Div,Z,Q>(); chk<R,Div,Z,R>();
    chk<ApN,Div,Q,D>(); chk<Q,Div,Q,I>(); chk<Q,Div,Q,Z>(); chk<Q,Div,Q,Q>(); chk<R,Div,Q,R>();
    chk<ApN,Div,R,D>(); chk<R,Div,R,I>(); chk<R,Div,R,Z>(); chk<R,Div,R,Q>(); chk<R,Div,R,R>();

    // Generic numbers
    chk<ExN,Div,ExN,ExN>(); chk<EfN,Div,ExN,EfN>(); chk<VaN,Div,ExN,VaN>(); chk<LoN,Div,ExN,UpN>(); chk<UpN,Div,ExN,LoN>(); chk<ApN,Div,ExN,ApN>();
    chk<EfN,Div,EfN,ExN>(); chk<EfN,Div,EfN,EfN>(); chk<VaN,Div,EfN,VaN>(); chk<LoN,Div,EfN,UpN>(); chk<UpN,Div,EfN,LoN>(); chk<ApN,Div,EfN,ApN>();
    chk<VaN,Div,VaN,ExN>(); chk<VaN,Div,VaN,EfN>(); chk<VaN,Div,VaN,VaN>(); chk<LoN,Div,VaN,UpN>(); chk<UpN,Div,VaN,LoN>(); chk<ApN,Div,VaN,ApN>();
    chk<UpN,Div,UpN,ExN>(); chk<UpN,Div,UpN,EfN>(); chk<UpN,Div,UpN,VaN>(); chk<ApN,Div,UpN,UpN>(); chk<UpN,Div,UpN,LoN>(); chk<ApN,Div,UpN,ApN>();
    chk<LoN,Div,LoN,ExN>(); chk<LoN,Div,LoN,EfN>(); chk<LoN,Div,LoN,VaN>(); chk<LoN,Div,LoN,UpN>(); chk<ApN,Div,LoN,LoN>(); chk<ApN,Div,LoN,ApN>();
    chk<ApN,Div,ApN,ExN>(); chk<ApN,Div,ApN,EfN>(); chk<ApN,Div,ApN,VaN>(); chk<ApN,Div,ApN,UpN>(); chk<ApN,Div,ApN,LoN>(); chk<ApN,Div,ApN,ApN>();

    // Mixed Generic # Concrete numbers
    chk<ApN,Div,ExN,D>(); chk<ExN,Div,ExN,I>(); chk<ExN,Div,ExN,Z>(); chk<ExN,Div,ExN,Q>(); chk<EfN,Div,ExN,R>();
    chk<ApN,Div,EfN,D>(); chk<EfN,Div,EfN,I>(); chk<EfN,Div,EfN,Z>(); chk<EfN,Div,EfN,Q>(); chk<EfN,Div,EfN,R>();
    chk<ApN,Div,VaN,D>(); chk<VaN,Div,VaN,I>(); chk<VaN,Div,VaN,Z>(); chk<VaN,Div,VaN,Q>(); chk<VaN,Div,VaN,R>();
    chk<ApN,Div,UpN,D>(); chk<UpN,Div,UpN,I>(); chk<UpN,Div,UpN,Z>(); chk<UpN,Div,UpN,Q>(); chk<UpN,Div,UpN,R>();
    chk<ApN,Div,LoN,D>(); chk<LoN,Div,LoN,I>(); chk<LoN,Div,LoN,Z>(); chk<LoN,Div,LoN,Q>(); chk<LoN,Div,LoN,R>();
    chk<ApN,Div,ApN,D>(); chk<ApN,Div,ApN,I>(); chk<ApN,Div,ApN,Z>(); chk<ApN,Div,ApN,Q>(); chk<ApN,Div,ApN,R>();

    // Mixed Concrete # Generic numbers
    chk<ApN,Div,D,ExN>(); chk<ApN,Div,D,EfN>(); chk<ApN,Div,D,VaN>(); chk<ApN,Div,D,UpN>(); chk<ApN,Div,D,LoN>(); chk<ApN,Div,D,ApN>();
    chk<ExN,Div,I,ExN>(); chk<EfN,Div,I,EfN>(); chk<VaN,Div,I,VaN>(); chk<LoN,Div,I,UpN>(); chk<UpN,Div,I,LoN>(); chk<ApN,Div,I,ApN>();
    chk<ExN,Div,Z,ExN>(); chk<EfN,Div,Z,EfN>(); chk<VaN,Div,Z,VaN>(); chk<LoN,Div,Z,UpN>(); chk<UpN,Div,Z,LoN>(); chk<ApN,Div,Z,ApN>();
    chk<ExN,Div,Q,ExN>(); chk<EfN,Div,Q,EfN>(); chk<VaN,Div,Q,VaN>(); chk<LoN,Div,Q,UpN>(); chk<UpN,Div,Q,LoN>(); chk<ApN,Div,Q,ApN>();
    chk<EfN,Div,R,ExN>(); chk<EfN,Div,R,EfN>(); chk<VaN,Div,R,VaN>(); chk<LoN,Div,R,UpN>(); chk<UpN,Div,R,LoN>(); chk<ApN,Div,R,ApN>();

    // Float numbers
    chk<VaF,Div,ExF,ExF>(); chk<MeF,Div,ExF,MeF>(); chk<BoF,Div,ExF,BoF>(); chk<LoF,Div,ExF,UpF>(); chk<UpF,Div,ExF,LoF>(); chk<ApF,Div,ExF,ApF>();
    chk<MeF,Div,MeF,ExF>(); chk<MeF,Div,MeF,MeF>(); chk<VaF,Div,MeF,BoF>(); chk<LoF,Div,MeF,UpF>(); chk<UpF,Div,MeF,LoF>(); chk<ApF,Div,MeF,ApF>();
    chk<BoF,Div,BoF,ExF>(); chk<VaF,Div,BoF,MeF>(); chk<BoF,Div,BoF,BoF>(); chk<LoF,Div,BoF,UpF>(); chk<UpF,Div,BoF,LoF>(); chk<ApF,Div,BoF,ApF>();
    chk<UpF,Div,UpF,ExF>(); chk<UpF,Div,UpF,MeF>(); chk<UpF,Div,UpF,BoF>(); chk<ApF,Div,UpF,UpF>(); chk<UpF,Div,UpF,LoF>(); chk<ApF,Div,UpF,ApF>();
    chk<LoF,Div,LoF,ExF>(); chk<LoF,Div,LoF,MeF>(); chk<LoF,Div,LoF,BoF>(); chk<LoF,Div,LoF,UpF>(); chk<ApF,Div,LoF,LoF>(); chk<ApF,Div,LoF,ApF>();
    chk<ApF,Div,ApF,ExF>(); chk<ApF,Div,ApF,MeF>(); chk<ApF,Div,ApF,BoF>(); chk<ApF,Div,ApF,UpF>(); chk<ApF,Div,ApF,LoF>(); chk<ApF,Div,ApF,ApF>();

    // Mixed Float # Generic and Generic # Float
    chk<VaF,Div,ExF,ExN>(); chk<VaF,Div,ExF,EfN>(); chk<VaF,Div,ExF,VaN>(); chk<LoF,Div,ExF,UpN>(); chk<UpF,Div,ExF,LoN>(); chk<ApF,Div,ExF,ApN>();
    chk<MeF,Div,MeF,ExN>(); chk<MeF,Div,MeF,EfN>(); chk<MeF,Div,MeF,VaN>(); chk<LoF,Div,MeF,UpN>(); chk<UpF,Div,MeF,LoN>(); chk<ApF,Div,MeF,ApN>();
    chk<BoF,Div,BoF,ExN>(); chk<BoF,Div,BoF,EfN>(); chk<BoF,Div,BoF,VaN>(); chk<LoF,Div,BoF,UpN>(); chk<UpF,Div,BoF,LoN>(); chk<ApF,Div,BoF,ApN>();
    chk<UpF,Div,UpF,ExN>(); chk<UpF,Div,UpF,EfN>(); chk<UpF,Div,UpF,VaN>(); chk<ApF,Div,UpF,UpN>(); chk<UpF,Div,UpF,LoN>(); chk<ApF,Div,UpF,ApN>();
    chk<LoF,Div,LoF,ExN>(); chk<LoF,Div,LoF,EfN>(); chk<LoF,Div,LoF,VaN>(); chk<LoF,Div,LoF,UpN>(); chk<ApF,Div,LoF,LoN>(); chk<ApF,Div,LoF,ApN>();
    chk<ApF,Div,ApF,ExN>(); chk<ApF,Div,ApF,EfN>(); chk<ApF,Div,ApF,VaN>(); chk<ApF,Div,ApF,UpN>(); chk<ApF,Div,ApF,LoN>(); chk<ApF,Div,ApF,ApN>();

    chk<VaF,Div,ExN,ExF>(); chk<MeF,Div,ExN,MeF>(); chk<BoF,Div,ExN,BoF>(); chk<LoF,Div,ExN,UpF>(); chk<UpF,Div,ExN,LoF>(); chk<ApF,Div,ExN,ApF>();
    chk<VaF,Div,EfN,ExF>(); chk<MeF,Div,EfN,MeF>(); chk<BoF,Div,EfN,BoF>(); chk<LoF,Div,EfN,UpF>(); chk<UpF,Div,EfN,LoF>(); chk<ApF,Div,EfN,ApF>();
    chk<VaF,Div,VaN,ExF>(); chk<MeF,Div,VaN,MeF>(); chk<BoF,Div,VaN,BoF>(); chk<LoF,Div,VaN,UpF>(); chk<UpF,Div,VaN,LoF>(); chk<ApF,Div,VaN,ApF>();
    chk<UpF,Div,UpN,ExF>(); chk<UpF,Div,UpN,MeF>(); chk<UpF,Div,UpN,BoF>(); chk<ApF,Div,UpN,UpF>(); chk<UpF,Div,UpN,LoF>(); chk<ApF,Div,UpN,ApF>();
    chk<LoF,Div,LoN,ExF>(); chk<LoF,Div,LoN,MeF>(); chk<LoF,Div,LoN,BoF>(); chk<LoF,Div,LoN,UpF>(); chk<ApF,Div,LoN,LoF>(); chk<ApF,Div,LoN,ApF>();
    chk<ApF,Div,ApN,ExF>(); chk<ApF,Div,ApN,MeF>(); chk<ApF,Div,ApN,BoF>(); chk<ApF,Div,ApN,UpF>(); chk<ApF,Div,ApN,LoF>(); chk<ApF,Div,ApN,ApF>();

    // Mixed Float # Concrete and Concrete # Float
    chk<ApF,Div,ExF,D>(); chk<VaF,Div,ExF,I>(); chk<VaF,Div,ExF,Z>(); chk<VaF,Div,ExF,Q>(); chk<VaF,Div,ExF,R>();
    chk<ApF,Div,MeF,D>(); chk<MeF,Div,MeF,I>(); chk<MeF,Div,MeF,Z>(); chk<MeF,Div,MeF,Q>(); chk<MeF,Div,MeF,R>();
    chk<ApF,Div,BoF,D>(); chk<BoF,Div,BoF,I>(); chk<BoF,Div,BoF,Z>(); chk<BoF,Div,BoF,Q>(); chk<BoF,Div,BoF,R>();
    chk<ApF,Div,UpF,D>(); chk<UpF,Div,UpF,I>(); chk<UpF,Div,UpF,Z>(); chk<UpF,Div,UpF,Q>(); chk<UpF,Div,UpF,R>();
    chk<ApF,Div,LoF,D>(); chk<LoF,Div,LoF,I>(); chk<LoF,Div,LoF,Z>(); chk<LoF,Div,LoF,Q>(); chk<LoF,Div,LoF,R>();
    chk<ApF,Div,ApF,D>(); chk<ApF,Div,ApF,I>(); chk<ApF,Div,ApF,Z>(); chk<ApF,Div,ApF,Q>(); chk<ApF,Div,ApF,R>();

    chk<ApF,Div,D,ExF>(); chk<ApF,Div,D,MeF>(); chk<ApF,Div,D,BoF>(); chk<ApF,Div,D,UpF>(); chk<ApF,Div,D,LoF>(); chk<ApF,Div,D,ApF>();
    chk<VaF,Div,I,ExF>(); chk<MeF,Div,I,MeF>(); chk<BoF,Div,I,BoF>(); chk<LoF,Div,I,UpF>(); chk<UpF,Div,I,LoF>(); chk<ApF,Div,I,ApF>();
    chk<VaF,Div,Z,ExF>(); chk<MeF,Div,Z,MeF>(); chk<BoF,Div,Z,BoF>(); chk<LoF,Div,Z,UpF>(); chk<UpF,Div,Z,LoF>(); chk<ApF,Div,Z,ApF>();
    chk<VaF,Div,Q,ExF>(); chk<MeF,Div,Q,MeF>(); chk<BoF,Div,Q,BoF>(); chk<LoF,Div,Q,UpF>(); chk<UpF,Div,Q,LoF>(); chk<ApF,Div,Q,ApF>();
    chk<VaF,Div,R,ExF>(); chk<MeF,Div,R,MeF>(); chk<BoF,Div,R,BoF>(); chk<LoF,Div,R,UpF>(); chk<UpF,Div,R,LoF>(); chk<ApF,Div,R,ApF>();
}


void CheckNumeric::check_equality() {
    // Test mixed equality with exact/effective concrete numbers

    if (IsEquivalent<SafeEqualsType<Rational,double>,ExactLogic>::value) {
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Integer,double>,ExactLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Rational,double>,ExactLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Real,double>,LowerLogic)

        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,double>,LowerLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,double>,LowerLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,double>,ExactLogic)

        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactFloat,double>,ExactLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<MetricFloat,double>,ValidatedLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<BoundFloat,double>,ValidatedLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperFloat,double>,LowerLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerFloat,double>,LowerLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateFloat,double>,ApproximateLogic)
    } else {
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Integer,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Rational,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Real,double>,ApproximateLogic)

        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,double>,ApproximateLogic)

        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactFloat,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<MetricFloat,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<BoundFloat,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperFloat,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerFloat,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateFloat,double>,ApproximateLogic)
    }

    // Concrete user classes
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Integer,int>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Rational,int>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Real,int>,LowerLogic)

    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Integer,Integer>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Integer,Rational>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Integer,Real>,LowerLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Rational,Integer>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Rational,Rational>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Rational,Real>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Real,Integer>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Real,Rational>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<Real,Real>,LowerLogic)

    // Generic classes
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,int>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,int>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,int>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,int>,ExactLogic)


    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,Integer>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,Rational>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,Real>,ApproximateLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,Integer>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,Rational>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,Real>,LowerLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,Integer>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,Rational>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,Real>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,Integer>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,Rational>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,Real>,LowerLogic)

    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,ApproximateNumber>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,ValidatedNumber>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,EffectiveNumber>,ApproximateLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,ExactNumber>,ApproximateLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,ApproximateNumber>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,ValidatedNumber>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,EffectiveNumber>,LowerLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,ExactNumber>,LowerLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,ApproximateNumber>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,ValidatedNumber>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,EffectiveNumber>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,ExactNumber>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,ApproximateNumber>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,ValidatedNumber>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,EffectiveNumber>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,ExactNumber>,ExactLogic)


    // Numerical Flt classes

    // Test mixed equality with int
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactFloat,int>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<MetricFloat,int>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<BoundFloat,int>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperFloat,int>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerFloat,int>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateFloat,int>,ApproximateLogic)

    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactFloat,Integer>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<MetricFloat,Integer>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<BoundFloat,Integer>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperFloat,Integer>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerFloat,Integer>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateFloat,Integer>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactFloat,Rational>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<MetricFloat,Rational>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<BoundFloat,Rational>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperFloat,Rational>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerFloat,Rational>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateFloat,Rational>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactFloat,Real>,LowerLogic)     // TODO: Should these return LowerLogic or ValidatedLogic
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<MetricFloat,Real>,LowerLogic) // Rationale for is that Real really is "Effective", so should never test for equality
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<BoundFloat,Real>,LowerLogic)    // Against is that Real is sufficiently "Concrete"
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperFloat,Real>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerFloat,Real>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateFloat,Real>,ApproximateLogic)

    // User numeric types
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactFloat,ExactFloat>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactFloat,MetricFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactFloat,BoundFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactFloat,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactFloat,LowerFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactFloat,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<MetricFloat,ExactFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<MetricFloat,MetricFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<MetricFloat,BoundFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<MetricFloat,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<MetricFloat,LowerFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<MetricFloat,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<BoundFloat,ExactFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<BoundFloat,MetricFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<BoundFloat,BoundFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<BoundFloat,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<BoundFloat,LowerFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<BoundFloat,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperFloat,ExactFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperFloat,MetricFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperFloat,BoundFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperFloat,UpperFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperFloat,LowerFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperFloat,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerFloat,ExactFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerFloat,MetricFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerFloat,BoundFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerFloat,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerFloat,LowerFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerFloat,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateFloat,ExactFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateFloat,MetricFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateFloat,BoundFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateFloat,UpperFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateFloat,LowerFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateFloat,ApproximateFloat>,ApproximateLogic)

    // Require reduction of concrete number paradigm
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,ExactFloat>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,ExactFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,ExactFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperNumber,ExactFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerNumber,ExactFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,ExactFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,MetricFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,MetricFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,MetricFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperNumber,MetricFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerNumber,MetricFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,MetricFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,BoundFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,BoundFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,BoundFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperNumber,BoundFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerNumber,BoundFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,BoundFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperNumber,UpperFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerNumber,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,UpperFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,LowerFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,LowerFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,LowerFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperNumber,LowerFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerNumber,LowerFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,LowerFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ExactNumber,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<EffectiveNumber,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ValidatedNumber,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<UpperNumber,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<LowerNumber,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeEqualsType<ApproximateNumber,ApproximateFloat>,ApproximateLogic)
}

void CheckNumeric::check_comparison() {
    // Test mixed comparison with double

    if (IsEquivalent<SafeLessType<Rational,double>,ExactLogic>::value) {
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Integer,double>,ExactLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Rational,double>,ExactLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Real,double>,EffectiveLogic)

        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerNumber,double>,LowerLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperNumber,double>,UpperLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,double>,ValidatedLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,double>,EffectiveLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,double>,ExactLogic)

        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactFloat,double>,ExactLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<MetricFloat,double>,ValidatedLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<BoundFloat,double>,ValidatedLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperFloat,double>,UpperLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerFloat,double>,LowerLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateFloat,double>,ApproximateLogic)
    } else {
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Integer,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Rational,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Real,double>,ApproximateLogic)

        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,double>,ApproximateLogic)

        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactFloat,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<MetricFloat,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<BoundFloat,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperFloat,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerFloat,double>,ApproximateLogic)
        ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateFloat,double>,ApproximateLogic)
    }

    // Test mixed equality with int
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Integer,int>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Rational,int>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Real,int>,EffectiveLogic)

    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,int>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,int>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,int>,EffectiveLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,int>,ExactLogic)

    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactFloat,int>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<MetricFloat,int>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<BoundFloat,int>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperFloat,int>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerFloat,int>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateFloat,int>,ApproximateLogic)

    // Concrete user classes
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Integer,Integer>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Integer,Rational>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Integer,Real>,EffectiveLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Rational,Integer>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Rational,Rational>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Rational,Real>,EffectiveLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Real,Integer>,EffectiveLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Real,Rational>,EffectiveLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<Real,Real>,EffectiveLogic)

    // Generic classes

    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,Integer>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,Rational>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,Real>,ApproximateLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,Integer>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,Rational>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,Real>,ValidatedLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,Integer>,EffectiveLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,Rational>,EffectiveLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,Real>,EffectiveLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,Integer>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,Rational>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,Real>,EffectiveLogic)

    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,ApproximateNumber>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,ValidatedNumber>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,EffectiveNumber>,ApproximateLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,ExactNumber>,ApproximateLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,ApproximateNumber>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,ValidatedNumber>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,EffectiveNumber>,ValidatedLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,ExactNumber>,ValidatedLogic);
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,ApproximateNumber>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,ValidatedNumber>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,EffectiveNumber>,EffectiveLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,ExactNumber>,EffectiveLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,ApproximateNumber>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,ValidatedNumber>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,EffectiveNumber>,EffectiveLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,ExactNumber>,ExactLogic)


    // Numerical Flt classes

    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactFloat,Integer>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<MetricFloat,Integer>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<BoundFloat,Integer>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperFloat,Integer>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerFloat,Integer>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateFloat,Integer>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactFloat,Rational>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<MetricFloat,Rational>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<BoundFloat,Rational>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperFloat,Rational>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerFloat,Rational>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateFloat,Rational>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactFloat,Real>,EffectiveLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<MetricFloat,Real>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<BoundFloat,Real>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperFloat,Real>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerFloat,Real>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateFloat,Real>,ApproximateLogic)

    // User numeric types
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactFloat,ExactFloat>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactFloat,MetricFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactFloat,BoundFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactFloat,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactFloat,LowerFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactFloat,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<MetricFloat,ExactFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<MetricFloat,MetricFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<MetricFloat,BoundFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<MetricFloat,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<MetricFloat,LowerFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<MetricFloat,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<BoundFloat,ExactFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<BoundFloat,MetricFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<BoundFloat,BoundFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<BoundFloat,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<BoundFloat,LowerFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<BoundFloat,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperFloat,ExactFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperFloat,MetricFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperFloat,BoundFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperFloat,UpperFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperFloat,LowerFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperFloat,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerFloat,ExactFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerFloat,MetricFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerFloat,BoundFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerFloat,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerFloat,LowerFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerFloat,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateFloat,ExactFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateFloat,MetricFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateFloat,BoundFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateFloat,UpperFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateFloat,LowerFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateFloat,ApproximateFloat>,ApproximateLogic)

    // Require reduction of concrete number paradigm
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,ExactFloat>,ExactLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,ExactFloat>,EffectiveLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,ExactFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperNumber,ExactFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerNumber,ExactFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,ExactFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,MetricFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,MetricFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,MetricFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperNumber,MetricFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerNumber,MetricFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,MetricFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,BoundFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,BoundFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,BoundFloat>,ValidatedLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperNumber,BoundFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerNumber,BoundFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,BoundFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperNumber,UpperFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerNumber,UpperFloat>,LowerLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,UpperFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,LowerFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,LowerFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,LowerFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperNumber,LowerFloat>,UpperLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerNumber,LowerFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,LowerFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ExactNumber,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<EffectiveNumber,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ValidatedNumber,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<UpperNumber,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<LowerNumber,ApproximateFloat>,ApproximateLogic)
    ARIADNE_TEST_EQUIVALENT_TYPE(SafeLessType<ApproximateNumber,ApproximateFloat>,ApproximateLogic)
}


int main() {
    CheckNumeric().check();
    return ARIADNE_TEST_FAILURES;
}
