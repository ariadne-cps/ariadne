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

class Plus;

namespace Detail {

template<class OP, class X1, class X2, class = Fallback> struct SafeTypedef {
    typedef Fallback Type; };
template<class OP, class X1, class X2> struct SafeTypedef<OP,X1,X2, EnableIf<IsConvertible<decltype(declval<OP>()(declval<X1>(),declval<X2>())),DontCare>,Fallback>> {
    typedef decltype(declval<OP>()(declval<X1>(),declval<X2>()) ) Type; };

/*
template<class OP, class X1, class X2> struct SafeTypedef { typedef decltype( OP()(declval<X1>(),declval<X2>()) ) Type; };
*/

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
        ++ARIADNE_TEST_FAILURES;
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

void output_check_explicitly_constructible_result(Bool p, String t, String f, Bool c) {
    std::cout << "Checking " << t << " is constructible from " << f << ":";
    if(p) {
        std::cout << " OK\n";
    } else {
        std::cout << " Failed!\n";
        if(c) {
            std::cerr << "ERROR: " << t << " is convertible from  " << f << ", but only explict construction should be allowed.\n";
        } else {
            std::cerr << "ERROR: " << t << " is not constructible from " << f << ".\n";
        }
    }
}

template<class T, class F> void check_explicitly_constructible() {
    Bool c=IsConvertible<F,T>::value;
    Bool p=IsConstructible<T,F>::value and not c;
    output_check_explicitly_constructible_result(p,class_name<T>(),class_name<F>(),c);
}


void output_check_convertible_result(Bool p, String t, String f) {
    std::cout << "Checking " << t << " is convertible from " << f << ":";
    if(p) {
        std::cout << " OK\n";
    } else {
        std::cout << " Failed!\n";
        std::cerr << "ERROR: " << t << " is not convertible from " << f << ".\n";
    }
}

template<class T, class F> void check_convertible() {
    Bool p=IsConvertible<F,T>::value;
    output_check_convertible_result(p,class_name<T>(),class_name<F>());
}

void output_check_not_constructible_result(Bool p, String t, String f) {
    std::cout << "Checking " << t << " is not constructible from " << f << ":";
    if(p) {
        std::cout << " OK\n";
    } else {
        std::cout << " Failed!\n";
        std::cerr << "ERROR: " << t << " is constructible from " << f << "\n";
    }
}

template<class T, class F> void check_not_constructible() {
    Bool p=not IsConstructible<T,F>::value;
    output_check_not_constructible_result(p,class_name<T>(),class_name<F>());
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
ARIADNE_CLASS_NAME(Logical<ValidatedLower>);
ARIADNE_CLASS_NAME(Logical<ValidatedUpper>);
ARIADNE_CLASS_NAME(Logical<Validated>);
ARIADNE_CLASS_NAME(Logical<EffectiveLower>);
ARIADNE_CLASS_NAME(Logical<EffectiveUpper>);
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

ARIADNE_CLASS_NAME(Number<Approximate>)
ARIADNE_CLASS_NAME(Number<Lower>)
ARIADNE_CLASS_NAME(Number<Upper>)
ARIADNE_CLASS_NAME(Number<Validated>)
ARIADNE_CLASS_NAME(Number<Effective>)
ARIADNE_CLASS_NAME(Number<Exact>)

ARIADNE_CLASS_NAME(ApproximateFloat64);
ARIADNE_CLASS_NAME(LowerFloat64);
ARIADNE_CLASS_NAME(UpperFloat64);
ARIADNE_CLASS_NAME(ValidatedFloat64);
ARIADNE_CLASS_NAME(ErrorFloat64);
ARIADNE_CLASS_NAME(ExactFloat64);

ARIADNE_CLASS_NAME(ApproximateFloatMP);
ARIADNE_CLASS_NAME(LowerFloatMP);
ARIADNE_CLASS_NAME(UpperFloatMP);
ARIADNE_CLASS_NAME(ValidatedFloatMP);
ARIADNE_CLASS_NAME(ErrorFloatMP);
ARIADNE_CLASS_NAME(ExactFloatMP);

#undef ARIADNE_CLASS_NAME

typedef int I; typedef double D;
typedef Integer Z; typedef Rational Q; typedef Real R;
typedef Number<Exact> ExN; typedef Number<Effective> EfN; typedef Number<Validated> VaN;
typedef Number<Upper> UpN; typedef Number<Lower> LoN; typedef Number<Approximate> ApN;
typedef ExactFloat ExF; typedef MetricFloat MeF; typedef BoundedFloat BoF;
typedef UpperFloat UpF; typedef LowerFloat LoF; typedef ApproximateFloat ApF;

typedef decltype(declval<ExF>() + declval<ExF>()) VaF;
typedef decltype(declval<UpF>() * declval<UpF>()) PrF;
typedef decltype(declval<double>() + declval<Number<Approximate>>()) ApD;
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
    ARIADNE_TEST_CALL(check_subtraction());
    ARIADNE_TEST_CALL(check_multiplication());
    ARIADNE_TEST_CALL(check_division());
    ARIADNE_TEST_WARN("Skipping check of less and equal functions.");
    return;
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
    ARIADNE_TEST_NOTIFY((String("Integer + double -> ")+class_name<SafeSumType<Integer,double>>()));
    ARIADNE_TEST_NOTIFY((String("ExactNumber + double -> ")+class_name<SafeSumType<ExactNumber,double>>()));
    ARIADNE_TEST_NOTIFY((String("ApproximateNumber + double -> ")+class_name<SafeSumType<ApproximateNumber,double>>()));
    ARIADNE_TEST_NOTIFY((String("ExactFloat + double -> ")+class_name<SafeSumType<ExactFloat,double>>()));
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
    chk< D ,Plus,D,D>(); chk<D,Plus,D,I>(); chk<ApD,Plus,D,Z>(); chk<ApD,Plus,D,Q>(); chk<ApD,Plus,D,R>();
    chk< D ,Plus,I,D>(); chk<I,Plus,I,I>(); chk<Z,Plus,I,Z>(); chk<Q,Plus,I,Q>(); chk<R,Plus,I,R>();
    chk<ApD,Plus,Z,D>(); chk<Z,Plus,Z,I>(); chk<Z,Plus,Z,Z>(); chk<Q,Plus,Z,Q>(); chk<R,Plus,Z,R>();
    chk<ApD,Plus,Q,D>(); chk<Q,Plus,Q,I>(); chk<Q,Plus,Q,Z>(); chk<Q,Plus,Q,Q>(); chk<R,Plus,Q,R>();
    chk<ApD,Plus,R,D>(); chk<R,Plus,R,I>(); chk<R,Plus,R,Z>(); chk<R,Plus,R,Q>(); chk<R,Plus,R,R>();

    // Generic numbers
    chk<ExN,Plus,ExN,ExN>(); chk<EfN,Plus,ExN,EfN>(); chk<VaN,Plus,ExN,VaN>(); chk<UpN,Plus,ExN,UpN>(); chk<LoN,Plus,ExN,LoN>(); chk<ApN,Plus,ExN,ApN>();
    chk<EfN,Plus,EfN,ExN>(); chk<EfN,Plus,EfN,EfN>(); chk<VaN,Plus,EfN,VaN>(); chk<UpN,Plus,EfN,UpN>(); chk<LoN,Plus,EfN,LoN>(); chk<ApN,Plus,EfN,ApN>();
    chk<VaN,Plus,VaN,ExN>(); chk<VaN,Plus,VaN,EfN>(); chk<VaN,Plus,VaN,VaN>(); chk<UpN,Plus,VaN,UpN>(); chk<LoN,Plus,VaN,LoN>(); chk<ApN,Plus,VaN,ApN>();
    chk<UpN,Plus,UpN,ExN>(); chk<UpN,Plus,UpN,EfN>(); chk<UpN,Plus,UpN,VaN>(); chk<UpN,Plus,UpN,UpN>(); chk<ApN,Plus,UpN,LoN>(); chk<ApN,Plus,UpN,ApN>();
    chk<LoN,Plus,LoN,ExN>(); chk<LoN,Plus,LoN,EfN>(); chk<LoN,Plus,LoN,VaN>(); chk<ApN,Plus,LoN,UpN>(); chk<LoN,Plus,LoN,LoN>(); chk<ApN,Plus,LoN,ApN>();
    chk<ApN,Plus,ApN,ExN>(); chk<ApN,Plus,ApN,EfN>(); chk<ApN,Plus,ApN,VaN>(); chk<ApN,Plus,ApN,UpN>(); chk<ApN,Plus,ApN,LoN>(); chk<ApN,Plus,ApN,ApN>();

    // Mixed Generic # Concrete numbers
    chk<ApD,Plus,ExN,D>(); chk<ExN,Plus,ExN,I>(); chk<ExN,Plus,ExN,Z>(); chk<ExN,Plus,ExN,Q>(); chk<EfN,Plus,ExN,R>();
    chk<ApD,Plus,EfN,D>(); chk<EfN,Plus,EfN,I>(); chk<EfN,Plus,EfN,Z>(); chk<EfN,Plus,EfN,Q>(); chk<EfN,Plus,EfN,R>();
    chk<ApD,Plus,VaN,D>(); chk<VaN,Plus,VaN,I>(); chk<VaN,Plus,VaN,Z>(); chk<VaN,Plus,VaN,Q>(); chk<VaN,Plus,VaN,R>();
    chk<ApD,Plus,UpN,D>(); chk<UpN,Plus,UpN,I>(); chk<UpN,Plus,UpN,Z>(); chk<UpN,Plus,UpN,Q>(); chk<UpN,Plus,UpN,R>();
    chk<ApD,Plus,LoN,D>(); chk<LoN,Plus,LoN,I>(); chk<LoN,Plus,LoN,Z>(); chk<LoN,Plus,LoN,Q>(); chk<LoN,Plus,LoN,R>();
    chk<ApD,Plus,ApN,D>(); chk<ApN,Plus,ApN,I>(); chk<ApN,Plus,ApN,Z>(); chk<ApN,Plus,ApN,Q>(); chk<ApN,Plus,ApN,R>();

    // Mixed Concrete # Generic numbers
    chk<ApD,Plus,D,ExN>(); chk<ApD,Plus,D,EfN>(); chk<ApD,Plus,D,VaN>(); chk<ApD,Plus,D,UpN>(); chk<ApD,Plus,D,LoN>(); chk<ApD,Plus,D,ApN>();
    chk<ExN,Plus,I,ExN>(); chk<EfN,Plus,I,EfN>(); chk<VaN,Plus,I,VaN>(); chk<UpN,Plus,I,UpN>(); chk<LoN,Plus,I,LoN>(); chk<ApN,Plus,I,ApN>();
    chk<ExN,Plus,Z,ExN>(); chk<EfN,Plus,Z,EfN>(); chk<VaN,Plus,Z,VaN>(); chk<UpN,Plus,Z,UpN>(); chk<LoN,Plus,Z,LoN>(); chk<ApN,Plus,Z,ApN>();
    chk<ExN,Plus,Q,ExN>(); chk<EfN,Plus,Q,EfN>(); chk<VaN,Plus,Q,VaN>(); chk<UpN,Plus,Q,UpN>(); chk<LoN,Plus,Q,LoN>(); chk<ApN,Plus,Q,ApN>();
    chk<EfN,Plus,R,ExN>(); chk<EfN,Plus,R,EfN>(); chk<VaN,Plus,R,VaN>(); chk<UpN,Plus,R,UpN>(); chk<LoN,Plus,R,LoN>(); chk<ApN,Plus,R,ApN>();

    // Float numbers
    chk<VaF,Plus,ExF,ExF>(); chk<MeF,Plus,ExF,MeF>(); chk<BoF,Plus,ExF,BoF>(); chk<UpF,Plus,ExF,UpF>(); chk<LoF,Plus,ExF,LoF>(); chk<ApF,Plus,ExF,ApF>();
    chk<MeF,Plus,MeF,ExF>(); chk<MeF,Plus,MeF,MeF>(); chk<VaF,Plus,MeF,BoF>(); chk<UpF,Plus,MeF,UpF>(); chk<LoF,Plus,MeF,LoF>(); chk<ApF,Plus,MeF,ApF>();
    chk<BoF,Plus,BoF,ExF>(); chk<VaF,Plus,BoF,MeF>(); chk<BoF,Plus,BoF,BoF>(); chk<UpF,Plus,BoF,UpF>(); chk<LoF,Plus,BoF,LoF>(); chk<ApF,Plus,BoF,ApF>();
    chk<UpF,Plus,UpF,ExF>(); chk<UpF,Plus,UpF,MeF>(); chk<UpF,Plus,UpF,BoF>(); chk<UpF,Plus,UpF,UpF>(); chk<ApF,Plus,UpF,LoF>(); chk<ApF,Plus,UpF,ApF>();
    chk<LoF,Plus,LoF,ExF>(); chk<LoF,Plus,LoF,MeF>(); chk<LoF,Plus,LoF,BoF>(); chk<ApF,Plus,LoF,UpF>(); chk<LoF,Plus,LoF,LoF>(); chk<ApF,Plus,LoF,ApF>();
    chk<ApF,Plus,ApF,ExF>(); chk<ApF,Plus,ApF,MeF>(); chk<ApF,Plus,ApF,BoF>(); chk<ApF,Plus,ApF,UpF>(); chk<ApF,Plus,ApF,LoF>(); chk<ApF,Plus,ApF,ApF>();

    // Mixed Float # Generic and Generic # Float
    chk<VaF,Plus,ExF,ExN>(); chk<VaF,Plus,ExF,EfN>(); chk<VaF,Plus,ExF,VaN>(); chk<UpF,Plus,ExF,UpN>(); chk<LoF,Plus,ExF,LoN>(); chk<ApF,Plus,ExF,ApN>();
    chk<MeF,Plus,MeF,ExN>(); chk<MeF,Plus,MeF,EfN>(); chk<MeF,Plus,MeF,VaN>(); chk<UpF,Plus,MeF,UpN>(); chk<LoF,Plus,MeF,LoN>(); chk<ApF,Plus,MeF,ApN>();
    chk<BoF,Plus,BoF,ExN>(); chk<BoF,Plus,BoF,EfN>(); chk<BoF,Plus,BoF,VaN>(); chk<UpF,Plus,BoF,UpN>(); chk<LoF,Plus,BoF,LoN>(); chk<ApF,Plus,BoF,ApN>();
    chk<UpF,Plus,UpF,ExN>(); chk<UpF,Plus,UpF,EfN>(); chk<UpF,Plus,UpF,VaN>(); chk<UpF,Plus,UpF,UpN>(); chk<ApF,Plus,UpF,LoN>(); chk<ApF,Plus,UpF,ApN>();
    chk<LoF,Plus,LoF,ExN>(); chk<LoF,Plus,LoF,EfN>(); chk<LoF,Plus,LoF,VaN>(); chk<ApF,Plus,LoF,UpN>(); chk<LoF,Plus,LoF,LoN>(); chk<ApF,Plus,LoF,ApN>();
    chk<ApF,Plus,ApF,ExN>(); chk<ApF,Plus,ApF,EfN>(); chk<ApF,Plus,ApF,VaN>(); chk<ApF,Plus,ApF,UpN>(); chk<ApF,Plus,ApF,LoN>(); chk<ApF,Plus,ApF,ApN>();

    chk<VaF,Plus,ExN,ExF>(); chk<MeF,Plus,ExN,MeF>(); chk<BoF,Plus,ExN,BoF>(); chk<UpF,Plus,ExN,UpF>(); chk<LoF,Plus,ExN,LoF>(); chk<ApF,Plus,ExN,ApF>();
    chk<VaF,Plus,EfN,ExF>(); chk<MeF,Plus,EfN,MeF>(); chk<BoF,Plus,EfN,BoF>(); chk<UpF,Plus,EfN,UpF>(); chk<LoF,Plus,EfN,LoF>(); chk<ApF,Plus,EfN,ApF>();
    chk<VaF,Plus,VaN,ExF>(); chk<MeF,Plus,VaN,MeF>(); chk<BoF,Plus,VaN,BoF>(); chk<UpF,Plus,VaN,UpF>(); chk<LoF,Plus,VaN,LoF>(); chk<ApF,Plus,VaN,ApF>();
    chk<UpF,Plus,UpN,ExF>(); chk<UpF,Plus,UpN,MeF>(); chk<UpF,Plus,UpN,BoF>(); chk<UpF,Plus,UpN,UpF>(); chk<ApF,Plus,UpN,LoF>(); chk<ApF,Plus,UpN,ApF>();
    chk<LoF,Plus,LoN,ExF>(); chk<LoF,Plus,LoN,MeF>(); chk<LoF,Plus,LoN,BoF>(); chk<ApF,Plus,LoN,UpF>(); chk<LoF,Plus,LoN,LoF>(); chk<ApF,Plus,LoN,ApF>();
    chk<ApF,Plus,ApN,ExF>(); chk<ApF,Plus,ApN,MeF>(); chk<ApF,Plus,ApN,BoF>(); chk<ApF,Plus,ApN,UpF>(); chk<ApF,Plus,ApN,LoF>(); chk<ApF,Plus,ApN,ApF>();

    // Mixed Float # Concrete and Concrete # Float
    chk<ApF,Plus,ExF,D>(); chk<VaF,Plus,ExF,I>(); chk<VaF,Plus,ExF,Z>(); chk<VaF,Plus,ExF,Q>(); chk<VaF,Plus,ExF,R>();
    chk<ApF,Plus,MeF,D>(); chk<MeF,Plus,MeF,I>(); chk<MeF,Plus,MeF,Z>(); chk<MeF,Plus,MeF,Q>(); chk<MeF,Plus,MeF,R>();
    chk<ApF,Plus,BoF,D>(); chk<BoF,Plus,BoF,I>(); chk<BoF,Plus,BoF,Z>(); chk<BoF,Plus,BoF,Q>(); chk<BoF,Plus,BoF,R>();
    chk<ApF,Plus,UpF,D>(); chk<UpF,Plus,UpF,I>(); chk<UpF,Plus,UpF,Z>(); chk<UpF,Plus,UpF,Q>(); chk<UpF,Plus,UpF,R>();
    chk<ApF,Plus,LoF,D>(); chk<LoF,Plus,LoF,I>(); chk<LoF,Plus,LoF,Z>(); chk<LoF,Plus,LoF,Q>(); chk<LoF,Plus,LoF,R>();
    chk<ApF,Plus,ApF,D>(); chk<ApF,Plus,ApF,I>(); chk<ApF,Plus,ApF,Z>(); chk<ApF,Plus,ApF,Q>(); chk<ApF,Plus,ApF,R>();

    chk<ApF,Plus,D,ExF>(); chk<ApF,Plus,D,MeF>(); chk<ApF,Plus,D,BoF>(); chk<ApF,Plus,D,UpF>(); chk<ApF,Plus,D,LoF>(); chk<ApF,Plus,D,ApF>();
    chk<VaF,Plus,I,ExF>(); chk<MeF,Plus,I,MeF>(); chk<BoF,Plus,I,BoF>(); chk<UpF,Plus,I,UpF>(); chk<LoF,Plus,I,LoF>(); chk<ApF,Plus,I,ApF>();
    chk<VaF,Plus,Z,ExF>(); chk<MeF,Plus,Z,MeF>(); chk<BoF,Plus,Z,BoF>(); chk<UpF,Plus,Z,UpF>(); chk<LoF,Plus,Z,LoF>(); chk<ApF,Plus,Z,ApF>();
    chk<VaF,Plus,Q,ExF>(); chk<MeF,Plus,Q,MeF>(); chk<BoF,Plus,Q,BoF>(); chk<UpF,Plus,Q,UpF>(); chk<LoF,Plus,Q,LoF>(); chk<ApF,Plus,Q,ApF>();
    chk<VaF,Plus,R,ExF>(); chk<MeF,Plus,R,MeF>(); chk<BoF,Plus,R,BoF>(); chk<UpF,Plus,R,UpF>(); chk<LoF,Plus,R,LoF>(); chk<ApF,Plus,R,ApF>();
}

void CheckNumeric::check_subtraction() {

    // Concrete numbers
    chk< D ,Minus,D,D>(); chk<D,Minus,D,I>(); chk<ApD,Minus,D,Z>(); chk<ApD,Minus,D,Q>(); chk<ApD,Minus,D,R>();
    chk< D ,Minus,I,D>(); chk<I,Minus,I,I>(); chk<Z,Minus,I,Z>(); chk<Q,Minus,I,Q>(); chk<R,Minus,I,R>();
    chk<ApD,Minus,Z,D>(); chk<Z,Minus,Z,I>(); chk<Z,Minus,Z,Z>(); chk<Q,Minus,Z,Q>(); chk<R,Minus,Z,R>();
    chk<ApD,Minus,Q,D>(); chk<Q,Minus,Q,I>(); chk<Q,Minus,Q,Z>(); chk<Q,Minus,Q,Q>(); chk<R,Minus,Q,R>();
    chk<ApD,Minus,R,D>(); chk<R,Minus,R,I>(); chk<R,Minus,R,Z>(); chk<R,Minus,R,Q>(); chk<R,Minus,R,R>();

    // Generic numbers
    chk<ExN,Minus,ExN,ExN>(); chk<EfN,Minus,ExN,EfN>(); chk<VaN,Minus,ExN,VaN>(); chk<LoN,Minus,ExN,UpN>(); chk<UpN,Minus,ExN,LoN>(); chk<ApN,Minus,ExN,ApN>();
    chk<EfN,Minus,EfN,ExN>(); chk<EfN,Minus,EfN,EfN>(); chk<VaN,Minus,EfN,VaN>(); chk<LoN,Minus,EfN,UpN>(); chk<UpN,Minus,EfN,LoN>(); chk<ApN,Minus,EfN,ApN>();
    chk<VaN,Minus,VaN,ExN>(); chk<VaN,Minus,VaN,EfN>(); chk<VaN,Minus,VaN,VaN>(); chk<LoN,Minus,VaN,UpN>(); chk<UpN,Minus,VaN,LoN>(); chk<ApN,Minus,VaN,ApN>();
    chk<UpN,Minus,UpN,ExN>(); chk<UpN,Minus,UpN,EfN>(); chk<UpN,Minus,UpN,VaN>(); chk<ApN,Minus,UpN,UpN>(); chk<UpN,Minus,UpN,LoN>(); chk<ApN,Minus,UpN,ApN>();
    chk<LoN,Minus,LoN,ExN>(); chk<LoN,Minus,LoN,EfN>(); chk<LoN,Minus,LoN,VaN>(); chk<LoN,Minus,LoN,UpN>(); chk<ApN,Minus,LoN,LoN>(); chk<ApN,Minus,LoN,ApN>();
    chk<ApN,Minus,ApN,ExN>(); chk<ApN,Minus,ApN,EfN>(); chk<ApN,Minus,ApN,VaN>(); chk<ApN,Minus,ApN,UpN>(); chk<ApN,Minus,ApN,LoN>(); chk<ApN,Minus,ApN,ApN>();

    // Mixed Generic # Concrete numbers
    chk<ApD,Minus,ExN,D>(); chk<ExN,Minus,ExN,I>(); chk<ExN,Minus,ExN,Z>(); chk<ExN,Minus,ExN,Q>(); chk<EfN,Minus,ExN,R>();
    chk<ApD,Minus,EfN,D>(); chk<EfN,Minus,EfN,I>(); chk<EfN,Minus,EfN,Z>(); chk<EfN,Minus,EfN,Q>(); chk<EfN,Minus,EfN,R>();
    chk<ApD,Minus,VaN,D>(); chk<VaN,Minus,VaN,I>(); chk<VaN,Minus,VaN,Z>(); chk<VaN,Minus,VaN,Q>(); chk<VaN,Minus,VaN,R>();
    chk<ApD,Minus,UpN,D>(); chk<UpN,Minus,UpN,I>(); chk<UpN,Minus,UpN,Z>(); chk<UpN,Minus,UpN,Q>(); chk<UpN,Minus,UpN,R>();
    chk<ApD,Minus,LoN,D>(); chk<LoN,Minus,LoN,I>(); chk<LoN,Minus,LoN,Z>(); chk<LoN,Minus,LoN,Q>(); chk<LoN,Minus,LoN,R>();
    chk<ApD,Minus,ApN,D>(); chk<ApN,Minus,ApN,I>(); chk<ApN,Minus,ApN,Z>(); chk<ApN,Minus,ApN,Q>(); chk<ApN,Minus,ApN,R>();

    // Mixed Concrete # Generic numbers
    chk<ApD,Minus,D,ExN>(); chk<ApD,Minus,D,EfN>(); chk<ApD,Minus,D,VaN>(); chk<ApD,Minus,D,UpN>(); chk<ApD,Minus,D,LoN>(); chk<ApD,Minus,D,ApN>();
    chk<ExN,Minus,I,ExN>(); chk<EfN,Minus,I,EfN>(); chk<VaN,Minus,I,VaN>(); chk<LoN,Minus,I,UpN>(); chk<UpN,Minus,I,LoN>(); chk<ApN,Minus,I,ApN>();
    chk<ExN,Minus,Z,ExN>(); chk<EfN,Minus,Z,EfN>(); chk<VaN,Minus,Z,VaN>(); chk<LoN,Minus,Z,UpN>(); chk<UpN,Minus,Z,LoN>(); chk<ApN,Minus,Z,ApN>();
    chk<ExN,Minus,Q,ExN>(); chk<EfN,Minus,Q,EfN>(); chk<VaN,Minus,Q,VaN>(); chk<LoN,Minus,Q,UpN>(); chk<UpN,Minus,Q,LoN>(); chk<ApN,Minus,Q,ApN>();
    chk<EfN,Minus,R,ExN>(); chk<EfN,Minus,R,EfN>(); chk<VaN,Minus,R,VaN>(); chk<LoN,Minus,R,UpN>(); chk<UpN,Minus,R,LoN>(); chk<ApN,Minus,R,ApN>();

    // Float numbers
    chk<VaF,Minus,ExF,ExF>(); chk<MeF,Minus,ExF,MeF>(); chk<BoF,Minus,ExF,BoF>(); chk<LoF,Minus,ExF,UpF>(); chk<UpF,Minus,ExF,LoF>(); chk<ApF,Minus,ExF,ApF>();
    chk<MeF,Minus,MeF,ExF>(); chk<MeF,Minus,MeF,MeF>(); chk<VaF,Minus,MeF,BoF>(); chk<LoF,Minus,MeF,UpF>(); chk<UpF,Minus,MeF,LoF>(); chk<ApF,Minus,MeF,ApF>();
    chk<BoF,Minus,BoF,ExF>(); chk<VaF,Minus,BoF,MeF>(); chk<BoF,Minus,BoF,BoF>(); chk<LoF,Minus,BoF,UpF>(); chk<UpF,Minus,BoF,LoF>(); chk<ApF,Minus,BoF,ApF>();
    chk<UpF,Minus,UpF,ExF>(); chk<UpF,Minus,UpF,MeF>(); chk<UpF,Minus,UpF,BoF>(); chk<ApF,Minus,UpF,UpF>(); chk<UpF,Minus,UpF,LoF>(); chk<ApF,Minus,UpF,ApF>();
    chk<LoF,Minus,LoF,ExF>(); chk<LoF,Minus,LoF,MeF>(); chk<LoF,Minus,LoF,BoF>(); chk<LoF,Minus,LoF,UpF>(); chk<ApF,Minus,LoF,LoF>(); chk<ApF,Minus,LoF,ApF>();
    chk<ApF,Minus,ApF,ExF>(); chk<ApF,Minus,ApF,MeF>(); chk<ApF,Minus,ApF,BoF>(); chk<ApF,Minus,ApF,UpF>(); chk<ApF,Minus,ApF,LoF>(); chk<ApF,Minus,ApF,ApF>();

    // Mixed Float # Generic and Generic # Float
    chk<VaF,Minus,ExF,ExN>(); chk<VaF,Minus,ExF,EfN>(); chk<VaF,Minus,ExF,VaN>(); chk<LoF,Minus,ExF,UpN>(); chk<UpF,Minus,ExF,LoN>(); chk<ApF,Minus,ExF,ApN>();
    chk<MeF,Minus,MeF,ExN>(); chk<MeF,Minus,MeF,EfN>(); chk<MeF,Minus,MeF,VaN>(); chk<LoF,Minus,MeF,UpN>(); chk<UpF,Minus,MeF,LoN>(); chk<ApF,Minus,MeF,ApN>();
    chk<BoF,Minus,BoF,ExN>(); chk<BoF,Minus,BoF,EfN>(); chk<BoF,Minus,BoF,VaN>(); chk<LoF,Minus,BoF,UpN>(); chk<UpF,Minus,BoF,LoN>(); chk<ApF,Minus,BoF,ApN>();
    chk<UpF,Minus,UpF,ExN>(); chk<UpF,Minus,UpF,EfN>(); chk<UpF,Minus,UpF,VaN>(); chk<ApF,Minus,UpF,UpN>(); chk<UpF,Minus,UpF,LoN>(); chk<ApF,Minus,UpF,ApN>();
    chk<LoF,Minus,LoF,ExN>(); chk<LoF,Minus,LoF,EfN>(); chk<LoF,Minus,LoF,VaN>(); chk<LoF,Minus,LoF,UpN>(); chk<ApF,Minus,LoF,LoN>(); chk<ApF,Minus,LoF,ApN>();
    chk<ApF,Minus,ApF,ExN>(); chk<ApF,Minus,ApF,EfN>(); chk<ApF,Minus,ApF,VaN>(); chk<ApF,Minus,ApF,UpN>(); chk<ApF,Minus,ApF,LoN>(); chk<ApF,Minus,ApF,ApN>();

    chk<VaF,Minus,ExN,ExF>(); chk<MeF,Minus,ExN,MeF>(); chk<BoF,Minus,ExN,BoF>(); chk<LoF,Minus,ExN,UpF>(); chk<UpF,Minus,ExN,LoF>(); chk<ApF,Minus,ExN,ApF>();
    chk<VaF,Minus,EfN,ExF>(); chk<MeF,Minus,EfN,MeF>(); chk<BoF,Minus,EfN,BoF>(); chk<LoF,Minus,EfN,UpF>(); chk<UpF,Minus,EfN,LoF>(); chk<ApF,Minus,EfN,ApF>();
    chk<VaF,Minus,VaN,ExF>(); chk<MeF,Minus,VaN,MeF>(); chk<BoF,Minus,VaN,BoF>(); chk<LoF,Minus,VaN,UpF>(); chk<UpF,Minus,VaN,LoF>(); chk<ApF,Minus,VaN,ApF>();
    chk<UpF,Minus,UpN,ExF>(); chk<UpF,Minus,UpN,MeF>(); chk<UpF,Minus,UpN,BoF>(); chk<ApF,Minus,UpN,UpF>(); chk<UpF,Minus,UpN,LoF>(); chk<ApF,Minus,UpN,ApF>();
    chk<LoF,Minus,LoN,ExF>(); chk<LoF,Minus,LoN,MeF>(); chk<LoF,Minus,LoN,BoF>(); chk<LoF,Minus,LoN,UpF>(); chk<ApF,Minus,LoN,LoF>(); chk<ApF,Minus,LoN,ApF>();
    chk<ApF,Minus,ApN,ExF>(); chk<ApF,Minus,ApN,MeF>(); chk<ApF,Minus,ApN,BoF>(); chk<ApF,Minus,ApN,UpF>(); chk<ApF,Minus,ApN,LoF>(); chk<ApF,Minus,ApN,ApF>();

    // Mixed Float # Concrete and Concrete # Float
    chk<ApF,Minus,ExF,D>(); chk<VaF,Minus,ExF,I>(); chk<VaF,Minus,ExF,Z>(); chk<VaF,Minus,ExF,Q>(); chk<VaF,Minus,ExF,R>();
    chk<ApF,Minus,MeF,D>(); chk<MeF,Minus,MeF,I>(); chk<MeF,Minus,MeF,Z>(); chk<MeF,Minus,MeF,Q>(); chk<MeF,Minus,MeF,R>();
    chk<ApF,Minus,BoF,D>(); chk<BoF,Minus,BoF,I>(); chk<BoF,Minus,BoF,Z>(); chk<BoF,Minus,BoF,Q>(); chk<BoF,Minus,BoF,R>();
    chk<ApF,Minus,UpF,D>(); chk<UpF,Minus,UpF,I>(); chk<UpF,Minus,UpF,Z>(); chk<UpF,Minus,UpF,Q>(); chk<UpF,Minus,UpF,R>();
    chk<ApF,Minus,LoF,D>(); chk<LoF,Minus,LoF,I>(); chk<LoF,Minus,LoF,Z>(); chk<LoF,Minus,LoF,Q>(); chk<LoF,Minus,LoF,R>();
    chk<ApF,Minus,ApF,D>(); chk<ApF,Minus,ApF,I>(); chk<ApF,Minus,ApF,Z>(); chk<ApF,Minus,ApF,Q>(); chk<ApF,Minus,ApF,R>();

    chk<ApF,Minus,D,ExF>(); chk<ApF,Minus,D,MeF>(); chk<ApF,Minus,D,BoF>(); chk<ApF,Minus,D,UpF>(); chk<ApF,Minus,D,LoF>(); chk<ApF,Minus,D,ApF>();
    chk<VaF,Minus,I,ExF>(); chk<MeF,Minus,I,MeF>(); chk<BoF,Minus,I,BoF>(); chk<LoF,Minus,I,UpF>(); chk<UpF,Minus,I,LoF>(); chk<ApF,Minus,I,ApF>();
    chk<VaF,Minus,Z,ExF>(); chk<MeF,Minus,Z,MeF>(); chk<BoF,Minus,Z,BoF>(); chk<LoF,Minus,Z,UpF>(); chk<UpF,Minus,Z,LoF>(); chk<ApF,Minus,Z,ApF>();
    chk<VaF,Minus,Q,ExF>(); chk<MeF,Minus,Q,MeF>(); chk<BoF,Minus,Q,BoF>(); chk<LoF,Minus,Q,UpF>(); chk<UpF,Minus,Q,LoF>(); chk<ApF,Minus,Q,ApF>();
    chk<VaF,Minus,R,ExF>(); chk<MeF,Minus,R,MeF>(); chk<BoF,Minus,R,BoF>(); chk<LoF,Minus,R,UpF>(); chk<UpF,Minus,R,LoF>(); chk<ApF,Minus,R,ApF>();
}

void CheckNumeric::check_multiplication() {

    // Concrete numbers
    chk< D ,Times,D,D>(); chk<D,Times,D,I>(); chk<ApD,Times,D,Z>(); chk<ApD,Times,D,Q>(); chk<ApD,Times,D,R>();
    chk< D ,Times,I,D>(); chk<I,Times,I,I>(); chk<Z,Times,I,Z>(); chk<Q,Times,I,Q>(); chk<R,Times,I,R>();
    chk<ApD,Times,Z,D>(); chk<Z,Times,Z,I>(); chk<Z,Times,Z,Z>(); chk<Q,Times,Z,Q>(); chk<R,Times,Z,R>();
    chk<ApD,Times,Q,D>(); chk<Q,Times,Q,I>(); chk<Q,Times,Q,Z>(); chk<Q,Times,Q,Q>(); chk<R,Times,Q,R>();
    chk<ApD,Times,R,D>(); chk<R,Times,R,I>(); chk<R,Times,R,Z>(); chk<R,Times,R,Q>(); chk<R,Times,R,R>();

    // Generic numbers
    chk<ExN,Times,ExN,ExN>(); chk<EfN,Times,ExN,EfN>(); chk<VaN,Times,ExN,VaN>(); chk<UpN,Times,ExN,UpN>(); chk<LoN,Times,ExN,LoN>(); chk<ApN,Times,ExN,ApN>();
    chk<EfN,Times,EfN,ExN>(); chk<EfN,Times,EfN,EfN>(); chk<VaN,Times,EfN,VaN>(); chk<UpN,Times,EfN,UpN>(); chk<LoN,Times,EfN,LoN>(); chk<ApN,Times,EfN,ApN>();
    chk<VaN,Times,VaN,ExN>(); chk<VaN,Times,VaN,EfN>(); chk<VaN,Times,VaN,VaN>(); chk<UpN,Times,VaN,UpN>(); chk<LoN,Times,VaN,LoN>(); chk<ApN,Times,VaN,ApN>();
    chk<UpN,Times,UpN,ExN>(); chk<UpN,Times,UpN,EfN>(); chk<UpN,Times,UpN,VaN>(); chk<UpN,Times,UpN,UpN>(); chk<ApN,Times,UpN,LoN>(); chk<ApN,Times,UpN,ApN>();
    chk<LoN,Times,LoN,ExN>(); chk<LoN,Times,LoN,EfN>(); chk<LoN,Times,LoN,VaN>(); chk<ApN,Times,LoN,UpN>(); chk<LoN,Times,LoN,LoN>(); chk<ApN,Times,LoN,ApN>();
    chk<ApN,Times,ApN,ExN>(); chk<ApN,Times,ApN,EfN>(); chk<ApN,Times,ApN,VaN>(); chk<ApN,Times,ApN,UpN>(); chk<ApN,Times,ApN,LoN>(); chk<ApN,Times,ApN,ApN>();

    // Mixed Generic # Concrete numbers
    chk<ApD,Times,ExN,D>(); chk<ExN,Times,ExN,I>(); chk<ExN,Times,ExN,Z>(); chk<ExN,Times,ExN,Q>(); chk<EfN,Times,ExN,R>();
    chk<ApD,Times,EfN,D>(); chk<EfN,Times,EfN,I>(); chk<EfN,Times,EfN,Z>(); chk<EfN,Times,EfN,Q>(); chk<EfN,Times,EfN,R>();
    chk<ApD,Times,VaN,D>(); chk<VaN,Times,VaN,I>(); chk<VaN,Times,VaN,Z>(); chk<VaN,Times,VaN,Q>(); chk<VaN,Times,VaN,R>();
    chk<ApD,Times,UpN,D>(); chk<UpN,Times,UpN,I>(); chk<UpN,Times,UpN,Z>(); chk<UpN,Times,UpN,Q>(); chk<UpN,Times,UpN,R>();
    chk<ApD,Times,LoN,D>(); chk<LoN,Times,LoN,I>(); chk<LoN,Times,LoN,Z>(); chk<LoN,Times,LoN,Q>(); chk<LoN,Times,LoN,R>();
    chk<ApD,Times,ApN,D>(); chk<ApN,Times,ApN,I>(); chk<ApN,Times,ApN,Z>(); chk<ApN,Times,ApN,Q>(); chk<ApN,Times,ApN,R>();

    // Mixed Concrete # Generic numbers
    chk<ApD,Times,D,ExN>(); chk<ApD,Times,D,EfN>(); chk<ApD,Times,D,VaN>(); chk<ApD,Times,D,UpN>(); chk<ApD,Times,D,LoN>(); chk<ApD,Times,D,ApN>();
    chk<ExN,Times,I,ExN>(); chk<EfN,Times,I,EfN>(); chk<VaN,Times,I,VaN>(); chk<UpN,Times,I,UpN>(); chk<LoN,Times,I,LoN>(); chk<ApN,Times,I,ApN>();
    chk<ExN,Times,Z,ExN>(); chk<EfN,Times,Z,EfN>(); chk<VaN,Times,Z,VaN>(); chk<UpN,Times,Z,UpN>(); chk<LoN,Times,Z,LoN>(); chk<ApN,Times,Z,ApN>();
    chk<ExN,Times,Q,ExN>(); chk<EfN,Times,Q,EfN>(); chk<VaN,Times,Q,VaN>(); chk<UpN,Times,Q,UpN>(); chk<LoN,Times,Q,LoN>(); chk<ApN,Times,Q,ApN>();
    chk<EfN,Times,R,ExN>(); chk<EfN,Times,R,EfN>(); chk<VaN,Times,R,VaN>(); chk<UpN,Times,R,UpN>(); chk<LoN,Times,R,LoN>(); chk<ApN,Times,R,ApN>();

    // Float numbers
    chk<VaF,Times,ExF,ExF>(); chk<MeF,Times,ExF,MeF>(); chk<BoF,Times,ExF,BoF>(); chk<UpF,Times,ExF,UpF>(); chk<LoF,Times,ExF,LoF>(); chk<ApF,Times,ExF,ApF>();
    chk<MeF,Times,MeF,ExF>(); chk<MeF,Times,MeF,MeF>(); chk<VaF,Times,MeF,BoF>(); chk<UpF,Times,MeF,UpF>(); chk<LoF,Times,MeF,LoF>(); chk<ApF,Times,MeF,ApF>();
    chk<BoF,Times,BoF,ExF>(); chk<VaF,Times,BoF,MeF>(); chk<BoF,Times,BoF,BoF>(); chk<UpF,Times,BoF,UpF>(); chk<LoF,Times,BoF,LoF>(); chk<ApF,Times,BoF,ApF>();
    chk<UpF,Times,UpF,ExF>(); chk<UpF,Times,UpF,MeF>(); chk<UpF,Times,UpF,BoF>(); chk<UpF,Times,UpF,UpF>(); chk<ApF,Times,UpF,LoF>(); chk<ApF,Times,UpF,ApF>();
    chk<LoF,Times,LoF,ExF>(); chk<LoF,Times,LoF,MeF>(); chk<LoF,Times,LoF,BoF>(); chk<ApF,Times,LoF,UpF>(); chk<LoF,Times,LoF,LoF>(); chk<ApF,Times,LoF,ApF>();
    chk<ApF,Times,ApF,ExF>(); chk<ApF,Times,ApF,MeF>(); chk<ApF,Times,ApF,BoF>(); chk<ApF,Times,ApF,UpF>(); chk<ApF,Times,ApF,LoF>(); chk<ApF,Times,ApF,ApF>();

    // Mixed Float # Generic and Generic # Float
    chk<VaF,Times,ExF,ExN>(); chk<VaF,Times,ExF,EfN>(); chk<VaF,Times,ExF,VaN>(); chk<UpF,Times,ExF,UpN>(); chk<LoF,Times,ExF,LoN>(); chk<ApF,Times,ExF,ApN>();
    chk<MeF,Times,MeF,ExN>(); chk<MeF,Times,MeF,EfN>(); chk<MeF,Times,MeF,VaN>(); chk<UpF,Times,MeF,UpN>(); chk<LoF,Times,MeF,LoN>(); chk<ApF,Times,MeF,ApN>();
    chk<BoF,Times,BoF,ExN>(); chk<BoF,Times,BoF,EfN>(); chk<BoF,Times,BoF,VaN>(); chk<UpF,Times,BoF,UpN>(); chk<LoF,Times,BoF,LoN>(); chk<ApF,Times,BoF,ApN>();
    chk<UpF,Times,UpF,ExN>(); chk<UpF,Times,UpF,EfN>(); chk<UpF,Times,UpF,VaN>(); chk<UpF,Times,UpF,UpN>(); chk<ApF,Times,UpF,LoN>(); chk<ApF,Times,UpF,ApN>();
    chk<LoF,Times,LoF,ExN>(); chk<LoF,Times,LoF,EfN>(); chk<LoF,Times,LoF,VaN>(); chk<ApF,Times,LoF,UpN>(); chk<LoF,Times,LoF,LoN>(); chk<ApF,Times,LoF,ApN>();
    chk<ApF,Times,ApF,ExN>(); chk<ApF,Times,ApF,EfN>(); chk<ApF,Times,ApF,VaN>(); chk<ApF,Times,ApF,UpN>(); chk<ApF,Times,ApF,LoN>(); chk<ApF,Times,ApF,ApN>();

    chk<VaF,Times,ExN,ExF>(); chk<MeF,Times,ExN,MeF>(); chk<BoF,Times,ExN,BoF>(); chk<UpF,Times,ExN,UpF>(); chk<LoF,Times,ExN,LoF>(); chk<ApF,Times,ExN,ApF>();
    chk<VaF,Times,EfN,ExF>(); chk<MeF,Times,EfN,MeF>(); chk<BoF,Times,EfN,BoF>(); chk<UpF,Times,EfN,UpF>(); chk<LoF,Times,EfN,LoF>(); chk<ApF,Times,EfN,ApF>();
    chk<VaF,Times,VaN,ExF>(); chk<MeF,Times,VaN,MeF>(); chk<BoF,Times,VaN,BoF>(); chk<UpF,Times,VaN,UpF>(); chk<LoF,Times,VaN,LoF>(); chk<ApF,Times,VaN,ApF>();
    chk<UpF,Times,UpN,ExF>(); chk<UpF,Times,UpN,MeF>(); chk<UpF,Times,UpN,BoF>(); chk<UpF,Times,UpN,UpF>(); chk<ApF,Times,UpN,LoF>(); chk<ApF,Times,UpN,ApF>();
    chk<LoF,Times,LoN,ExF>(); chk<LoF,Times,LoN,MeF>(); chk<LoF,Times,LoN,BoF>(); chk<ApF,Times,LoN,UpF>(); chk<LoF,Times,LoN,LoF>(); chk<ApF,Times,LoN,ApF>();
    chk<ApF,Times,ApN,ExF>(); chk<ApF,Times,ApN,MeF>(); chk<ApF,Times,ApN,BoF>(); chk<ApF,Times,ApN,UpF>(); chk<ApF,Times,ApN,LoF>(); chk<ApF,Times,ApN,ApF>();

    // Mixed Float # Concrete and Concrete # Float
    chk<ApF,Times,ExF,D>(); chk<VaF,Times,ExF,I>(); chk<VaF,Times,ExF,Z>(); chk<VaF,Times,ExF,Q>(); chk<VaF,Times,ExF,R>();
    chk<ApF,Times,MeF,D>(); chk<MeF,Times,MeF,I>(); chk<MeF,Times,MeF,Z>(); chk<MeF,Times,MeF,Q>(); chk<MeF,Times,MeF,R>();
    chk<ApF,Times,BoF,D>(); chk<BoF,Times,BoF,I>(); chk<BoF,Times,BoF,Z>(); chk<BoF,Times,BoF,Q>(); chk<BoF,Times,BoF,R>();
    chk<ApF,Times,UpF,D>(); chk<UpF,Times,UpF,I>(); chk<UpF,Times,UpF,Z>(); chk<UpF,Times,UpF,Q>(); chk<UpF,Times,UpF,R>();
    chk<ApF,Times,LoF,D>(); chk<LoF,Times,LoF,I>(); chk<LoF,Times,LoF,Z>(); chk<LoF,Times,LoF,Q>(); chk<LoF,Times,LoF,R>();
    chk<ApF,Times,ApF,D>(); chk<ApF,Times,ApF,I>(); chk<ApF,Times,ApF,Z>(); chk<ApF,Times,ApF,Q>(); chk<ApF,Times,ApF,R>();

    chk<ApF,Times,D,ExF>(); chk<ApF,Times,D,MeF>(); chk<ApF,Times,D,BoF>(); chk<ApF,Times,D,UpF>(); chk<ApF,Times,D,LoF>(); chk<ApF,Times,D,ApF>();
    chk<VaF,Times,I,ExF>(); chk<MeF,Times,I,MeF>(); chk<BoF,Times,I,BoF>(); chk<UpF,Times,I,UpF>(); chk<LoF,Times,I,LoF>(); chk<ApF,Times,I,ApF>();
    chk<VaF,Times,Z,ExF>(); chk<MeF,Times,Z,MeF>(); chk<BoF,Times,Z,BoF>(); chk<UpF,Times,Z,UpF>(); chk<LoF,Times,Z,LoF>(); chk<ApF,Times,Z,ApF>();
    chk<VaF,Times,Q,ExF>(); chk<MeF,Times,Q,MeF>(); chk<BoF,Times,Q,BoF>(); chk<UpF,Times,Q,UpF>(); chk<LoF,Times,Q,LoF>(); chk<ApF,Times,Q,ApF>();
    chk<VaF,Times,R,ExF>(); chk<MeF,Times,R,MeF>(); chk<BoF,Times,R,BoF>(); chk<UpF,Times,R,UpF>(); chk<LoF,Times,R,LoF>(); chk<ApF,Times,R,ApF>();

}

void CheckNumeric::check_division() {

    // Concrete numbers
    chk< D ,Divides,D,D>(); chk<D,Divides,D,I>(); chk<ApD,Divides,D,Z>(); chk<ApD,Divides,D,Q>(); chk<ApD,Divides,D,R>();
    chk< D ,Divides,I,D>(); chk<I,Divides,I,I>(); chk<Q,Divides,I,Z>(); chk<Q,Divides,I,Q>(); chk<R,Divides,I,R>();
    chk<ApD,Divides,Z,D>(); chk<Q,Divides,Z,I>(); chk<Q,Divides,Z,Z>(); chk<Q,Divides,Z,Q>(); chk<R,Divides,Z,R>();
    chk<ApD,Divides,Q,D>(); chk<Q,Divides,Q,I>(); chk<Q,Divides,Q,Z>(); chk<Q,Divides,Q,Q>(); chk<R,Divides,Q,R>();
    chk<ApD,Divides,R,D>(); chk<R,Divides,R,I>(); chk<R,Divides,R,Z>(); chk<R,Divides,R,Q>(); chk<R,Divides,R,R>();

    // Generic numbers
    chk<ExN,Divides,ExN,ExN>(); chk<EfN,Divides,ExN,EfN>(); chk<VaN,Divides,ExN,VaN>(); chk<LoN,Divides,ExN,UpN>(); chk<UpN,Divides,ExN,LoN>(); chk<ApN,Divides,ExN,ApN>();
    chk<EfN,Divides,EfN,ExN>(); chk<EfN,Divides,EfN,EfN>(); chk<VaN,Divides,EfN,VaN>(); chk<LoN,Divides,EfN,UpN>(); chk<UpN,Divides,EfN,LoN>(); chk<ApN,Divides,EfN,ApN>();
    chk<VaN,Divides,VaN,ExN>(); chk<VaN,Divides,VaN,EfN>(); chk<VaN,Divides,VaN,VaN>(); chk<LoN,Divides,VaN,UpN>(); chk<UpN,Divides,VaN,LoN>(); chk<ApN,Divides,VaN,ApN>();
    chk<UpN,Divides,UpN,ExN>(); chk<UpN,Divides,UpN,EfN>(); chk<UpN,Divides,UpN,VaN>(); chk<ApN,Divides,UpN,UpN>(); chk<UpN,Divides,UpN,LoN>(); chk<ApN,Divides,UpN,ApN>();
    chk<LoN,Divides,LoN,ExN>(); chk<LoN,Divides,LoN,EfN>(); chk<LoN,Divides,LoN,VaN>(); chk<LoN,Divides,LoN,UpN>(); chk<ApN,Divides,LoN,LoN>(); chk<ApN,Divides,LoN,ApN>();
    chk<ApN,Divides,ApN,ExN>(); chk<ApN,Divides,ApN,EfN>(); chk<ApN,Divides,ApN,VaN>(); chk<ApN,Divides,ApN,UpN>(); chk<ApN,Divides,ApN,LoN>(); chk<ApN,Divides,ApN,ApN>();

    // Mixed Generic # Concrete numbers
    chk<ApD,Divides,ExN,D>(); chk<ExN,Divides,ExN,I>(); chk<ExN,Divides,ExN,Z>(); chk<ExN,Divides,ExN,Q>(); chk<EfN,Divides,ExN,R>();
    chk<ApD,Divides,EfN,D>(); chk<EfN,Divides,EfN,I>(); chk<EfN,Divides,EfN,Z>(); chk<EfN,Divides,EfN,Q>(); chk<EfN,Divides,EfN,R>();
    chk<ApD,Divides,VaN,D>(); chk<VaN,Divides,VaN,I>(); chk<VaN,Divides,VaN,Z>(); chk<VaN,Divides,VaN,Q>(); chk<VaN,Divides,VaN,R>();
    chk<ApD,Divides,UpN,D>(); chk<UpN,Divides,UpN,I>(); chk<UpN,Divides,UpN,Z>(); chk<UpN,Divides,UpN,Q>(); chk<UpN,Divides,UpN,R>();
    chk<ApD,Divides,LoN,D>(); chk<LoN,Divides,LoN,I>(); chk<LoN,Divides,LoN,Z>(); chk<LoN,Divides,LoN,Q>(); chk<LoN,Divides,LoN,R>();
    chk<ApD,Divides,ApN,D>(); chk<ApN,Divides,ApN,I>(); chk<ApN,Divides,ApN,Z>(); chk<ApN,Divides,ApN,Q>(); chk<ApN,Divides,ApN,R>();

    // Mixed Concrete # Generic numbers
    chk<ApD,Divides,D,ExN>(); chk<ApD,Divides,D,EfN>(); chk<ApD,Divides,D,VaN>(); chk<ApD,Divides,D,UpN>(); chk<ApD,Divides,D,LoN>(); chk<ApD,Divides,D,ApN>();
    chk<ExN,Divides,I,ExN>(); chk<EfN,Divides,I,EfN>(); chk<VaN,Divides,I,VaN>(); chk<LoN,Divides,I,UpN>(); chk<UpN,Divides,I,LoN>(); chk<ApN,Divides,I,ApN>();
    chk<ExN,Divides,Z,ExN>(); chk<EfN,Divides,Z,EfN>(); chk<VaN,Divides,Z,VaN>(); chk<LoN,Divides,Z,UpN>(); chk<UpN,Divides,Z,LoN>(); chk<ApN,Divides,Z,ApN>();
    chk<ExN,Divides,Q,ExN>(); chk<EfN,Divides,Q,EfN>(); chk<VaN,Divides,Q,VaN>(); chk<LoN,Divides,Q,UpN>(); chk<UpN,Divides,Q,LoN>(); chk<ApN,Divides,Q,ApN>();
    chk<EfN,Divides,R,ExN>(); chk<EfN,Divides,R,EfN>(); chk<VaN,Divides,R,VaN>(); chk<LoN,Divides,R,UpN>(); chk<UpN,Divides,R,LoN>(); chk<ApN,Divides,R,ApN>();

    // Float numbers
    chk<VaF,Divides,ExF,ExF>(); chk<MeF,Divides,ExF,MeF>(); chk<BoF,Divides,ExF,BoF>(); chk<LoF,Divides,ExF,UpF>(); chk<UpF,Divides,ExF,LoF>(); chk<ApF,Divides,ExF,ApF>();
    chk<MeF,Divides,MeF,ExF>(); chk<MeF,Divides,MeF,MeF>(); chk<VaF,Divides,MeF,BoF>(); chk<LoF,Divides,MeF,UpF>(); chk<UpF,Divides,MeF,LoF>(); chk<ApF,Divides,MeF,ApF>();
    chk<BoF,Divides,BoF,ExF>(); chk<VaF,Divides,BoF,MeF>(); chk<BoF,Divides,BoF,BoF>(); chk<LoF,Divides,BoF,UpF>(); chk<UpF,Divides,BoF,LoF>(); chk<ApF,Divides,BoF,ApF>();
    chk<UpF,Divides,UpF,ExF>(); chk<UpF,Divides,UpF,MeF>(); chk<UpF,Divides,UpF,BoF>(); chk<ApF,Divides,UpF,UpF>(); chk<UpF,Divides,UpF,LoF>(); chk<ApF,Divides,UpF,ApF>();
    chk<LoF,Divides,LoF,ExF>(); chk<LoF,Divides,LoF,MeF>(); chk<LoF,Divides,LoF,BoF>(); chk<LoF,Divides,LoF,UpF>(); chk<ApF,Divides,LoF,LoF>(); chk<ApF,Divides,LoF,ApF>();
    chk<ApF,Divides,ApF,ExF>(); chk<ApF,Divides,ApF,MeF>(); chk<ApF,Divides,ApF,BoF>(); chk<ApF,Divides,ApF,UpF>(); chk<ApF,Divides,ApF,LoF>(); chk<ApF,Divides,ApF,ApF>();

    // Mixed Float # Generic and Generic # Float
    chk<VaF,Divides,ExF,ExN>(); chk<VaF,Divides,ExF,EfN>(); chk<VaF,Divides,ExF,VaN>(); chk<LoF,Divides,ExF,UpN>(); chk<UpF,Divides,ExF,LoN>(); chk<ApF,Divides,ExF,ApN>();
    chk<MeF,Divides,MeF,ExN>(); chk<MeF,Divides,MeF,EfN>(); chk<MeF,Divides,MeF,VaN>(); chk<LoF,Divides,MeF,UpN>(); chk<UpF,Divides,MeF,LoN>(); chk<ApF,Divides,MeF,ApN>();
    chk<BoF,Divides,BoF,ExN>(); chk<BoF,Divides,BoF,EfN>(); chk<BoF,Divides,BoF,VaN>(); chk<LoF,Divides,BoF,UpN>(); chk<UpF,Divides,BoF,LoN>(); chk<ApF,Divides,BoF,ApN>();
    chk<UpF,Divides,UpF,ExN>(); chk<UpF,Divides,UpF,EfN>(); chk<UpF,Divides,UpF,VaN>(); chk<ApF,Divides,UpF,UpN>(); chk<UpF,Divides,UpF,LoN>(); chk<ApF,Divides,UpF,ApN>();
    chk<LoF,Divides,LoF,ExN>(); chk<LoF,Divides,LoF,EfN>(); chk<LoF,Divides,LoF,VaN>(); chk<LoF,Divides,LoF,UpN>(); chk<ApF,Divides,LoF,LoN>(); chk<ApF,Divides,LoF,ApN>();
    chk<ApF,Divides,ApF,ExN>(); chk<ApF,Divides,ApF,EfN>(); chk<ApF,Divides,ApF,VaN>(); chk<ApF,Divides,ApF,UpN>(); chk<ApF,Divides,ApF,LoN>(); chk<ApF,Divides,ApF,ApN>();

    chk<VaF,Divides,ExN,ExF>(); chk<MeF,Divides,ExN,MeF>(); chk<BoF,Divides,ExN,BoF>(); chk<LoF,Divides,ExN,UpF>(); chk<UpF,Divides,ExN,LoF>(); chk<ApF,Divides,ExN,ApF>();
    chk<VaF,Divides,EfN,ExF>(); chk<MeF,Divides,EfN,MeF>(); chk<BoF,Divides,EfN,BoF>(); chk<LoF,Divides,EfN,UpF>(); chk<UpF,Divides,EfN,LoF>(); chk<ApF,Divides,EfN,ApF>();
    chk<VaF,Divides,VaN,ExF>(); chk<MeF,Divides,VaN,MeF>(); chk<BoF,Divides,VaN,BoF>(); chk<LoF,Divides,VaN,UpF>(); chk<UpF,Divides,VaN,LoF>(); chk<ApF,Divides,VaN,ApF>();
    chk<UpF,Divides,UpN,ExF>(); chk<UpF,Divides,UpN,MeF>(); chk<UpF,Divides,UpN,BoF>(); chk<ApF,Divides,UpN,UpF>(); chk<UpF,Divides,UpN,LoF>(); chk<ApF,Divides,UpN,ApF>();
    chk<LoF,Divides,LoN,ExF>(); chk<LoF,Divides,LoN,MeF>(); chk<LoF,Divides,LoN,BoF>(); chk<LoF,Divides,LoN,UpF>(); chk<ApF,Divides,LoN,LoF>(); chk<ApF,Divides,LoN,ApF>();
    chk<ApF,Divides,ApN,ExF>(); chk<ApF,Divides,ApN,MeF>(); chk<ApF,Divides,ApN,BoF>(); chk<ApF,Divides,ApN,UpF>(); chk<ApF,Divides,ApN,LoF>(); chk<ApF,Divides,ApN,ApF>();

    // Mixed Float # Concrete and Concrete # Float
    chk<ApF,Divides,ExF,D>(); chk<VaF,Divides,ExF,I>(); chk<VaF,Divides,ExF,Z>(); chk<VaF,Divides,ExF,Q>(); chk<VaF,Divides,ExF,R>();
    chk<ApF,Divides,MeF,D>(); chk<MeF,Divides,MeF,I>(); chk<MeF,Divides,MeF,Z>(); chk<MeF,Divides,MeF,Q>(); chk<MeF,Divides,MeF,R>();
    chk<ApF,Divides,BoF,D>(); chk<BoF,Divides,BoF,I>(); chk<BoF,Divides,BoF,Z>(); chk<BoF,Divides,BoF,Q>(); chk<BoF,Divides,BoF,R>();
    chk<ApF,Divides,UpF,D>(); chk<UpF,Divides,UpF,I>(); chk<UpF,Divides,UpF,Z>(); chk<UpF,Divides,UpF,Q>(); chk<UpF,Divides,UpF,R>();
    chk<ApF,Divides,LoF,D>(); chk<LoF,Divides,LoF,I>(); chk<LoF,Divides,LoF,Z>(); chk<LoF,Divides,LoF,Q>(); chk<LoF,Divides,LoF,R>();
    chk<ApF,Divides,ApF,D>(); chk<ApF,Divides,ApF,I>(); chk<ApF,Divides,ApF,Z>(); chk<ApF,Divides,ApF,Q>(); chk<ApF,Divides,ApF,R>();

    chk<ApF,Divides,D,ExF>(); chk<ApF,Divides,D,MeF>(); chk<ApF,Divides,D,BoF>(); chk<ApF,Divides,D,UpF>(); chk<ApF,Divides,D,LoF>(); chk<ApF,Divides,D,ApF>();
    chk<VaF,Divides,I,ExF>(); chk<MeF,Divides,I,MeF>(); chk<BoF,Divides,I,BoF>(); chk<LoF,Divides,I,UpF>(); chk<UpF,Divides,I,LoF>(); chk<ApF,Divides,I,ApF>();
    chk<VaF,Divides,Z,ExF>(); chk<MeF,Divides,Z,MeF>(); chk<BoF,Divides,Z,BoF>(); chk<LoF,Divides,Z,UpF>(); chk<UpF,Divides,Z,LoF>(); chk<ApF,Divides,Z,ApF>();
    chk<VaF,Divides,Q,ExF>(); chk<MeF,Divides,Q,MeF>(); chk<BoF,Divides,Q,BoF>(); chk<LoF,Divides,Q,UpF>(); chk<UpF,Divides,Q,LoF>(); chk<ApF,Divides,Q,ApF>();
    chk<VaF,Divides,R,ExF>(); chk<MeF,Divides,R,MeF>(); chk<BoF,Divides,R,BoF>(); chk<LoF,Divides,R,UpF>(); chk<UpF,Divides,R,LoF>(); chk<ApF,Divides,R,ApF>();
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
