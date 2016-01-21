#include <type_traits>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <memory>

using std::declval;
static const bool dummy = true;
template<class P, class T=bool> using EnableIf = typename std::enable_if<P::value,T>::type;
template<class P, class T=bool> using DisableIf = typename std::enable_if<not P::value,T>::type;
using True = std::true_type;
using False = std::false_type;
template<class F, class T> using IsConvertible = std::is_convertible<F,T>;
template<class T, class... F> using IsConstructible = std::is_constructible<T,F...>;
template<class T1, class T2> using IsSame = std::is_same<T1,T2>;
template<class X> using InitializerList = std::initializer_list<X>;
using Void = void;
using SizeType = std::size_t;

using OutputStream = std::ostream;

#define TEMPLATE_FLOAT
#ifndef TEMPLATE_FLOAT
#define FRIEND_FLOAT
#endif

auto const RESET   = "\033[0m";
auto const BLACK   = "\033[30m";      /* Black */
auto const RED     = "\033[31m";      /* Red */
auto const GREEN   = "\033[32m";      /* Green */
auto const YELLOW  = "\033[33m";      /* Yellow */
auto const BLUE    = "\033[34m";      /* Blue */
auto const MAGENTA = "\033[35m";      /* Magenta */
auto const CYAN    = "\033[36m";      /* Cyan */
auto const WHITE   = "\033[37m";      /* White */

// ---------------- Indexing --------------------------------------------------------------------------------------- //

class Index { SizeType _i; public: Index(SizeType i) : _i(i) { } operator SizeType() const { return _i; } };


// ---------------- Operators -------------------------------------------------------------------------------------- //

struct Cnst { };
struct Var { };

struct And {
    template<class A1, class A2> auto operator() (A1&& a1, A2&& a2) const -> decltype(conjunction(a1,a2)) {
        return conjunction(std::forward<A1>(a1),std::forward<A2>(a2)); }
};
struct Or {
    template<class A1, class A2> auto operator() (A1&& a1, A2&& a2) const -> decltype(disjunction(a1,a2)) {
        return disjunction(std::forward<A1>(a1),std::forward<A2>(a2)); }
};
struct Not { template<class A> auto operator() (A&& a) const -> decltype(negation(a)) { return negation(std::forward<A>(a)); } };


struct Max {
    template<class A1, class A2> auto operator() (A1&& a1, A2&& a2) const -> decltype(max(a1,a2)) {
        return max(std::forward<A1>(a1),std::forward<A2>(a2)); }
};
struct Add {
    template<class A1, class A2> auto operator() (A1&& a1, A2&& a2) const -> decltype(add(a1,a2)) {
        return add(std::forward<A1>(a1),std::forward<A2>(a2)); }
};
struct Mul {
    template<class A1, class A2> auto operator() (A1&& a1, A2&& a2) const -> decltype(mul(a1,a2)) {
        return mul(std::forward<A1>(a1),std::forward<A2>(a2)); }
};

struct Neg { template<class A> auto operator() (A&& a) const -> decltype(neg(a)) { return neg(std::forward<A>(a)); } };
struct Rec { template<class A> auto operator() (A&& a) const -> decltype(rec(a)) { return rec(std::forward<A>(a)); } };
struct Exp { template<class A> auto operator() (A&& a) const -> decltype(exp(a)) { return exp(std::forward<A>(a)); } };
struct Log { template<class A> auto operator() (A&& a) const -> decltype(log(a)) { return log(std::forward<A>(a)); } };

// ---------------- Traits ---------------------------------------------------------------------------------------- //

template<class T> struct Self { typedef T Type; };
template<class T> using SelfType = typename Self<T>::Type;
template<class T> using GenericType = typename T::GenericType;
template<class T> using NumericType = typename T::NumericType;

template<class T> struct IndexTraits { typedef Void Type; };
template<class T> using IndexType = typename IndexTraits<T>::Type;

//template<class T> struct NumericTypedef { typedef typename T::NumericType Type; };
//template<class T> using NumericType = typename NumericTypedef<T>::Type;

template<class T> using GenericNumericType = GenericType<NumericType<T>>;

template<class T1, class T2=Void> struct ArithmeticTraits {
    typedef decltype(add(declval<T1>(),declval<T2>())) AddType;
    typedef decltype(mul(declval<T1>(),declval<T2>())) MulType;
};

template<class T> struct ArithmeticTraits<T> {
    typedef decltype(add(declval<T>(),declval<T>())) AddType;
    typedef decltype(mul(declval<T>(),declval<T>())) MulType;
    typedef decltype(neg(declval<T>())) NegType;
    typedef decltype(rec(declval<T>())) RecType;
    typedef decltype(exp(declval<T>())) ExpType;
    typedef decltype(log(declval<T>())) LogType;
};

template<class T1, class T2=T1> using AddType = typename ArithmeticTraits<T1,T2>::AddType;
template<class T1, class T2=T1> using MulType = typename ArithmeticTraits<T1,T2>::MulType;
template<class T> using NegType = typename ArithmeticTraits<T>::NegType;
template<class T> using RecType = typename ArithmeticTraits<T>::RecType;
template<class T> using ExpType = typename ArithmeticTraits<T>::ExpType;
template<class T> using LogType = typename ArithmeticTraits<T>::LogType;

// ---------------- Declarators ------------------------------------------------------------------------------------ //

template<class X> struct DeclareRingOperators {
    friend X add(X,X); friend X mul(X,X); friend X neg(X);
};
template<class X> struct DeclareFieldOperators : DeclareRingOperators<X> {
   friend X rec(X);
};
template<class X> struct DeclareOperators : DeclareFieldOperators<X> {
   friend X exp(X);friend X log(X);
   friend X max(X,X);
};

template<class A, class X=typename A::NumericType> struct DeclareAlgbraOperators : DeclareOperators<A> {
    friend A add(A,X); friend A add(X,A); friend A mul(A,X); friend A mul(X,A);
};

template<class A, class X=typename A::NumericType> struct DeclareStaticAlgbraOperators  {
    static A _add(A,A); static A _mul(A,A);
    static A _neg(A); static A _rec(A); static A _exp(A); static A _log(A);
    static A _add(A,X); static A _add(X,A); static A _mul(A,X); static A _mul(X,A);
};

template<class A> struct Operations {
    typedef typename A::NumericType X;
    static A _add(A,A); static A _mul(A,A);
    static A _neg(A); static A _rec(A); static A _exp(A); static A _log(A);
    static A _add(A,X); static A _add(X,A); static A _mul(A,X); static A _mul(X,A);
};

template<class A> struct StaticDispatchNumericOperators {
    friend A add(A a1, A a2) { return Operations<A>::_add(a1,a2); }
    friend A mul(A a1, A a2) { return Operations<A>::_mul(a1,a2); }
    friend A neg(A a) { return Operations<A>::_neg(a); }
    friend A rec(A a) { return Operations<A>::_rec(a); }
    friend A exp(A a) { return Operations<A>::_exp(a); }
    friend A log(A a) { return Operations<A>::_log(a); }
};

template<class A, class X=typename A::NumericType> struct StaticDispatchAlgbraOperators : StaticDispatchNumericOperators<A> {
    friend A add(A a1, X x2) { return Operations<A>::_add(a1,x2); }
    friend A add(X x1, A a2) { return Operations<A>::_add(x1,a2); }
    friend A mul(A a1, X x2) { return Operations<A>::_mul(a1,x2); }
    friend A mul(X x1, A a2) { return Operations<A>::_mul(x1,a2); }
};

template<template<typename>class F> struct DeclareTemplateOperators {
    template<class X> friend F<X> add(F<X>,F<X>);
    template<class X> friend F<X> mul(F<X>,F<X>);
    template<class X> friend F<X> neg(F<X>);
    template<class X> friend F<X> rec(F<X>);
    template<class X> friend F<X> exp(F<X>);
    template<class X> friend F<X> log(F<X>);
};

template<template<typename>class F> struct DeclareTemplateMixedSelfOperators {
    template<class X> friend F<X> add(F<X>,SelfType<F<X>>);
    template<class X> friend F<X> add(SelfType<F<X>>,F<X>);
    template<class X> friend F<X> mul(F<X>,SelfType<F<X>>);
    template<class X> friend F<X> mul(SelfType<F<X>>,F<X>);
};

template<template<typename>class F> struct DeclareTemplateAlgebraOperators {
    template<class X> friend F<X> add(F<X>,F<X>);
    template<class X> friend F<X> mul(F<X>,F<X>);
    template<class X> friend F<X> neg(F<X>);
    template<class X> friend F<X> rec(F<X>);
    template<class X> friend F<X> exp(F<X>);
    template<class X> friend F<X> log(F<X>);
    template<class X> friend F<X> add(F<X>,NumericType<F<X>>);
    template<class X> friend F<X> add(NumericType<F<X>>,F<X>);
    template<class X> friend F<X> mul(F<X>,NumericType<F<X>>);
    template<class X> friend F<X> mul(NumericType<F<X>>,F<X>);
};

template<class X> struct DeclareLogicalOperators {
   friend X conjunction(X,X); friend X disjunction(X,X);
   friend X negation(X);
};

template<template<typename>class F> struct DeclareTemplateLogicalOperators {
    template<class X> friend F<X> conjunction(F<X>,F<X>);
    template<class X> friend F<X> disjunction(F<X>,F<X>);
    template<class X> friend F<X> negation(F<X>);
};

template<template<typename>class F> struct DefineTemplateMixedSelfOperators {
    template<class X> friend F<X> add(F<X> x1, SelfType<F<X>> x2) { return add(x1,x2); }
    template<class X> friend F<X> add(SelfType<F<X>> x1, F<X> x2) { return add(x1,x2); }
    template<class X> friend F<X> mul(F<X> x1, SelfType<F<X>> x2) { return mul(x1,x2); }
    template<class X> friend F<X> mul(SelfType<F<X>> x1, F<X> x2) { return mul(x1,x2); }
};

struct Fooby { };
//Fooby conjunction(Fooby,Fooby);
//Fooby disjunction(Fooby,Fooby);
//Fooby negation(Fooby);


// ---------------- Effort ----------------------------------------------------------------------------------------- //

struct Effort { };

// ---------------- Logical ---------------------------------------------------------------------------------------- //

class Boolean;

class Logical {
  public:
    enum class Value : char { FALSE=-2, UNLIKELY=-1, INDETERMINATE=0, LIKELY=+1, TRUE=2 };
  private:
    Value _v;
  public:
    Logical(Value v) : _v(v) { }
    explicit Logical(bool b) : _v(b?Value::TRUE:Value::FALSE) { }
    explicit Logical(Boolean b);
    friend Logical conjunction(Logical l1, Logical l2) { return l1._v<=l2._v?l1:l2; }
    friend Logical disjunction(Logical l1, Logical l2) { return l1._v>=l2._v?l1:l2; }
    friend Logical negation(Logical l) { return Logical(static_cast<Value>(-static_cast<char>(l._v))); }
    friend bool definitely(Logical l) { return l._v==Value::TRUE; }
    friend bool possibly(Logical l) { return l._v!=Value::FALSE; }
    friend bool likely(Logical l) { return l._v>=Value::LIKELY; }
    friend bool decide(Logical l) { return l._v>=Value::LIKELY; }
    friend OutputStream& operator<<(OutputStream&, Logical);
};

class Boolean {
    bool _b;
  public:
    Boolean(bool b) : _b(b) { }
    operator bool() const { return this->_b; }
    Logical value() const { return Logical(this->_b); }
    explicit operator Logical() const { return this->value(); }
    friend Boolean conjunction(Boolean l1, Boolean l2);
    friend Boolean disjunction(Boolean l1, Boolean l2);
    friend Boolean negation(Boolean l);
    friend bool decide(Boolean l) { return l._b; }
    friend OutputStream& operator<<(OutputStream&, Boolean);
};

inline Logical::Logical(Boolean l) : Logical(static_cast<bool>(l)) { }

class ValidatedKleenean {
    Logical _l;
  public:
    ValidatedKleenean(bool b) : _l(b) { }
    explicit ValidatedKleenean(Logical l) : _l(l) { }
    explicit operator Logical() const { return _l; }
    Logical check(Effort e) const;
    friend Logical check(ValidatedKleenean l, Effort e);
    friend ValidatedKleenean conjunction(ValidatedKleenean l1, ValidatedKleenean l2);
    friend ValidatedKleenean disjunction(ValidatedKleenean l1, ValidatedKleenean l2);
    friend ValidatedKleenean negation(ValidatedKleenean l);
    friend bool decide(ValidatedKleenean l, Effort e);
    friend OutputStream& operator<<(OutputStream&, ValidatedKleenean);
};

class ValidatedNegatedSierpinskian;

class ValidatedSierpinskian {
    Logical _l;
  public:
    ValidatedSierpinskian(bool b) : _l(b) { }
    ValidatedSierpinskian(ValidatedKleenean l) : _l(reinterpret_cast<Logical const&>(l)) { }
    explicit ValidatedSierpinskian(Logical l) : _l(l) { }
    explicit operator Logical() const { return _l; }
    Logical check(Effort e) const;
    friend Logical check(ValidatedSierpinskian l, Effort e);
    friend ValidatedSierpinskian conjunction(ValidatedSierpinskian l1, ValidatedSierpinskian l2);
    friend ValidatedSierpinskian disjunction(ValidatedSierpinskian l1, ValidatedSierpinskian l2);
    friend ValidatedNegatedSierpinskian negation(ValidatedSierpinskian l);
    friend bool decide(ValidatedSierpinskian l, Effort e);
    friend OutputStream& operator<<(OutputStream&, ValidatedSierpinskian);
};

class ValidatedNegatedSierpinskian {
    Logical _l;
  public:
    explicit ValidatedNegatedSierpinskian(Logical l) : _l(l) { }
    explicit operator Logical() const { return _l; }
    friend ValidatedSierpinskian negation(ValidatedNegatedSierpinskian l);
    friend OutputStream& operator<<(OutputStream&, ValidatedNegatedSierpinskian);
};

template<class L> using IsLogical = IsConstructible<Logical,L>;


template<class T1, class T2=Void> struct LogicalTraits {
    typedef decltype(conjunction(declval<T1>(),declval<T2>())) AndType;
    typedef decltype(disjunction(declval<T1>(),declval<T2>())) OrType;
};

template<class T> struct LogicalTraits<T> {
    typedef decltype(conjunction(declval<T>(),declval<T>())) AndType;
    typedef decltype(disjunction(declval<T>(),declval<T>())) OrType;
    typedef decltype(negation(declval<T>())) NotType;
};
template<class T1, class T2=T1> using AndType = typename LogicalTraits<T1,T2>::AndType;
template<class T1, class T2=T1> using OrType = typename LogicalTraits<T1,T2>::OrType;
template<class T> using NotType = typename LogicalTraits<T>::NotType;

// ---------------- Integer --------------------------------------------------------------------------------------- //

class Integer : DeclareRingOperators<Integer> {
  private: public:
    int _n;
  public:
    Integer(double) = delete;
    Integer(int n);
    int get_si() const { return _n; }
    friend Integer abs(Integer);
    friend Integer quot(Integer,Integer);
    friend Integer rem(Integer,Integer);
    friend Integer gcd(Integer,Integer);
    friend OutputStream& operator<<(OutputStream&,Integer const&);
};

// ---------------- Rational --------------------------------------------------------------------------------------- //

class Rational : DeclareFieldOperators<Rational> {
  private: public:
    Integer _num; Integer _den;
    Rational(double) = delete;
 public:
    typedef Rational GenericType;
    Rational(int n) : Rational(Integer(n)) { }
    Rational(Integer n, Integer d=1);
    double get_d() const;
    friend OutputStream& operator<<(OutputStream&,Rational const&);
};

// ---------------- Real ------------------------------------------------------------------------------------------- //

template<class F> class Bounds;
class Precision64;
class Float64;
typedef Bounds<Float64> Float64Bounds;

class RealInterface;

struct Real : DeclareOperators<Real> {
    Real(double) = delete;
  private: public:
    std::shared_ptr<RealInterface> _ptr;
    Real(RealInterface* p);
  public:
    typedef Real GenericType;
    Real(int n);
    Real(Integer z);
    Real(Rational q);
    Rational lower_bound(Effort e) const;
    Rational upper_bound(Effort e) const;
    Float64Bounds get(Precision64) const;

    friend OutputStream& operator<<(OutputStream&,Real const&);

    // FIXME: Should not be necessary...
    Real create_constant(Real) const;
};
class ValidatedReal : DeclareOperators<ValidatedReal> {
  private: public:
    Real _r;
  public:
    typedef ValidatedReal GenericType;
    template<class T, EnableIf<IsConvertible<T,Real>> =dummy> ValidatedReal(T t) : ValidatedReal(Real(t)) { }
    ValidatedReal(Real);
    Float64Bounds get(Precision64) const;

    friend OutputStream& operator<<(OutputStream&,ValidatedReal const&);
    ValidatedReal create_constant(ValidatedReal) const;
};

// ---------------- Float ------------------------------------------------------------------------------------------ //

#define TEMPLATE_FLOAT

struct Precision64 { };

class Float64 {
  private: public:
    double _dbl;
  public:
    typedef Precision64 PrecisionType;
    explicit Float64(double x=0.0) : _dbl(x) { }
    double get_d() { return _dbl; }
    friend OutputStream& operator<<(OutputStream&, Float64 const&);
};

template<class F> class Bounds
#ifdef TEMPLATE_FLOAT
    : DeclareTemplateOperators<Bounds>
#else
    : DeclareOperators<Bounds<F>>
#endif
{
  private: public:
    F _l; F _u;
  public:
    typedef F DataType;
    typedef typename F::PrecisionType PrecisionType;
    typedef ValidatedReal GenericType;

    Bounds<F>(double l, double u);
    Bounds<F>(F l, F u);
    Bounds<F>(PrecisionType);
    Bounds<F>(Integer z, PrecisionType pr) : Bounds(Rational(z),pr) { };
    Bounds<F>(Rational, PrecisionType);
    Bounds<F>(GenericType, PrecisionType);
    operator GenericType() const;
    Bounds<F> create(GenericType) const;
    Bounds<F> create_constant(GenericType) const;

    F lower_bound() const;
    F upper_bound() const;
    PrecisionType precision() const;
    F error_bound() const;
    friend OutputStream& operator<<(OutputStream&,Bounds<F> const&);
};

template<class F> Bounds<F> add(Bounds<F> x1, Bounds<F> x2);
template<class F> Bounds<F> mul(Bounds<F> x1, Bounds<F> x2);

template<class X> X add(X const& x, GenericType<X> const& y) { return add(x,x.create(y)); }
template<class X> X add(GenericType<X> const& y, X const& x) { return add(x.create(y),x); }
template<class X> X mul(X const& x, GenericType<X> const& y) { return mul(x,x.create(y)); }
template<class X> X mul(GenericType<X> const& y, X const& x) { return mul(x.create(y),x); }

typedef Bounds<Float64> Float64Bounds;


// ---------------- Field ------------------------------------------------------------------------------------------ //

class Field {
    typedef Field SelfType;

    friend SelfType add(SelfType,SelfType);
    friend SelfType mul(SelfType,SelfType);

    friend SelfType neg(SelfType);
    friend SelfType rec(SelfType);
    friend SelfType exp(SelfType);
    friend SelfType log(SelfType);
};

// ---------------- Scalar ----------------------------------------------------------------------------------------- //

template<class X> using Scalar = X;

// ---------------- Vector ----------------------------------------------------------------------------------------- //

template<class X> class Vector {
    std::vector<X> _ary;
  public:
    Vector(SizeType n, X x) : _ary(n,x) { }
    Vector(InitializerList<X> lst) : _ary(lst) { }
    Vector(std::vector<X> ary) : _ary(ary) { }
    X create_constant(X c) const { return c; }
    template<class Y, DisableIf<IsConvertible<Y,X>> =dummy> auto
        create_constant(Y const& c) const -> decltype(declval<X>().create_constant(c)) {
            const X& x0=_ary[0]; return x0.create_constant(c); }
    SizeType size() const { return _ary.size(); }
    X const& operator[] (Index i) const { return _ary[i]; }
    X& operator[] (Index i) { return _ary[i]; }
    friend OutputStream& operator<<(OutputStream& os, Vector<X> const& v) {
        for(SizeType i=0; i!=v.size(); ++i) { os << (i==0?"[":",") << v[i]; } return os << "]"; }
};
typedef Vector<Real> RealVector;
typedef Vector<ValidatedReal> ValidatedRealVector;

// ---------------- Covector --------------------------------------------------------------------------------------- //


template<class X> class Covector {
    std::vector<X> _ary;
  public:
    typedef Index IndexType;
    Covector(SizeType n, X x) : _ary(n,x) { }
    Covector(InitializerList<X> lst) : _ary(lst) { }
    Covector(std::vector<X> ary) : _ary(ary) { }
    X create_constant(X c) const { return c; }
    SizeType size() const { return _ary.size(); }
    X const& operator[] (Index j) const { return _ary[j]; }
    X& operator[] (Index j) { return _ary[j]; }
    friend OutputStream& operator<<(OutputStream& os, Covector<X> const& v) {
        for(SizeType j=0; j!=v.size(); ++j) { os << (j==0?"[":",") << v[j]; } return os << "]"; }
    friend Covector<AddType<X,X>> add(Covector<X> const&, Covector<X> const&);
    friend Covector<MulType<X>> mul(Covector<X> const&, Scalar<X> const&);
    friend Covector<MulType<X>> mul(Scalar<X> const&, Covector<X> const);
};
typedef Vector<Real> RealVector;
typedef Vector<ValidatedReal> ValidatedRealVector;


// ---------------- Algebra ---------------------------------------------------------------------------------------- //


template<class X> class Algebra;

template<class X> class AlgebraInterface {
  public:
    virtual ~AlgebraInterface<X> () = default;
  private: public:
    friend class Algebra<X>;
    virtual AlgebraInterface<X>* _create_constant(X c) const = 0;
    virtual AlgebraInterface<X>* _add(AlgebraInterface<X> const*) const = 0;
    virtual AlgebraInterface<X>* _mul(AlgebraInterface<X> const*) const = 0;
    virtual AlgebraInterface<X>* _neg() const = 0;
    virtual AlgebraInterface<X>* _rec() const = 0;
    virtual AlgebraInterface<X>* _exp() const = 0;
    virtual AlgebraInterface<X>* _log() const = 0;
    virtual AlgebraInterface<X>* _add(X const&) const = 0;
    virtual AlgebraInterface<X>* _radd(X const&) const = 0;
    virtual AlgebraInterface<X>* _mul(X const&) const = 0;
    virtual AlgebraInterface<X>* _rmul(X const&) const = 0;
    virtual OutputStream& _write(OutputStream&) const = 0;
};

template<class X> class Algebra
// : public Field
{
    std::shared_ptr<AlgebraInterface<X>> _ptr;
  private: public:
    explicit Algebra(AlgebraInterface<X>* p) : _ptr(p) { }
  public:
    typedef Algebra<X> SelfType;
    typedef X NumericType;
    typedef GenericType<X> GenericTypeNumericType;
    template<class A> A extract() const;
    friend OutputStream& operator<<(OutputStream& os, Algebra<X> const& a) { return a._ptr->_write(os); }

    SelfType create_constant(NumericType c) const { return Algebra<X>(this->_ptr->_create_constant(c)); }
    SelfType& operator=(NumericType c) { return (*this) = this->create_constant(c); }

    friend SelfType add(SelfType a1, SelfType a2) { return Algebra<X>(a1._ptr->_add(a2._ptr.operator->())); }
    friend SelfType mul(SelfType a1, SelfType a2) { return Algebra<X>(a1._ptr->_mul(a2._ptr.operator->())); }

    friend SelfType add(SelfType a1, NumericType c2) { return Algebra<X>(a1._ptr->_add(c2)); }
    friend SelfType add(NumericType c1, SelfType a2) { return Algebra<X>(a2._ptr->_radd(c1)); }
    friend SelfType mul(SelfType a1, NumericType c2) { return Algebra<X>(a1._ptr->_mul(c2)); }
    friend SelfType mul(NumericType c1, SelfType a2) { return Algebra<X>(a2._ptr->_rmul(c1)); }

    friend SelfType neg(SelfType a) { return Algebra<X>(a._ptr->_neg()); }
    friend SelfType rec(SelfType a) { return Algebra<X>(a._ptr->_rec()); }
    friend SelfType exp(SelfType a) { return Algebra<X>(a._ptr->_exp()); }
    friend SelfType log(SelfType a) { return Algebra<X>(a._ptr->_log()); }
};
template<class X> Algebra<X> make_handle(AlgebraInterface<X>* aip) { return Algebra<X>(aip); }

template<class A, class X> using IfAlgebra = EnableIf<IsConstructible<Algebra<X>,A>,A>;

using RealAlgebra = Algebra<Real>;
using ValidatedRealAlgebra = Algebra<ValidatedReal>;

// ---------------- Creation --------------------------------------------------------------------------------------- //

template<class X> struct Creator {
    struct Foo {}; typedef Foo P;

    X create(GenericType<X>) const;
    X create_constant(NumericType<X>) const;

    static X concrete(P, GenericType<X>);
    static X constant(P, NumericType<X>);
    static X coordinate(P, IndexType<X>);
    static X variable(P, NumericType<X>, IndexType<X>);
    static X identity(P);
};

// ---------------- Expression ------------------------------------------------------------------------------------- //

class Identifier : public std::string {
  public:
    Identifier(const char* s) : std::string(s) { }
    Identifier(std::string s) : std::string(s) { }
};

template<class X> class Expression;
template<class X> class Valuation;

template<class X> class Variable
    : StaticDispatchNumericOperators<Expression<X>>
{
    Identifier _name;
  public:
    explicit Variable(std::string name) : _name(name) { }
    Identifier name() const { return _name; }
    template<class Y> Expression<Y> create_constant(Y const& c) { return Expression<Y>(c); }
    // template<class Y> Y create_constant(Y const& c) { return c; }
    friend bool operator<(Variable<X> v1, Variable<X> v2) { return v1._name < v2._name; }
    friend OutputStream& operator<<(OutputStream& os, Variable<X> const& v) { return os << v._name; }
  private:
    friend X evaluate(Expression<X>, Valuation<X>);
};

template<class X> class Valuation {
    std::map<Identifier,X> _vals;
    typedef typename std::map<Identifier,X>::const_iterator ConstIterator;
  public:
    Valuation() : _vals() { };
    Valuation(std::map<Identifier,X> m) : _vals(m) { };
    Valuation(std::map<Variable<Real>,X> m) : _vals() {
        for(auto v:m) { this->insert(v.first.name(),v.second); } }
    std::map<Identifier,X> map() const { return this->_vals; }
    ConstIterator begin() const { return _vals.begin(); }
    ConstIterator end() const { return _vals.end(); }
    X create_constant(X const& c) const { return c; }
    template<class Y, DisableIf<IsConvertible<Y,X>> =dummy> auto
        create_constant(Y const& c) const -> decltype(declval<X>().create_constant(c)) {
            const X& x0=_vals.begin()->second; return x0.create_constant(c); }
    X const& operator[] (Identifier s) const { return _vals.find(s)->second; }
    X& operator[] (Identifier s) { return _vals[s]; }
    X const& operator[] (Variable<Real> v) const { return (*this)[v.name()]; }
    X& operator[] (Variable<Real> v) { return (*this)[v.name()]; }
    Void insert(Identifier s, X c) { this->_vals.insert(std::make_pair(s,c)); }
    Void insert(Variable<Real> v, X c) { this->insert(v.name(),c); }
    friend OutputStream& operator<<(OutputStream& os, Valuation<X> const& m) {
        bool first = true; for(auto v : m) { os << (first?"{":",") << v.first << ":" << v.second; first=false; } return os << "}"; }
};

template<class X> class ExpressionInterface;

template<class X> class Expression
    : StaticDispatchNumericOperators<Expression<X>>
{
  private: public:
    std::shared_ptr<ExpressionInterface<X>> _ptr;
  private: public:
    explicit Expression<X>(ExpressionInterface<X>* p);
  public:
    typedef X NumericType;
    typedef Algebra<X> AlgebraType;
    typedef Variable<Real> VariableType;

    template<class SX, EnableIf<IsConvertible<SX,X>> =dummy> Expression<X>(Expression<SX> e);
    template<class SX, EnableIf<IsConvertible<SX,X>> =dummy> Expression<X>(SX c)
        : Expression(NumericType(c)) { }
    explicit operator Algebra<NumericType> () const;
    explicit Expression(Algebra<NumericType> const& a);

    Expression(NumericType);
    Expression(VariableType);
    static Expression constant(NumericType c);
    static Expression variable(VariableType v);
    Expression create_constant(NumericType c) const;

    NumericType operator() (Valuation<NumericType>) const;
    AlgebraType operator() (Valuation<AlgebraType> a) const;
    template<class A> IfAlgebra<A,X> operator() (Valuation<A> const& a);

    friend NumericType evaluate(Expression, Valuation<NumericType>);
    friend AlgebraType evaluate(Expression, Valuation<AlgebraType>);
    friend OutputStream& operator<<(OutputStream& os, Expression<X> const& e);
};

template<class X, class A> IfAlgebra<A,X> evaluate_algebra(Expression<X> e, Valuation<A> x) {
    Valuation<Algebra<X>> ax;
    for(auto xv : x) { ax.insert(xv.first,Algebra<X>(xv.second)); }
    Algebra<X> eax=evaluate(e,ax);
    return eax.template extract<A>();
}

template<class X> template<class A> auto
Expression<X>::operator() (Valuation<A> const& x) -> IfAlgebra<A,X> {
    return evaluate_algebra(*this,x);
}

typedef Variable<Real> RealVariable;
typedef Expression<Real> RealExpression;
typedef Valuation<Real> RealValuation;

typedef Expression<ValidatedReal> ValidatedRealExpression;
typedef Valuation<ValidatedReal> ValidatedRealValuation;

// ---------------- Differential ----------------------------------------------------------------------------------- //

/*
template<class A> struct AlgebraOperations {
    typedef typename A::NumericType X;
    static A _add(A,A); static A _add(A,X); static A _add(X,A);
    static A _mul(A,A); static A _mul(A,X); static A _mul(X,A);
    static A _neg(A); static A _rec(A); static A _exp(A); static A _log(A);
    //static A _sub(A,A); static A _sub(A,X); static A _sub(X,A);
};
*/

template<class X, template<typename>class D> class Differential
    : StaticDispatchAlgbraOperators<Differential<X,D>,GenericType<X>>
{
    typedef GenericType<X> Y;
    typedef IndexType<D<X>> I;
  private: public:
    X _vx; D<X> _dx;
  public:
    typedef I IndexType;
    typedef X NumericType;
    typedef Y GenericNumericType;
    Differential<X,D>(X vx, D<X> dx);
    static Differential<X,D> constant(NumericType c);
    static Differential<X,D> coordinate(NumericType c);
//    static Differential<X,D> coordinate(NumericType c, IndexType i);
    Differential<X,D> create_constant(GenericNumericType c) const;
    explicit operator Algebra<Y> () const;
    explicit Differential<X,D>(Algebra<Y>);
    Differential<X,D>& operator=(Y);
  public:
    friend OutputStream& operator<<(OutputStream& os, Differential<X,D> const& dx) {
        return os << "(" << dx._vx << ";" << dx._dx << ")"; }
};


template<class D> struct IsDifferential : False { };
template<class X, template<typename>class D> struct IsDifferential<Differential<X,D>> : True { };

template<class D, EnableIf<IsDifferential<D>> =dummy> D add(D,D);


// ---------------- Polynomial ------------------------------------------------------------------------------------- //

template<class X> class Polynomial
    : StaticDispatchAlgbraOperators<Polynomial<X>,GenericType<X>>
{
    typedef GenericType<X> Y;
  public:
    typedef Index IndexType;
    typedef X NumericType;
    typedef Y GenericNumericType;
    static Polynomial<X> constant(NumericType c);
    static Polynomial<X> coordinate(NumericType c, IndexType j);
    Polynomial<X> create_constant(GenericNumericType c) const;
    explicit operator Algebra<Y> () const;
    explicit Polynomial(Algebra<Y>);
    Polynomial<X>& operator=(Y);
  public:
    friend OutputStream& operator<<(OutputStream& os, Polynomial<X> const& p);
};

// ---------------- Function --------------------------------------------------------------------------------------- //

class Coordinate {
    SizeType _ind;
  public:
    explicit Coordinate(SizeType ind) : _ind(ind) { }
    SizeType index() const { return _ind; }
};

struct EuclideanSpace {
    SizeType _dim;
  public:
    explicit EuclideanSpace(SizeType dim) : _dim(dim) { }
    SizeType dimension() const { return _dim; }
    friend EuclideanSpace intersection(EuclideanSpace spc1, EuclideanSpace spc2);
};

template<class X> class FormulaInterface;

template<class X> struct Formula {
    std::shared_ptr<FormulaInterface<X>> _ptr;
    Formula(FormulaInterface<X>* p) : _ptr(p) { }
    Formula(std::shared_ptr<FormulaInterface<X>> p) : _ptr(p) { }
    static Formula<X> constant(X const& c);
};

template<class X> class Function
    : DeclareTemplateAlgebraOperators<Function>
    , DefineTemplateMixedSelfOperators<Function>
//    : StaticDispatchAlgbraOperators<Function<X>,X>
{
  public:
    typedef EuclideanSpace DomainType;
    typedef X NumericType;
    typedef Algebra<X> AlgebraType;
    typedef SizeType IndexType;
    typedef Coordinate CoordinateType;
    typedef Formula<X> FormulaType;
  private:
    DomainType _dom;
    FormulaType _frml;
  public:
    explicit Function<X>(DomainType dom, Formula<X> frml);
    template<class SX, EnableIf<IsConvertible<SX,X>> =dummy> Function<X>(Function<SX>);

    explicit operator Algebra<NumericType> () const;
    explicit Function(Algebra<NumericType> const& a);

    Function(Vector<RealVariable> a, Expression<X> e);

    Function(DomainType d, NumericType);
    Function(DomainType d, CoordinateType);
    static Function constant(DomainType d, NumericType c);
    static Function coordinate(DomainType d, IndexType j);
    Function create_constant(NumericType c) const;
    Function create_coordinate(IndexType j) const;
    Function& operator=(NumericType c);

    DomainType domain() const;
    FormulaType formula() const;

    Float64Bounds operator() (Vector<Float64Bounds>) const;
    NumericType operator() (Vector<NumericType>) const;
    AlgebraType operator() (Vector<AlgebraType>) const;
    template<class A> IfAlgebra<A,X> operator() (Vector<A> const& a) const;

    friend NumericType evaluate(Function, Vector<NumericType>);
    friend AlgebraType evaluate(Function, Vector<AlgebraType>);
    friend Float64Bounds evaluate(Function, Vector<Float64Bounds>);

    friend OutputStream& operator<<(OutputStream& os, Function<X> const& e);
};

template<class X, class A> IfAlgebra<A,X> evaluate_algebra(Function<X> const& f, Vector<A> const& x) {
    std::vector<Algebra<X>> va; for(SizeType i=0; i!=x.size(); ++i) { va.push_back(Algebra<X>(x[i])); } Vector<Algebra<X>> ax(va);
    Algebra<X> fax=f(ax); return fax.template extract<A>();
}

template<class X> template<class A> IfAlgebra<A,X> Function<X>::operator() (Vector<A> const& x) const {
    return evaluate_algebra(*this,x);
}

typedef Function<Real> RealFunction;
typedef Function<ValidatedReal> ValidatedRealFunction;

Real evaluate(RealFunction, Vector<Real>);
ValidatedReal evaluate(ValidatedRealFunction, Vector<ValidatedReal>);

