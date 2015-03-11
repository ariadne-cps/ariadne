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
template<class F, class T> using IsConvertible = std::is_convertible<F,T>;
template<class T, class... F> using IsConstructible = std::is_constructible<T,F...>;
template<class T1, class T2> using IsSame = std::is_same<T1,T2>;
template<class X> using InitializerList = std::initializer_list<X>;
using Void = void;
using SizeType = std::size_t;

using OutputStream = std::ostream;

template<class T> struct Self { typedef T Type; };
template<class T> using SelfType = typename Self<T>::Type;
template<class T> using GenericType = typename T::GenericType;
//template<class T> using NumericType = typename T::NumericType;

template<class T> struct NumericTypedef { typedef typename T::NumericType Type; };
template<class T> using NumericType = typename NumericTypedef<T>::Type;

struct Cnst { };
struct Var { };

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

template<class T1, class T2=Void> struct ArithmeticTraits {
    typedef decltype(add(declval<T1>(),declval<T2>())) AddType;
    typedef decltype(mul(declval<T1>(),declval<T2>())) MulType;
};

template<class T> struct ArithmeticTraits<T,Void> {
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

struct Effort { };

// ---------------- Integer --------------------------------------------------------------------------------------- //

class Integer {
    int _n;
  public:
    Integer(double) = delete;
    Integer(int n);
    int get_si() const { return _n; }
    friend Integer add(Integer,Integer);
    friend Integer mul(Integer,Integer);
    friend Integer neg(Integer);
    friend Integer abs(Integer);
    friend Integer quot(Integer,Integer);
    friend Integer rem(Integer,Integer);
    friend Integer gcd(Integer,Integer);
    friend OutputStream& operator<<(OutputStream&,Integer const&);
};

// ---------------- Rational --------------------------------------------------------------------------------------- //

class Rational {
    Integer _num; Integer _den;
    Rational(double) = delete;
 public:
    Rational(int n) : Rational(Integer(n)) { }
    Rational(Integer n, Integer d=1);
    double get_d() const;
    friend Rational add(Rational,Rational);
    friend Rational mul(Rational,Rational);
    friend Rational neg(Rational);
    friend Rational rec(Rational);
    friend OutputStream& operator<<(OutputStream&,Rational const&);
};

// ---------------- Real ------------------------------------------------------------------------------------------- //

template<class F> class Bounds;
class Precision64;
class Float64;
typedef Bounds<Float64> Float64Bounds;

class RealInterface;
class Real {
    Real(double) = delete;
    std::shared_ptr<RealInterface> _ptr;
  public:
    Real(RealInterface* p);
  public:
    Real(int n);
    Real(Integer z);
    Real(Rational q);
    Rational lower_bound(Effort e) const;
    Rational upper_bound(Effort e) const;
    Float64Bounds get(Precision64) const;

    friend Real add(Real,Real);
    friend Real mul(Real,Real);

    friend Real neg(Real);
    friend Real rec(Real);

    friend Real exp(Real);
    friend Real log(Real);
    friend OutputStream& operator<<(OutputStream&,Real const&);

    // FIXME: Should not be necessary...
    Real create_constant(Real);
};

class ValidatedReal {
    Real _r;
  public:
    template<class T, EnableIf<IsConvertible<T,Real>> =dummy> ValidatedReal(T t) : ValidatedReal(Real(t)) { }
    ValidatedReal(Real);
    Float64Bounds get(Precision64) const;

    friend ValidatedReal add(ValidatedReal,ValidatedReal);
    friend ValidatedReal mul(ValidatedReal,ValidatedReal);

    friend ValidatedReal neg(ValidatedReal);
    friend ValidatedReal rec(ValidatedReal);

    friend ValidatedReal exp(ValidatedReal);
    friend ValidatedReal log(ValidatedReal);
    friend OutputStream& operator<<(OutputStream&,ValidatedReal const&);
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

template<class F> class Bounds {
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

    F lower_bound() const;
    F upper_bound() const;
    PrecisionType precision() const;
    friend OutputStream& operator<<(OutputStream&,Bounds<F> const&);

#ifndef TEMPLATE_FLOAT
    friend Bounds<F> add(Bounds<F>,Bounds<F>);
    friend Bounds<F> mul(Bounds<F>,Bounds<F>);

    friend Bounds<F> neg(Bounds<F>);
    friend Bounds<F> rec(Bounds<F>);

    friend Bounds<F> exp(Bounds<F>);
    friend Bounds<F> log(Bounds<F>);
#endif
};

#ifdef TEMPLATE_FLOAT
template<class F> Bounds<F> add(Bounds<F>,Bounds<F>);
template<class F> Bounds<F> mul(Bounds<F>,Bounds<F>);
template<class F> Bounds<F> neg(Bounds<F>);
template<class F> Bounds<F> rec(Bounds<F>);
template<class F> Bounds<F> exp(Bounds<F>);
template<class F> Bounds<F> log(Bounds<F>);
#endif

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

template<class A, class X=NumericType<A>> class AlgebraWrapper
    : public virtual AlgebraInterface<X>, public A
{
    friend A; friend class Algebra<X>;
  public:
    AlgebraWrapper(A&& a) : A(std::move(a)) { }
    AlgebraWrapper(A const& a) : A(a) { }
    virtual ~AlgebraWrapper<A,X> () = default;
  private:
    //virtual AlgebraInterface<X>* _create_constant(X c) const { AlgebraWrapper<A,X> r(*this); r=c; return r; }
    virtual AlgebraInterface<X>* _create_constant(X c) const { return new AlgebraWrapper<A,X>(this->create_constant(c)); }
    virtual AlgebraInterface<X>* _add(AlgebraInterface<X> const* other) const { return this->_apply(Add(),other); }
    virtual AlgebraInterface<X>* _mul(AlgebraInterface<X> const* other) const { return this->_apply(Mul(),other); }
    virtual AlgebraInterface<X>* _neg() const { return new AlgebraWrapper<A,X>(neg(*this)); }
    virtual AlgebraInterface<X>* _rec() const { return new AlgebraWrapper<A,X>(rec(*this)); }
    virtual AlgebraInterface<X>* _exp() const { return new AlgebraWrapper<A,X>(exp(*this)); }
    virtual AlgebraInterface<X>* _log() const { return new AlgebraWrapper<A,X>(log(*this)); }
    virtual AlgebraInterface<X>* _add(X const& c) const { return new AlgebraWrapper<A,X>(add(*this,c)); }
    virtual AlgebraInterface<X>* _radd(X const& c) const { return new AlgebraWrapper<A,X>(add(c,*this)); }
    virtual AlgebraInterface<X>* _mul(X const& c) const { return new AlgebraWrapper<A,X>(mul(*this,c)); }
    virtual AlgebraInterface<X>* _rmul(X const& c) const { return new AlgebraWrapper<A,X>(mul(c,*this)); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<A const&>(*this); }
  private:
    template<class OP> AlgebraInterface<X>* _apply(OP op, AlgebraInterface<X> const* ai2) const;
    template<class OP> AlgebraInterface<X>* _apply(OP op) const;
};

template<class A, class X> template<class OP> AlgebraInterface<X>*
AlgebraWrapper<A,X>::_apply(OP op, AlgebraInterface<X> const* ai2) const {
    AlgebraWrapper<A,X> const* aw1=this; AlgebraWrapper<A,X> const* aw2=dynamic_cast<AlgebraWrapper<A,X>const*>(ai2);
    if(not aw2) {
        std::cerr << "aw1="; aw1->_write(std::cerr);
        std::cerr << ", ai2="; ai2->_write(std::cerr);std::cerr << "\n";
        throw std::bad_cast(); }
    return new AlgebraWrapper<A,X>(op(static_cast<A const&>(*aw1),static_cast<A const&>(*aw2)));
}

template<class A, class X> template<class OP>
AlgebraInterface<X>* AlgebraWrapper<A,X>::_apply(OP op) const {
       return new AlgebraWrapper<A,X>(op(static_cast<A const&>(*this)));
}

template<class X> class Algebra
// : public Field
{
    std::shared_ptr<AlgebraInterface<X>> _ptr;
  private:
    explicit Algebra(AlgebraInterface<X>* p) : _ptr(p) { }
  public:
    typedef Algebra<X> SelfType;
    typedef X NumericType;
    template<class XX> friend Algebra<XX> make_handle(AlgebraInterface<XX>*);
    template<class A> A const& extract() const { return dynamic_cast<AlgebraWrapper<A,X>const&>(*this->_ptr); }
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

// ---------------- Expression ------------------------------------------------------------------------------------- //

class Identifier : public std::string {
  public:
    Identifier(const char* s) : std::string(s) { }
    Identifier(std::string s) : std::string(s) { }
};

template<class X> class Expression;
template<class X> class Valuation;

template<class X> class Variable {
    Identifier _name;
  public:
    explicit Variable(std::string name) : _name(name) { }
    Identifier name() const { return _name; }
    template<class Y> Y create_constant(Y const& c) { return Expression<Y>(c); }
    friend bool operator<(Variable<X> v1, Variable<X> v2) { return v1._name < v2._name; }
    friend OutputStream& operator<<(OutputStream& os, Variable<X> const& v) { return os << v._name; }
  private:
    friend Expression<X> add(Expression<X>,Expression<X>);
    friend Expression<X> mul(Expression<X>,Expression<X>);
    friend Expression<X> neg(Expression<X>);
    friend Expression<X> rec(Expression<X>);
    friend Expression<X> exp(Expression<X>);
    friend Expression<X> log(Expression<X>);
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

template<class X> class Expression {
    std::shared_ptr<ExpressionInterface<X>> _ptr;
  private: public:
    explicit Expression<X>(ExpressionInterface<X>* p);
  public:
    typedef X NumericType;
    typedef Algebra<X> AlgebraType;
    typedef Variable<X> VariableType;

//     template<class SX, EnableIf<IsConvertible<SX,X>> =dummy> Expression<X>(Expression<SX>);
    template<class SX, EnableIf<IsConvertible<SX,X>> =dummy> Expression<X>(SX c) : Expression(NumericType(c)) { }
    explicit operator Algebra<NumericType> () const;
    explicit Expression<X>(Algebra<X> const& a);

    Expression(NumericType);
    Expression(VariableType);
    static Expression constant(NumericType c);
    static Expression variable(VariableType v);
    Expression create_constant(NumericType c) const;

    friend Expression add(Expression,Expression);
    friend Expression mul(Expression,Expression);

    friend Expression neg(Expression);
    friend Expression rec(Expression);
    friend Expression exp(Expression);
    friend Expression log(Expression);

    NumericType operator() (Valuation<NumericType>) const;
    AlgebraType operator() (Valuation<AlgebraType> a) const;
    template<class A> IfAlgebra<A,X> operator() (Valuation<A> const& a);

    friend NumericType evaluate(Expression, Valuation<NumericType>);
    friend AlgebraType evaluate(Expression, Valuation<AlgebraType>);
    template<class A> friend IfAlgebra<A,NumericType>  evaluate_algebra(Expression, Valuation<A>);
    friend OutputStream& operator<<(OutputStream& os, Expression<X> const& e);
};

template<class X> template<class A> auto
Expression<X>::operator() (Valuation<A> const& x) -> IfAlgebra<A,X> {
    Valuation<Algebra<X>> ax;
    for(auto xv : x) { ax.insert(xv.first,Algebra<X>(xv.second)); }
    Algebra<X> eax=evaluate(*this,ax);
    return eax.template extract<A>();
}

template<class A> IfAlgebra<A,Real> evaluate_algebra(Expression<Real> e, Valuation<A> x) {
    return e(x); }

typedef Variable<Real> RealVariable;
typedef Expression<Real> RealExpression;
typedef Valuation<Real> RealValuation;

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
    X const& operator[] (SizeType i) const { return _ary[i]; }
    X& operator[] (SizeType i) { return _ary[i]; }
    friend OutputStream& operator<<(OutputStream& os, Vector<X> const& v) {
        for(SizeType i=0; i!=v.size(); ++i) { os << (i==0?"[":",") << v[i]; } return os << "]"; }
};
typedef Vector<Real> RealVector;
typedef Vector<ValidatedReal> ValidatedRealVector;


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
    Formula<X> constant(X const& c);
    friend X evaluate(Formula<X>, Vector<X>);
};

template<class X> class Function {
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

    DomainType domain() const;
    FormulaType formula() const;

    friend Function add(Function,Function);
    friend Function mul(Function,Function);

    friend Function neg(Function);
    friend Function rec(Function);
    friend Function exp(Function);
    friend Function log(Function);

    friend Function add(Function,NumericType);
    friend Function add(NumericType,Function);
    friend Function mul(Function,NumericType);
    friend Function mul(NumericType,Function);

    NumericType operator() (Vector<NumericType>) const;
    AlgebraType operator() (Vector<AlgebraType>) const;
//    template<class A> IfAlgebra<A,X> operator() (Vector<A> const& a);

    friend NumericType evaluate(Function, Vector<NumericType>);
    friend NumericType evaluate(Function, Vector<AlgebraType>);

    friend OutputStream& operator<<(OutputStream& os, Function<X> const& e);
};

template<class X, class A> IfAlgebra<A,X> evaluate_algebra(Function<X> const& f, Vector<A> const& x) {
    std::vector<Algebra<X>> va; for(SizeType i=0; i!=x.size(); ++i) { va.push_back(Algebra<X>(x[i])); } Vector<Algebra<X>> ax(va);
    Algebra<X> fax=f(ax); return fax.template extract<A>(); }


typedef Function<Real> RealFunction;
typedef Function<ValidatedReal> ValidatedRealFunction;
