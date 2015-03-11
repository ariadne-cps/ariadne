#include "simple.h"

// ---------------- Misc --------------------------------------------------------------------------------------- //

template<class K, class V> OutputStream& operator<<(OutputStream& os, std::map<K,V> const& m) {
    bool first = true; for(auto v : m) { os << (first?"{":",") << v.first << ":" << v.second; first=false; } return os << "}"; }

// ---------------- Symbolic --------------------------------------------------------------------------------------- //

template<class OP, class... AS> struct Symbolic;

template<class X> struct Symbolic<Cnst,X> {
    X _val;
};

template<class I> struct Symbolic<Var,I> {
    I _ind;
};

template<class OP, class A> struct Symbolic<OP,A> {
    A _arg;
};

template<class OP, class A1, class A2> struct Symbolic<OP,A1,A2> {
    A1 _arg1; A1 _arg2;
};

template<class OP, class... AS> Symbolic<OP,AS...> make_symbolic(OP op, AS... as) {
    return Symbolic<OP,AS...>{as...}; }

template<class X, class C, class V> X evaluate(Symbolic<Cnst,C> const& e, V const& v) {
    return v.create_constant(e._val); }
template<class X, class I, class V> X evaluate(Symbolic<Var,I> const& e, V const& v) {
    return v[e._ind]; }
template<class X, class OP, class A, class V> X evaluate(Symbolic<OP,A> const& e, V const& v) {
    OP _op; return _op(evaluate(e._arg,v)); }
template<class X, class OP, class A1, class A2, class V> X evaluate(Symbolic<OP,A1,A2> const& e, V const& v) {
    OP _op; return _op(evaluate(e._arg1,v),evaluate(e._arg2,v)); }

template<class X> OutputStream& operator<<(OutputStream& os, Symbolic<Cnst,X> e) {
    return os << e._val; }
OutputStream& operator<<(OutputStream& os, Symbolic<Var,SizeType> e) {
    return os << "[" << e._ind << "]"; }
OutputStream& operator<<(OutputStream& os, Symbolic<Var,Identifier> e) {
    return os << e._ind; }
template<class OP, class A> OutputStream& operator<<(OutputStream& os, Symbolic<OP,A> e) {
    return os << OP() << "(" << e._arg << ")"; }
template<class OP, class A1, class A2> OutputStream& operator<<(OutputStream& os, Symbolic<OP,A1,A2> e) {
    return os << OP() << "(" << e._arg1 << "," << e._arg2 << ")"; }

OutputStream& operator<<(OutputStream& os, Add op) { return os << "add"; }
OutputStream& operator<<(OutputStream& os, Mul op) { return os << "mul"; }
OutputStream& operator<<(OutputStream& os, Neg op) { return os << "neg"; }
OutputStream& operator<<(OutputStream& os, Rec op) { return os << "rec"; }
OutputStream& operator<<(OutputStream& os, Exp op) { return os << "exp"; }
OutputStream& operator<<(OutputStream& os, Log op) { return os << "log"; }

// ---------------- Integer --------------------------------------------------------------------------------------- //

Integer::Integer(int n) : _n(n) { }
Integer add(Integer z1, Integer z2) { return Integer(z1._n+z2._n); }
Integer mul(Integer z1, Integer z2) { return Integer(z1._n*z2._n); }
Integer neg(Integer z) { return Integer(-z._n); }
Integer abs(Integer z) { return Integer(std::max(z._n,-z._n)); }
Integer quot(Integer z1, Integer z2) { return Integer(z1._n/z2._n); }
Integer rem(Integer z1, Integer z2) { return Integer(z1._n%z2._n); }
OutputStream& operator<<(OutputStream& os, Integer const& z) { return os << z._n; }

Integer gcd_pos(Integer z1, Integer z2) {
    if(z2.get_si()==0) { return z1; }
    return gcd_pos(z2,rem(z1,z2));
}
Integer gcd(Integer z1, Integer z2) {
    Integer az1=abs(z1); Integer az2=abs(z2);
    return (az1.get_si()>=az2.get_si()) ? gcd_pos(az1,az2) : gcd_pos(az2,az1);
}

// ---------------- Rational --------------------------------------------------------------------------------------- //

double Rational::get_d() const {
    return double(_num.get_si()) / _den.get_si();
}

Rational::Rational(Integer num, Integer den) : _num(num), _den(den) {
    if(den.get_si()<0) { _num=neg(_num); _den=neg(_den); }
    Integer _gcd=gcd(num,den); _num=quot(_num,_gcd); _den=quot(_den,_gcd);
}

Rational add(Rational q1, Rational q2) {
    return Rational(add(mul(q1._num,q2._den),mul(q1._den,q2._num)),mul(q1._den,q2._den));
}

Rational mul(Rational q1, Rational q2) {
    return Rational(mul(q1._num,q1._num),mul(q1._den,q2._den));
}

Rational neg(Rational q) {
    return Rational(neg(q._num),q._den);
}

Rational rec(Rational q) {
    return Rational(q._den,q._num);
}

OutputStream& operator<<(OutputStream& os, Rational const& q) {
    return os << q._num << "/" << q._den;
}

// ---------------- Real ------------------------------------------------------------------------------------------- //

class RealInterface {
  public:
    virtual ~RealInterface() = default;
    virtual Float64Bounds _get(Precision64 pr) const = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

template<class OP, class... AS> struct RealWrapper : public virtual RealInterface, public Symbolic<OP,AS...> {
    RealWrapper(AS... as) : Symbolic<OP,AS...>{as...} { }
    virtual Float64Bounds _get(Precision64 pr) const { return evaluate<Float64Bounds> (this->_arg1.get(pr),this->_arg2.get(pr)); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<Symbolic<OP,AS...>const&>(*this); }
};


template<class X> struct RealWrapper<Cnst,X> : public virtual RealInterface, public Symbolic<Cnst,X> {
    RealWrapper(X const& c) : Symbolic<Cnst,X>{c} { }
    Float64Bounds _get(Precision64 pr) const { return Float64Bounds(this->_val,pr); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<Symbolic<Cnst,X>const&>(*this); }
};
template<class OP> struct RealWrapper<OP,Real> : public virtual RealInterface, public Symbolic<OP,Real> {
    RealWrapper(Real a) : Symbolic<OP,Real>{a} { }
    virtual Float64Bounds _get(Precision64 pr) const { OP _op; return _op(this->_arg.get(pr)); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<Symbolic<OP,Real>const&>(*this); }
};
template<class OP> struct RealWrapper<OP,Real,Real> : public virtual RealInterface, public Symbolic<OP,Real,Real> {
    RealWrapper(Real a1, Real a2) : Symbolic<OP,Real,Real>{a1,a2} { }
    virtual Float64Bounds _get(Precision64 pr) const { OP _op; return _op(this->_arg1.get(pr),this->_arg2.get(pr)); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<Symbolic<OP,Real,Real>const&>(*this); }
};

template<class OP, class... AS> Real make_real(AS... as) { return Real(new RealWrapper<OP,AS...>{as...}); }
template<> Real make_real<Cnst,Float64Bounds>(Float64Bounds c) { return Real(new RealWrapper<Cnst,Float64Bounds>{c}); }

Real::Real(RealInterface* p) : _ptr(p) { }
Real::Real(int n) : Real(Integer(n)) { }
Real::Real(Integer z) : _ptr(new RealWrapper<Cnst,Integer>{z}) { }
Real::Real(Rational q) : _ptr(new RealWrapper<Cnst,Rational>{q}) { }
Float64Bounds Real::get(Precision64 pr) const { return this->_ptr->_get(pr); }
Real add(Real r1, Real r2) { return make_real<Add>(r1,r2); }
Real mul(Real r1, Real r2) { return make_real<Mul>(r1,r2); }
Real neg(Real r) { return make_real<Neg>(r); }
Real rec(Real r) { return make_real<Rec>(r); }
Real exp(Real r) { return make_real<Exp>(r); }
Real log(Real r) { return make_real<Log>(r); }

//OutputStream& operator<<(OutputStream& os, Real const& r) { os << "Real("; r._ptr->_write(os); return os << ")"; }
OutputStream& operator<<(OutputStream& os, Real const& r) { r._ptr->_write(os); return os; }

ValidatedReal::ValidatedReal(Real r) : _r(r) { }
Float64Bounds ValidatedReal::get(Precision64 pr) const { return this->_r.get(pr); }
OutputStream& operator<<(OutputStream& os, ValidatedReal const& vr) { return os << "Validated" << vr._r; }


// ---------------- Algebra ------------------------------------------------------------------------------------- //

// ---------------- Expression ------------------------------------------------------------------------------------- //

template<class X> class ExpressionInterface {
  public:
    ~ExpressionInterface<X> () = default;
    virtual X _evaluate(Valuation<X> const&) = 0;
    virtual Algebra<X> _evaluate(Valuation<Algebra<X>> const&) = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

template<class X> Float64Bounds evaluate(Symbolic<Cnst,X> const& e, Precision64 pr) { return Float64Bounds(e._val,pr); }

template<class X, class OP, class... AS> struct ExpressionWrapper
    : virtual ExpressionInterface<X>, Symbolic<OP,AS...>
{
    ExpressionWrapper(AS... as)
        : Symbolic<OP,AS...>{as...} { }
    virtual X _evaluate(Valuation<X> const& v) {
        return evaluate<X>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual Algebra<X> _evaluate(Valuation<Algebra<X>> const& v) {
        return evaluate<Algebra<X>>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual OutputStream& _write(OutputStream& os) const {
        return os << static_cast<Symbolic<OP,AS...>const&>(*this); }
};

template<class X, class OP, class... AS> Expression<X> make_expression(AS... as) {
    return Expression<X>(new ExpressionWrapper<X,OP,AS...>(as...)); }

template<class X> Expression<X>::Expression(ExpressionInterface<X>* p) : _ptr(p) { }
template<class X> Expression<X>::operator Algebra<NumericType> () const { return make_handle(new AlgebraWrapper<Expression<X>,X>(*this)); }
template<class X> Expression<X>::Expression(Algebra<X> const& a) : Expression(a.template extract<Expression<X>>()) { }

template<class X> Expression<X>::Expression(X c) : Expression(new ExpressionWrapper<X,Cnst,X>{c}) { }
template<class X> Expression<X>::Expression(Variable<X> v) : Expression(new ExpressionWrapper<X,Var,Identifier>{v.name()}) { }
template<class X> Expression<X> Expression<X>::create_constant(X c) const { return Expression<X>(c); }

template class Expression<Real>;

template<class X> Expression<X> add(Expression<X> e1,Expression<X> e2) { return make_expression<X,Add>(e1,e2); }
template<class X> Expression<X> mul(Expression<X> e1,Expression<X> e2) { return make_expression<X,Mul>(e1,e2); }
template<class X> Expression<X> neg(Expression<X> e) { return make_expression<X,Neg>(e); }
template<class X> Expression<X> rec(Expression<X> e) { return make_expression<X,Rec>(e); }
template<class X> Expression<X> exp(Expression<X> e) { return make_expression<X,Exp>(e); }
template<class X> Expression<X> log(Expression<X> e) { return make_expression<X,Log>(e); }
template Expression<Real> exp(Expression<Real>);

RealExpression add(RealExpression e1, RealExpression e2) { return make_expression<Real,Add>(e1,e2); }
RealExpression mul(RealExpression e1, RealExpression e2) { return make_expression<Real,Mul>(e1,e2); }
RealExpression neg(RealExpression e) { return make_expression<Real,Neg>(e); }
RealExpression rec(RealExpression e) { return make_expression<Real,Rec>(e); }
RealExpression exp(RealExpression e) { return make_expression<Real,Exp>(e); }
RealExpression log(RealExpression e) { return make_expression<Real,Log>(e); }

template<class X> X Expression<X>::operator()(Valuation<X> v) const { return evaluate(*this,v); }
Real evaluate(RealExpression e, RealValuation v) { return e._ptr->_evaluate(v); }
RealAlgebra evaluate(RealExpression e, Valuation<RealAlgebra> v) { return e._ptr->_evaluate(v); }
//OutputStream& operator<<(OutputStream& os, RealExpression const& e) { os << "RealExpression("; e._ptr->_write(os); return os << ")"; }
OutputStream& operator<<(OutputStream& os, RealExpression const& e) { e._ptr->_write(os); return os; }

// ---------------- Domain ----------------------------------------------------------------------------------------- //

EuclideanSpace intersection(EuclideanSpace spc1, EuclideanSpace spc2) {
    if(spc1.dimension() != spc2.dimension()) { throw std::runtime_error("Dimension mismatch"); }
    return spc1;
}

// ---------------- Function --------------------------------------------------------------------------------------- //

template<class X> class FormulaInterface {
  public:
    ~FormulaInterface () = default;
    virtual X _evaluate(Vector<X> const&) = 0;
    virtual Algebra<X> _evaluate(Vector<Algebra<X>> const&) = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
};
template<class X> OutputStream& operator<<(OutputStream& os, Formula<X> const& f) { f._ptr->_write(os); return os; }
template<class X> X evaluate(Formula<X> const& f, Vector<X> const& v) { return f._ptr->_evaluate(v); }
template<class X> Algebra<X> evaluate(Formula<X> const& f, Vector<Algebra<X>> const& v) { return f._ptr->_evaluate(v); }

template<class X, class OP, class... AS> struct FormulaWrapper
    : virtual FormulaInterface<X>, Symbolic<OP,AS...>
{
    FormulaWrapper(AS... as) : Symbolic<OP,AS...>{as...} { }
    virtual X _evaluate(Vector<X> const& v) { return evaluate<X>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual Algebra<X> _evaluate(Vector<Algebra<X>> const& v) { return evaluate<Algebra<X>>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<Symbolic<OP,AS...>const&>(*this); }
};

template<class X, class OP, class... AS> Function<X> make_function(typename Function<X>::DomainType dom, AS... as) {
    return Function<X>(dom, Formula<X>(new FormulaWrapper<X,OP,AS...>(as...))); }

template<class X> Function<X>::Function(EuclideanSpace dom, Formula<X> frml) : _dom(dom), _frml(frml) { }
template<class X> Function<X>::operator Algebra<NumericType> () const { return make_handle(new AlgebraWrapper<Function<X>,X>(*this)); }
template<class X> Function<X>::Function(Algebra<X> const& a) : Function(a.template extract<Function<X>>()) { }

template<class X> Function<X>::Function(DomainType dom, NumericType c) : Function(dom,new FormulaWrapper<X,Cnst,NumericType>{c}) { }
template<class X> Function<X>::Function(DomainType dom, CoordinateType j) : Function(dom,new FormulaWrapper<X,Var,IndexType>{j.index()}) { }
template<class X> Function<X> Function<X>::create_constant(NumericType c) const { return Function<X>(this->domain(),c); }
template<class X> Function<X> Function<X>::create_coordinate(IndexType j) const { return Function<X>(this->domain(),Coordinate(j)); }
template<class X> Function<X> Function<X>::constant(DomainType dom, NumericType c) { return Function<X>(dom,c); }
template<class X> Function<X> Function<X>::coordinate(DomainType dom, IndexType j) { return Function<X>(dom,Coordinate(j)); }

template<class X> auto Function<X>::domain() const -> DomainType { return this->_dom; }
template<class X> auto Function<X>::formula() const -> FormulaType { return this->_frml; }

template<class X> Function<X>::Function(Vector<RealVariable> a, Expression<X> e) : Function(make_function(a,e)) { }

template<class X> auto Function<X>::operator() (Vector<AlgebraType> v) const -> AlgebraType {
    return evaluate(this->formula(),v);
}

static_assert(IsSame<IfAlgebra<Function<Real>,Real>,Function<Real>>::value,"");

template<class X> Function<X> make_function(Vector<RealVariable> a, Expression<X> e) {
    EuclideanSpace dom(a.size());
    Valuation<Function<X>> v;
    for(SizeType i=0; i!=a.size(); ++i) { v.insert(a[i],Function<X>::coordinate(dom,i)); }
    std::cout << "val="<<v.map() << "\n";
    return evaluate_algebra(e,v);
}

template class Function<Real>;
typedef Formula<Real> RealFormula;

template<class X> Function<X> add(Function<X> f1,Function<X> f2) {
    return make_function<X,Add>(intersection(f1.domain(),f2.domain()),f1.formula(),f2.formula()); }
template<class X> Function<X> mul(Function<X> f1,Function<X> f2) {
    return make_function<X,Mul>(intersection(f1.domain(),f2.domain()),f1.formula(),f2.formula()); }
template<class X> Function<X> neg(Function<X> f) { return make_function<X,Neg>(f.domain(),f.formula()); }
template<class X> Function<X> rec(Function<X> f) { return make_function<X,Rec>(f.domain(),f.formula()); }
template<class X> Function<X> exp(Function<X> f) { return make_function<X,Exp>(f.domain(),f.formula()); }
template<class X> Function<X> log(Function<X> f) { return make_function<X,Log>(f.domain(),f.formula()); }
template Function<Real> exp(Function<Real>);

RealFunction add(RealFunction f1, RealFunction f2) {
    return make_function<Real,Add>(intersection(f1.domain(),f2.domain()),f1.formula(),f2.formula()); }
RealFunction mul(RealFunction f1, RealFunction f2) {
    return make_function<Real,Mul>(intersection(f1.domain(),f2.domain()),f1.formula(),f2.formula()); }
RealFunction neg(RealFunction f) { return make_function<Real,Neg>(f.domain(),f.formula()); }
RealFunction rec(RealFunction f) { return make_function<Real,Rec>(f.domain(),f.formula()); }
RealFunction exp(RealFunction f) { return make_function<Real,Exp>(f.domain(),f.formula()); }
RealFunction log(RealFunction f) { return make_function<Real,Log>(f.domain(),f.formula()); }

RealFunction add(RealFunction f1, Real c2) { return add(f1,f1.create_constant(c2)); }
RealFunction add(Real c1, RealFunction f2) { return add(f2.create_constant(c1),f2); }
RealFunction mul(RealFunction f1, Real c2) { return mul(f1,f1.create_constant(c2)); }
RealFunction mul(Real c1, RealFunction f2) { return mul(f2.create_constant(c1),f2); }


Real evaluate(RealFunction f, RealVector v) { return Real(0); }
Real evaluate(RealFormula f, RealVector v) { return Real(0); }
//OutputStream& operator<<(OutputStream& os, RealFunction const& e) { os << "RealFunction("; e._ptr->_write(os); return os << ")"; }
OutputStream& operator<<(OutputStream& os, RealFunction const& f) { f._frml._ptr->_write(os); return os; }
