#include "simple.h"
#include <type_traits>
#include <tuple>
#include <cassert>

// ---------------- Misc --------------------------------------------------------------------------------------- //

template<SizeType POS> struct SizeTag { };

template <class TUP, SizeType POS>
OutputStream& print_tuple(OutputStream& os, const TUP& t, SizeTag<POS> ) {
    print_tuple(os, t, SizeTag<POS-1>());
    os << ',' << std::get<POS-1>(t);
    return os; }

template <class TUP>
OutputStream& print_tuple(OutputStream& os, const TUP& t, SizeTag<1> ) {
    return os << std::get<0>(t); }

template <class... AS>
OutputStream& operator<<(OutputStream& os, const std::tuple<AS...>& t) {
    os << '('; print_tuple(os, t, SizeTag<sizeof...(AS)>()); return os << ')'; }

template<class T> OutputStream& operator<<(OutputStream& os, std::vector<T> const& v) {
    bool first = true; for(auto t : v) { os << (first?"[":",") << t; first=false; } return os << "]"; }

template<class K, class V> OutputStream& operator<<(OutputStream& os, std::map<K,V> const& m) {
    bool first = true; for(auto kv : m) { os << (first?"{":",") << kv.first << ":" << kv.second; first=false; } return os << "}"; }

bool initialize() { return true; }
const bool initialized = initialize();

// ---------------- Symbolic --------------------------------------------------------------------------------------- //

template<class OP, class... AS> struct Symbolic;

template<class X> struct Symbolic<Cnst,X> { X _val; };

template<class I> struct Symbolic<Var,I> { I _ind; };

template<class OP, class A> struct Symbolic<OP,A> { A _arg; };

template<class OP, class A1, class A2> struct Symbolic<OP,A1,A2> { A1 _arg1; A2 _arg2; };

template<class OP, class A1, class A2, class A3> struct Symbolic<OP,A1,A2,A3> { A1 _arg1; A2 _arg2; A3 _arg3; };

template<class OP, class... AS> Symbolic<OP,AS...> make_symbolic(OP op, AS... as) { return Symbolic<OP,AS...>{as...}; }

template<class X, class C, class V> X evaluate_as(Symbolic<Cnst,C> const& e, V const& v) {
    return v.create_constant(e._val); }
template<class X, class I, class V> X evaluate_as(Symbolic<Var,I> const& e, V const& v) {
    return v[e._ind]; }
template<class X, class OP, class A, class V> X evaluate_as(Symbolic<OP,A> const& e, V const& v) {
    OP _op; return _op(evaluate(e._arg,v)); }
template<class X, class OP, class A1, class A2, class V> X evaluate_as(Symbolic<OP,A1,A2> const& e, V const& v) {
    OP _op; return _op(evaluate(e._arg1,v),evaluate(e._arg2,v)); }
template<class X, class OP, class A1, class A2, class A3, class V> X evaluate_as(Symbolic<OP,A1,A2,A3> const& e, V const& v) {
    OP _op; return _op(evaluate(e._arg1,v),evaluate(e._arg2,v),evaluate(e._arg3,v)); }

template<class X> OutputStream& operator<<(OutputStream& os, Symbolic<Cnst,X> const& e) {
    return os << e._val; }
OutputStream& operator<<(OutputStream& os, Symbolic<Var,SizeType> const& e) {
    return os << "[" << e._ind << "]"; }
OutputStream& operator<<(OutputStream& os, Symbolic<Var,Identifier> const& e) {
    return os << e._ind; }
template<class OP, class A> OutputStream& operator<<(OutputStream& os, Symbolic<OP,A> const& e) {
    return os << OP() << "(" << e._arg << ")"; }
template<class OP, class A1, class A2> OutputStream& operator<<(OutputStream& os, Symbolic<OP,A1,A2> const& e) {
    return os << OP() << "(" << e._arg1 << "," << e._arg2 << ")"; }

OutputStream& operator<<(OutputStream& os, Cnst op) { return os << "cnst"; }
OutputStream& operator<<(OutputStream& os, Var op) { return os << "var"; }
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

template<class X, class C> X evaluate_as(Symbolic<Cnst,C> const& e, Precision64 const& pr) {
    return Float64Bounds(e._val,pr); }

class RealInterface {
  public:
    virtual ~RealInterface() = default;
    virtual Float64Bounds _get(Precision64 pr) const = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

template<class X> Float64Bounds evaluate(Symbolic<Cnst,X> const& e, Precision64 pr) { return Float64Bounds(e._val,pr); }

template<class OP, class... AS> struct RealWrapper : public virtual RealInterface, public Symbolic<OP,AS...> {
    RealWrapper(AS... as) : Symbolic<OP,AS...>{as...} { }
    virtual Float64Bounds _get(Precision64 pr) const { return evaluate_as<Float64Bounds>(static_cast<Symbolic<OP,AS...>const&>(*this),pr); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<Symbolic<OP,AS...>const&>(*this); }
};

template<class OP, class... AS> Real make_real(AS... as) { return Real(new RealWrapper<OP,AS...>{as...}); }


inline Float64Bounds evaluate(Real const& r, Precision64 pr) { return r._ptr->_get(pr); }

Real::Real(RealInterface* p) : _ptr(p) { }
Real::Real(int n) : Real(Integer(n)) { }
Real::Real(Integer z) : _ptr(new RealWrapper<Cnst,Integer>{z}) { }
Real::Real(Rational q) : _ptr(new RealWrapper<Cnst,Rational>{q}) { }
Real Real::create_constant(Real r) const { return r; }
Float64Bounds Real::get(Precision64 pr) const { return this->_ptr->_get(pr); }
Real add(Real r1, Real r2) { return make_real<Add>(r1,r2); }
Real mul(Real r1, Real r2) { return make_real<Mul>(r1,r2); }
Real neg(Real r) { return make_real<Neg>(r); }
Real rec(Real r) { return make_real<Rec>(r); }
Real exp(Real r) { return make_real<Exp>(r); }
Real log(Real r) { return make_real<Log>(r); }
OutputStream& operator<<(OutputStream& os, Real const& r) { os << MAGENTA; r._ptr->_write(os); os << RESET; return os; }


inline Float64Bounds evaluate(ValidatedReal const& vr, Precision64 pr) { return vr.get(pr); }

ValidatedReal::ValidatedReal(Real r) : _r(r) { }
Float64Bounds ValidatedReal::get(Precision64 pr) const { return this->_r.get(pr); }
ValidatedReal add(ValidatedReal r1, ValidatedReal r2) { return make_real<Add>(r1,r2); }
ValidatedReal mul(ValidatedReal r1, ValidatedReal r2) { return make_real<Mul>(r1,r2); }
ValidatedReal neg(ValidatedReal r) { return make_real<Neg>(r); }
ValidatedReal rec(ValidatedReal r) { return make_real<Rec>(r); }
ValidatedReal exp(ValidatedReal r) { return make_real<Exp>(r); }
ValidatedReal log(ValidatedReal r) { return make_real<Log>(r); }
OutputStream& operator<<(OutputStream& os, ValidatedReal const& vr) { os << BLUE; vr._r._ptr->_write(os); os << RESET; return os; }

// ---------------- Real/Float ------------------------------------------------------------------------------------- //

template<class F> Bounds<F>::operator GenericType() const { return make_real<Cnst,Bounds<F>>(*this); }
template Bounds<Float64>::operator GenericType() const;

// ---------------- Algebra ---------------------------------------------------------------------------------------- //

// ---------------- Differential ----------------------------------------------------------------------------------- //

template<class X> Differential<X>::Differential(X vx, X dx)
    : _vx(vx), _dx(dx) { }
template<class X> Differential<X> Differential<X>::constant(NumericType c) {
    return Differential<X>(c,c.create_constant(1)); }
template<class X> Differential<X> Differential<X>::coordinate(NumericType c) {
    return Differential<X>(c,c.create_constant(1)); }
template<class X> Differential<X> Differential<X>::create_constant(GenericNumericType c) const {
    return Differential<X>::constant(this->_vx.create_constant(c)); }

template<class X> template<class Z,EnableIf<IsConvertible<X,Z>>> Differential<X>::operator Algebra<Z> () const {
    return Algebra<Z>(new AlgebraWrapper<Differential<X>,Z>(*this));
}
template<class X> template<class Z,EnableIf<IsConvertible<X,Z>>> Differential<X>::Differential(Algebra<Z> a)
    : Differential(a.template extract<Differential<X>>()) { }

template class Differential<Real>;
template class Differential<Float64Bounds>;

template<class X> Differential<X> Differential<X>::_add(Differential<X> dx1, Differential<X> dx2) {
    return Differential<X>(add(dx1._vx,dx2._vx),add(dx1._dx,dx2._dx)); }
template<class X> Differential<X> Differential<X>::_mul(Differential<X> dx1, Differential<X> dx2) {
    return Differential<X>(mul(dx1._vx,dx2._vx),add(mul(dx1._dx,dx2._vx),mul(dx1._vx,dx2._dx))); }
template<class X> Differential<X> Differential<X>::_neg(Differential<X> dx) {
    return Differential<X>(neg(dx._vx),neg(dx._dx)); }
template<class X> Differential<X> Differential<X>::_rec(Differential<X> dx) {
    return Differential<X>(rec(dx._vx),mul(rec(neg(mul(dx._vx,dx._vx))),dx._dx)); }
template<class X> Differential<X> Differential<X>::_exp(Differential<X> dx) {
    return Differential<X>(exp(dx._vx),mul(exp(dx._vx),dx._dx)); }
template<class X> Differential<X> Differential<X>::_log(Differential<X> dx) {
    return Differential<X>(log(dx._vx),mul(rec(dx._vx),dx._dx)); }

template<class X> Differential<X> Differential<X>::_add(Differential<X> dx, Y c) {
    return Differential<X>(add(dx._vx,c),dx._dx); }
template<class X> Differential<X> Differential<X>::_add(Y c, Differential<X> dx) {
    return Differential<X>(add(c,dx._vx),dx._dx); }
template<class X> Differential<X> Differential<X>::_mul(Differential<X> dx, Y c) {
    return Differential<X>(mul(dx._vx,c),mul(dx._dx,c)); }
template<class X> Differential<X> Differential<X>::_mul(Y c, Differential<X> dx) {
    return Differential<X>(mul(dx._vx,c),mul(dx._dx,c)); }

// FIXME: Should not need to explicitly instantiate these
template Differential<Real>::operator Algebra<Real>() const;
template Differential<Float64Bounds>::operator Algebra<ValidatedReal>() const;

// ---------------- Expression ------------------------------------------------------------------------------------- //

template<class X> class ExpressionInterface {
  public:
    ~ExpressionInterface<X> () = default;
    virtual X _evaluate(Valuation<X> const&) = 0;
    virtual Algebra<X> _evaluate(Valuation<Algebra<X>> const&) = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

template<> class ExpressionInterface<Real> : public virtual ExpressionInterface<ValidatedReal> {
    typedef Real X;
  public:
    virtual X _evaluate(Valuation<X> const&) = 0;
    virtual Algebra<X> _evaluate(Valuation<Algebra<X>> const&) = 0;
};

template<class X, class OP, class... AS> struct ExpressionWrapper
    : virtual ExpressionInterface<X>, Symbolic<OP,AS...>
{
    ExpressionWrapper(AS... as)
        : Symbolic<OP,AS...>{as...} { }
    virtual X _evaluate(Valuation<X> const& v) {
        return evaluate_as<X>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual Algebra<X> _evaluate(Valuation<Algebra<X>> const& v) {
        return evaluate_as<Algebra<X>>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual OutputStream& _write(OutputStream& os) const {
        os << static_cast<Symbolic<OP,AS...>const&>(*this); return os; }
};

template<class OP, class... AS> struct ExpressionWrapper<Real,OP,AS...>
    : virtual ExpressionInterface<Real>, Symbolic<OP,AS...>
{
    typedef Real X;
    typedef ValidatedReal VX;
    ExpressionWrapper(AS... as)
        : Symbolic<OP,AS...>{as...} { }
    virtual VX _evaluate(Valuation<VX> const& v) {
        return evaluate_as<VX>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual Algebra<VX> _evaluate(Valuation<Algebra<VX>> const& v) {
        return evaluate_as<Algebra<VX>>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual X _evaluate(Valuation<X> const& v) {
        return evaluate_as<X>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual Algebra<X> _evaluate(Valuation<Algebra<X>> const& v) {
        return evaluate_as<Algebra<X>>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual OutputStream& _write(OutputStream& os) const {
        os << static_cast<Symbolic<OP,AS...>const&>(*this); return os; }
};

template<class X, class OP, class... AS> Expression<X> make_expression(AS... as) {
    return Expression<X>(new ExpressionWrapper<X,OP,AS...>(as...)); }


template<class X> template<class SX, EnableIf<IsConvertible<SX,X>>> Expression<X>::Expression(Expression<SX> e)
    : _ptr(std::dynamic_pointer_cast<ExpressionInterface<X>>(e._ptr)) { if(!_ptr) { std::cerr << "ERROR: e="<<e<<"\n"; throw std::bad_cast(); } }
template Expression<ValidatedReal>::Expression(Expression<Real>);

template<class X> Expression<X>::Expression(ExpressionInterface<X>* p) : _ptr(p) { }
template<class X> Expression<X>::operator Algebra<NumericType> () const { return make_handle(new AlgebraWrapper<Expression<X>,X>(*this)); }
template<class X> Expression<X>::Expression(Algebra<X> const& a) : Expression(a.template extract<Expression<X>>()) { }

template<class X> Expression<X>::Expression(NumericType c) : Expression(new ExpressionWrapper<X,Cnst,X>{c}) { }
template<class X> Expression<X>::Expression(VariableType v) : Expression(new ExpressionWrapper<X,Var,Identifier>{v.name()}) { }
template<class X> Expression<X> Expression<X>::constant(X c) { return Expression<X>(c); }
template<class X> Expression<X> Expression<X>::create_constant(X c) const { return Expression<X>(c); }

template<class X> X Expression<X>::operator()(Valuation<X> v) const { return evaluate(*this,v); }

template class Expression<Real>;
template class Expression<ValidatedReal>;

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

Real evaluate(RealExpression e, RealValuation v) { return e._ptr->_evaluate(v); }
RealAlgebra evaluate(RealExpression e, Valuation<RealAlgebra> v) { return e._ptr->_evaluate(v); }
//OutputStream& operator<<(OutputStream& os, RealExpression const& e) { os << "RealExpression("; e._ptr->_write(os); return os << ")"; }
OutputStream& operator<<(OutputStream& os, RealExpression const& e) { os << MAGENTA; e._ptr->_write(os); os << RESET; return os; }

ValidatedRealExpression add(ValidatedRealExpression e1, ValidatedRealExpression e2) { return make_expression<ValidatedReal,Add>(e1,e2); }
ValidatedRealExpression mul(ValidatedRealExpression e1, ValidatedRealExpression e2) { return make_expression<ValidatedReal,Mul>(e1,e2); }
ValidatedRealExpression neg(ValidatedRealExpression e) { return make_expression<ValidatedReal,Neg>(e); }
ValidatedRealExpression rec(ValidatedRealExpression e) { return make_expression<ValidatedReal,Rec>(e); }
ValidatedRealExpression exp(ValidatedRealExpression e) { return make_expression<ValidatedReal,Exp>(e); }
ValidatedRealExpression log(ValidatedRealExpression e) { return make_expression<ValidatedReal,Log>(e); }
ValidatedReal evaluate(ValidatedRealExpression e, Valuation<ValidatedReal> v) { return e._ptr->_evaluate(v); }
ValidatedRealAlgebra evaluate(ValidatedRealExpression e, Valuation<ValidatedRealAlgebra> v) { return e._ptr->_evaluate(v); }
OutputStream& operator<<(OutputStream& os, ValidatedRealExpression const& e) { os << BLUE; e._ptr->_write(os); os << RESET; return os; }

// ---------------- Domain ----------------------------------------------------------------------------------------- //

EuclideanSpace intersection(EuclideanSpace spc1, EuclideanSpace spc2) {
    if(spc1.dimension() != spc2.dimension()) { throw std::runtime_error("Dimension mismatch"); }
    return spc1;
}

// ---------------- Function --------------------------------------------------------------------------------------- //

template<class X> class FormulaInterface {
  public:
    ~FormulaInterface () = default;
    virtual Float64Bounds _evaluate(Vector<Float64Bounds> const&) = 0;
    virtual X _evaluate(Vector<X> const&) = 0;
    virtual Algebra<X> _evaluate(Vector<Algebra<X>> const&) = 0;
    virtual OutputStream& _write(OutputStream& os) const = 0;
};
template<> class FormulaInterface<Real> : public virtual FormulaInterface<ValidatedReal> {
    typedef Real X;
  public:
    using FormulaInterface<ValidatedReal>::_evaluate;
    virtual X _evaluate(Vector<X> const&) = 0;
    virtual Algebra<X> _evaluate(Vector<Algebra<X>> const&) = 0;
};

template<class X> OutputStream& operator<<(OutputStream& os, Formula<X> const& f) { f._ptr->_write(os); return os; }
template<class X, class Y> Y evaluate(Formula<X> const& f, Vector<Y> const& v) { return f._ptr->_evaluate(v); }
template<class X, class Y> Algebra<Y> evaluate(Formula<X> const& f, Vector<Algebra<Y>> const& v) { return f._ptr->_evaluate(v); }

template<class X, class OP, class... AS> struct FormulaWrapper
    : virtual FormulaInterface<X>, Symbolic<OP,AS...>
{
    FormulaWrapper(AS... as) : Symbolic<OP,AS...>{as...} { }
    virtual Float64Bounds _evaluate(Vector<Float64Bounds> const& v) {
        return evaluate_as<Float64Bounds>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual X _evaluate(Vector<X> const& v) {
        return evaluate_as<X>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual Algebra<X> _evaluate(Vector<Algebra<X>> const& v) {
        return evaluate_as<Algebra<X>>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<Symbolic<OP,AS...>const&>(*this); }
};

template<class OP, class... AS> struct FormulaWrapper<Real,OP,AS...>
    : virtual FormulaInterface<Real>, Symbolic<OP,AS...>
{
    typedef Real X;
    typedef ValidatedReal VX;
    FormulaWrapper(AS... as) : Symbolic<OP,AS...>{as...} { }
    virtual Float64Bounds _evaluate(Vector<Float64Bounds> const& v) {
        return evaluate_as<Float64Bounds>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual VX _evaluate(Vector<VX> const& v) {
        return evaluate_as<VX>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual Algebra<VX> _evaluate(Vector<Algebra<VX>> const& v) {
        return evaluate_as<Algebra<VX>>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual X _evaluate(Vector<X> const& v) {
        return evaluate_as<X>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual Algebra<X> _evaluate(Vector<Algebra<X>> const& v) {
        return evaluate_as<Algebra<X>>(static_cast<Symbolic<OP,AS...>const&>(*this),v); }
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<Symbolic<OP,AS...>const&>(*this); }
};

template<class X, class OP, class... AS> Function<X> make_function(typename Function<X>::DomainType dom, AS... as) {
    return Function<X>(dom, Formula<X>(new FormulaWrapper<X,OP,AS...>(as...))); }

template<class X> template<class SX, EnableIf<IsConvertible<SX,X>>> Function<X>::Function(Function<SX> f)
    : _dom(f.domain()), _frml(std::dynamic_pointer_cast<FormulaInterface<X>>(f.formula()._ptr))
{
    if(!this->_frml._ptr) { std::cerr << "ERROR: f="<<f<<"\n"; std::cerr<<typeid(*f.formula()._ptr).name()<<"\n"; throw std::bad_cast(); }
}

template Function<ValidatedReal>::Function(Function<Real>);

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

template<class X> auto Function<X>::operator() (Vector<NumericType> v) const -> NumericType {
    return evaluate(this->formula(),v);
}
template<class X> auto Function<X>::operator() (Vector<AlgebraType> v) const -> AlgebraType {
    return evaluate(this->formula(),v);
}
template<class X> auto Function<X>::operator() (Vector<Float64Bounds> v) const -> Float64Bounds {
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

template<class X> Function<X> add(Function<X> f1,Function<X> f2) {
    return make_function<X,Add>(intersection(f1.domain(),f2.domain()),f1.formula(),f2.formula()); }
template<class X> Function<X> mul(Function<X> f1,Function<X> f2) {
    return make_function<X,Mul>(intersection(f1.domain(),f2.domain()),f1.formula(),f2.formula()); }
template<class X> Function<X> neg(Function<X> f) { return make_function<X,Neg>(f.domain(),f.formula()); }
template<class X> Function<X> rec(Function<X> f) { return make_function<X,Rec>(f.domain(),f.formula()); }
template<class X> Function<X> exp(Function<X> f) { return make_function<X,Exp>(f.domain(),f.formula()); }
template<class X> Function<X> log(Function<X> f) { return make_function<X,Log>(f.domain(),f.formula()); }


template class Function<Real>;
typedef Formula<Real> RealFormula;

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

Float64Bounds evaluate(RealFunction f, Vector<Float64Bounds> v) { return evaluate(f.formula(),v); }
Real evaluate(RealFunction f, RealVector v) { return evaluate(f.formula(),v); }
Real evaluate(RealFormula f, RealVector v) { return f._ptr->_evaluate(v); }
//OutputStream& operator<<(OutputStream& os, RealFunction const& e) { os << "RealFunction("; e._ptr->_write(os); return os << ")"; }
OutputStream& operator<<(OutputStream& os, RealFunction const& f) { f._frml._ptr->_write(os); return os; }


template class Function<ValidatedReal>;
typedef Formula<ValidatedReal> ValidatedRealFormula;

ValidatedRealFunction add(ValidatedRealFunction f1, ValidatedRealFunction f2) {
    return make_function<ValidatedReal,Add>(intersection(f1.domain(),f2.domain()),f1.formula(),f2.formula()); }
ValidatedRealFunction mul(ValidatedRealFunction f1, ValidatedRealFunction f2) {
    return make_function<ValidatedReal,Mul>(intersection(f1.domain(),f2.domain()),f1.formula(),f2.formula()); }
ValidatedRealFunction neg(ValidatedRealFunction f) { return make_function<ValidatedReal,Neg>(f.domain(),f.formula()); }
ValidatedRealFunction rec(ValidatedRealFunction f) { return make_function<ValidatedReal,Rec>(f.domain(),f.formula()); }
ValidatedRealFunction exp(ValidatedRealFunction f) { return make_function<ValidatedReal,Exp>(f.domain(),f.formula()); }
ValidatedRealFunction log(ValidatedRealFunction f) { return make_function<ValidatedReal,Log>(f.domain(),f.formula()); }

ValidatedRealFunction add(ValidatedRealFunction f1, ValidatedReal c2) { return add(f1,f1.create_constant(c2)); }
ValidatedRealFunction add(ValidatedReal c1, ValidatedRealFunction f2) { return add(f2.create_constant(c1),f2); }
ValidatedRealFunction mul(ValidatedRealFunction f1, ValidatedReal c2) { return mul(f1,f1.create_constant(c2)); }
ValidatedRealFunction mul(ValidatedReal c1, ValidatedRealFunction f2) { return mul(f2.create_constant(c1),f2); }

ValidatedReal evaluate(ValidatedRealFunction f, ValidatedRealVector v) { return evaluate(f.formula(),v); }
ValidatedReal evaluate(ValidatedRealFormula f, ValidatedRealVector v) { return f._ptr->_evaluate(v); }
//OutputStream& operator<<(OutputStream& os, ValidatedRealFunction const& e) { os << "ValidatedRealFunction("; e._ptr->_write(os); return os << ")"; }
OutputStream& operator<<(OutputStream& os, ValidatedRealFunction const& f) { f._frml._ptr->_write(os); return os; }

template<class X> SizeType instantiate_operators() {
    auto add_ptr=(X(*)(X,X)) &add;
    auto mul_ptr=(X(*)(X,X)) &mul;
    auto neg_ptr=(X(*)(X)) &neg;
    auto rec_ptr=(X(*)(X)) &rec;
    auto exp_ptr=(X(*)(X)) &exp;
    auto log_ptr=(X(*)(X)) &log;

    return (SizeType)add_ptr + (SizeType)mul_ptr
        + (SizeType)neg_ptr + (SizeType)rec_ptr
        + (SizeType)exp_ptr + (SizeType)log_ptr;
}
template SizeType instantiate_operators<RealFunction>();
template SizeType instantiate_operators<ValidatedRealFunction>();
