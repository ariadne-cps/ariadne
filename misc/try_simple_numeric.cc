#include <utility>
#include <iostream>
#include <type_traits>

using OutputStream = std::ostream;
using std::declval;

class Integer {
    long long int _n;
  public:
    Integer(int n) : _n(n) { }
    Integer(double d) = delete;
};

enum Signum  { LESS=-1, EQUAL=0, GREATER=+1, UNKNOWN=0 };

class Rational {
    Integer _num; Integer _den;
  public:
    Rational(Integer num, Integer den=1);
    friend Rational max(Rational,Rational);
    friend Rational add(Rational,Rational);
    friend Rational mul(Rational,Rational);
    friend Rational neg(Rational);
    friend Rational rec(Rational);
    friend Rational abs(Rational);
    friend Signum cmp(Rational, Rational);
};

Rational max(Rational,Rational);
Rational add(Rational,Rational);
Rational mul(Rational,Rational);
Rational neg(Rational);
Rational rec(Rational);
Rational abs(Rational);

struct Effort {
    int _bits;
    Effort(int bits) : _bits(bits) { }
    Effort& operator++() { ++_bits; return *this; }
};

class Comparison {
    Signum check(Effort e);
};

class LowerReal;
class UpperReal;
class Real;
class PositiveLowerReal;
class PositiveUpperReal;
class PositiveReal;

class LowerReal {
  public:
    LowerReal();
    Rational get_lower(Effort e);
    friend LowerReal max(LowerReal,LowerReal);
    friend LowerReal add(LowerReal,LowerReal);
    friend LowerReal neg(UpperReal);
    friend UpperReal neg(LowerReal);
};

class UpperReal {
  public:
    UpperReal();
    UpperReal(Rational q);
    Rational get_upper(Effort e);
    friend LowerReal max(LowerReal,LowerReal);
    friend LowerReal add(LowerReal,LowerReal);
    friend LowerReal neg(UpperReal);
    friend UpperReal neg(LowerReal);
};

class Real {
  public:
    operator LowerReal() const;
    operator UpperReal() const;

    Real();
    Real(Rational q);
    Rational get_lower(Effort e);
    Rational get_upper(Effort e);

    friend Real max(Real,Real);
    friend Real add(Real,Real);
    friend Real mul(Real,Real);
    friend Real neg(Real);
    friend Real rec(Real);
    friend PositiveReal abs(Real);

    friend PositiveReal exp(Real);

    friend PositiveReal dist(Real,Real);
    Comparison cmp(Real, Real);
};

class PositiveLowerReal : LowerReal {
  public:
    PositiveLowerReal();
    PositiveLowerReal(Rational q);
    friend PositiveLowerReal max(PositiveLowerReal,PositiveLowerReal);
    friend PositiveLowerReal add(PositiveLowerReal,PositiveLowerReal);
    friend PositiveLowerReal mul(PositiveLowerReal,PositiveLowerReal);
    friend PositiveLowerReal div(PositiveLowerReal,PositiveUpperReal);
    friend PositiveLowerReal rec(PositiveUpperReal);
    friend PositiveUpperReal rec(PositiveLowerReal);
};

class PositiveUpperReal : UpperReal {
  public:
    PositiveUpperReal();
    PositiveUpperReal(Rational q);
    friend PositiveUpperReal max(PositiveUpperReal,PositiveUpperReal);
    friend PositiveUpperReal add(PositiveUpperReal,PositiveUpperReal);
    friend PositiveUpperReal mul(PositiveUpperReal,PositiveUpperReal);
    friend PositiveUpperReal div(PositiveUpperReal,PositiveLowerReal);
    friend PositiveUpperReal rec(PositiveLowerReal);
    friend PositiveLowerReal rec(PositiveUpperReal);
};

class PositiveReal : Real {
  public:
    PositiveReal();
    PositiveReal(Rational q);
};

struct ApproximateTag { };
struct ValidatedTag { };
struct EffectiveTag { };
struct ExactTag { };

struct LowerTag { };
struct UpperTag { };
struct RangeTag { };
struct MetricTag { };


template<class A1, class A2> inline auto min(A1 const& a1, A2 const& a2) -> decltype(neg(max(neg(a1),neg(a2)))) {
    return neg(max(neg(a1),neg(a2))); }
template<class A1, class A2> inline auto sub(A1&& x1, A2&& x2) -> decltype(add(x1,neg(x2))) {
    return add(std::forward<A1>(x1),neg(std::forward<A2>(x2))); }
template<class A1, class A2> inline auto div(A1&& x1, A2&& x2) -> decltype(mul(x1,rec(x2))) {
    return mul(std::forward<A1>(x1),rec(std::forward<A2>(x2))); }

template<class A1, class A2> inline auto generic_max(A1 a1, A2 a2) -> decltype(max(a1,a2)) { return max(a1,a2); }
template<class A1, class A2, class... ARGS> inline auto generic_max(A1 a1, A2 a2, ARGS... as) -> decltype(generic_max(max(a1,a2),as...)) {
    return generic_max(max(a1,a2),as...); }


struct DefineArithmeticOperators {
template<class X1, class X2> friend auto operator==(X1 const& x1, X2 const& x2) -> decltype(cmp(x1,x2)==0) { return cmp(x1,x2)==0; }
template<class X1, class X2> friend auto operator!=(X1 const& x1, X2 const& x2) -> decltype(cmp(x1,x2)!=0) { return cmp(x1,x2)!=0; }
template<class X1, class X2> friend auto operator<=(X1 const& x1, X2 const& x2) -> decltype(cmp(x1,x2)<=0) { return cmp(x1,x2)<=0; }
template<class X1, class X2> friend auto operator>=(X1 const& x1, X2 const& x2) -> decltype(cmp(x1,x2)>=0) { return cmp(x1,x2)>=0; }
template<class X1, class X2> friend auto operator< (X1 const& x1, X2 const& x2) -> decltype(cmp(x1,x2)< 0) { return cmp(x1,x2)< 0; }
template<class X1, class X2> friend auto operator> (X1 const& x1, X2 const& x2) -> decltype(cmp(x1,x2)> 0) { return cmp(x1,x2)> 0; }

template<class X> friend auto operator+(X const& x) -> decltype(neg(neg(x))) { return x; }
template<class X> friend auto operator-(X const& x) -> decltype(neg(x)) { return neg(x); }

template<class X1, class X2> friend auto operator+(X1 const& x1, X2 const& x2) -> decltype(add(x1,x2)) { return add(x1,x2); }
template<class X1, class X2> friend auto operator-(X1 const& x1, X2 const& x2) -> decltype(sub(x1,x2)) { return sub(x1,x2); }
template<class X1, class X2> friend auto operator*(X1 const& x1, X2 const& x2) -> decltype(mul(x1,x2)) { return mul(x1,x2); }
template<class X1, class X2> friend auto operator/(X1 const& x1, X2 const& x2) -> decltype(div(x1,x2)) { return div(x1,x2); }

template<class X1, class X2> friend auto operator+=(X1& x1, X2 const& x2) -> decltype(x1=add(x1,x2)) { return x1=add(x1,x2); }
template<class X1, class X2> friend auto operator-=(X1& x1, X2 const& x2) -> decltype(x1=sub(x1,x2)) { return x1=sub(x1,x2); }
template<class X1, class X2> friend auto operator*=(X1& x1, X2 const& x2) -> decltype(x1=mul(x1,x2)) { return x1=mul(x1,x2); }
template<class X1, class X2> friend auto operator/=(X1& x1, X2 const& x2) -> decltype(x1=div(x1,x2)) { return x1=div(x1,x2); }
};

class Float64 : DefineArithmeticOperators {
    volatile double _dbl;
  public:
    Float64() : _dbl() { }
    Float64(int n) : _dbl(n) { }
    Float64(float f) : _dbl(f) { }
    Float64(double d) : _dbl(d) { }
    explicit operator Rational() const;
  public:
    static Float64 inf() { return 1.0/0.0; }
  public:
    friend Float64 neg_exact(Float64 x) { return -x._dbl; }
    friend Float64 max_exact(Float64 x1, Float64 x2) { return x1._dbl>=x2._dbl ? x1._dbl : x2._dbl; }
    friend Float64 add_down(Float64 x1, Float64 x2) { return x1._dbl+x2._dbl; }
    friend Float64 add_up(Float64 x1,Float64 x2) { return x1._dbl+x2._dbl; }
    friend Float64 add_near(Float64 x1, Float64 x2) { return x1._dbl+x2._dbl; }
    friend Float64 mul_down(Float64 x1, Float64 x2) { return x1._dbl*x2._dbl; }
    friend Float64 mul_up(Float64 x1, Float64 x2) { return x1._dbl*x2._dbl; }
    friend Float64 mul_near(Float64 x1, Float64 x2) { return x1._dbl*x2._dbl; }
    friend Float64 rec_down(Float64 x) { return 1.0/x._dbl; }
    friend Float64 rec_up(Float64 x) { return 1.0/x._dbl; }
    friend Float64 rec_near(Float64 x) { return 1.0/x._dbl; }
    friend Float64 abs_exact(Float64 x) { return x._dbl>=0.0 ? x._dbl : -x._dbl; }

    friend Float64 max(Float64 x1, Float64 x2) { return max_exact(x1,x2); }
    friend Float64 abs(Float64 x) { return abs_exact(x); }
    friend Float64 neg(Float64 x) { return neg_exact(x); }

    friend Signum cmp(Float64 x1, Float64 x2) { return x1._dbl==x2._dbl ? EQUAL : x1._dbl>x2._dbl ? GREATER : LESS; }
    friend OutputStream& operator<<(OutputStream& os, Float64 x) { return os << x._dbl; }
};

class Float32 : DefineArithmeticOperators {
    volatile float _flt;
  public:
    Float32() : _flt() { }
    Float32(int n) : _flt(n) { }
    Float32(float d) : _flt(d) { }
    explicit operator Rational() const;
  public:
    static Float32 inf() { return 1.0f/0.0f; }
  public:
    friend Float32 neg_exact(Float32 x) { return -x._flt; }
    friend Float32 max_exact(Float32 x1, Float32 x2) { return x1._flt>=x2._flt ? x1._flt : x2._flt; }
    friend Float32 add_down(Float32 x1, Float32 x2) { return x1._flt+x2._flt; }
    friend Float32 add_up(Float32 x1,Float32 x2) { return x1._flt+x2._flt; }
    friend Float32 add_near(Float32 x1, Float32 x2) { return x1._flt+x2._flt; }
    friend Float32 mul_down(Float32 x1, Float32 x2) { return x1._flt*x2._flt; }
    friend Float32 mul_up(Float32 x1, Float32 x2) { return x1._flt*x2._flt; }
    friend Float32 mul_near(Float32 x1, Float32 x2) { return x1._flt*x2._flt; }
    friend Float32 rec_down(Float32 x) { return 1.0f/x._flt; }
    friend Float32 rec_up(Float32 x) { return 1.0f/x._flt; }
    friend Float32 rec_near(Float32 x) { return 1.0f/x._flt; }
    friend Float32 abs_exact(Float32 x) { return x._flt>=0.0f ? x._flt : -x._flt; }

    friend Float32 max(Float32 x1, Float32 x2) { return max_exact(x1,x2); }
    friend Float32 abs(Float32 x) { return abs_exact(x); }
    friend Float32 neg(Float32 x) { return neg_exact(x); }

    friend Signum cmp(Float32 x1, Float32 x2) { return x1._flt==x2._flt ? EQUAL : x1._flt>x2._flt ? GREATER : LESS; }
    friend OutputStream& operator<<(OutputStream& os, Float32 x) { return os << x._flt; }
};

template<class F> class Approximate;
template<class F> class Lower;
template<class F> class Upper;
template<class F> class Range;
template<class F, class EF> class Ball;

template<class X> using NegType = decltype(neg(declval<X>()));
template<class X> using RecType = decltype(rec(declval<X>()));
template<class X1, class X2=X1> using MaxType = decltype(max(declval<X1>(),declval<X2>()));
template<class X1, class X2=X1> using AddType = decltype(add(declval<X1>(),declval<X2>()));
template<class X1, class X2=X1> using MulType = decltype(mul(declval<X1>(),declval<X2>()));

template<class F, class R> struct ConcreteTypeTraits;
template<class F> struct ConcreteTypeTraits<F,LowerReal> { typedef Lower<F> Type; };
template<class F> struct ConcreteTypeTraits<F,UpperReal>  { typedef Upper<F> Type; };
template<class F> struct ConcreteTypeTraits<F,Real>  { typedef Range<F> Type; };
template<class F, class R> using ConcreteType = typename ConcreteTypeTraits<F,R>::Type;

template<class X> struct Traits;
template<class F> struct Traits<Range<F>> { typedef F RawType; typedef Real ModelType; };
template<class F> struct Traits<Lower<F>> { typedef F RawType; typedef LowerReal ModelType; };
template<class F> struct Traits<Upper<F>> { typedef F RawType; typedef UpperReal ModelType; };

template<class T> using ModelType = typename Traits<T>::ModelType;

template<class X> struct DeclareDirectedOperations {
    typedef typename Traits<X>::RawType F;
    typedef ConcreteType<F,NegType<ModelType<X>>> NEGX;
    friend auto neg(X) -> ConcreteType<F,NegType<ModelType<X>>>;
    friend auto max(X,X) -> ConcreteType<F,MaxType<ModelType<X>>>;
    friend auto add(X,X) -> ConcreteType<F,AddType<ModelType<X>>>;
//    friend auto sub(X,NEGX) -> ConcreteType<F,AddType<ModelType<X>>>;
};

template<class X> struct DeclarePositiveDirectedOperations {
    typedef typename Traits<X>::RawType F;
    typedef ConcreteType<F,RecType<ModelType<X>>> RECX;
    friend auto rec(X) -> ConcreteType<F,RecType<ModelType<X>>>;
    friend auto max(X,X) -> ConcreteType<F,MaxType<ModelType<X>>>;
    friend auto add(X,X) -> ConcreteType<F,AddType<ModelType<X>>>;
    friend auto mul(X,X) -> ConcreteType<F,MulType<ModelType<X>>>;
//    friend auto div(X,RECX) -> ConcreteType<F,MulType<ModelType<X>>>;
};

template<class X> struct DeclareOperations : DeclareDirectedOperations<X>, DeclarePositiveDirectedOperations<X> {
};

template<class F> class Approximate : DefineArithmeticOperators, DeclareOperations<Approximate<F>> {
  protected:
    F _a;
  public:
    Approximate<F>(F a) : _a(a) { }
    F raw() { return _a; }
    friend Approximate<F> max(Approximate<F> x1, Approximate<F> x2) { return Approximate<F>(max_exact(x1._a, x2._a)); }
    friend Approximate<F> add(Approximate<F> x1, Approximate<F> x2) { return Approximate<F>(add_near(x1._a, x2._a)); }
    friend Approximate<F> mul(Approximate<F> x1, Approximate<F> x2) { return Approximate<F>(mul_near(x1._a, x2._a)); }
    friend Approximate<F> neg(Approximate<F> x) { return Approximate<F>(neg_exact(x._a)); }
    friend Approximate<F> rec(Approximate<F> x) { return Approximate<F>(rec_near(x._a)); }
};

template<class F> class Lower : public DefineArithmeticOperators, DeclareDirectedOperations<Lower<F>> {
  protected:
    F _l;
  public:
    typedef LowerReal ModelType;
    typedef ValidatedTag InformationTag;
    Lower<F>(F l) : _l(l) { }
    F lower_raw() { return _l; }
    friend Lower<F> max(Lower<F> x1, Lower<F> x2) { return Lower<F>(max(x1._l, x2._l)); }
    friend Lower<F> add(Lower<F> x1, Lower<F> x2) { return Lower<F>(add_down(x1._l, x2._l)); }
    friend Lower<F> neg(Upper<F> x);
    friend Upper<F> neg(Lower<F> x);
    friend OutputStream& operator<<(OutputStream& os, Lower<F> x) { return os << x._l << ":"; }
};

template<class F> class Upper : public DefineArithmeticOperators, DeclareDirectedOperations<Lower<F>> {
  protected:
    F _u;
  public:
    typedef UpperReal ModelType;
    typedef ValidatedTag InformationTag;
    Upper<F>(F u) : _u(u) { }
    F upper_raw() { return _u; }
    friend Upper<F> max(Upper<F> x1, Upper<F> x2) { return Upper<F>(max(x1._u, x2._u)); }
    friend Upper<F> add(Upper<F> x1, Upper<F> x2) { return Upper<F>(add_up(x1._u, x2._u)); }
    friend Upper<F> neg(Lower<F> x) { return Upper<F>(neg_exact(x._l)); }
    friend Lower<F> neg(Upper<F> x) { return Lower<F>(neg_exact(x._u)); }
    friend OutputStream& operator<<(OutputStream& os, Upper<F> x) { return os << ":" << x._u; }
};

template<class F> class PositiveRange;

template<class F> class Range;
template<class F> Range<F> generic_max(Range<F> const& x1, Range<F> const& x2);
template<class F> Range<F> generic_add(Range<F> const& x1, Range<F> const& x2);
template<class F> Range<F> generic_sub(Range<F> const& x1, Range<F> const& x2);
template<class F> Range<F> generic_mul(Range<F> const& x1, Range<F> const& x2);
template<class F> Range<F> generic_neg(Range<F> const& x);
template<class F> Range<F> generic_rec(Range<F> const& x);

template<class F> class Range : public Lower<F>, public Upper<F>, public DefineArithmeticOperators, public DeclareOperations<Range<F>> {
    template<class FF> friend Range<FF> generic_max(Range<FF> const& x1, Range<FF> const& x2);
    template<class FF> friend Range<FF> generic_add(Range<FF> const& x1, Range<FF> const& x2);
    template<class FF> friend Range<FF> generic_sub(Range<FF> const& x1, Range<FF> const& x2);
    template<class FF> friend Range<FF> generic_mul(Range<FF> const& x1, Range<FF> const& x2);
    template<class FF> friend Range<FF> generic_neg(Range<FF> const& x);
    template<class FF> friend Range<FF> generic_rec(Range<FF> const& x);
  public:
    typedef Real ModelType;
    typedef ValidatedTag InformationTag;
    Range<F>(Lower<F> l, Upper<F> u) : Lower<F>(l), Upper<F>(u) { }
    Range<F>(F l, F u) : Lower<F>(l), Upper<F>(u) { }

    friend Range<F> max(Range<F> x1, Range<F> x2) { return generic_max<F>(x1,x2); }
    friend Range<F> max(Upper<F> x1, Range<F> x2) { return Range<F>(x2._l,max(x1._u,x2._u)); }
    friend Range<F> max(Range<F> x1, Upper<F> x2) { return Range<F>(x1._l,max(x1._u,x2._u)); }
//    friend Range<F> min(Range<F> x1, Range<F> x2) { return neg(max(neg(x1),neg(x2))); }
    friend Range<F> add(Range<F> x1, Range<F> x2) { return generic_add<F>(x1,x2); }
//    friend Range<F> sub(Range<F> x1, Range<F> x2) { return generic_sub<F>(x1,x2); }
    friend Range<F> mul(Range<F> x1, Range<F> x2) { return generic_mul<F>(x1,x2); }
//    friend Range<F> div(Range<F> x1, Range<F> x2) { return mul(x1,rec(x2)); }
    friend Range<F> neg(Range<F> x) { return Range<F>(neg_exact(x._u),neg_exact(x._l)); }
    friend Range<F> rec(Range<F> x) { return generic_rec<F>(x); }
    friend PositiveRange<F> abs(Range<F> x);

    friend PositiveRange<F> dist(Range<F> x1, Range<F> x2) { return abs(sub(x1,x2)); }
    friend Signum cmp(Range<F> x1, Range<F> x2) { return x1._l>x2._u ? GREATER : x1._u<x2._l ? LESS : UNKNOWN; }

    friend OutputStream& operator<<(OutputStream& os, Range<F> x) { return os << x._l << ":" << x._u; }
};

template<class F> class PositiveLower;
template<class F> class PositiveUpper;

template<class F> class PositiveLower : public Lower<F> {
  public:
    PositiveLower(Lower<F> x);
    PositiveLower(F u) : Lower<F>(u) { }

    friend PositiveLower<F> add(PositiveLower<F> x1, PositiveLower<F> x2) {  return PositiveLower<F>(add_down(x1._u,x2._u)); }
    friend PositiveLower<F> mul(PositiveLower<F> x1, PositiveLower<F> x2) {  return PositiveLower<F>(mul_down(x1._u,x2._u)); }
    friend PositiveUpper<F> rec(PositiveLower<F> x);
    friend PositiveLower<F> rec(PositiveUpper<F> x);
};

template<class F> class PositiveUpper : public Upper<F> {
  public:
    PositiveUpper(Upper<F> x);
    PositiveUpper(F u) : Upper<F>(u) { }

    friend PositiveUpper<F> add(PositiveUpper<F> x1, PositiveUpper<F> x2) {  return PositiveUpper<F>(add_up(x1._u,x2._u)); }
    friend PositiveUpper<F> mul(PositiveUpper<F> x1, PositiveUpper<F> x2) {  return PositiveUpper<F>(mul_up(x1._u,x2._u)); }
    friend PositiveUpper<F> rec(PositiveLower<F> x) {  return PositiveUpper<F>(rec_up(x._l)); }
    friend PositiveLower<F> rec(PositiveUpper<F> x) {  return PositiveLower<F>(rec_down(x._u)); }
};

template<class F> class PositiveRange : public Range<F> {
  public:
    PositiveRange(Range<F> x);
    PositiveRange(F l, F u) : Range<F>(l,u) { }

    friend PositiveRange<F> abs(Range<F> x) { return PositiveRange(max(max(x._l,neg(x._u)),0.0f),max(neg(x._l),x._u)); }
};


template<class F> Range<F> generic_max(Range<F> const& x1, Range<F> const& x2) {
    return Range<F>(max(x1._l,x2._l),max(x1._u,x2._u)); }

template<class F> Range<F> generic_add(Range<F> const& x1, Range<F> const& x2) {
    return Range<F>(add_down(x1._l,x2._l),add_up(x1._u,x2._u)); }

template<class F> Range<F> generic_sub(Range<F> const& x1, Range<F> const& x2) {
    return Range<F>(add_down(x1._l,neg_exact(x2._u)),add_up(x1._u,neg_exact(x2._l))); }

template<class F> Range<F> generic_mul(Range<F> const& x1, Range<F> const& x2) {
    F rl=min(min(mul_down(x1._l,x2._l),mul_down(x1._l,x2._u)),min(mul_down(x1._u,x2._l),mul_down(x1._u,x2._u)));
    F ru=max(max(mul_up(x1._l,x2._l),mul_up(x1._l,x2._u)),max(mul_up(x1._u,x2._l),mul_up(x1._u,x2._u)));
    return Range<F>(rl,ru); }

template<class F> Range<F> generic_neg(Range<F> const& x) {
    return Range<F>(neg_exact(x._u),neg_exact(x._l)); }

template<class F> Range<F> generic_rec(Range<F> const& x) {
    if(x._l>0 || x._u<0) { return Range<F>(rec_down(x._u),rec_up(x._l)); }
    else { return Range<F>(neg(F::inf()),F::inf()); } }

#define PRINT(expression) std::cout << #expression << ": " << (expression) << "\n";

template<class F> void test_float() {
    Range<F> x(-0.5f,1.5f);
    Range<F> y(1.0f,3.5f);

    PRINT(x)
    PRINT(y)
    PRINT(max(x,y))
    PRINT(min(x,y))
    PRINT(add(x,y))
    PRINT(sub(x,y))
    PRINT(mul(x,y))
    PRINT(div(x,y))
    PRINT(neg(x))
    PRINT(rec(y))
    PRINT(abs(x))

    PRINT(dist(x,y));
    PRINT(cmp(x,y));

    Upper<F> ux(-1.5f);
    Lower<F> lx(-3.25f);
    PRINT(neg(lx))
    PRINT(add(lx,lx))
    PRINT(add(ux,ux))
    PRINT(sub(lx,ux))
    PRINT(sub(ux,lx))

    PositiveUpper<F> pux(1.5f);
    PRINT(pux);
    PRINT(neg(pux));
    PRINT(add(pux,pux));
    PRINT(mul(pux,pux));
    PRINT(div(pux,rec(pux)));
    PRINT(rec(pux));

    +x; -x; x+x; x-x; x*x; x/x; x+=x; x-=x; x*=x; x/=x;
    auto r=(x==x) xor (x!=x) xor x<=x xor x>=x xor x<x xor x>x;
}

int main() {
    test_float<Float32>();
    std::cout << std::endl;
    test_float<Float64>();
}

Real::Real() { }
LowerReal::LowerReal() { }
UpperReal::UpperReal() { }
PositiveReal::PositiveReal() { }
PositiveLowerReal::PositiveLowerReal() { }
PositiveUpperReal::PositiveUpperReal() { }
PositiveUpperReal rec(PositiveLowerReal plr) { return PositiveUpperReal(); }
PositiveLowerReal rec(PositiveUpperReal pur) { return PositiveLowerReal(); }
