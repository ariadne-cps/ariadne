#include "simple.h"

#include <limits>
#include <cmath>

// ---------------- Float ------------------------------------------------------------------------------------------ //

#if ( defined __SSE_MATH__ &&  defined __SSE2__ )
    #define ARIADNE_SSE_ROUNDING
#else
    #define ARIADNE_GCC_ROUNDING
#endif

#if defined ARIADNE_SSE_ROUNDING
#include <xmmintrin.h>
typedef unsigned short rounding_mode_t;
inline void set_rounding_to_nearest() { _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);  }
inline void set_rounding_downward() { _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);  }
inline void set_rounding_upward() { _MM_SET_ROUNDING_MODE(_MM_ROUND_UP);  }
inline void set_rounding_toward_zero() { _MM_SET_ROUNDING_MODE(_MM_ROUND_TOWARD_ZERO);  }
inline void set_rounding_mode(rounding_mode_t rnd) { _MM_SET_ROUNDING_MODE(rnd); }
inline rounding_mode_t get_rounding_mode() { return _MM_GET_ROUNDING_MODE(); }
#else
typedef unsigned short rounding_mode_t;
const rounding_mode_t ROUND_NEAR = 895;
const rounding_mode_t ROUND_DOWN = 895+1024;
const rounding_mode_t ROUND_UP = 895+2048;
inline void set_rounding_to_nearest() { asm volatile ("fldcw %0" : : "m" (ROUND_NEAR) ); }
inline void set_rounding_downward() { asm volatile ("fldcw %0" : : "m" (ROUND_DOWN) ); }
inline void set_rounding_upward() { asm volatile ("fldcw %0" : : "m" (ROUND_UP) ); }
inline void set_rounding_mode(rounding_mode_t rnd) { asm volatile ("fldcw %0" : : "m" (rnd) ); }
inline rounding_mode_t get_rounding_mode() { rounding_mode_t rnd; asm volatile ("fstcw %0" : "=m" (rnd) ); return rnd; }
#endif

bool initialize_rounding_mode() { set_rounding_upward(); }
static const bool rounding_mode_initialized = initialize_rounding_mode();

static const double inf = std::numeric_limits<double>::infinity();

Float64 max_exact(Float64 x1, Float64 x2) { return Float64(std::max(x1._dbl,x2._dbl)); }
Float64 min_exact(Float64 x1, Float64 x2) { return Float64(std::min(x1._dbl,x2._dbl)); }
Float64 add_down(Float64 x1, Float64 x2) { volatile double negx1=-x1._dbl; return Float64(-(negx1-x2._dbl)); }
Float64 add_up(Float64 x1, Float64 x2) { return Float64(x1._dbl+x2._dbl); }
Float64 mul_down(Float64 x1, Float64 x2) { volatile double negx1=-x1._dbl; return Float64(-(negx1*x2._dbl)); }
Float64 mul_up(Float64 x1, Float64 x2) { return Float64(x1._dbl*x2._dbl); }
Float64 neg_exact(Float64 x) { return Float64(-x._dbl); }
Float64 rec_up(Float64 x) { return Float64(1.0/x._dbl); }
Float64 rec_down(Float64 x) { volatile double negx=-x._dbl; return Float64(-(1.0/negx)); }
Float64 half_exact(Float64 x) { return Float64(x._dbl/2); }
bool operator<(Float64 x1, Float64 x2) { return x1._dbl < x2._dbl; }
bool operator>(Float64 x1, Float64 x2) { return x1._dbl > x2._dbl; }
OutputStream& operator<<(OutputStream& os, Float64 const& x) { return os << x._dbl; }

bool operator<(Float64 x1, double x2) { return x1._dbl < x2; }
bool operator>(Float64 x1, double x2) { return x1._dbl > x2; }

template<class F> Bounds<F> generic_mul(Bounds<F> const& x1, Bounds<F> const& x2) {
    F rl=min_exact(min_exact(mul_down(x1._l,x2._l),mul_down(x1._l,x2._u)),min_exact(mul_down(x1._u,x2._l),mul_down(x1._u,x2._u)));
    F ru=max_exact(max_exact(mul_up(x1._l,x2._l),mul_up(x1._l,x2._u)),max_exact(mul_up(x1._u,x2._l),mul_up(x1._u,x2._u)));
    return Bounds<F>(rl,ru);
}

template<class F> Bounds<F> generic_rec(Bounds<F> const& x) {
    if(x._l>0.0 || x._u<0.0) { return Bounds<F>(rec_down(x._u),rec_up(x._l)); }
    else { return Bounds<F>(-inf,+inf); } }


template<class F> Bounds<F>::Bounds(double l, double u) : _l(l), _u(u) { }
template<class F> Bounds<F>::Bounds(F l, F u) : _l(l), _u(u) { }
template<class F> Bounds<F>::Bounds(PrecisionType) : _l(0.0), _u(0.0) { }
template<class F> Bounds<F>::Bounds(Rational y, PrecisionType pr) : Bounds<F>(y.get_d()-2e-16,y.get_d()+2e-16) { }
template<class F> Bounds<F>::Bounds(GenericType y, PrecisionType pr) : Bounds<F>(y.get(this->precision())) { }
//template<class F> Bounds<F>::operator GenericType() const { return make_real<Cnst,Bounds<F>>(*this); }
template<class F> Bounds<F> Bounds<F>::create(GenericType y) const { return Bounds<F>(y,this->precision()); }
template<class F> Bounds<F> Bounds<F>::create_constant(GenericType y) const { return Bounds<F>(y,this->precision()); }
template<class F> auto Bounds<F>::lower_bound() const -> DataType { return this->_l; }
template<class F> auto Bounds<F>::upper_bound() const -> DataType { return this->_u; }
template<class F> auto Bounds<F>::error_bound() const -> DataType { return half_exact(add_up(this->_u,neg_exact(this->_l))); }
template<class F> auto Bounds<F>::precision() const -> PrecisionType { return PrecisionType(); }
template class Bounds<Float64>;

OutputStream& operator<<(OutputStream& os, Bounds<Float64> const& x) { return os << BLUE << "{" << x._l << ":" << x._u << "}" << RESET; }
#ifndef TEMPLATE_FLOAT
Bounds<Float64> add(Bounds<Float64> x1, Bounds<Float64> x2) { return Bounds<Float64>(x1._l._dbl+x2._l._dbl,x1._u._dbl+x2._u._dbl); }
Bounds<Float64> mul(Bounds<Float64> x1, Bounds<Float64> x2) { return Bounds<Float64>(-inf,+inf); }
#else
template<class F> Bounds<F> add(Bounds<F> x1, Bounds<F> x2) { return Bounds<F>(x1._l._dbl+x2._l._dbl,x1._u._dbl+x2._u._dbl); }
template<class F> Bounds<F> mul(Bounds<F> x1, Bounds<F> x2) { return generic_mul(x1,x2); }
template<class F> Bounds<F> neg(Bounds<F> x) { return Bounds<F>(-x._u._dbl,x._l._dbl); }
template<class F> Bounds<F> rec(Bounds<F> x) { return generic_rec(x); }
template<class F> Bounds<F> exp(Bounds<F> x) { return Bounds<F>(std::exp(x._l._dbl),std::exp(x._u._dbl)); }
template<class F> Bounds<F> log(Bounds<F> x) { return Bounds<F>(std::log(x._l._dbl),std::log(x._u._dbl)); }
template<class F> OutputStream& operator<<(OutputStream& os, Bounds<F>const& x) { return os << "{" << x._l << ":" << x._u << "}"; }

template<class F> SizeType instantiate_floats() {
    auto add_ptr=(Bounds<F>(*)(Bounds<F>,Bounds<F>)) &add;
    auto mul_ptr=(Bounds<F>(*)(Bounds<F>,Bounds<F>)) &mul;
    auto neg_ptr=(Bounds<F>(*)(Bounds<F>)) &neg;
    auto rec_ptr=(Bounds<F>(*)(Bounds<F>)) &rec;
    auto exp_ptr=(Bounds<F>(*)(Bounds<F>)) &exp;
    auto log_ptr=(Bounds<F>(*)(Bounds<F>)) &log;
    auto write_ptr=(OutputStream&(*)(OutputStream&,Bounds<F>const&)) &operator<<;
    return (SizeType)add_ptr + (SizeType)mul_ptr
        + (SizeType)neg_ptr + (SizeType)rec_ptr
        + (SizeType)exp_ptr + (SizeType)log_ptr
        + (SizeType)write_ptr;
}

template SizeType instantiate_floats<Float64>();
template Bounds<Float64> add(Bounds<Float64>,Bounds<Float64>);
template Bounds<Float64> mul(Bounds<Float64>,Bounds<Float64>);
#endif


