#include <limits>
#include <cmath>
#include "simple.h"

// ---------------- Float ------------------------------------------------------------------------------------------ //

static const double inf = std::numeric_limits<double>::infinity();

OutputStream& operator<<(OutputStream& os, Float64 const& x) { return os << x._dbl; }

template<class OP, class... AS> Real make_real(AS... as);
template<> Real make_real<Cnst,Float64Bounds>(Float64Bounds);

template<class F> Bounds<F>::Bounds(double l, double u) : _l(l), _u(u) { }
template<class F> Bounds<F>::Bounds(F l, F u) : _l(l), _u(u) { }
template<class F> Bounds<F>::Bounds(PrecisionType) : _l(0.0), _u(0.0) { }
template<class F> Bounds<F>::Bounds(Rational y, PrecisionType pr) : Bounds<F>(y.get_d()-2e-16,y.get_d()+2e-16) { }
template<class F> Bounds<F>::Bounds(GenericType y, PrecisionType pr) : Bounds<F>(y.get(this->precision())) { }
template<class F> Bounds<F>::operator GenericType() const { return make_real<Cnst,Bounds<F>>(*this); }
template<class F> Bounds<F> Bounds<F>::create(GenericType y) const { return Bounds<F>(y,this->precision()); }
template<class F> auto Bounds<F>::lower_bound() const -> DataType { return this->_l; }
template<class F> auto Bounds<F>::upper_bound() const -> DataType { return this->_u; }
template<class F> auto Bounds<F>::precision() const -> PrecisionType { return PrecisionType(); }
template class Bounds<Float64>;

OutputStream& operator<<(OutputStream& os, Bounds<Float64> const& x) { return os << "{" << x._l << ":" << x._u << "}"; }
#ifndef TEMPLATE_FLOAT
Bounds<Float64> add(Bounds<Float64> x1, Bounds<Float64> x2) { return Bounds<Float64>(x1._l._dbl+x2._l._dbl,x1._u._dbl+x2._u._dbl); }
Bounds<Float64> mul(Bounds<Float64> x1, Bounds<Float64> x2) { return Bounds<Float64>(-inf,+inf); }
#else
template<class F> Bounds<F> add(Bounds<F> x1, Bounds<F> x2) { return Bounds<F>(x1._l._dbl+x2._l._dbl,x1._u._dbl+x2._u._dbl); }
template<class F> Bounds<F> mul(Bounds<F> x1, Bounds<F> x2) { return Bounds<F>(-inf,+inf); }
template<class F> Bounds<F> neg(Bounds<F> x) { return Bounds<F>(-x._u._dbl,x._l._dbl); }
template<class F> Bounds<F> rec(Bounds<F> x) { return Bounds<F>(-inf,+inf); }
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


