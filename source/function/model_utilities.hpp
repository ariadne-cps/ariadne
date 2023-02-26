/***************************************************************************
 *            model_utilities.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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


namespace Ariadne {

namespace {

template<class FLT> struct ValidatedApproximation {
    Bounds<FLT> _v; Approximation<FLT> _a;
    ValidatedApproximation(Bounds<FLT>const& x) : _v(x), _a(x) { }
    operator Bounds<FLT> const& () const { return _v; }
    LowerBound<FLT> lower() const { return _v.lower(); }
    Approximation<FLT> middle() const { return _a; }
    UpperBound<FLT> upper() const { return _v.upper(); }
    FLT const& lower_raw() const { return _v.lower_raw(); }
    FLT const& middle_raw() const { return _a.raw(); }
    FLT const& upper_raw() const { return _v.upper_raw(); }
    friend OutputStream& operator<<(OutputStream& os, ValidatedApproximation<FLT> const& x) {
        return os << "{"<<x._v.lower()<<":"<<x._a<<":"<<x._v.upper()<<"}"; }
};
template<class FLT> Approximation<FLT> const& make_validated_approximation(Approximation<FLT> const& x) { return x; }
template<class FLT> ValidatedApproximation<FLT> make_validated_approximation(Bounds<FLT> const& x) { return ValidatedApproximation<FLT>(x); }


template<ARawFloat FLT> Rounded<FLT> const& cast_rounded(FLT const& x) { return reinterpret_cast<Rounded<FLT>const&>(x); }
template<ARawFloat FLT> Rounded<FLT>& cast_rounded(Error<FLT>& x) { return reinterpret_cast<Rounded<FLT>&>(x); }

template<class FLT> inline Void acc_err(FLT const& ml, FLT const& u, FLT& e) {
    e=add(rounded,e,hlf(add(rounded,ml,u)));
}

template<class FLT, class PRE> Ball<FLT,RawFloatType<PRE>> add(FLT const& x1, FLT const& x2, PRE pre) {
    FLT mx1=-x1;
    FLT::set_rounding_to_nearest();
    FLT r(x1.raw() + x2.raw());
    FLT::set_rounding_upward();
    FLT u=x1.raw()+x2.raw();
    FLT ml=mx1.raw()-x2.raw();
    Error<RawFloatType<PRE>> e(pre);
    e += max(u-r,ml+r);
    return Ball(r,e);
}

template<class FLT, class PRE> Ball<FLT,RawFloatType<PRE>> mul(FLT const& x1, FLT const& x2, PRE pre) {
    FLT mx1=-x1;
    FLT::set_rounding_to_nearest();
    FLT r(x1.raw() * x2.raw());
    FLT::set_rounding_upward();
    FLT u=x1.raw()*x2.raw();
    FLT ml=mx1.raw()*x2.raw();
    Error<RawFloatType<PRE>> e(pre);
    e += max(u-r,ml+r);
    return Ball(r,e);
}



template<class FLT> FLT add_err(FLT const& x1, FLT const& x2, Error<FLT>& e) {
    FLT mx1=-x1;
    FLT::set_rounding_to_nearest();
    FLT r(add(rounded,x1,x2));
    FLT::set_rounding_upward();
    FLT u=add(rounded,x1,x2);
    FLT ml=sub(rounded,mx1,x2);
    acc_err(ml,u,e.raw());
    return r;
}

template<class FLT> FLT add_err(FLT const& x, ValidatedApproximation<FLT> const& c, Error<FLT>& e) {
    FLT const& xv=x;
    FLT const& cl=c.lower_raw();
    FLT const& cm=c.middle_raw();
    FLT const& cu=c.upper_raw();
    FLT& re=e.raw();
    FLT::set_rounding_to_nearest();
    FLT rv=add(rounded,xv,cm);
    FLT::set_rounding_upward();
    FLT u=add(rounded,xv,cu);
    FLT ml=add(rounded,(-xv),cl);
    acc_err(ml,u,re);
    return FLT(rv);
}

template<class FLT> FLT add_err(FLT const& x, Bounds<FLT> const& c, Error<FLT>& e) {
    return add_err(x,ValidatedApproximation<FLT>(c),e);
}

template<class FLT> Bounds<FLT> add_err(Bounds<FLT> const& x1, Bounds<FLT> const& x2, Error<FLT>& e) {
    return add(x1,x2);
}

template<class FLT> UpperInterval<FLT> add_err(UpperInterval<FLT> const& x1, UpperInterval<FLT> const& x2, Error<FLT>& e) {
    return add(x1,x2);
}

template<class FLT> UpperInterval<FLT> add_err(UpperInterval<FLT> const& x1, Nat n2, Error<FLT>& e) {
    return add(x1,n2);
}

template<class FLT> Approximation<FLT> add_err(Approximation<FLT> const& x1, Approximation<FLT> const& x2, UnknownError<FLT>& e) {
    return add(x1,x2);
}

template<class FLT> FLT sub_err(FLT const& x1, FLT const& x2, Error<FLT>& e) {
    FLT mx1=-x1;
    FLT::set_rounding_to_nearest();
    FLT r(sub(rounded,x1.raw(),x2.raw()));
    FLT::set_rounding_upward();
    FLT u=sub(rounded,x1.raw(),x2.raw());
    FLT ml=add(rounded,mx1.raw(),x2.raw());
    acc_err(ml,u,e.raw());
    return r;
}

template<class FLT> Bounds<FLT> sub_err(Bounds<FLT> const& x1, Bounds<FLT> const& x2, Error<FLT>& e) {
    return sub(x1,x2);
}

template<class FLT> UpperInterval<FLT> sub_err(UpperInterval<FLT> const& x1, UpperInterval<FLT> const& x2, Error<FLT>& e) {
    return sub(x1,x2);
}

template<class FLT> UpperInterval<FLT> sub_err(UpperInterval<FLT> const& x1, Nat n2, Error<FLT>& e) {
    return sub(x1,n2);
}

template<class FLT> Approximation<FLT> sub_err(Approximation<FLT> const& x1, Approximation<FLT> const& x2, UnknownError<FLT>& e) {
    return sub(x1,x2);
}

template<class FLT> FLT mul_no_err(FLT const& x1, FLT const& x2) {
    FLT::set_rounding_to_nearest();
    FLT r(x1.raw() * x2.raw());
    FLT::set_rounding_upward();
    return r;
}

template<class FLT> Approximation<FLT> mul_no_err(Approximation<FLT> const& x1, Approximation<FLT> const& x2) {
    return x1*x2;
}

template<class FLT> FLT mul_err(FLT const& x1, FLT const& x2, Error<FLT>& e) {
    FLT mx1=-x1;
    FLT::set_rounding_to_nearest();
    FLT r(mul(rounded,x1,x2));
    FLT::set_rounding_upward();
    FLT u=mul(rounded,x1,x2);
    FLT ml=mul(rounded,mx1,x2);
    acc_err(ml,u,e.raw());
    return r;
}

template<class FLT> FLT mul_err(FLT const& x, ValidatedApproximation<FLT> const& c, Error<FLT>& e) {
    FLT const& xv=x;
    FLT const& cu=c.upper_raw();
    FLT const& cm=c.middle_raw();
    FLT const& cl=c.lower_raw();
    FLT& re=e.raw();
    FLT::set_rounding_to_nearest();
    FLT rv=mul(rounded,xv,cm);
    FLT::set_rounding_upward();
    if(xv>=0) {
        FLT mcl=-cl;
        FLT u=mul(rounded,xv,cu);
        FLT ml=mul(rounded,xv,mcl);
        acc_err(ml,u,re);
    } else {
        FLT mcu=-cu;
        FLT u=mul(rounded,xv,cl);
        FLT ml=mul(rounded,xv,mcu);
        acc_err(ml,u,re);
    }
    return FLT(rv);
}

template<class FLT> FLT mul_err(FLT const& x, Bounds<FLT> const& c, Error<FLT>& e) {
    return mul_err(x,ValidatedApproximation<FLT>(c),e);
}

template<class FLT> FLT mul_err(FLT const& x1, Nat n2, Error<FLT>& e) {
    return mul_err(x1,FLT(n2,x1.precision()),e);
}

template<class FLT> Bounds<FLT> mul_err(Bounds<FLT> const& x1, Bounds<FLT> const& x2, Error<FLT>& e) {
    return mul(x1,x2);
}

template<class FLT> Bounds<FLT> mul_err(Bounds<FLT> const& x1, ValidatedApproximation<FLT> const& x2, Error<FLT>& e) {
    return mul(x1,static_cast<Bounds<FLT>>(x2));
}

template<class FLT> Bounds<FLT> mul_err(Bounds<FLT> const& x1, Nat const& n2, Error<FLT>& e) {
    return mul(x1,n2);
}

template<class FLT> UpperInterval<FLT> mul_err(UpperInterval<FLT> const& x1, UpperInterval<FLT> const& x2, Error<FLT>& e) {
    return mul(x1,x2);
}

template<class FLT> UpperInterval<FLT> mul_err(UpperInterval<FLT> const& x1, Nat n2, Error<FLT>& e) {
    return mul(x1,n2);
}

template<class FLT> Approximation<FLT> mul_err(Approximation<FLT> const& x1, Approximation<FLT> const& x2, UnknownError<FLT>& e) {
    return mul(x1,x2);
}

template<class FLT> Approximation<FLT> mul_err(Approximation<FLT> const& x1, Nat n2, UnknownError<FLT>& e) {
    return mul_err(x1,Approximation<FLT>(n2,x1.precision()),e);
}

template<class FLT> FLT div_err(FLT const& x1, FLT const& x2, Error<FLT>& e) {
    FLT mx1=-x1;
    FLT::set_rounding_to_nearest();
    FLT r(div(rounded,x1.raw(),x2.raw()));
    FLT::set_rounding_upward();
    FLT u=div(rounded,x1.raw(),x2.raw());
    FLT ml=div(rounded,mx1.raw(),x2.raw());
    acc_err(ml,u,e.raw());
    return r;
}

template<class FLT> FLT div_err(FLT const& x1, Nat n2, Error<FLT>& e) {
    return div_err(x1,FLT(n2,x1.precision()),e);
}

template<class FLT> Bounds<FLT> div_err(Bounds<FLT> const& x1, Bounds<FLT> const& x2, Error<FLT>& e) {
    return div(x1,x2);
}

template<class FLT> Bounds<FLT> div_err(Bounds<FLT> const& x1, Nat n2, Error<FLT>& e) {
    return div(x1,n2);
}

template<class FLT> UpperInterval<FLT> div_err(UpperInterval<FLT> const& x1, UpperInterval<FLT> const& x2, Error<FLT>& e) {
    return div(x1,x2);
}

template<class FLT> UpperInterval<FLT> div_err(UpperInterval<FLT> const& x1, Nat n2, Error<FLT>& e) {
    return div(x1,n2);
}

template<class FLT> Approximation<FLT> div_err(Approximation<FLT> const& x1, Approximation<FLT> const& x2, UnknownError<FLT>& e) {
    return div(x1,x2);
}

template<class FLT> Approximation<FLT> div_err(Approximation<FLT> const& x1, Nat n2, UnknownError<FLT>& e) {
    return div_err(x1,Approximation<FLT>(n2,x1.precision()),e);
}




template<class FLT> FLT fma_err(FLT const& x, FLT const& y, FLT z, Error<FLT>& e) {
    FLT const& xv=x.raw();
    FLT const& yv=y.raw();
    FLT const& zv=z.raw();
    FLT& re=e.raw();
    FLT::set_rounding_to_nearest();
    FLT rv=fma(rounded,xv,yv,zv);
    FLT::set_rounding_upward();
    FLT myv=-yv;
    FLT mzv=-zv;
    FLT u=fma(rounded,xv,yv,zv);
    FLT ml=fma(rounded,xv,myv,mzv);
    acc_err(ml,u,re);
    return FLT(rv);
}

template<class FLT> FLT fma_err(ValidatedApproximation<FLT> const& c, FLT const& x, FLT y, Error<FLT>& e) {
    FLT const& xv=x.raw();
    FLT const& cu=c.upper_raw();
    FLT const& cm=c.middle_raw();
    FLT const& cl=c.lower_raw();
    FLT const& yv=y.raw();
    FLT& re=e.raw();
    FLT::set_rounding_to_nearest();
    FLT rv=fma(rounded,xv,cm,yv);
    FLT::set_rounding_upward();
    FLT u,ml;
    if(xv>=0) {
        FLT mcl=-cl;
        FLT myv=-yv;
        u=fma(xv,cu,yv);
        ml=fma(xv,mcl,myv);
    } else {
        FLT mcu=-cu;
        FLT myv=-yv;
        u=fma(xv,cl,yv);
        ml=fma(xv,mcu,myv);
    }
    acc_err(ml,u,e);
    return FLT(rv);
}

template<class FLT> FLT fma_err(Bounds<FLT> const& c, FLT const& x, FLT y, Error<FLT>& e) {
    return fma_err(ValidatedApproximation<FLT>(c),x,y,e);
}

template<class FLT> Bounds<FLT> fma_err(Bounds<FLT> const& x, Bounds<FLT> const& y, Bounds<FLT> z, Error<FLT>& e) {
    return fma(x,y,z);
}

template<class FLT> UpperInterval<FLT> fma_err(UpperInterval<FLT> const& x, UpperInterval<FLT> const& y, UpperInterval<FLT> z, Error<FLT>& e) {
    return fma(x,y,z);
}

template<class FLT> Bounds<FLT> fma_err(ValidatedApproximation<FLT> const& x, Bounds<FLT> const& y, Bounds<FLT> z, Error<FLT>& e) {
    return fma(static_cast<Bounds<FLT>>(x),y,z);
}

template<class FLT> Approximation<FLT> fma_err(Approximation<FLT> const& x, Approximation<FLT> const& y, Approximation<FLT> z, UnknownError<FLT>& e) {
    return fma(x,y,z);
}

// Returns an approximation to a1*b2+a2*b2, adding error to e
template<class FLT> FLT lin_err(FLT const& a1, FLT const& b1, FLT const& a2, FLT const& b2, Error<FLT>& e) {
    FLT::set_rounding_to_nearest();
    FLT r=cast_exact(cast_rounded(a1)*cast_rounded(b2)+cast_rounded(b1)*cast_rounded(a2));
    FLT::set_rounding_upward();
    FLT mb1=-b1;
    FLT mb2=-b2;
    Rounded<FLT> u=cast_rounded(a1)*cast_rounded(b2)+cast_rounded(b1)*cast_rounded(a2);
    Rounded<FLT> ml=cast_rounded(a1)*cast_rounded(mb2)+cast_rounded(mb1)*cast_rounded(a2);
    cast_rounded(e) += hlf(ml+u);
    return r;
}

} // namespace

} //namespace Ariadne


