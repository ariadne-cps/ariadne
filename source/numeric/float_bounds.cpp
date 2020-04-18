/***************************************************************************
 *            numeric/float_bounds.cpp
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

#include "float_bounds.hpp"
#include "float_bounds.inl.hpp"
#include "float_bounds.tpl.hpp"

#include "floatdp.hpp"
#include "floatmp.hpp"

namespace Ariadne {

Value<FloatDP> midpoint(Bounds<FloatDP> const& x) { return x.value(); } // DEPRECATED

template<> auto Operations<FloatBounds<DoublePrecision>>::_write(OutputStream& os, const FloatBounds<DoublePrecision>& x) -> OutputStream&;
template<> auto Operations<FloatBounds<MultiplePrecision>>::_write(OutputStream& os, const FloatBounds<MultiplePrecision>& x) -> OutputStream&;

template class Bounds<FloatDP>;
template class Operations<Bounds<FloatDP>>;
template class Bounds<FloatMP>;
template class Operations<Bounds<FloatMP>>;

template<> String class_name<Bounds<FloatDP>>() { return "FloatDPBounds"; }
template<> String class_name<Bounds<FloatMP>>() { return "FloatMPBounds"; }


template<> OutputStream& Operations<FloatBounds<MultiplePrecision>>::_write(OutputStream& os, const FloatBounds<MultiplePrecision>& x)
{
    static const double log2ten = 3.3219280948873621817;
    using std::max; using std::min;
    FloatMP const& l=x.lower_raw();
    FloatMP const& u=x.upper_raw();
    double ldbl=l.get_d();
    double udbl=u.get_d();
    if(ldbl==0.0 && udbl==0.0) { return os << "0.0[:]"; }

    int errplc=static_cast<int>(FloatError<MultiplePrecision>::output_places);
    //int bndplc=FloatBounds<MultiplePrecision>::output_places;
    int precplc=x.precision()/log2ten;
    int log10wdth=abslog10floor(udbl-ldbl);
    int log10mag=abslog10floor(max(-ldbl,udbl));
    int dgtswdth=errplc-(log10wdth+1); // Digits appropriate given width of interval
    //int dgtsbnd=bndplc-(log10mag+1); // Digits appropriate given asked-for precision of bounded objects
    int dgtsprec=precplc-(log10mag+1); // Digits appropriate given precision of objects
    Nat dgts=static_cast<Nat>(max(min(dgtswdth,dgtsprec),1));
    DecimalPlaces plcs{dgts}
;
    String lstr=print(l,plcs,MPFR_RNDD);
    String ustr=print(u,plcs,MPFR_RNDU);
    auto lcstr=lstr.c_str();
    auto ucstr=ustr.c_str();
    size_t cpl=0;
    if(ldbl*udbl>=0 && abslog10floor(ldbl)==abslog10floor(udbl)) {
        while(lcstr[cpl]!='\0' && lcstr[cpl]==ustr[cpl]) { ++cpl; }

    }

    char ocstr[1024];
    ocstr[0]='\0';
    strncat(ocstr,lcstr,cpl);
    strcat(ocstr,"[");
    strcat(ocstr,lcstr+cpl);
    strcat(ocstr,":");
    strcat(ocstr,ucstr+cpl);
    strcat(ocstr,"]");
    return os << ocstr;
}


template<> OutputStream& Operations<FloatBounds<DoublePrecision>>::_write(OutputStream& os, const FloatBounds<DoublePrecision>& x)
{
    MultiplePrecision prec(64);
    return os << FloatBounds<MultiplePrecision>(FloatMP(x.lower_raw(),prec),FloatMP(x.upper_raw(),prec));
}


} // namespace Ariadne
