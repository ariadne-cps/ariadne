/***************************************************************************
 *            geometry.h
 *
 *  Copyright 2008  Pieter Collins
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

/*! \file geometry.h
 *  \brief Geometric operations on abstract sets.
 */
#ifndef ARIADNE_GEOMETRY_H
#define ARIADNE_GEOMETRY_H

#include <vector>

#include "geometry/box.h"
#include "function/function_interface.h"

namespace Ariadne {

enum class SplitPart : char;


inline
SizeType irmax(const ExactBoxType& bx) {
    Float64 dmax=0.0;
    Nat imax=0;
    for(Nat i=0; i!=bx.size(); ++i) {
        Float64 d=bx[i].width().upper().raw();
        if(d>dmax) {
            imax=i;
            dmax=d;
        }
    }
    return imax;
}

inline Float64Value mid(Float64Value l, Float64Value u) {
    return Float64Value(med_approx(l.raw(),u.raw()));
}

inline
ExactBoxType split(const ExactBoxType& bx, SizeType i, SplitPart lr) {
    ExactBoxType result(bx);
    ExactIntervalType& ivl=result[i];
    const Float64Value& l=ivl.lower();
    const Float64Value& u=ivl.upper();
    Float64Value c=mid(l,u);
    if(lr==SplitPart::MIDDLE) {
        ivl.set(mid(l,c),mid(c,u));
    } else {
        if(lr==SplitPart::LOWER) { ivl.set_upper(c); }
        else { ivl.set_lower(c); }
    }
    return result;
}

inline
Pair<ExactBoxType,ExactBoxType> split(const ExactBoxType& bx, SizeType i)
{
    Pair<ExactBoxType,ExactBoxType> result(bx,bx);
    Float64Value c=mid(bx[i].lower(),bx[i].upper());
    result.first[i].set_upper(c);
    result.second[i].set_lower(c);
    return result;
}

inline
ExactBoxType split(const ExactBoxType& bx, SplitPart lr)
{
    SizeType i=irmax(bx);
    return split(bx,i,lr);
}

inline
Pair<ExactBoxType,ExactBoxType> split(const ExactBoxType& bx) {
    return split(bx,irmax(bx));
}


template<class F>
Kleenean
image_separated(const ExactBoxType& d, const F& f, const ExactBoxType& b, const RawFloat64& eps)
{

    ExactBoxType fd=f.evaluate(d);
    ExactBoxType fc=f.evaluate(ExactBoxType(midpoint(d)));

    //cout << "called with " << d << ", having fd=" << fd << " and fc=" << fc << endl;
    if(disjoint(fd,b)) {
        //cout << "evaluation is disjoint\n";
        return true;
    } else if(definitely(inside(fc,b))) {
        //cout << "evaluation of the midpoint is inside\n";
        return false;
    } else if(d.radius().upper()<cast_exact(eps)) {
        //cout << "radius limit reached\n";
        return indeterminate;
    } else {
        Nat i=irmax(d);
        //cout << "splitting\n";
        return separated(split(d,i,SplitPart::LOWER),f,b,eps) && separated(split(d,i,SplitPart::UPPER),f,b,eps);
    }
}


template<class F>
Kleenean
image_inside(const ExactBoxType& d, const F& f, const ExactBoxType& b, const RawFloat64& eps)
{

    ExactBoxType fd=f.evaluate(d);
    ExactBoxType fc=f.evaluate(ExactBoxType(midpoint(d)));

    //cout << "called with " << d << ", having fd=" << fd << " and fc=" << fc << endl;
    if(disjoint(fc,b)) {
        //cout << "evaluation of the midpoint is disjoint\n";
        return false;
    } else if(definitely(inside(fd,b))) {
        //cout << "evaluation is definitely inside\n";
        return true;
    } else if(d.radius().upper().raw()<eps) {
        //cout << "radius limit reached\n";
        return indeterminate;
    } else {
        Nat i=irmax(d);
        //cout << "splitting\n";
        return inside(split(d,i,SplitPart::LOWER),f,b,eps) && inside(split(d,i,SplitPart::UPPER),f,b,eps);
    }
}

template<class DS>
DS remove_subsets(const DS& ls)
{
    DS result;
    for(Nat i=0; i!=ls.size(); ++i) {
        for(Nat j=0; j!=ls.size(); ++j) {
            if(inside(ls[i],ls[j])) {
                break;
            }
        }
        result.adjoin(ls[i]);
    }
}

template<class DS>
DS remove_supersets(const DS& ls)
{
    DS result;
    for(Nat i=0; i!=ls.size(); ++i) {
        for(Nat j=0; j!=ls.size(); ++j) {
            if(inside(ls[j],ls[i])) {
                break;
            }
        }
        result.adjoin(ls[i]);
    }
}


//! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
Kleenean overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float64& eps);

//! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
Kleenean inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float64& eps);

//! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
Kleenean separated(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float64& eps);


//! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
Sierpinski overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops, const Float64& eps);

//! \brief Tests if \a cps is a inside of \a ops, to a tolerance of \a eps.
Sierpinski inside(const CompactSetInterface& cps, const OpenSetInterface& ops, const Float64& eps);

//! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
Sierpinski separated(const CompactSetInterface& cps, const ClosedSetInterface& cls, const Float64& eps);




//! \brief Tests if the intersection of \a ls and \a bx overlaps \a rs, to a tolerance of \a eps.
Kleenean overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const ExactBoxType& bx, const Float64& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
Kleenean inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const ExactBoxType& bx, const Float64& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
Kleenean separated(const LocatedSetInterface& ls, const RegularSetInterface& rs, const ExactBoxType& bx, const Float64& eps);


//! \brief Tests if the intersection of \a ls and \a bx overlaps \a rs, to a tolerance of \a eps.
Sierpinski intersection_overlap(const OvertSetInterface& ls, const OpenSetInterface& rs, const ExactBoxType& bx, const Float64& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
Sierpinski intersection_inside(const ClosedSetInterface& ls, const OpenSetInterface& rs, const ExactBoxType& bx, const Float64& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
Sierpinski intersection_separated(const ClosedSetInterface& ls, const ClosedSetInterface& rs, const ExactBoxType& bx, const Float64& eps);


} // namespace Ariadne

#endif // ARIADNE_GEOMETRY_H
