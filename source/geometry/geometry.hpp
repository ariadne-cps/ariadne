/***************************************************************************
 *            geometry/geometry.hpp
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

/*! \file geometry/geometry.hpp
 *  \brief Geometric operations on abstract sets.
 */
#ifndef ARIADNE_GEOMETRY_HPP
#define ARIADNE_GEOMETRY_HPP

#include <vector>

#include "../geometry/box.hpp"
#include "../function/function_interface.hpp"

namespace Ariadne {

enum class SplitPart : char;


inline
SizeType irmax(const ExactBoxType& bx) {
    FloatDP dmax=0.0;
    Nat imax=0;
    for(Nat i=0; i!=bx.size(); ++i) {
        FloatDP d=bx[i].width().upper().raw();
        if(d>dmax) {
            imax=i;
            dmax=d;
        }
    }
    return imax;
}

inline FloatDPValue mid(FloatDPValue l, FloatDPValue u) {
    return FloatDPValue(med(approx,l.raw(),u.raw()));
}

inline
ExactBoxType split(const ExactBoxType& bx, SizeType i, SplitPart lr) {
    ExactBoxType result(bx);
    ExactIntervalType& ivl=result[i];
    const FloatDPValue& l=ivl.lower();
    const FloatDPValue& u=ivl.upper();
    FloatDPValue c=mid(l,u);
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
    FloatDPValue c=mid(bx[i].lower(),bx[i].upper());
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
ValidatedKleenean
image_separated(const ExactBoxType& d, const F& f, const ExactBoxType& b, const RawFloatDP& eps)
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
    } else if(d.radius().upper().raw()<eps) {
        //cout << "radius limit reached\n";
        return indeterminate;
    } else {
        Nat i=irmax(d);
        //cout << "splitting\n";
        return separated(split(d,i,SplitPart::LOWER),f,b,eps) && separated(split(d,i,SplitPart::UPPER),f,b,eps);
    }
}


template<class F>
ValidatedKleenean
image_inside(const ExactBoxType& d, const F& f, const ExactBoxType& b, const RawFloatDP& eps)
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
ValidatedKleenean overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const FloatDP& eps);

//! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
ValidatedKleenean inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const FloatDP& eps);

//! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
ValidatedKleenean separated(const LocatedSetInterface& ls, const RegularSetInterface& rs, const FloatDP& eps);


//! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
ValidatedLowerKleenean overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops, const FloatDP& eps);

//! \brief Tests if \a cps is a inside of \a ops, to a tolerance of \a eps.
ValidatedLowerKleenean inside(const CompactSetInterface& cps, const OpenSetInterface& ops, const FloatDP& eps);

//! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
ValidatedLowerKleenean separated(const CompactSetInterface& cps, const ClosedSetInterface& cls, const FloatDP& eps);




//! \brief Tests if the intersection of \a ls and \a bx overlaps \a rs, to a tolerance of \a eps.
ValidatedKleenean overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const ExactBoxType& bx, const FloatDP& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
ValidatedKleenean inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const ExactBoxType& bx, const FloatDP& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
ValidatedKleenean separated(const LocatedSetInterface& ls, const RegularSetInterface& rs, const ExactBoxType& bx, const FloatDP& eps);


//! \brief Tests if the intersection of \a ls and \a bx overlaps \a rs, to a tolerance of \a eps.
ValidatedLowerKleenean intersection_overlap(const OvertSetInterface& ls, const OpenSetInterface& rs, const ExactBoxType& bx, const FloatDP& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
ValidatedSierpinskian intersection_inside(const ClosedSetInterface& ls, const OpenSetInterface& rs, const ExactBoxType& bx, const FloatDP& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
ValidatedSierpinskian intersection_separated(const ClosedSetInterface& ls, const ClosedSetInterface& rs, const ExactBoxType& bx, const FloatDP& eps);


} // namespace Ariadne

#endif // ARIADNE_GEOMETRY_HPP
