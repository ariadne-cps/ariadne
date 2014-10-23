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

enum Piece { left=0, right=1, middle=2 };

inline
uint irmax(const ExactBox& bx) {
    Float dmax=0.0;
    uint imax=0;
    for(uint i=0; i!=bx.size(); ++i) {
        Float d=bx[i].width().raw();
        if(d>dmax) {
            imax=i;
            dmax=d;
        }
    }
    return imax;
}

inline ExactFloat mid(ExactFloat l, ExactFloat u) {
    return ExactFloat(med_approx(l.raw(),u.raw()));
}

inline
ExactBox split(const ExactBox& bx, uint i, Piece lr) {
    ExactBox result(bx);
    ExactInterval& ivl=result[i];
    const ExactFloat& l=ivl.lower();
    const ExactFloat& u=ivl.upper();
    ExactFloat c=mid(l,u);
    if(lr==middle) {
        ivl.set(mid(l,c),mid(c,u));
    } else {
        if(lr==left) { ivl.set_upper(c); }
        else { ivl.set_lower(c); }
    }
    return result;
}

inline
std::pair<ExactBox,ExactBox> split(const ExactBox& bx, uint i)
{
    std::pair<ExactBox,ExactBox> result(bx,bx);
    ExactFloat c=mid(bx[i].lower(),bx[i].upper());
    result.first[i].set_upper(c);
    result.second[i].set_lower(c);
    return result;
}

inline
ExactBox split(const ExactBox& bx, Piece lr)
{
    uint i=irmax(bx);
    return split(bx,i,lr);
}

inline
std::pair<ExactBox,ExactBox> split(const ExactBox& bx) {
    return split(bx,irmax(bx));
}


template<class F>
Tribool
separated(const ExactBox& d, const F& f, const ExactBox& b, const RawFloat& eps)
{

    ExactBox fd=f.evaluate(d);
    ExactBox fc=f.evaluate(ExactBox(midpoint(d)));

    //cout << "called with " << d << ", having fd=" << fd << " and fc=" << fc << endl;
    if(disjoint(fd,b)) {
    	//cout << "evaluation is disjoint\n";
        return true;
    } else if(definitely(inside(fc,b))) {
    	//cout << "evaluation of the midpoint is inside\n";
        return false;
    } else if(d.radius().raw()<eps) {
    	//cout << "radius limit reached\n";
        return indeterminate;
    } else {
        uint i=irmax(d);
        //cout << "splitting\n";
        return separated(split(d,i,left),f,b,eps) && separated(split(d,i,right),f,b,eps);
    }
}


template<class F>
Tribool
inside(const ExactBox& d, const F& f, const ExactBox& b, const RawFloat& eps)
{

    ExactBox fd=f.evaluate(d);
    ExactBox fc=f.evaluate(ExactBox(midpoint(d)));

    //cout << "called with " << d << ", having fd=" << fd << " and fc=" << fc << endl;
    if(disjoint(fc,b)) {
    	//cout << "evaluation of the midpoint is disjoint\n";
        return false;
    } else if(definitely(inside(fd,b))) {
    	//cout << "evaluation is definitely inside\n";
        return true;
    } else if(d.radius().raw()<eps) {
    	//cout << "radius limit reached\n";
        return indeterminate;
    } else {
        uint i=irmax(d);
        //cout << "splitting\n";
        return inside(split(d,i,left),f,b,eps) && inside(split(d,i,right),f,b,eps);
    }
}

template<class DS>
DS remove_subsets(const DS& ls)
{
    DS result;
    for(uint i=0; i!=ls.size(); ++i) {
        for(uint j=0; j!=ls.size(); ++j) {
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
    for(uint i=0; i!=ls.size(); ++i) {
        for(uint j=0; j!=ls.size(); ++j) {
            if(inside(ls[j],ls[i])) {
                break;
            }
        }
        result.adjoin(ls[i]);
    }
}


//! \brief Tests if \a ls overlaps \a rs, to a tolerance of \a eps.
Tribool overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);

//! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
Tribool inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);

//! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
Tribool separated(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);


//! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
Tribool overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops, const Float& eps);

//! \brief Tests if \a cps is a inside of \a ops, to a tolerance of \a eps.
Tribool inside(const CompactSetInterface& cps, const OpenSetInterface& ops, const Float& eps);

//! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
Tribool separated(const CompactSetInterface& cps, const ClosedSetInterface& cls, const Float& eps);




//! \brief Tests if the intersection of \a ls and \a bx overlaps \a rs, to a tolerance of \a eps.
Tribool overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const ExactBox& bx, const Float& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
Tribool inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const ExactBox& bx, const Float& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
Tribool separated(const LocatedSetInterface& ls, const RegularSetInterface& rs, const ExactBox& bx, const Float& eps);


//! \brief Tests if the intersection of \a ls and \a bx overlaps \a rs, to a tolerance of \a eps.
Tribool overlap(const OvertSetInterface& ls, const OpenSetInterface& rs, const ExactBox& bx, const Float& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
Tribool inside(const ClosedSetInterface& ls, const OpenSetInterface& rs, const ExactBox& bx, const Float& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
Tribool separated(const ClosedSetInterface& ls, const ClosedSetInterface& rs, const ExactBox& bx, const Float& eps);


} // namespace Ariadne

#endif // ARIADNE_GEOMETRY_H
