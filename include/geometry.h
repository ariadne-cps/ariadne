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
#define ARIANDE_GEOMETRY_H

#include <vector>

#include "box.h"
#include "function_interface.h"

namespace Ariadne {


inline
uint irmax(const Box& bx) {
    Float dmax=0.0;
    uint imax=0;
    for(uint i=0; i!=bx.size(); ++i) {
        Float d=bx[i].width();
        if(d>dmax) {
            imax=i;
            dmax=d;
        }
    }
    return imax;
}
 
inline
Box split(const Box& bx, uint i, bool lr) {
    Box result(bx);
    Float c=med_approx(bx[i].lower(),bx[i].upper());
    if(lr==false) { result[i].u=c; }
    else { result[i].l=c; }
    return result;
}

inline
std::pair<Box,Box> split(const Box& bx, uint i)
{
    std::pair<Box,Box> result(bx,bx);
    Float c=med_approx(bx[i].lower(),bx[i].upper());
    result.first[i].u=c;
    result.second[i].l=c;
    return result;
}

inline
Box split(const Box& bx, bool lr)
{
    uint i=irmax(bx);
    return split(bx,i,lr);
}

inline
std::pair<Box,Box> split(const Box& bx) {
    return split(bx,irmax(bx));
}


template<class F>
tribool 
disjoint(const Box& d, const F& f, const Box& b, const Float& eps)
{
    Box fd=f.evaluate(d);
    Box fc=f.evaluate(Box(midpoint(d)));
    if(disjoint(fd,b)) { 
        return true;
    } else if(inside(fc,b)) {
        return false;
    } else if(d.radius()<eps) {
        return indeterminate;
    } else {
        uint i=irmax(d);
        return disjoint(split(d,i,0),f,b,eps) || disjoint(split(d,i,0),f,b,eps);
    }
}

template<class F>
tribool 
inside(const Box& d, const F& f, const Box& b, const Float& eps)
{
    Box fd=f(d);
    Box fc=f(Box(midpoint(d)));
    if(inside(fc,b)) { 
        return true;
    } else if(disjoint(fd,b)) {
        return false;
    } else if(d.radius()<eps) {
        return indeterminate;
    } else {
        uint i=irmax(d);
        return inside(split(d,i,0),f,b,eps) && inside(split(d,i,0),f,b,eps);
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
tribool overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);

//! \brief Tests if \a ls is a inside of \a rs, to a tolerance of \a eps.
tribool inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);

//! \brief Tests if \a ls is disjoint from \a rs, to a tolerance of \a eps.
tribool disjoint(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps);


//! \brief Tests if \a ovs overlaps \a ops, to a tolerance of \a eps.
tribool overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops, const Float& eps);

//! \brief Tests if \a cps is a inside of \a ops, to a tolerance of \a eps.
tribool inside(const CompactSetInterface& cps, const OpenSetInterface& ops, const Float& eps);

//! \brief Tests if \a cps is disjoint from \a cls, to a tolerance of \a eps.
tribool disjoint(const CompactSetInterface& cps, const ClosedSetInterface& cls, const Float& eps);




//! \brief Tests if the intersection of \a ls and \a bx overlaps \a rs, to a tolerance of \a eps.
tribool overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Box& bx, const Float& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
tribool inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Box& bx, const Float& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
tribool disjoint(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Box& bx, const Float& eps);


//! \brief Tests if the intersection of \a ls and \a bx overlaps \a rs, to a tolerance of \a eps.
tribool overlap(const OvertSetInterface& ls, const OpenSetInterface& rs, const Box& bx, const Float& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
tribool inside(const ClosedSetInterface& ls, const OpenSetInterface& rs, const Box& bx, const Float& eps);

//! \brief Tests if the intersection of \a ls and \a bx is a inside of \a rs, to a tolerance of \a eps.
tribool disjoint(const ClosedSetInterface& ls, const ClosedSetInterface& rs, const Box& bx, const Float& eps);


} // namespace Ariadne

#endif // ARIADNE_GEOMETRY_H
