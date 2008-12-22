/***************************************************************************
 *            geometry.cc
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
 
#include "geometry.h"

namespace Ariadne {


tribool 
disjoint(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps)
{
    Box bb=ls.bounding_box();
    if(bb.empty()) { return true; }
    return disjoint(ls,rs,Box(ls.bounding_box()),eps);
}


tribool 
overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps)
{
    Box bb=ls.bounding_box();
    if(bb.empty()) { return false; }
    return overlap(ls,rs,ls.bounding_box(),eps);
}


tribool 
inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Float& eps)
{
    Box bb=ls.bounding_box();
    if(bb.empty()) { return true; }
    return inside(ls,rs,ls.bounding_box(),eps);
}


tribool 
overlap(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Box& bx, const Float& eps)
{
    if(ls.disjoint(bx)) { 
        return false; 
    }
    if(rs.disjoint(bx)) { 
        return false; 
    }
    else if(rs.covers(bx)) {
        return true; 
    }
    else if(bx.radius()<eps) {
        return indeterminate;
    } else {
        Box bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(ls.disjoint(bx1)) {
            return overlap(ls,rs,bx2,eps);
        } else if(ls.disjoint(bx2)) {
            return overlap(ls,rs,bx1,eps);
        } else {
            return overlap(ls,rs,bx1,eps) || overlap(ls,rs,bx2,eps);
        }
    }
}
    
  
tribool 
inside(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Box& bx, const Float& eps)
{
    if(ls.disjoint(bx) || rs.covers(bx)) { 
        return true; 
    } else if(bx.radius()<eps) {
        return indeterminate;
    } else {
        Box bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(ls.disjoint(bx1)) {
            return inside(ls,rs,bx2,eps);
        } else if(ls.disjoint(bx2)) {
            return inside(ls,rs,bx1,eps);
        } else {
            return inside(ls,rs,bx1,eps) && inside(ls,rs,bx2,eps);
        }
    }
}
 
   
tribool 
disjoint(const LocatedSetInterface& ls, const RegularSetInterface& rs, const Box& bx, const Float& eps)
{
    if(ls.disjoint(bx) || rs.disjoint(bx)) { 
        return true; 
    } else if(bx.radius()<eps) {
        return indeterminate;
    } else {
        Box bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(ls.disjoint(bx1)) {
            return disjoint(ls,rs,bx2,eps);
        } else if(ls.disjoint(bx2)) {
            return disjoint(ls,rs,bx1,eps);
        } else {
            return disjoint(ls,rs,bx1,eps) && disjoint(ls,rs,bx2,eps);
        }
    }
}
    
  


tribool 
overlap(const OvertSetInterface& ovs, const OpenSetInterface& ops, const Box& bx, const Float& eps)
{
    if(ovs.overlaps(bx)) {
        if(ops.covers(bx)) { 
            return true; 
        } else if(bx.radius()<eps) {
            return indeterminate;
        } else {
            Box bx1,bx2;
            make_lpair(bx1,bx2)=split(bx);
            if(overlap(ovs,ops,bx1,eps)) {
                return true;
            } else {
                return overlap(ovs,ops,bx2,eps);
            }
        }
    } else {
        return indeterminate;
    }
}
    
  
tribool 
inside(const ClosedSetInterface& cls, const OpenSetInterface& ops, const Box& bx, const Float& eps)
{
    if(cls.disjoint(bx) || ops.covers(bx)) { 
        return true; 
    } else if(bx.radius()<eps) {
        return indeterminate;
    } else {
        Box bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(inside(cls,ops,bx1,eps)) {
            return inside(cls,ops,bx2,eps);
        } else {
            return indeterminate;
        }
    }
}

    
tribool 
disjoint(const ClosedSetInterface& cls1, const ClosedSetInterface& cls2, const Box& bx, const Float& eps)
{
    if(cls1.disjoint(bx) || cls2.disjoint(bx)) {
        return true;
    } else if(bx.radius()<eps) {
        return indeterminate;
    } else {
        Box bx1,bx2;
        make_lpair(bx1,bx2)=split(bx);
        if(disjoint(cls1,cls2,bx1,eps)) {
            return disjoint(cls1,cls2,bx2,eps);
        } else {
            return indeterminate;
        }
    }
}

    


} // namespace Ariadne

