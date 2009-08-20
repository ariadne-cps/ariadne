/***************************************************************************
 *            function_set.cc
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
 
#include "macros.h"
#include "function.h"
#include "function_set.h"
#include "geometry.h"

namespace Ariadne {


ImageSet::ImageSet()
    : _domain(), _function_ptr(new IdentityFunction(0))
{ 
}

ImageSet::ImageSet(const Vector<Interval>& dom)
    //  : _domain(dom), _function_ptr(new IdentityFunction(dom.size()))
    : _domain(Vector<Interval>(dom.size(),Interval(-1,1))), 
      _function_ptr(new ScalingFunction(dom))
{ 
}


ImageSet::ImageSet(const Vector<Interval>& dom, const FunctionInterface& fn)
    : _domain(dom), _function_ptr(fn.clone()) 
{ 
    ARIADNE_ASSERT(dom.size()==fn.argument_size()); 
}


ImageSet*
ImageSet::clone() const 
{ 
    return new ImageSet(*this);
}


uint
ImageSet::dimension() const 
{ 
    return this->_function_ptr->result_size();
}


tribool
ImageSet::empty() const 
{ 
    return this->_domain.empty();
}


tribool
ImageSet::disjoint(const Box& bx) const 
{ 
    if(dynamic_cast<const IdentityFunction*>(&*this->_function_ptr)) {
        return Ariadne::disjoint(bx,this->_domain);
    } else if(dynamic_cast<const ScalingFunction*>(&*this->_function_ptr)) {
        return Ariadne::disjoint(bx,this->_function_ptr->evaluate(this->_domain));
    } else {
        static const int MAX_SUBDIVISIONS=8;
        return Ariadne::disjoint(this->domain(),*this->_function_ptr,bx,radius(bx)/MAX_SUBDIVISIONS);
    }
}


tribool
ImageSet::overlaps(const Box& bx) const 
{ 
    static const int MAX_SUBDIVISIONS=8;
    return !Ariadne::disjoint(this->domain(),*this->_function_ptr,bx,radius(bx)/MAX_SUBDIVISIONS);
}


tribool
ImageSet::inside(const Box& bx) const 
{ 
    if(Ariadne::inside(this->bounding_box(),bx)) {
        return true;
    } else {
        return indeterminate;
    }
}


Box
ImageSet::bounding_box() const 
{ 
    return this->_function_ptr->evaluate(this->_domain);
}


void
ImageSet::draw(CanvasInterface& os) const
{ 
    return this->bounding_box().draw(os);
}


std::ostream&
ImageSet::write(std::ostream& os) const 
{ 
    ARIADNE_NOT_IMPLEMENTED;
}





ConstraintSet::ConstraintSet(const Vector<Interval>& codom, const FunctionInterface& fn) 
    : _codomain(codom), _function_ptr(fn.clone()) 
{ 
    ARIADNE_ASSERT(codom.size()==fn.result_size());
}

ConstraintSet*
ConstraintSet::clone() const 
{ 
    return new ConstraintSet(*this);
}


uint
ConstraintSet::dimension() const 
{ 
    return this->_function_ptr->argument_size();
}


tribool
ConstraintSet::disjoint(const Box& bx) const 
{ 
    ARIADNE_NOT_IMPLEMENTED;
}


tribool
ConstraintSet::overlaps(const Box& bx) const 
{ 
    ARIADNE_NOT_IMPLEMENTED;
}


tribool
ConstraintSet::covers(const Box& bx) const 
{ 
    ARIADNE_NOT_IMPLEMENTED;
}


std::ostream&
ConstraintSet::write(std::ostream& os) const 
{ 
    ARIADNE_NOT_IMPLEMENTED;
}



} // namespace Ariadne
