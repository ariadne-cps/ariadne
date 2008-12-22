/***************************************************************************
 *            curve.code.h
 *
 *  Copyright  2007  Pieter Collins
 *  pieter.collins@cwi.nl
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
 
#include <cassert>

#include "vector.h"
#include "matrix.h"
#include "point.h"

#include "curve.h"


namespace Ariadne {


Curve::~Curve() 
{
}



Curve::Curve(const FunctionInterface& f) 
    : _function_ptr(f.clone())
{
    assert(this->_function_ptr->argument_size()==1);
}


Curve::Curve(const Curve& c) 
    : _function_ptr(c._function_ptr->clone())
{
}


Curve* 
Curve::clone() const 
{
    return new Curve(*this);
}


uint 
Curve::dimension() const 
{
    return this->_function_ptr->result_size();
}


ushort 
Curve::smoothness() const 
{
    return this->_function_ptr->smoothness();
}



Point 
Curve::value(const Float& s) const 
{
    Vector<Float> v(1,&s);
    return Point(this->_function_ptr->evaluate(v));
}


Vector< Float > 
Curve::tangent(const Float& s) const 
{
    Vector<Float> v(1,&s);
    return column(this->_function_ptr->jacobian(v),0);
}



std::ostream& 
Curve::write(std::ostream& os) const 
{
    return os << "Curve( function=" << *this->_function_ptr << " )";
}

 
}
