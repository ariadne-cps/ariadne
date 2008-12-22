/***************************************************************************
 *            box.cc
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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
 
#include <sstream>
#include <string>
#include <vector>

#include <cstdarg>

#include "box.h"
#include "stlio.h"

namespace Ariadne {

Box::Box(uint d, const Float& x0l, const Float& x0u, ...)
    : Vector<Interval>(d)
{
    assert(d>=1); 
    va_list args; 
    va_start(args,x0u);
    (*this)[0]=Interval(x0l,x0u);
    for(uint i=1; i!=d; ++i) { 
        Float xil=va_arg(args,Float);
        Float xiu=va_arg(args,Float);
        (*this)[i]=Interval(xil,xiu); 
    } 
    va_end(args);
}

Box::Box(const std::string& str) 
{
    *this=make_box(str);
}


Box make_box(const std::string& str)
{
    // Representation as a literal 
    //   "[a1,b1]x[a2,b2]x...x[an,bn]" 

    std::stringstream ss(str);
    std::vector<Interval> vec;
    Interval ivl; 
    char c;

    c='x';
    while(c=='x') {
        ss >> ivl;
        vec.push_back(ivl);
        c=' ';
        while( ss && c==' ') {
            ss >> c;
        }
    }
    if(ss) {
        ss.putback(c);
    }
    return Box(vec.size(),&vec[0]);
}

} //namespace Ariadne
