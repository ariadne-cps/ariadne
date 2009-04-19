/***************************************************************************
 *            test_formula.cc
 *
 *  Copyright 2009  Pieter Collins
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

#include <iostream>
#include "formula.h"

using namespace Ariadne;

int main() {


    Variable x=Variable("x");
    Variable y=Variable("y");
    Variable z=Variable("z");
    Variable w=x;

    Space spc; spc,x,y;
    std::cout<<"spc.dimension()="<<spc.dimension()<<std::endl;
    std::cout<<spc<<std::endl;

    std::cout << eval<Polynomial<Interval> >(spc,x+y*y)<< "\n";
    std::cerr << (x+y*y) << "\n";
    std::cerr << (w=(x+y*y)) << "\n";
    std::cerr << (dot(x)=(x+y*y)) << "\n";
    std::cerr << (y>(x+y*y)) << "\n";
    std::cerr << (y<x) << "\n";
    std::cerr << (y+1.0) <<std::endl;
    Reset r1=(Reset(spc,spc),(next(x)=x+y),(next(y)=x-y*y));
    std::cerr<<r1<<std::endl;
    //Reset r3={Reset(2,2),(next(x)=x+y),(next(y)=x-y*y)};
    //Reset r4(2,2),(next(x)=x+y),(next(y)=x-y*y);

    Space spc3; spc3,x,y,z;
    std::cerr<<spc3[0]<<spc3[1]<<spc3[2]<<"\n";
    std::cerr<<spc3.index(x)<<spc3.index(y)<<spc3.index(z)<<"\n";

    Dynamic d=(Dynamic(spc),(dot(x)=x+y),(dot(y)=x-y*y+1));
    std::cerr<<d<<std::endl;

    Guard g1=(Guard(spc),(x*x>y));
    std::cerr<<g1<<std::endl;

    Guard g2=(Guard(spc),(x*x<y));
    std::cerr<<g2<<std::endl;

    Guard g3=(Guard(spc3),(x*x<y));
    std::cerr<<g3<<std::endl;

    Polynomial<Interval> p=eval<Polynomial<Interval> >(spc3,1.5+0.375*y-x*x+z-2*x*y*z-3*x*x*z);
    std::cerr<<p<<" "<<spc3<<std::endl;
    std::cerr<<std::make_pair(p,spc3)<<std::endl;
};