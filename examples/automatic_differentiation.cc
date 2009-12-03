/***************************************************************************
 *            automatic_differentiation.cc
 *
 *  Copyright  2009  Pieter Collins
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

#include <cstdarg>
#include "ariadne.h"

using namespace Ariadne;

int main() {

    const uint dimension = 3;
    const uint degree = 2;
    // Construct identity function of the given degree
    VectorFunction x=VectorFunction::identity(3);

    // Construct differential object to compute derivatives at zero in three variables
    Vector< Differential<Float> > dx=Differential<Float>::variables(degree,Vector<Float>(dimension,0.0));
    std::cout << "dx=" << dx << "\n";

    // Construct a vector-valued function componentwise
    VectorFunction f((x[0]+2*x[1]+3*x[2]+5*x[0]*x[0],7.0+x[0]+x[1]*x[2]));
    std::cout << "f=" << f << "\n";

    // Compute the derivatives by calling the function on the differential vector.
    Vector< Differential<Float> > df=f(dx);
    std::cout << "df=" << df << "\n";

    // Use multi-index construnctor MultiIndex(n,a0,a1,...,an-1)
    std::cout << "df[0][0,0,0]=" << df[0][MultiIndex(3,0,0,0)] << "\n";
    std::cout << "df[0][1,0,0]=" << df[0][MultiIndex(3,1,0,0)] << "\n";
    std::cout << "df[0][0,1,0]=" << df[0][MultiIndex(3,0,1,0)] << "\n";
    std::cout << "df[0][0,0,1]=" << df[0][MultiIndex(3,0,0,1)] << "\n";
    std::cout << "df[0][2,0,0]=" << df[0][MultiIndex(3,2,0,0)] << "\n";
    std::cout << "df[0][1,1,0]=" << df[0][MultiIndex(3,1,1,0)] << "\n";
    std::cout << "df[1][0,0,0]=" << df[0][MultiIndex(3,0,0,0)] << "\n";
    std::cout << std::endl;

    // Construct a scalar-valued function componentwise
    ScalarFunction g(x[1]*x[1]);
    // Compute the derivatives by calling the function on the differential vector.
    Differential<Float> dg=g(dx);
    std::cout << "g=" << g << "\n";
    std::cout << "dg=" << dg << "\n";
    std::cout << "dg[0][0,2,0]=" << dg[MultiIndex(3,0,2,0)] << "\n";

}