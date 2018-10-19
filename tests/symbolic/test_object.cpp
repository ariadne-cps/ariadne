/***************************************************************************
 *            symbolic/test_object.cpp
 *
 *  Copyright 2013-18  Pieter Collins
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

#warning "Test object functionality should move when complete"

#include "symbolic/object.hpp"


using namespace Ariadne;

void test() {
    Symbolic<Real> sr(3);
    sr=add(sr,sr);
    sr=sqrt(sr);
    RealInterface const& ri=sr;
    Real r=sr;

    auto sf = SymbolicFunction<Real(Real)>::identity();
    sf(r);
    std::cout << "sr=" << sr << "\n";
    std::cout << "r=" << r << "\n";
    std::cout << "sf(r)=" << sf(r) << "\n";

}


int main() {
    test();
    return 0;
}
