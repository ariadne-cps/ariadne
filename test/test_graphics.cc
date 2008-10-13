/***************************************************************************
 *            test_graphics.cc
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
 
#include "graphics.h"
#include "box.h"

using namespace Ariadne;

int main(int argc, char **argv) 
{

    Box bx1(2); bx1[0]=Interval(0.1,0.3); bx1[1]=Interval(0.05,0.15);
    Box bx2(2); bx2[0]=Interval(0.2,0.4); bx2[1]=Interval(0.10,0.25);
    Box bx3(2); bx3[0]=Interval(0.25,0.5); bx3[1]=Interval(0.20,0.50);
    Box bx4(2); bx4[0]=Interval(0.4,0.8); bx4[1]=Interval(0.40,1.1);
    Graphic g;
    g.plot(bx1);
    g.plot(bx2);
    g.plot(bx3);
    g.plot(bx4);
    g.write("test_graphics");
    g.display();
    g.clear();
    g.plot(bx1);
    g.plot(bx4);
    g.display();

}

