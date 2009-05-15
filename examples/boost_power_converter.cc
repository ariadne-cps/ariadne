/***************************************************************************
 *            boost_power_converter.cc
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

/*! \file boost_power_converter.cc */

#include "ariadne.h"

using namespace Ariadne;

HybridAutomaton create_boost_power_converter()
{
    // Variables V, I
    // Parameters L,R,C,E
    double E,L,R,C;
    double RC=R*C;

    // Open, I>0, dot(V)=(I-V/R)/C; dot(I)=(E-V)/L
    AffineFunction open_dyamic(2, 0.0,E/L, 2, -1/RC,1/C, -1/L,0.0);     
    // Open, V<E  dot(V)=-V/RC; I=0;
    AffineFunction saturated_dyamic(1, 0.0, 1, -1/RC);
    // Closed,   dot(V)=-V/RC; dot(I)=E/L
    AffineFunction closed_dyamic(2, 0.0,E/L, 2, -1/RC,0.0, 0.0,0.0);     

    

}

int main() {
}

