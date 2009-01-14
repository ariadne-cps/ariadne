/***************************************************************************
 *            sucessor.h
 *
 *  Copyright  2006-8  Pieter Collins
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
 
/*! \file sucessor.h
 *  \brief Methods for computing the evolution of a system on a grid, suitable for model checking.
 */

#ifndef ARIADNE_SUCESSOR_H
#define ARIADNE_SUCESSOR_H

#include <bitset>
#include <vector>

namespace Ariadne {

typedef double Float;
class Grid;
class GridCell;
class VectorField;
template<class Sys, class BS> class DiscretiserInterface;

struct Cell {
    unsigned char height;
    unsigned char depth;
    std::bitset<112> word;
};

enum EvolutionKind { REACH, EVOLVE };

GridCell make_grid_cell(const Cell&, const Grid&);
Cell make_cell(const GridCell&);
std::ostream& operator<<(std::ostream& os, const Cell& c);
std::vector<Cell> successor(const DiscretiserInterface<VectorField,GridCell>&, const Grid&, const VectorField&, const Cell&, Float, EvolutionKind);

} // namespace Ariadne



#endif /* ARIADNE_SUCESSOR_H */
