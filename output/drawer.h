/***************************************************************************
 *            drawer.h
 *
 *  Copyright  2011-12  Pieter Collins
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

/*! \file drawer.h
 *  \brief Class for drawing affine and nonlinear sets.
 */

#ifndef ARIADNE_DRAWER_H
#define ARIADNE_DRAWER_H

#include <iosfwd>
#include "utility/declarations.h"
#include "output/drawer_interface.h"

namespace Ariadne {

typedef Int Int;

struct Depth { Int _d; explicit Depth(Int d) : _d(d) { } operator Int() const { return _d; } };

//! \brief Draw a box over-approximation to the set.
class BoxDrawer : public DrawerInterface
{
    Void draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set);
};

//! \brief Draw an affine over-approximation to the set.
class AffineDrawer : public DrawerInterface
{
    Int _depth;
  public:
    AffineDrawer(Int depth) : _depth(depth) { }
    Void draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set);
};

//! \brief Subdivide the set and draw affine approximations to small pieces.
class SubdivisionDrawer : public DrawerInterface
{
    Void draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set);
};

//! \brief Pave the set and draw the computed cells.
class GridDrawer : public DrawerInterface
{
    Void draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set);
};


} //namespace Ariadne

#endif /* ARIADNE_DRAWER_H */
