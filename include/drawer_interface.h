/***************************************************************************
 *            drawer_interface.h
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

/*! \file drawer_interface.h
 *  \brief Interface for drawing affine and nonlinear sets.
 */

#ifndef ARIADNE_DRAWER_INTERFACE_H
#define ARIADNE_DRAWER_INTERFACE_H

#include <iosfwd>

namespace Ariadne {

typedef void Void;
class CanvasInterface;
class PlanarProjectionMap;
typedef PlanarProjectionMap Projection2d;
class IntervalConstrainedImageSet;

//! \brief A class for computing outer approximations to sets defined by functions.
class DrawerInterface
{
    virtual Void draw(CanvasInterface& cnvs, const Projection2d& proj, const IntervalConstrainedImageSet& set) = 0;
};


} //namespace Ariadne

#endif /* ARIADNE_DRAWER_INTERFACE_H */
