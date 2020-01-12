/***************************************************************************
 *            output/drawer_interface.hpp
 *
 *  Copyright  2011-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file output/drawer_interface.hpp
 *  \brief Interface for drawing affine and nonlinear sets.
 */

#ifndef ARIADNE_DRAWER_INTERFACE_HPP
#define ARIADNE_DRAWER_INTERFACE_HPP

#include <iosfwd>
#include "../utility/writable.hpp"

namespace Ariadne {

typedef Void Void;
class CanvasInterface;
struct PlanarProjectionMap;
typedef PlanarProjectionMap Projection2d;
class ValidatedConstrainedImageSet;

//! \related TaylorConstrainedImageSet \brief The possible types of method used to draw a nonlinear set.
enum class DrawingMethod : std::uint8_t { CURVE, BOX, AFFINE, GRID };
//! \related TaylorConstrainedImageSet \brief The type of method currently used to draw a set.
//! HACK: May be replaced by more advanced functionality in the future.
extern DrawingMethod DRAWING_METHOD;
//! \related TaylorConstrainedImageSet \brief The accuracy used to draw a set.
//! HACK: May be replaced by more advanced functionality in the future.
extern uint DRAWING_ACCURACY;

//! \brief A class for computing outer approximations to sets defined by functions.
class DrawerInterface : public WritableInterface
{
  public:
    virtual Void draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set) const = 0;
};

//! \brief A class for computing outer approximations to sets defined by functions.
class Drawer
{
    SharedPointer<const DrawerInterface> _ptr;
  public:
    Drawer(SharedPointer<const DrawerInterface> ptr) : _ptr(ptr) { }
    Drawer(const DrawerInterface* ptr) : _ptr(ptr) { }
    inline Void draw(CanvasInterface& cnvs, const Projection2d& proj, const ValidatedConstrainedImageSet& set) const {
        return this->_ptr->draw(cnvs,proj,set); }
    friend inline OutputStream& operator<<(OutputStream& os, Drawer const& drw) {
        return drw._ptr->_write(os); }
};


} //namespace Ariadne

#endif /* ARIADNE_DRAWER_INTERFACE_HPP */
