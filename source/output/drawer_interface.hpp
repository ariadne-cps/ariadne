/***************************************************************************
 *            drawer_interface.hpp
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

/*! \file drawer_interface.hpp
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
