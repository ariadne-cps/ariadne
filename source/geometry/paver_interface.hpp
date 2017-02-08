/***************************************************************************
 *            paver_interface.hpp
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

/*! \file paver_interface.hpp
 *  \brief Interface for class to compute approximations to function sets.
 */

#ifndef ARIADNE_PAVER_INTERFACE_HPP
#define ARIADNE_PAVER_INTERFACE_HPP

#include "utility/declarations.hpp"

namespace Ariadne {

class PavingInterface;
class ValidatedConstrainedImageSet;

//! \brief A class for computing outer approximations to sets defined by functions.
class PaverInterface
{
  public:
    typedef ValidatedConstrainedImageSet SetType;
  public:
    virtual Void  adjoin_outer_approximation(PavingInterface& paving, SetType const& set, Int depth) const = 0;
};

} //namespace Ariadne

#endif /* ARIADNE_PAVER_INTERFACE_HPP */
