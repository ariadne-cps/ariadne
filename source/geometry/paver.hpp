/***************************************************************************
 *            paver.hpp
 *
 *  Copyright  2011-12  Pieter Collins
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

/*! \file paver.hpp
 *  \brief Class for computing outer approximations to nonlinear sets.
 */

#ifndef ARIADNE_PAVER_HPP
#define ARIADNE_PAVER_HPP

#include "../geometry/paver_interface.hpp"

namespace Ariadne {

class Paver {
    SharedPointer<const PaverInterface> _ptr;
  public:
    using SetType = PaverInterface::SetType;
    Paver(SharedPointer<const PaverInterface> ptr) : _ptr(ptr) { }
    Paver(const PaverInterface* ptr) : _ptr(ptr) { }
    Void adjoin_outer_approximation(PavingInterface& paving, const SetType& set, Nat depth) const {
        this->_ptr->adjoin_outer_approximation(paving,set,depth); }
    friend OutputStream& operator<<(OutputStream& os, Paver const& pv) { return pv._ptr->_write(os); }
};

//! \brief A class for computing outer approximations to sets defined by functions.
class AffinePaver : public PaverInterface
{
  public:
    Void adjoin_outer_approximation(PavingInterface& paving, const SetType& set, Nat depth) const;
    OutputStream& _write(OutputStream& os) const;
};

//! \brief A class for computing outer approximations to sets defined by functions.
class SubdivisionPaver : public PaverInterface
{
  public:
    Void adjoin_outer_approximation(PavingInterface& paving, const SetType& set, Nat depth) const;
    Void adjoin_outer_approximation_recursion(PavingInterface& paving, ValidatedConstrainedImageSet const& set, Nat depth, const Vector<FloatDPValue>& max_errors) const;
    OutputStream& _write(OutputStream& os) const;
};

//! \brief A class for computing outer approximations to sets defined by functions.
class ReducePaver : public PaverInterface
{
  public:
    Void adjoin_outer_approximation(PavingInterface& paving, const SetType& set, Nat depth) const;
    OutputStream& _write(OutputStream& os) const;
};

//! \brief A class for computing outer approximations to sets defined by functions.
class ConstraintPaver : public PaverInterface
{
  public:
    Void adjoin_outer_approximation(PavingInterface& paving, const SetType& set, Nat depth) const;
    OutputStream& _write(OutputStream& os) const;
};

//! \brief A class for computing outer approximations to sets defined by functions.
class OptimalConstraintPaver : public PaverInterface
{
  public:
    Void adjoin_outer_approximation(PavingInterface& paving, const SetType& set, Nat depth) const;
    OutputStream& _write(OutputStream& os) const;
};

} //namespace Ariadne

#endif /* ARIADNE_PAVER_HPP */
