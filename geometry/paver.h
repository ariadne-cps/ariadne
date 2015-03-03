/***************************************************************************
 *            paver.h
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

/*! \file paver.h
 *  \brief Class for computing outer approximations to nonlinear sets.
 */

#ifndef ARIADNE_PAVER_H
#define ARIADNE_PAVER_H

#include "geometry/paver_interface.h"

namespace Ariadne {

//! \brief A class for computing outer approximations to sets defined by functions.
class AffinePaver : public PaverInterface
{
  public:
    virtual Void
    adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunction& space_function,
                               const ValidatedVectorFunction& constraint_function, const ExactBoxType& constraint_bounds, Int depth) const;
};

//! \brief A class for computing outer approximations to sets defined by functions.
class SubdivisionPaver : public PaverInterface
{
  public:
    virtual Void
    adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunction& space_function,
                               const ValidatedVectorFunction& constraint_function, const ExactBoxType& constraint_bounds, Int depth) const;
};

//! \brief A class for computing outer approximations to sets defined by functions.
class ConstraintPaver : public PaverInterface
{
  public:
    virtual Void
    adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunction& space_function,
                               const ValidatedVectorFunction& constraint_function, const ExactBoxType& constraint_bounds, Int depth) const;
};

//! \brief A class for computing outer approximations to sets defined by functions.
class OptimalConstraintPaver : public PaverInterface
{
  public:
    virtual Void
    adjoin_outer_approximation(PavingInterface& paving, const DomainType& domain, const ValidatedVectorFunction& space_function,
                               const ValidatedVectorFunction& constraint_function, const ExactBoxType& constraint_bounds, Int depth) const;
};

} //namespace Ariadne

#endif /* ARIADNE_PAVER_H */
