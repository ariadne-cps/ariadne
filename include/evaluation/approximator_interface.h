/***************************************************************************
 *            approximator_interface.h
 *
 *  Copyright  2007  Pieter Collins
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
 
/*! \file approximator_interface.h
 *  \brief Interface for approximating basic sets
 */

#ifndef ARIADNE_APPROXIMATOR_INTERFACE_H
#define ARIADNE_APPROXIMATOR_INTERFACE_H

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"

#include "geometry/list_set.h"

namespace Ariadne {
  
  
    /*! \brief Interface for approximating enclosure sets on a paving of space.
     *  \ingroup EvaluatorInterfaces \ingroup Approximators
     */
    template<class Aprx, class ES> 
    class ApproximatorInterface 
    { 
      typedef typename Aprx::BasicSet BasicSet;
      typedef typename Aprx::Paving Paving;
      typedef typename Aprx::CoverListSet CoverListSet;
      typedef typename Aprx::PartitionListSet PartitionListSet;
      typedef typename Aprx::PartitionTreeSet PartitionTreeSet;
      typedef ES EnclosureSet;
      typedef ListSet<ES> EnclosureSetList;
     public:
      /*! \brief Virtual destructor. */
      virtual ~ApproximatorInterface() { }

      /*! \brief Create a dynamically-allocated copy. */
      virtual ApproximatorInterface<Aprx,ES>* clone() const = 0;

      /*! \brief Computes an over-approximation of a set from a rectangle. */
      virtual EnclosureSet enclosure_set(const BasicSet& bs) const = 0;

      /*! \brief Computes a bounding box for a set. */
      virtual BasicSet bounding_box(const EnclosureSet& es) const = 0;

      /*! \brief Computes an outer-approximation of a set on a grid. */
      virtual CoverListSet lower_approximation(const EnclosureSet& es) const = 0;

      /*! \brief Computes an outer-approximation of a set. */
      virtual PartitionListSet inner_approximation(const EnclosureSet& es) const = 0;

      /*! \brief Computes an outer-approximation of a set. */
      virtual PartitionListSet outer_approximation(const EnclosureSet& es) const = 0;

      /*! \brief Computes an outer-approximation of a set on a grid. */
      virtual PartitionListSet inner_approximation(const EnclosureSet& es, const Paving& pv) const = 0;

      /*! \brief Computes an outer-approximation of a set on a grid. */
      virtual PartitionListSet outer_approximation(const EnclosureSet& es, const Paving& pv) const = 0;


      /*! \brief Computes a bounding box for a set. */
      virtual BasicSet bounding_box(const EnclosureSetList& esl) const = 0;

      /*! \brief Computes an lower-approximation of a set. */
      virtual CoverListSet lower_approximation(const EnclosureSetList& esl) const = 0;

      /*! \brief Computes an outer-approximation of a set. */
      virtual PartitionListSet outer_approximation(const EnclosureSetList& esl) const = 0;

      /*! \brief Computes an outer-approximation of a set. */
      virtual PartitionListSet outer_approximation(const EnclosureSetList& esl, const Paving& pv) const = 0;


      /*! \brief Adjoins an outer approximation to a basic set to a grid mask set. */
      virtual void adjoin_outer_approximation(PartitionTreeSet& pts, const EnclosureSet& es) const = 0;

      /*! \brief Computes and over-approximation of a set from a rectangle. */
      virtual void adjoin_outer_approximation(PartitionTreeSet& pts, const EnclosureSetList& esl) const = 0;

      /*! \brief Write to an output stream. */
      virtual std::ostream& write(std::ostream& os) const = 0;
    };

    template<class Aprx, class ES>
    std::ostream& operator<<(std::ostream& os, const ApproximatorInterface<Aprx,ES>& ai) {
      return ai.write(os); 
    }




  
} // namespace Ariadne


#endif /* ARIADNE_APPROXIMATOR_INTERFACE_H */
