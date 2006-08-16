/***************************************************************************
 *            apply.h
 *
 *  17 January 2006
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file apply.h
 *  \brief Methods for computing the images of sets under maps.
 */

#ifndef _ARIADNE_APPLY_H
#define _ARIADNE_APPLY_H

#include "../declarations.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for computing the image of a set under a map. */
    template<typename R, template<typename> class BS>
    class Applicator {
     public:
      /*! \brief Compute the image of a basic set under a continuous function. */
      virtual ~Applicator();
      
      /*! \brief Compute the image of a basic set under a continuous function. */
      virtual BS<R> apply(const System::Map<R>& f, const BS<R>& s) const = 0;

      /*! \brief Compute the image of a list set under a map. */
      virtual 
      Geometry::ListSet<R,BS> 
      apply(const System::Map<R>& f, const Geometry::ListSet<R,BS>& ds) const;
       
      /*! \brief Compute the image of \a map starting in \a initial_set while remaining in \a bounding_set. */
      virtual
      Geometry::GridMaskSet<R> 
      apply(const System::Map<R>& map, 
            const Geometry::GridMaskSet<R>& initial_set,
            const Geometry::GridMaskSet<R>& bounding_set) const;

      /*! \brief Compute the chain-reachable set of \a map starting in \a initial_set while staying within \a bounding_set. */
      virtual
      Geometry::GridMaskSet<R> 
      chainreach(const System::Map<R>& map, 
                 const Geometry::GridMaskSet<R>& initial_set, 
                 const Geometry::GridMaskSet<R>& bounding_set) const;
    };
    
    /*! \brief A class for computing the image of a continuous map on a rectangle. */
    template<typename R>
    class C0Applicator
      : public Applicator<R,Geometry::Rectangle> 
    {
     public:
      /*! \brief Compute the image of a rectangle under a continuous function. */
      virtual Geometry::Rectangle<R> apply(const System::Map<R>& f, const Geometry::Rectangle<R>& r) const;
    };

    /*! \brief A class for computing the image of a differentiable map on a parallelotope. */
    template<typename R>
    class C1Applicator
      : public Applicator<R,Geometry::Parallelotope> 
    {
     public:
      /*! \brief Compute the image of a parallelotope under a continuous function. */
      virtual Geometry::Parallelotope<R> apply(const System::Map<R>& f, const Geometry::Parallelotope<R>& p) const;
    };

    
    
    
    /*! \brief Compute the image of a rectangle under a continuous function. */
    template<typename R>
    Geometry::Rectangle<R> 
    apply(const System::Map<R>& f, const Geometry::Rectangle<R>& s) {
      return C0Applicator<R>().apply(f,s);
    }
    
    /*! \brief Compute the image of a parallelotope under a differentiable function. */
    template<typename R>
    inline
    Geometry::Parallelotope<R> 
    apply(const System::Map<R>& f, const Geometry::Parallelotope<R>& s) {
      return C1Applicator<R>().apply(f,s);
    }
    
    /*! \brief Compute the image of a parallelotope under a differentiable function. */
    template<typename R>
    inline
    Geometry::ListSet<R,Geometry::Parallelotope>
    apply(const System::Map<R>& f, const Geometry::ListSet<R,Geometry::Parallelotope>& s) {
      return C1Applicator<R>().Applicator<R,Geometry::Parallelotope>::apply(f,s);
    }
    
    /*! \brief Compute the chain-reachable set of \a map starting in \a initial_set on the grid \a grid while staying within \a bounds. */
    template<typename R>
    inline
    Geometry::GridMaskSet<R> 
    apply(const System::Map<R>& map, 
          const Geometry::GridMaskSet<R>& initial_set, 
          const Geometry::GridMaskSet<R>& bounding_set) 
    {
      return C1Applicator<R>().Applicator<R,Geometry::Parallelotope>::apply(map,initial_set,bounding_set);
    }

    /*! \brief Compute the chain-reachable set of \a map starting in \a initial_set on the grid \a grid while staying within \a bounds. */
    template<typename R>
    inline
    Geometry::GridMaskSet<R> 
    chainreach(const System::Map<R>& map, 
               const Geometry::GridMaskSet<R>& initial_set, 
               const Geometry::GridMaskSet<R>& bounding_set) 
    {
      return C1Applicator<R>().chainreach(map,initial_set,bounding_set);
    }

  }
}

#endif /* _ARIADNE_APPLY_H */
