/***************************************************************************
 *            applicator.h
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
 
/*! \file applicator.h
 *  \brief Methods for computing the images of sets under maps.
 */

#ifndef _ARIADNE_APPLICATOR_H
#define _ARIADNE_APPLICATOR_H

#include "../declarations.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief A class for computing the image of a set under a map. 
     *  \ingroup Applicators
     */
    template<class R>
    class Applicator {
     public:
      /*! \brief Compute the image of a basic set under a continuous function. */
      virtual ~Applicator();
      
      /*! \brief Compute the image of a rectangle under a continuous function. */
      virtual Geometry::Rectangle<R> image(const System::Map<R>& f, const Geometry::Rectangle<R>& s) const;

      /*! \brief Compute the image of a zonotope under a continuous function. */
      virtual Geometry::Parallelotope<R> image(const System::Map<R>& f, const Geometry::Parallelotope<R>& s) const;

      /*! \brief Compute the image of a parallelotope under a continuous function. */
      virtual Geometry::Zonotope<R> image(const System::Map<R>& f, const Geometry::Zonotope<R>& s) const;

      /*! \brief Compute the image of an interval parallelotope under a continuous function. */
      virtual Geometry::Parallelotope< Interval<R> > image(const System::Map<R>& f, const Geometry::Parallelotope< Interval<R> >& s) const;

      /*! \brief Compute the image of an interval zonotope under a continuous function. */
      virtual Geometry::Zonotope< Interval<R> > image(const System::Map<R>& f, const Geometry::Zonotope< Interval<R> >& s) const;
     
     protected:
      /*! \brief Template for integrating a list set. */
      template<class Rl,template<class> class BS>
      Geometry::ListSet<Rl,BS> 
      image_list_set(const System::Map<R>& f, 
                     const Geometry::ListSet<Rl,BS>& initial_set) const;

      
      /*! \brief Template for integrating a basic set. */
      template<template<class> class BS>
      BS<R>
      image_basic_set(const System::Map<R>& f, 
                      const BS<R>& initial_set) const;
     public:
      /*! \brief Compute the image of a list set under a map. */
      virtual 
      Geometry::ListSet<R,Geometry::Rectangle> 
      image(const System::Map<R>& f, const Geometry::ListSet<R,Geometry::Rectangle>& ds) const;
       
      virtual 
      Geometry::ListSet<R,Geometry::Parallelotope> 
      image(const System::Map<R>& f, const Geometry::ListSet<R,Geometry::Parallelotope>& ds) const;
       
      virtual 
      Geometry::ListSet<R,Geometry::Zonotope> 
      image(const System::Map<R>& f, const Geometry::ListSet<R,Geometry::Zonotope>& ds) const;
      
      
      /*! \brief Compute the image of a list set under a map. */
      virtual 
      Geometry::ListSet<Interval<R>,Geometry::Parallelotope> 
      image(const System::Map<R>& f, const Geometry::ListSet<Interval<R>,Geometry::Parallelotope>& ds) const;
      
      virtual 
      Geometry::ListSet<Interval<R>,Geometry::Zonotope> 
      image(const System::Map<R>& f, const Geometry::ListSet<Interval<R>,Geometry::Zonotope>& ds) const;
      
      
      /*! \brief Compute the image of \a map starting in \a initial_set computing the result on \a grid. */
      virtual
      Geometry::GridCellListSet<R> 
      image(const System::Map<R>& map, 
            const Geometry::GridCellListSet<R>& initial_set,
            const Geometry::Grid<R>& grid) const;

      /*! \brief Compute the image of \a map starting in \a initial_set computing the result on \a grid. */
      virtual
      Geometry::GridMaskSet<R> 
      image(const System::Map<R>& map, 
            const Geometry::GridMaskSet<R>& initial_set,
            const Geometry::FiniteGrid<R>& grid) const;

      /*! \brief Compute the image of \a map starting in \a initial_set while remaining in \a bounding_set. */
      virtual
      Geometry::GridMaskSet<R> 
      image(const System::Map<R>& map, 
            const Geometry::GridMaskSet<R>& initial_set,
            const Geometry::GridMaskSet<R>& bounding_set) const;

      /*! \brief Compute the reachable set of \a map starting in \a initial_set. */
      virtual
      Geometry::GridMaskSet<R> 
      reach(const System::Map<R>& map, 
            const Geometry::GridMaskSet<R>& initial_set) const;

      
      /*! \brief Compute the chain-reachable set of \a map starting in \a initial_set while staying within \a bounding_set. */
      virtual
      Geometry::GridMaskSet<R> 
      chainreach(const System::Map<R>& map, 
                 const Geometry::GridMaskSet<R>& initial_set, 
                 const Geometry::GridMaskSet<R>& bounding_set) const;
    
      /*! \brief Attempt to verify that the reachable set of \a map starting in \a initial_set remains in \a safe_set. */
      virtual
      bool
      verify(const System::Map<R>& map, 
             const Geometry::GridMaskSet<R>& initial_set, 
             const Geometry::GridMaskSet<R>& safe_set) const;
    };
    
    
    
    
    /*! \brief Compute the image of a rectangle under a continuous function. 
     *  \ingroup Apply
     */
    template<class R>
    Geometry::Rectangle<R> 
    image(const System::Map<R>& f, const Geometry::Rectangle<R>& s) {
      return Applicator<R>().image(f,s);
    }
    
    /*! \brief Compute the image of a parallelotope under a differentiable function. 
     *  \ingroup Apply
     */
    template<class R>
    inline
    Geometry::Parallelotope<R> 
    image(const System::Map<R>& f, const Geometry::Parallelotope<R>& s) {
      return Applicator<R>().image(f,s);
    }
    
    /*! \brief Compute the image of a parallelotope under a differentiable function.  
     *  \ingroup Apply
     */
    template<class R>
    inline
    Geometry::ListSet<R,Geometry::Parallelotope>
    image(const System::Map<R>& f, const Geometry::ListSet<R,Geometry::Parallelotope>& s) {
      return Applicator<R>().image(f,s);
    }
    
    /*! \brief Compute the chain-reachable set of \a map starting in \a initial_set on the grid \a grid while staying within \a bounds.  
     *  \ingroup Apply
     */
    template<class R>
    inline
    Geometry::GridMaskSet<R> 
    image(const System::Map<R>& map, 
          const Geometry::GridMaskSet<R>& initial_set, 
          const Geometry::GridMaskSet<R>& bounding_set) 
    {
      return Applicator<R>().image(map,initial_set,bounding_set);
    }

    /*! \brief Compute the chain-reachable set of \a map starting in \a initial_set on the grid \a grid while staying within \a bounds.  
     *  \ingroup Apply
     */
    template<class R>
    inline
    Geometry::GridMaskSet<R> 
    chainreach(const System::Map<R>& map, 
               const Geometry::GridMaskSet<R>& initial_set, 
               const Geometry::GridMaskSet<R>& bounding_set) 
    {
      return Applicator<R>().chainreach(map,initial_set,bounding_set);
    }

    /*! \brief Attempt to verify that the reachable set of \a map starting in \a initial_set remains in \a safe_set. 
     *  \ingroup Apply
     */
    template<class R>
    inline
    bool
    verify(const System::Map<R>& map, 
               const Geometry::GridMaskSet<R>& initial_set, 
               const Geometry::GridMaskSet<R>& safe_set) 
    {
      return Applicator<R>().verify(map,initial_set,safe_set);
    }

  }
}

#endif /* _ARIADNE_APPLY_H */
