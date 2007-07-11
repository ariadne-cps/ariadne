/***************************************************************************
 *            system/declarations.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file system/declarations.h
 *  \brief Forward declarations of classes.
 */

#ifndef ARIADNE_SYSTEM_DECLARATIONS_H
#define ARIADNE_SYSTEM_DECLARATIONS_H

namespace Ariadne { 
  namespace System {
    template<class R> class Function;
      
    template<class R> class MapInterface;
    template<class R> class AffineMap;
    template<class R, template<class> class BS > class AffineMultiMap;
    template<class R> class GridMultiMap;
    template<class R> class Monomial;
    template<class R> class Polynomial;
    template<class R> class PolynomialMap;
    template<class R> class PolynomialMatrix;
    template<class R> class DiscreteTimeSystem;
    template<class R> class VectorFieldInterface;
    template<class R> class AffineVectorField;

    template<class R> class DiscreteMode;
    template<class R> class DiscreteTransition;
    template<class R> class HybridAutomaton;

    template<class R> class ConstraintDiscreteMode;
    template<class R> class ConstraintDiscreteTransition;
    template<class R> class ConstraintHybridAutomaton;
  }
}

#endif /* ARIADNE_SYSTEM_DECLARATIONS_H */
