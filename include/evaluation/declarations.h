/***************************************************************************
 *            evaluation/declarations.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file evaluation/declarations.h
 *  \brief Forward declarations of classes in the Evaluation module.
 */

#ifndef ARIADNE_EVALUATION_DECLARATIONS_H
#define ARIADNE_EVALUATION_DECLARATIONS_H

namespace Ariadne { 
  namespace Evaluation {

    template<class R> class SolverInterface;
  
    template<class R> class BounderInterface;
    template<class R> class DetectorInterface;

    template<class BS> class ApproximatorInterface;
    template<class BS> class ApplicatorInterface;
    template<class BS> class IntegratorInterface;
    template<class BS> class DifferentiableIntegratorInterface;
    template<class BS> class MapOrbiterInterface;
    template<class BS> class VectorFieldOrbiterInterface;

    template<class R> class EvolutionParameters;

    template<class R> class ModelChecker;
    template<class R> class MapEvolver;
    template<class R> class VectorFieldEvolver;
    template<class R> class SetBasedHybridEvolver;
    template<class R> class ConstraintBasedHybridEvolver;

  }
}

#endif /* ARIADNE_EVALUATION_DECLARATIONS_H */
