/***************************************************************************
 *            function/declarations.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
/*! \file function/declarations.h
 *  \brief Forward declarations of classes in the function module.
 */

#ifndef ARIADNE_FUNCTION_DECLARATIONS_H
#define ARIADNE_FUNCTION_DECLARATIONS_H

namespace Ariadne { 
  namespace Function {
    class MultiIndex;
    class SortedIndex;
    class PositionIndex;

    template<class R> class FunctionInterface;
    template<class R> class InterpretedFunction;
     
    template<class R> class ConstantFunction;
    template<class R> class AffineFunction;
    template<class R> class PolynomialFunction;

    template<class R> class AffineVariable;
    template<class R> class AffineDerivative;
    template<class R> class AffineModel;

    template<class R> class TaylorSeries;
    template<class R> class TaylorVariable;
    template<class R> class TaylorDerivative;
    template<class R> class TaylorModel;

  }
}

#endif /* ARIADNE_FUNCTION_DECLARATIONS_H */
