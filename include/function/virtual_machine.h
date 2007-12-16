/***************************************************************************
 *            virtual_machine.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
/*! \file function/virtual_machine.h
 *  \brief Virtual machine for function evaluation.
 */
 
#ifndef ARIADNE_FUNCTION_VIRTUAL_MACHINE_H
#define ARIADNE_FUNCTION_VIRTUAL_MACHINE_H

#include <vector>
#include <string>
#include <exception>
#include <stdexcept>

#include "base/types.h"
#include "base/array.h"
#include "numeric/traits.h"
#include "linear_algebra/declarations.h"
#include "function/interpreted_function.h"

namespace Ariadne {
  namespace Function {

    /*!\ingroup Function
     * \brief A virtual machine for evaluating functions.
     */
    class VirtualMachine
    {
     public:
      /*! \brief The operation codes of the virtual machine. */
      enum Operation { PUSH, PULL, CONST, MAX, MIN, ABS, POS, NEG, ADD, SUB, MUL, DIV, POW, EXP, LOG, SIN, COS, TAN, ASIN, ACOS, ATAN };
    
      /*! \brief An index of a value. */
      struct Index { short i,j; short& operator[](int n) { return (&i)[n]; } };

      /*! \brief A byte code of a virtual machine; either an operation, an index/pointer, or a value.. */
      union ByteCode { Operation op; Index ind; int val; };

      /*! \brief Evaluate the virtual machine on a program. */
      template<class X> void evaluate(const array<ByteCode>& program, X** arguments) const;
     private:
      array<ByteCode> _operations;
    };

   
  }
}  



#endif /* ARIADNE_FUNCTION_VIRTUAL_MACHINE_H */
