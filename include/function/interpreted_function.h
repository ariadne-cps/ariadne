/***************************************************************************
 *            interpreted_function.h
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
 
/*! \file system/interpreted_function.h
 *  \brief Functions running on a virtual machine.
 */
 
#ifndef ARIADNE_SYSTEM_INTERPRETED_FUNCTION_H
#define ARIADNE_SYSTEM_INTERPRETED_FUNCTION_H

#include <vector>
#include <string>
#include <exception>
#include <stdexcept>

#include "base/types.h"
#include "base/array.h"
#include "numeric/traits.h"
#include "linear_algebra/declarations.h"
#include "function/function_interface.h"
#include "function/virtual_machine.h"
#include "function/multi_index.h"
#include "function/variable.h"

namespace Ariadne {
  namespace Function {
      
    /*!\ingroup Function
     * \brief A function running on a virtual machine.
     */
    template<class R>
    class InterpretedFunction
      : public FunctionInterface<R>
    {
      typedef typename Numeric::traits<R>::arithmetic_type A; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      
      /*! \brief  destructor. */
      virtual ~InterpretedFunction() { }
     
      /*! \brief Default constructor. */
      InterpretedFunction();
     
      /*! \brief Construct from a string literal. */
      InterpretedFunction(const std::string& filename);
     
      /*! \brief Construct by reading from an input stream. */
      InterpretedFunction( std::istream& is);
     
      /*! \brief Construct from a list of variables and an algorithm. */
      InterpretedFunction(const std::string& name,
                          const std::vector<FunctionVariable>& variables,
                          const std::vector<Numeric::Rational>& constants,
                          const std::vector<VirtualMachine::ByteCode>& algorithm);
     
      /*! \brief Make a copy (clone) of the vector field. */
      virtual InterpretedFunction<R>* clone() const;
     
      /*! \brief A bound for the vector field over a set of vectors. */
      virtual LinearAlgebra::Vector<A> evaluate(const LinearAlgebra::Vector<A>& x) const;
 
      /*! \brief A bound for the vector field over a set of vectors. */
      virtual LinearAlgebra::Matrix<A> jacobian(const LinearAlgebra::Vector<A>& x) const;

      /*! \brief A bound for the vector field over a set of vectors. */
      virtual Function::TaylorDerivative<A> derivative(const LinearAlgebra::Vector<A>& x, const smoothness_type& s) const;

      /*! \brief The degree of differentiability of the function. */
      virtual smoothness_type smoothness() const;
      /*! \brief The size of the function argument. */
      virtual size_type argument_size() const;
      /*! \brief The size of the function result. */
      virtual size_type result_size() const;

      /*! \brief The name of the function. */
      std::string name() const;
      
      /*! \brief Read from the file \a filename. */
      void read(const std::string& filename);
      
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
      void _initialise();
     private:
      std::string _name;
      array<FunctionVariable> _variables;
      array<VirtualMachine::ByteCode> _operations;
      array<A> _constants;
      size_type _argument_size;
      size_type _result_size;
      size_type _intermediates_size;
      size_type _constants_size;
    };
   

    template<class R> inline
    std::istream& operator>>(std::istream& is, InterpretedFunction<R>& f) {
      return f.read(is);
    }

    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const  InterpretedFunction<R>& f) {
      return f.write(os);
    };



  }
}

#endif /* ARIADNE_SYSTEM_INTERPRETED_FUNCTION_H */
