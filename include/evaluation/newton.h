/***************************************************************************
 *            newton.h
 *
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
 
/*! \file newton.h
 *  \brief Newton and interval Newton methods.
 */

#ifndef _ARIADNE_NEWTON_H
#define _ARIADNE_NEWTON_H

#include <exception>
#include <stdexcept>
#include <string>

#include "../declarations.h"

#include "linear_algebra/interval_vector.h"
#include "linear_algebra/interval_matrix.h"
#include "geometry/rectangle.h"
#include "system/map.h"
#include "system/vector_field.h"

namespace Ariadne {
  namespace Evaluation {
   
    class EvaluationException
      : public std::exception 
    {
     public:
       EvaluationException(const std::string& s) : _what(s) { }
      ~EvaluationException() throw () { }
      const char* what() const throw () { return this->_what.c_str(); }
     private:
      std::string _what;
    };
    
    template<typename Real>
    Geometry::Rectangle<Real>
    interval_newton(const System::VectorField<Real>& f, 
                    const Geometry::Rectangle<Real>& r, 
                    const Real& e,
                    uint max_steps=64);

  }
}


namespace Ariadne {
  namespace System {

    template<typename R> class DifferenceMap
      : public VectorField<R>
    {
     public:
      DifferenceMap(const Map<R>& f) : _base(f) { assert(f.argument_dimension()==f.result_dimension()); }
      DifferenceMap(const HenonMap<R>& f) : _base(f) { }
      virtual dimension_type dimension() const { return _base.argument_dimension(); }
      virtual LinearAlgebra::IntervalVector<R> operator() (const Geometry::Rectangle<R>& r) const {
        return _base(r)-r; }
      virtual LinearAlgebra::IntervalMatrix<R> derivative(const Geometry::Rectangle<R>& r) const {
        LinearAlgebra::IntervalMatrix<R> d=_base.derivative(r);
        LinearAlgebra::IntervalMatrix<R> i=LinearAlgebra::IntervalMatrix<R>::identity(this->dimension());
        return d-i; }
      virtual std::string name() const { return "DifferenceMap"; }
     private:
      const Map<R>& _base;
    };
  }
}


#endif /* _ARIADNE_NEWTON_H */
