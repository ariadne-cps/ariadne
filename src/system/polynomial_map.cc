/***************************************************************************
 *            polynomial_map.cc
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

#include "numeric/float.h"

#include "system/polynomial_map.h"
#include "system/polynomial_map.code.h"

namespace Ariadne {
  namespace System {
    using namespace Numeric;
    
#ifdef ENABLE_FLOAT64
    template class Monomial<Float64>;
    template class Polynomial<Float64>;
    template class PolynomialMap<Float64>;
    template class PolynomialMatrix<Float64>;

    template bool operator<(const Monomial<Float64>&, const Monomial<Float64>&);
    
    template std::ostream& operator<<(std::ostream&, const Monomial<Float64>&);
    template std::ostream& operator<<(std::ostream&, const Polynomial<Float64>&);
    template std::ostream& operator<<(std::ostream&, const PolynomialMap<Float64>&);
    template std::ostream& operator<<(std::ostream&, const PolynomialMatrix<Float64>&);

    template std::istream& operator>>(std::istream&, Monomial<Float64>&);
    template std::istream& operator>>(std::istream&, Polynomial<Float64>&);
    template std::istream& operator>>(std::istream&, PolynomialMap<Float64>&);
#endif
    
#ifdef ENABLE_FLOATMP
    template class Monomial<FloatMP>;
    template class Polynomial<FloatMP>;
    template class PolynomialMap<FloatMP>;
    template class PolynomialMatrix<FloatMP>;

    template bool operator<(const Monomial<FloatMP>&, const Monomial<FloatMP>&);
    
    template std::ostream& operator<<(std::ostream&, const Monomial<FloatMP>&);
    template std::ostream& operator<<(std::ostream&, const Polynomial<FloatMP>&);
    template std::ostream& operator<<(std::ostream&, const PolynomialMap<FloatMP>&);
    template std::ostream& operator<<(std::ostream&, const PolynomialMatrix<FloatMP>&);

    template std::istream& operator>>(std::istream&, Monomial<FloatMP>&);
    template std::istream& operator>>(std::istream&, Polynomial<FloatMP>&);
    template std::istream& operator>>(std::istream&, PolynomialMap<FloatMP>&);
#endif

  }
}
