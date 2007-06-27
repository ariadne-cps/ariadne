/***************************************************************************
 *      lptst.template.h
 *
 * Copyright 2006 Alberto Casagrande, Pieter Collins
 * casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

namespace Ariadne {
  namespace LinearAlgebra {
    
    
    template<class R, class AP>
    tribool lptst(const Matrix<R>& A, const Vector<R>& b, Vector<AP>& x) {
      tribool ans;
      bool indeterminate_state = false;
      
      if (verbosity > 2)
        std::clog << "lptst() \n A:" << A << ", b:" << b << ", x:" << x << std::endl;
      
      ARIADNE_CHECK_SIZE(b, A.number_of_rows(), __PRETTY_FUNCTION__)
      ARIADNE_CHECK_SIZE(x, A.number_of_columns(), __PRETTY_FUNCTION__)
      
      for (size_type i = 0; i != b.size(); i++) {
        AP sum = 0;
        for (size_type j = 0; j != x.size(); j++)
          sum += A(i, j) * x(j);
        
        if (verbosity > 3) std::clog << "sum: " << sum << ", b(i): " << b(i) << std::endl;
        
        ans = sum <= b(i);
        if (!ans) return ans;
        if (indeterminate(ans))  indeterminate_state = true;
      }
      if (indeterminate_state)
        return tribool(indeterminate);
      else
        return true;
    }
    
    
    template<>
    inline tribool lptst<Numeric::Float64, Numeric::Float64>(const Matrix<Numeric::Float64>& A, const Vector<Numeric::Float64>& b, const Vector<Numeric::Float64>& x) {
      Vector < Numeric::Interval < Numeric::Float64 > > _x(x);
      return lptst(A, b, _x);
    }
    
    
    template<class R, class AP>
    tribool lptstopt(const Matrix<R>& A, const Vector<R>& b, const Vector<R>& c, const Vector<AP>& x, const Vector<AP>& y) {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    
  } // namespace LinearAlgebra
} // namespace Ariadne
