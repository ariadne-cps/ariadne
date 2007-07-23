/***************************************************************************
 *            ellipsoid.inline.h
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

namespace Ariadne {
  namespace Geometry {
     
    template<class R> inline
    Ellipsoid<R>::Ellipsoid(const Ellipsoid<R>& original)
      : _centre(original._centre), _bilinear_form(original._bilinear_form)
    {
    }
      
    template<class R> inline
    Ellipsoid<R>& 
    Ellipsoid<R>::operator=(const Ellipsoid<R>& original)
    {
      if(this != &original) {
        this->_centre = original._centre;
        this->_bilinear_form = original._bilinear_form;
      }
      return *this;
    }

    template<class R> inline
    bool 
    Ellipsoid<R>::operator==(const Ellipsoid<R>& other) const
    {
      return this->_centre==other._centre && this->_bilinear_form==other._bilinear_form;
    }
      
    template<class R> inline
    bool 
    Ellipsoid<R>::operator!=(const Ellipsoid<R>& other) const
    {
      return !(*this == other);
    }

    template<class R> inline
    const Point<R>& 
    Ellipsoid<R>::centre() const
    {
      return this->_centre;
    }
      
    template<class R> inline
    const LinearAlgebra::Matrix<R>&
    Ellipsoid<R>::bilinear_form() const
    {
      return this->_bilinear_form;
    }



    template<class R> inline
    size_type 
    Ellipsoid<R>::dimension() const
    {
      return this->_centre.dimension();
    }
      
    template<class R> inline
    bool 
    Ellipsoid<R>::empty() const
    {
      return false;
    }
      
    template<class R> inline
    bool 
    Ellipsoid<R>::empty_interior() const
    {
      return LinearAlgebra::singular(this->_bilinear_form);
    }
      



    template<class R> inline 
    tribool 
    disjoint(const Ellipsoid<R>& A, const Ellipsoid<R>& B) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    template<class R> inline 
    tribool 
    disjoint(const Ellipsoid<R>& A, const Rectangle<R>& B) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    template<class R> inline 
    tribool 
    disjoint(const Rectangle<R>& A, const Ellipsoid<R>& B) 
    {
      return disjoint(B,A);
    }
    
    
    
    template<class R> inline
    tribool 
    subset(const Ellipsoid<R>& A, 
           const Ellipsoid<R>& B) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    template<class R> inline
    tribool 
    subset(const Ellipsoid<R>& A, 
           const Rectangle<R>& B) 
    {
      return subset(A.bounding_box(),B);
    }
    
    template<class R> inline
    tribool 
    subset(const Rectangle<R>& A, 
           const Ellipsoid<R>& B) 
    {
      array< Point<R> > vertices=A.vertices();
      for(typename Rectangle<R>::vertex_iterator vertex_iter=vertices.begin(); vertex_iter!=vertices.end(); ++vertex_iter) {
        if(! B.contains(*vertex_iter) ) {
          return false;
        }
      }
      return true;
    }

    
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Ellipsoid<R>& e) {
      return e.write(os);
    }

    
    template<class R> inline
    std::istream& operator>>(std::istream& is, Ellipsoid<R>& e) {
      return e.read(is);
    }

  }
}
