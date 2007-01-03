/***************************************************************************
 *            polyhedron.inline.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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

namespace Ariadne {  
  namespace Geometry {

    template<class R> template<class Rl> inline 
    Polyhedron<R>::Polyhedron(const Polyhedron<Rl>& original)
      : _A(original.A()), _b(original.b()) 
    { 
    }
        
    template<class R> inline 
    const LinearAlgebra::Matrix<R>& 
    Polyhedron<R>::A() const 
    { 
      return this->_A; 
    }
    
    template<class R> inline 
    const LinearAlgebra::Vector<R>& 
    Polyhedron<R>::b() const 
    { 
      return this->_b; 
    }
    
    template<class R> inline 
    size_type 
    Polyhedron<R>::number_of_constraints() const 
    { 
      return _A.number_of_rows(); 
    }
    
    
    
    
    
    template<class R> template<class Rl1, class Rl2> inline
    Polyhedron< Interval<R> >::Polyhedron(const LinearAlgebra::Matrix<Rl1> A, const LinearAlgebra::Vector<Rl2> b) 
      : _A(A), _b(b) 
    { 
    }
    
    
    template<class R> inline
    dimension_type 
    Polyhedron< Interval<R> >::dimension() const
    { 
      return this->_A.number_of_columns(); 
    }
    
    
    template<class R> inline
    size_type 
    Polyhedron< Interval<R> >::number_of_constraints() const
    { 
      return this->_A.number_of_rows(); 
    }
    
    
    template<class R> inline
    const LinearAlgebra::Matrix< Interval<R> >& 
    Polyhedron< Interval<R> >::A() const 
    { 
      return this->_A; 
    }
    
    
    template<class R> inline
    const LinearAlgebra::Vector< Interval<R> >& 
    Polyhedron< Interval<R> >::b() const 
    { 
      return this->_b; 
    }
    
    
    
    
    
    
    template<class R> inline
    Constraint<R>::Constraint(const dimension_type d, const R* a, const R& b)
      : _d(d), _a(a), _b(&b) 
    { 
    }
    
    
    template<class R1> template<class R2> inline 
    tribool Constraint<R1>::satisfied_by(const Point<R2>& pt) const
    {
      typedef typename Numeric::traits<R1,R2>::arithmetic_type F;
      F prod=0;
      for(dimension_type i=0; i!=pt.dimension(); ++i) {
        prod+=this->_a[i]*pt[i]; 
      }
      return prod<=*this->_b;
    }
    
    
    
    
    
    
    template<class R>
    class PolyhedronConstraintsIterator
      : public boost::iterator_facade<PolyhedronConstraintsIterator<R>,
                                      Constraint<R>,
                                      boost::forward_traversal_tag,
                                      const Constraint<R>&,
                                      const Constraint<R>*
                                     >
    {
     public:
      PolyhedronConstraintsIterator(const Polyhedron<R>& ply, const size_type& n)
        : _c(ply.dimension(),ply.A().begin()+n*ply.dimension(),ply.b().begin()[n]) { }
      bool equal(const PolyhedronConstraintsIterator<R>& other) const { 
        return this->_c._a==other._c._a; }
      const Constraint<R>& dereference() const { return _c; }
      void increment() { _c._a+=_c._d; _c._b+=1; }
     private:
      Constraint<R> _c;
    };
    
    
    
    
    
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Constraint<R>& c) 
    {
      return c.write(os); 
    }
    
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Polyhedron<R>& p) 
    {
      return p.write(os); 
    }
    
    
    template<class R> inline
    std::istream& operator>>(std::istream& os, Polyhedron<R>& p) 
    {
      return p.read(os); 
    }
   



  }
}
