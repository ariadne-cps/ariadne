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
    Polyhedron<R>::Polyhedron(const Polytope<Rl>& p)
    { 
      (*this)=polyhedron(p);
    }
        
    template<class R> template<class Rl> inline 
    Polyhedron<R>::Polyhedron(const Polyhedron<Rl>& p)
      : _dimension(p.dimension()), 
        _number_of_constraints(p.number_of_constraints()), 
        _data(p.data())
    { 
    }
        


    template<class R> inline 
    array<R>& 
    Polyhedron<R>::data() 
    { 
      return this->_data; 
    }
    
    template<class R> inline 
    const array<R>& 
    Polyhedron<R>::data() const
    { 
      return this->_data; 
    }
    
     
    template<class R> inline 
    R* 
    Polyhedron<R>::begin() 
    {
      return this->_data.begin(); 
    }
    

    template<class R> inline 
    const R* 
    Polyhedron<R>::begin() const 
    {
      return this->_data.begin(); 
    }



    template<class R> inline 
    const LinearAlgebra::MatrixSlice<R> 
    Polyhedron<R>::constraints() const 
    { 
      return LinearAlgebra::MatrixSlice<R>(this->_number_of_constraints,
                                           this->_dimension+1u,
                                           const_cast<R*>(this->begin()),
                                           this->_dimension+1u,
                                           1u);
    }
    

    template<class R> inline 
    LinearAlgebra::Matrix<R> 
    Polyhedron<R>::A() const 
    { 
      return -LinearAlgebra::MatrixSlice<R>(this->number_of_constraints(),
                                            this->dimension(),
                                            const_cast<R*>(this->begin()),
                                            this->dimension()+1,1u);
    }
    
    template<class R> inline 
    LinearAlgebra::Vector<R> 
    Polyhedron<R>::b() const 
    { 
      return LinearAlgebra::VectorSlice<R>(this->number_of_constraints(),
                                           const_cast<R*>(this->begin()+this->dimension()),
                                           this->dimension()+1u);
    }
    

    template<class R> inline 
    size_type 
    Polyhedron<R>::number_of_constraints() const 
    { 
      return this->_number_of_constraints;
    }
    
    

    
    
    template<class R> inline
    PolyhedralConstraint<R>::PolyhedralConstraint(const dimension_type d, const R* a)
      : _d(d), _a(a) 
    { 
    }
    
    
    template<class R> inline
    dimension_type
    PolyhedralConstraint<R>::dimension() const
    { 
      return this->_d;
    }
    
    
    template<class R1> template<class R2> inline 
    tribool PolyhedralConstraint<R1>::satisfied_by(const Point<R2>& pt) const
    {
      typedef typename Numeric::traits<R1,R2>::arithmetic_type F;
      F prod=0;
      for(dimension_type i=0; i!=pt.dimension(); ++i) {
        prod+=this->_a[i]*pt[i]; 
      }
      prod+=this->_a[pt.dimension()]; 
      return prod>=0;
    }
    
    
    
    
    
    
    template<class R>
    class PolyhedronConstraintsIterator
      : public boost::iterator_facade<PolyhedronConstraintsIterator<R>,
                                      PolyhedralConstraint<R>,
                                      boost::forward_traversal_tag,
                                      const PolyhedralConstraint<R>&,
                                      const PolyhedralConstraint<R>*
                                     >
    {
     public:
      PolyhedronConstraintsIterator(const Polyhedron<R>& ply, const size_type& n)
        : _c(ply.dimension(),ply.constraints().begin()+n*(ply.dimension()+1)) { }
      bool equal(const PolyhedronConstraintsIterator<R>& other) const { 
        return this->_c._a==other._c._a; }
      const PolyhedralConstraint<R>& dereference() const { return _c; }
      void increment() { _c._a+=_c._d+1u; ; }
     private:
      PolyhedralConstraint<R> _c;
    };
    
    
    
    
    
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const PolyhedralConstraint<R>& c) 
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
