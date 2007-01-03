/***************************************************************************
 *            zonotope.inline.h
 *
 *  6 February 2006
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
    Zonotope<R>::Zonotope()
      : _centre(0), _generators(0,0)
    { 
    }
     
    
    template<class R> inline
    Zonotope<R>::Zonotope(dimension_type d)
      : _centre(d), _generators(d,0)
    { 
    }
     
      
    template<class R> inline
    Zonotope<R>::Zonotope(dimension_type d, size_type m)
      : _centre(d), _generators(d,m)
    { 
    }
     

    template<class R> template<class R1, class R2> inline
      Zonotope<R>::Zonotope(const Point<R1>& c, const LinearAlgebra::Matrix<R2>& g)
        : _centre(c), _generators(g)
      {
        if(c.dimension()!=g.number_of_rows()) { 
//          throw InvalidGenerators("Zonotope<R>::Zonotope(Point<R1>,Matrix<R2>): "
//                                  "The matrix of principal directions does not have the same number of rows as the point dimension.");
          throw InvalidGenerators(__PRETTY_FUNCTION__);
        }
        this->minimize_generators();
      }

 
    template<class R> template<class R1, class R2, class R3> inline
    Zonotope<R>::Zonotope(const Point<R1>& c, const LinearAlgebra::Matrix<R2>& g1, const LinearAlgebra::Vector<R3>& g2)
      : _centre(c), _generators(c.dimension(),g1.number_of_columns()+1) 
    { 
      if(c.dimension()!=g1.number_of_rows() || c.dimension()!=g2.size()) { 
//          throw InvalidGenerators("Zonotope<R>::Zonotope(Point<R1>,Matrix<R2>,Vector<R3>): "
//                                  "The principal directions do not all have the same size as the point dimension.");
        throw InvalidGenerators(__PRETTY_FUNCTION__);
      }
      for(size_type i=0; i!=this->dimension();++i) {
        for(size_type j1=0; j1!=g1.number_of_columns(); ++j1) {
          this->_generators(i,j1)=g1(i,j1);
        }
        this->_generators(i,g1.number_of_columns())=g2(i);
      }
    }
    
    
    template<class R> template<class R1, class R2, class R3> inline
    Zonotope<R>::Zonotope(const Point<R1>& c, const LinearAlgebra::Matrix<R2>& g1, const LinearAlgebra::Matrix<R3>& g2)
      : _centre(c), _generators(c.dimension(),g1.number_of_columns()+g2.number_of_columns()) 
    { 
      if(c.dimension()!=g1.number_of_rows() || c.dimension()!=g2.number_of_rows()) { 
//          throw InvalidGenerators("Zonotope<R>::Zonotope(Point<R1>,Matrix<R2>,Matrix<R3>): "
//                                  "The principal directions do not all have the same size as the point dimension.");
        throw InvalidGenerators(__PRETTY_FUNCTION__);
      }
      for(size_type i=0; i!=this->dimension();++i) {
        for(size_type j1=0; j1!=g1.number_of_columns(); ++j1) {
          this->_generators(i,j1)=g1(i,j1);
        }
        for(size_type j2=0; j2!=g2.number_of_columns(); ++j2) {
          this->_generators(i,g1.number_of_columns()+j2)=g2(i,j2);
        }
      }
    }


    template<class R> template<class Rl> inline
    Zonotope<R>::Zonotope(const Zonotope<Rl>& original)
      : _centre(original.centre()),
        _generators(original.generators())
    { 
    }
    
    
    template<class R> inline
    Zonotope<R>& 
    Zonotope<R>::operator=(const Rectangle<R>& r) {
      this->_centre=r.centre(); this->_generators.resize(r.dimension(),r.dimension());
      for(size_type i=0; i!=r.dimension(); ++i) {
        for(size_type j=0; j!=r.dimension(); ++j) {
          this->_generators(i,j)=0;
        }
        this->_generators(i,i)=div_up(sub_up(r.upper_bound(i),r.lower_bound(i)),R(2));
      }
      return *this;
    }
    
    
    template<class R> inline
    Zonotope<R>& 
    Zonotope<R>::operator=(const Zonotope<R>& original) {
      if(this != &original) {
        this->_centre = original._centre;
        this->_generators = original._generators;
      }
      return *this;
    }



    template<class R> inline
    Point<R> 
    Zonotope<R>::centre() const 
    { 
      return this->_centre; 
    }


    template<class R> inline
    const LinearAlgebra::Matrix<R>& 
    Zonotope<R>::generators() const 
    {
      return this->_generators;
    }
   

    template<class R> inline
    size_type 
    Zonotope<R>::number_of_generators() const 
    {
      return this->_generators.number_of_columns();
    }


    template<class R> inline
    LinearAlgebra::Vector<R> 
    Zonotope<R>::generator(size_type n) const
    {
      return this->_generators.column(n);
    }


    
    template<class R> inline
    dimension_type 
    Zonotope<R>::dimension() const 
    {
      return this->_centre.dimension();
    }
    
    
    template<class R> inline
    tribool 
    Zonotope<R>::empty() const 
    { 
      return false; 
    }
    
    
    template<class R> inline
    tribool 
    Zonotope<R>::bounded() const 
    { 
      return true; 
    }
    
    
    template<class R> inline
    R
    Zonotope<R>::radius() const 
    {
      return this->bounding_box().radius();
    }
    
    
      
      
      
    template<class R> inline
    Zonotope< Interval<R> >::Zonotope(dimension_type d)
      : _centre(d), _generators(d,d) 
    { 
    }
    
      
    template<class R> template<class R1, class R2> inline
    Zonotope< Interval<R> >::Zonotope(const Point<R1>& c, const LinearAlgebra::Matrix<R2>& g)
      : _centre(c), _generators(g) 
    {
    }
      
      
    template<class R> template<class R1, class R2, class R3> inline
    Zonotope< Interval<R> >::Zonotope(const Point<R1>& c, const LinearAlgebra::Matrix<R2>& g1, const LinearAlgebra::Vector<R3>& g2)
      : _centre(c), _generators(c.dimension(),g1.number_of_columns()+1) 
    { 
      for(size_type i=0; i!=this->dimension();++i) {
        for(size_type j1=0; j1!=g1.number_of_columns(); ++j1) {
          this->_generators(i,j1)=g1(i,j1);
        }
        this->_generators(i,g1.number_of_columns())=g2(i);
      }
    }
    
    
    template<class R> template<class R1, class R2, class R3> inline
    Zonotope< Interval<R> >::Zonotope(const Point<R1>& c, const LinearAlgebra::Matrix<R2>& g1, const LinearAlgebra::Matrix<R3>& g2)
      : _centre(c), _generators(c.dimension(),g1.number_of_columns()+g2.number_of_columns()) 
    { 
      for(size_type i=0; i!=this->dimension();++i) {
        for(size_type j1=0; j1!=g1.number_of_columns(); ++j1) {
          this->_generators(i,j1)=g1(i,j1);
        }
        for(size_type j2=0; j2!=g2.number_of_columns(); ++j2) {
          this->_generators(i,g1.number_of_columns()+j2)=g2(i,j2);
        }
      }
    }
    
    
    template<class R> inline
    Zonotope< Interval<R> >::Zonotope(const Zonotope<R>& z)
      : _centre(z.centre()), _generators(z.generators()) 
    {
    }
    
    
    template<class R> inline
    dimension_type 
    Zonotope< Interval<R> >::dimension() const 
    { 
      return this->_centre.dimension(); 
    }
    
    
    template<class R> inline
    size_type 
    Zonotope< Interval<R> >::number_of_generators() const 
    { 
      return this->_generators.number_of_columns();
    }
    
    
    template<class R> inline
    const Point< Interval<R> >& 
    Zonotope< Interval<R> >::centre() const 
    { 
      return _centre; 
    }
    
    
    template<class R> inline
    const LinearAlgebra::Matrix< Interval<R> >& 
    Zonotope< Interval<R> >::generators() const 
    { 
      return _generators; 
    }
    
    
    
    
    
    template<class R> inline
    Zonotope<R> 
    over_approximation(const Zonotope<R>& z) 
    {
      return z;
    }

    template<class R> inline
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_sum(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return Geometry::minkowski_sum(Zonotope<R>(A),B);
    }

    
    template<class R> inline
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_sum(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Geometry::minkowski_sum(A,Zonotope<R>(B));
    }

    
    template<class R> inline
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_difference(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return Geometry::minkowski_difference(Zonotope<R>(A),B);
    }

    
    template<class R> inline
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_difference(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Geometry::minkowski_difference(A,Zonotope<R>(B));
    }
    
    
    
    template<class R> 
    inline
    Zonotope<R> 
    operator+(const Zonotope<R>& z, const LinearAlgebra::Matrix<R>& A) 
    {
      return Zonotope<R>(z.centre(),concatenate_columns(z.generators(),A));
    }
    
    
    template<class R> 
    inline
    Zonotope<typename Zonotope<R>::F> 
    scale(const Zonotope<R>& z, const R& sf) {
      return Zonotope<R>::scale(z,sf);
    }
    
    
    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const Zonotope<R>& z) 
    {
      return z.write(os);
    }
    
    template<class R> inline
    std::istream& operator>>(std::istream& is, Zonotope<R>& z) 
    {
      return z.read(is);
    }





    template<class R>
    class ZonotopeVerticesIterator 
      : public boost::iterator_facade<ZonotopeVerticesIterator<R>,
                                      Point<typename Numeric::traits<R>::arithmetic_type>,
                                      boost::forward_traversal_tag,
                                      Point<typename Numeric::traits<R>::arithmetic_type> const&,
                                      Point<typename Numeric::traits<R>::arithmetic_type> const*
                                     >
    {
      friend class Zonotope<R>;
      typedef typename Numeric::traits<R>::arithmetic_type F;
      const Zonotope<R>* _z; long unsigned int _i; bool _parity; Point<F> _pt;
     public:
      ZonotopeVerticesIterator(const Zonotope<R>& z, bool end) 
        : _z(&z), _i(end ? (1u<<(z.number_of_generators()-1))*3 : 0), _parity(0), _pt(z.centre())
      {
        if(end) { return; }
        for(dimension_type i=0; i!=z.dimension(); ++i) {
          for(dimension_type j=0; j!=z.number_of_generators(); ++j) {
            this->_pt[i]-=z.generators()(i,j); } }
      }
      bool equal(const ZonotopeVerticesIterator<R>& other) const {
        //std::cerr << "ZonotopeVerticesIterator<R>::equal" << std::endl;
        return this->_i==other._i && this->_z==other._z; }
      const Point<F>& dereference() const { 
        //std::cerr << "ZonotopeVerticesIterator<R>::dereference" << std::endl;
        return this->_pt; }
      void increment() { 
        //std::cerr << "ZonotopeVerticesIterator<R>::increment" << std::endl;
        uint j=0; uint m=1; if(this->_parity) { while(!(m&(this->_i))) { ++j; m*=2u; } ++j; m*=2u; }
        this->_parity=!this->_parity;
        if(j==this->_z->number_of_generators()) { this->_i+=m; return; }
        if(m&(this->_i)) { this->_pt=this->_pt-R(2)*this->_z->generator(j); this->_i-=m; }
        else { this->_pt=this->_pt+R(2)*this->_z->generator(j); this->_i+=m; }
      }

      std::ostream& write(std::ostream& os) const {
        return os << _z << " " << _i << " " << _parity << std::endl; }
    }; 

    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const ZonotopeVerticesIterator<R>& iter) {
      return iter.write(os);
    }


    
  

  }
}
