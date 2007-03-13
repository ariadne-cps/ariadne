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
 
#include "../geometry/rectangle_expression.h"
#include "../geometry/rectangle.h"

namespace Ariadne {
  namespace Geometry {
    
    template<class R> inline
    void over_approximation(Zonotope<R>& z, const Rectangle<R>& r) 
    {
      dimension_type d=r.dimension();
      R* cptr=z.begin();
      R* gptr=cptr+d;
      const R* rptr=r.begin();
      for(size_type i=0; i!=d; ++i) {
        cptr[i]=med_approx(r.lower_bound(i),r.upper_bound(i));
        for(size_type j=0; j!=d; ++j) {
          gptr[i*d+j]=0;
        }
        gptr[(d+1)*i]=div_up(sub_up(r.upper_bound(i),r.lower_bound(i)),2);
      }
    }
    
    template<class R> inline
    Zonotope<R>::Zonotope()
      : _dimension(0),_number_of_generators(0), _data(0)
    { 
    }
     
    
    template<class R> inline
    Zonotope<R>::Zonotope(dimension_type d)
      : _dimension(d),_number_of_generators(0), _data(d,static_cast<R>(0))
    { 
    }
     
      
    template<class R> inline
    Zonotope<R>::Zonotope(dimension_type d, size_type m)
      : _dimension(d), _number_of_generators(m), _data(d*(m+1),static_cast<R>(0))
    { 
    }
     


 

    template<class R> template<class E> inline
    Zonotope<R>::Zonotope(const RectangleExpression<E>& re)
      : _dimension(re().dimension()), 
        _number_of_generators(re().dimension()), 
        _data(re().dimension()*(re().dimension()+1u))
    {
      (*this)=re;
    }
    
    
    template<class R> template<class Rl> inline
    Zonotope<R>::Zonotope(const Zonotope<Rl>& z)
      : _dimension(z.dimension()), 
        _number_of_generators(z.number_of_generators()), 
        _data(z.data())
    { 
    }
    
    template<class R> inline
    Zonotope<R>& 
    Zonotope<R>::operator=(const Zonotope<R>& z) {
      if(this != &z) {
        this->_dimension = z._dimension;
        this->_number_of_generators = z._number_of_generators;
        this->_data = z._data;
      }
      return *this;
    }

    
    
    
    template<class R, class E> inline
    void assign_zonotope(LinearAlgebra::VectorSlice<R> c, 
                         LinearAlgebra::MatrixSlice<R> g, 
                         const RectangleExpression<E>& re) 
    {
      typedef typename Numeric::traits<R>::number_type N;
      R zero=static_cast<R>(0);
      N two=static_cast<N>(2);
      const E& r(re());

      g=zero;
      for(size_type i=0; i!=r.dimension(); ++i) {
        c(i)=med_approx(r.lower_bound(i),r.upper_bound(i));
        g(i,i)=div_up(sub_up(r.upper_bound(i),r.lower_bound(i)),two);
      }
    }
      
    template<class R>
    void assign_zonotope(LinearAlgebra::VectorSlice< Interval<R> > c, 
                         LinearAlgebra::MatrixSlice< Interval<R> > g, 
                         const Rectangle<R>& r) 
    {
      c=r.centre().position_vector(); 
      g=static_cast<R>(0);
      for(size_type i=0; i!=r.dimension(); ++i) {
        g(i,i)=(r.upper_bound(i)-r.lower_bound(i))/static_cast<R>(2);
      }
    }
      
    template<class R> template<class E> inline      
    Zonotope<R>& 
    Zonotope<R>::operator=(const RectangleExpression<E>& re) 
    {
      this->resize(re().dimension(),re().dimension());
      assign_zonotope(this->_centre(),this->_generators(),re);
      return *this;
    }
    
    
    template<class R> template<class Rl> inline
    Zonotope<R>& 
    Zonotope<R>::operator=(const Zonotope<Rl>& z) {
      this->_dimension = z.dimension();
      this->_number_of_generators = z.number_of_generators();
      this->_data = z.data();
      return *this;
    }



    template<class R> inline
    const array<R>&
    Zonotope<R>::data() const 
    { 
      return this->_data;
    }

    template<class R> inline
    array<R>&
    Zonotope<R>::data()  
    { 
      return this->_data;
    }

    template<class R> inline
    const R& 
    Zonotope<R>::centre(size_type i) const
    {
      return this->data()[i];
    }

    template<class R> inline
    const R& 
    Zonotope<R>::generators(size_type i, size_type j) const
    {
      return this->data()[(j+1)*this->dimension()+i];

    }  
    

    template<class R> inline
    Point<R> 
    Zonotope<R>::centre() const 
    { 
      return Point<R>(this->dimension(),this->data().begin());
    }

    template<class R> inline
    size_type 
    Zonotope<R>::data_size() const 
    {
      return this->_dimension*(this->_number_of_generators+1); 
    }
      
    template<class R> inline
    void 
    Zonotope<R>::resize(dimension_type d, size_type m) 
    { 
      if(this->_data.size()!=d*(m+1)) { this->_data.resize(d*(m+1)); }
      this->_dimension=d;
      this->_number_of_generators=m;
    }



    template<class R> inline
    LinearAlgebra::VectorSlice<R> 
    Zonotope<R>::_centre()  
    { 
      return LinearAlgebra::VectorSlice<R>(this->dimension(),this->data().begin());
    }


    template<class R> inline
    LinearAlgebra::MatrixSlice<R> 
    Zonotope<R>::_generators() 
    {
      return LinearAlgebra::MatrixSlice<R>(this->dimension(),this->number_of_generators(),
                                           this->data().begin()+this->dimension(),
                                           1u,this->dimension());
    }
   
    template<class R> inline
    const LinearAlgebra::MatrixSlice<R> 
    Zonotope<R>::generators() const 
    {
      return LinearAlgebra::MatrixSlice<R>(this->dimension(),this->number_of_generators(),
                                           const_cast<R*>(this->data().begin()+this->dimension()),
                                           1u,this->dimension());
    }
   

    template<class R> inline
    size_type 
    Zonotope<R>::number_of_generators() const 
    {
      return this->_number_of_generators;
    }


    template<class R> inline
    LinearAlgebra::Vector<R> 
    Zonotope<R>::generator(size_type n) const
    {
      return LinearAlgebra::VectorSlice<R>(this->dimension(),const_cast<R*>(this->_data.begin()+(n+1)*this->dimension()),1);
    }


    
    template<class R> inline
    dimension_type 
    Zonotope<R>::dimension() const 
    {
      return this->_dimension;
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
    ListSet< Zonotope<R> >
    Zonotope<R>::divide() const 
    {
      return Geometry::divide(*this);
    }
    
    template<class R> inline
    ListSet< Zonotope<R> >
    Zonotope<R>::subdivide() const 
    {
      return Geometry::subdivide(*this);
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

      std::ostream& write(std::ostream& os) const;
    };

    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const ZonotopeVerticesIterator<R>& iter) {
      return iter.write(os);
    }


    
  

  }
}
