/***************************************************************************
 *            orthotope.code.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
 
#include <iostream>
#include <vector>
#include <algorithm>

#include "orthotope.h"

#include "../base/array.h"
#include "exceptions.h"
#include "../numeric/conversion.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/lu_matrix.h"
#include "../linear_algebra/qr_matrix.h"

#include "../linear_programming/linear_program.h"

#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"
#include "../geometry/polytope.h"
#include "../geometry/polyhedron.h"
#include "../geometry/parallelotope.h"
#include "../geometry/list_set.h"

#include "../output/logging.h"



namespace {
  
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;


template<class R> inline
void 
accumulate(R& value, R& error, uint n, const R* aptr, const R* bptr) {
  Interval<R> v=value;
  for(uint i=0; i!=n; ++i) {
    v+=aptr[i]*bptr[i];
  }
  value=v.midpoint();
  error=add_up(error,v.radius());
}


template<class R> inline
tribool norm_grtr(const LinearAlgebra::Vector<R>& v1, const LinearAlgebra::Vector<R>& v2) 
{
  return LinearAlgebra::norm(v1)>LinearAlgebra::norm(v2);
}


template<class R>
tribool
disjoint_approx(const Orthotope< Interval<R>, Interval<R> >& z, const Rectangle<R>& r);

template<class R>
tribool
disjoint_approx(const Orthotope<Interval<R>,R>& z, const Rectangle<R>& r);

template<class R>
tribool
disjoint_approx(const Orthotope<R,R>& z, const Rectangle<R>& r);


tribool
disjoint_exact(const Orthotope<Rational,Rational>& z, const Rectangle<Rational>& r);

tribool
disjoint_exact(const Orthotope<Rational,Rational>& z, const Orthotope<Rational,Rational>& r);



template<class R> inline
tribool
contains_approx(const Orthotope<R,R>& z, const Point<R>& r);

template<class R> inline
tribool
contains_approx(const Orthotope<Interval<R>,R>& z, const Point<R>& r);

template<class R> inline
tribool
contains_approx(const Orthotope< Interval<R>,Interval<R> >& z, const Point<R>& r);


tribool
contains_exact(const Orthotope<Rational>& z, const Point<Rational>& r);





template<class R> inline
tribool
superset_approx(const Orthotope<R,R>& z, const Rectangle<R>& r);

template<class R> inline
tribool
superset_approx(const Orthotope<Interval<R>,R>& z, const Rectangle<R>& r);

template<class R> inline
tribool
superset_approx(const Orthotope< Interval<R>,Interval<R> >& z, const Rectangle<R>& r);


tribool
superset_exact(const Orthotope<Rational,Rational>& z, const Rectangle<Rational>& r);




template<class R> inline
void
adjoin_subdivision(ListSet< Orthotope<R,R> >&, const Orthotope<R,R>& z);

template<class R> inline
void
adjoin_subdivision(ListSet< Orthotope<Interval<R>,R> >&, const Orthotope<Interval<R>,R>& z);

template<class R> inline
void
adjoin_subdivision(ListSet< Orthotope< Interval<R>,Interval<R> > >&, const Orthotope< Interval<R>,Interval<R> >& z);


template<class R> inline
void
adjoin_division(ListSet< Orthotope<R,R> >&, const Orthotope<R,R>& z);

template<class R> inline
void
adjoin_division(ListSet< Orthotope<Interval<R>,R> >&, const Orthotope<Interval<R>,R>& z);

template<class R> inline
void
adjoin_division(ListSet< Orthotope< Interval<R>,Interval<R> > >&, const Orthotope< Interval<R>,Interval<R> >& z);



template<class R> inline
void
convert(Rectangle<Rational>& qr, const Rectangle<R>& r) {
  qr=r;
}

template<class R> inline
void
convert(Rectangle<Rational>& qr, const Rectangle< Interval<R> >& r) {
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R> inline
void
convert(Orthotope<Rational>& qz, const Orthotope<R,R>& z) {
  qz=z;
}

template<class R> inline
void
convert(Orthotope<Rational>& qz, const Orthotope<Interval<R>,R>& z) {
  Point< Interval<Rational> > ic=z.centre();
  Point< Rational > c=midpoint(ic);
  Matrix< Rational > G=z.generators();
  Matrix< Rational > E(z.dimension(),z.dimension());
  for(dimension_type i=0; i!=z.dimension(); ++i) {
    E(i,i)=ic[i].radius();
  }
  qz=Orthotope<Rational>(c,G,E);
}

template<class R> inline
void
convert(Orthotope<Rational>& qz, const Orthotope< Interval<R>, Interval<R> >& z) {
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R,class RC,class RG> inline
void
convert(Orthotope<RC,RG>& z, const Rectangle<R>& r) {
  z=r;
}

template<class R,class RC,class RG> inline 
void
convert(Polyhedron<R>& p, const Orthotope<RC,RG>& z) {
  p=Polyhedron<R>(Polytope<R>(z.vertices()));
}

template<class R,class RC,class RG> inline 
void
convert(Polytope<R>& p, const Orthotope<RC,RG>& z) {
  p=Polytope<R>(z.vertices());
}

template<class R,class RC,class RG> inline 
void
convert(Polyhedron< Interval<R> >& p, const Orthotope<RC,RG>& z) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R,class RC,class RG> inline 
void
convert(Polytope< Interval<R> >& p, const Orthotope<RC,RG>& z) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R> void instantiate_orthotope();
template<> void instantiate_orthotope<Numeric::Rational>();

} // namespace 






namespace Ariadne {

extern int Geometry::verbosity;


template<class RC,class RG>
Orthotope<RC,RG>::Orthotope(const std::string& s)
  : _centre(), _generators()
{
  std::stringstream ss(s);
  ss >> *this;
}


template<class RC,class RG>
typename Geometry::Orthotope<RC,RG>::vertices_const_iterator
Geometry::Orthotope<RC,RG>::vertices_begin() const 
{
  return OrthotopeVerticesIterator<RC,RG>(*this,false);
}


template<class RC,class RG>
typename Geometry::Orthotope<RC,RG>::vertices_const_iterator
Geometry::Orthotope<RC,RG>::vertices_end() const 
{
  return OrthotopeVerticesIterator<RC,RG>(*this,true);
}


template<class RC,class RG>
Geometry::PointList<typename Numeric::traits<RC,RG>::arithmetic_type>
Geometry::Orthotope<RC,RG>::vertices() const
{
  //std::clog << "Orthotope<RC,RG>::vertices()" << std::endl;
  PointList<typename Numeric::traits<R>::arithmetic_type> v(this->dimension());
  for(typename Orthotope<RC,RG>::vertices_const_iterator vi=this->vertices_begin();
      vi!=this->vertices_end(); ++vi)
    {
      v.push_back(*vi);
    }
  return v;
}



template<class RC,class RG>
Geometry::Orthotope<RC,RG>::operator Polytope<typename Numeric::traits<RC,RG>::arithmetic_type> () const
{
  return Geometry::polytope(*this);
}

template<class RC,class RG>
Geometry::Orthotope<RC,RG>::operator Polyhedron<typename Numeric::traits<RC,RG>::arithmetic_type> () const
{
  return Geometry::polyhedron(*this);
}


template<class RC,class RG>
Geometry::Polytope<typename Numeric::traits<RC,RG>::arithmetic_type> 
Geometry::polytope(const Orthotope<RC,RG>& z) 
{
  typedef typename Numeric::traits<RC,RG>::arithmetic_type A;
  return Geometry::Polytope<A>(z.vertices());
}

template<class RC,class RG>
Geometry::Polyhedron<typename Numeric::traits<RC,RG>::arithmetic_type> 
Geometry::polyhedron(const Orthotope<RC,RG>& z) 
{
  typedef typename Numeric::traits<RC,RG>::arithmetic_type A;
  return Geometry::Polyhedron<A>(Geometry::Polytope<A>(z.vertices()));
}



template<class RC,class RG>
tribool 
Geometry::equal(const Orthotope<RC,RG>& A, const Orthotope<RC,RG>& B)
{
  return subset(A,B) && subset(B,A);
}


template<class RC,class RG>
Geometry::Rectangle<typename Orthotope<RC,RG>::real_type>
Geometry::bounding_box(const Orthotope<RC,RG>& z)
{
  typedef typename Orthotope<RC,RG>::real_type R;
  using Numeric::Interval;
  LinearAlgebra::Vector< Interval<R> > v(z.number_of_generators(),Interval<R>(-1,1));
  LinearAlgebra::Vector< Interval<R> > b=z.centre().position_vector()+z.generators()*v;
  return Rectangle<R>(b);
}



template<class RC,class RG>
Geometry::Rectangle<typename Geometry::Orthotope<RC,RG>::R>
Geometry::Orthotope<RC,RG>::bounding_box() const
{
  return Geometry::bounding_box(*this);
}



template<class RC, class RG>
Geometry::ListSet< Geometry::Orthotope<RC,RG> >
Geometry::subdivide(const Orthotope<RC,RG>& z) 
{
  ListSet< Orthotope<RC,RG> > result(z.dimension());
  ::adjoin_subdivision(result,z);
  return result;
}


template<class RC, class RG>
Geometry::ListSet< Geometry::Orthotope<RC,RG> >
Geometry::divide(const Orthotope<RC,RG>& z) 
{
  ListSet< Orthotope<RC,RG> > result(z.dimension());
  ::adjoin_division(result,z);
  return result;
}


template<class RC,class RG>
Geometry::ListSet< Geometry::Orthotope<RC,RG> >
Geometry::Orthotope<RC,RG>::divide() const 
{
  return Geometry::divide(*this);
}

template<class RC,class RG>
Geometry::ListSet< Geometry::Orthotope<RC,RG> >
Geometry::Orthotope<RC,RG>::subdivide() const 
{
  return Geometry::subdivide(*this);
}









template<class RC,class RG>
void 
Geometry::Orthotope<RC,RG>::minimize_generators(void) 
{
  return;
}



template<class RC,class RG>
void 
Geometry::Orthotope<RC,RG>::sort_generators(void)
{
  std::vector< LinearAlgebra::Vector<RG> > generator_vectors;
  for(size_type j=0; j!=this->number_of_generators(); ++j) {
    generator_vectors.push_back(this->generators().column(j));
  }
  
  std::stable_sort(generator_vectors.begin(),generator_vectors.end(),&norm_grtr<RG>);
  
  for(size_type j=0; j!=this->number_of_generators(); ++j) {
    for(size_type i=0; i!=this->dimension(); ++i) {
      this->_generators(i,j)=generator_vectors[j](i);
    }
  }
}



template<class RC,class RG> template<class R0>
tribool 
Geometry::Orthotope<RC,RG>::contains(const Point<R0>& pt) const 
{
  return ::contains_approx(*this,pt);
}


template<class R,class RC,class RG>
tribool 
Geometry::contains(const Orthotope<RC,RG>& z, const Point<R>& pt)
{
  return ::contains_approx(z,pt);
}


template<class R,class RC,class RG>
tribool
Geometry::disjoint(const Orthotope<RC,RG>& z, const Rectangle<R>& r)
{
  return disjoint_approx(z,r);
}

template<class R,class RC,class RG>
tribool
Geometry::disjoint(const Rectangle<R>& r, const Orthotope<RC,RG>& z)
{
  return disjoint_approx(z,r);
}

template<class RC,class RG>
tribool
Geometry::disjoint(const Orthotope<RC,RG>& z1, const Orthotope<RC,RG>& z2)
{
  ARIADNE_LOG(6,"disjoint(Orthotope z1, Orthotope z2)\n");
  ARIADNE_LOG(7,"z1="<<z1<<", z2="<<z2<<"\n");
  Geometry::Orthotope<Numeric::Rational> qz1;
  Geometry::Orthotope<Numeric::Rational> qz2;
  convert(qz1,z1);
  convert(qz2,z2);
  
  return disjoint_exact(qz1,qz2);
}


template<class R,class RC,class RG>
tribool 
Geometry::subset(const Rectangle<R>& r, const Orthotope<RC,RG>& z) 
{
  return ::superset_approx(z,r);
}


template<class R,class RC,class RG>
tribool 
Geometry::subset(const Orthotope<RC,RG>& z, const Rectangle<R>& r) 
{
  return Geometry::subset(z.bounding_box(),r);
}


template<class RC,class RG>
tribool 
Geometry::subset(const Orthotope<RC,RG>& z1, const Orthotope<RC,RG>& z2) 
{
  typedef typename Numeric::traits<RC,RG>::arithmetic_type F;
  return Geometry::subset(z1.operator Polytope<F>(),z2.operator Polyhedron<F>());
} 


template<class R,class RC,class RG>
void 
Geometry::over_approximation(Orthotope<RC,RG>& z, const Rectangle<R>& r) 
{
  dimension_type d=r.dimension();
  RC* cptr=z.begin();
  RG* gptr=cptr+d;
  const R* rptr=r.begin();
  for(size_type i=0; i!=d; ++i) {
    cptr[i]=med_approx(r.lower_bound(i),r.upper_bound(i));
    for(size_type j=0; j!=d; ++j) {
      gptr[i*d+j]=0;
    }
    gptr[(d+1)*i]=div_up(sub_up(r.upper_bound(i),r.lower_bound(i)),2);
  }
}


template<class RC,class RG> 
Geometry::Orthotope<typename Numeric::traits<RC>::arithmetic_type,RG>
Geometry::minkowski_sum(const Orthotope<RC,RG>& z1, const Orthotope<RC,RG>& z2)
{
  typedef typename Numeric::traits<RC>::arithmetic_type F;
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(z1,z2,"Orthotope minkowski_sum(Orthotope z1, Orthotope z2)");
  
  return Orthotope<F,RG>(Geometry::minkowski_sum(z1.centre(),z2.centre()),
                        LinearAlgebra::concatenate_columns(z1.generators(),z2.generators()));
}


template<class RC,class RG> 
Geometry::Orthotope<typename Numeric::traits<RC>::arithmetic_type,RG> 
Geometry::minkowski_difference(const Orthotope<RC,RG>& z1, const Orthotope<RC,RG>& z2)
{
  typedef typename Numeric::traits<RC>::arithmetic_type F;
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(z1,z2,"Orthotope minkowski_difference(Orthotope z1, Orthotope z2)");
  
  return Orthotope<F,RG>(Geometry::minkowski_difference(z1.centre(),z2.centre()),
                        LinearAlgebra::concatenate_columns(z1.generators(),z2.generators()));
}



template<class R,class RC,class RG>
Geometry::Orthotope<typename Numeric::traits<R,RC>::arithmetic_type,typename Numeric::traits<R,RG>::arithmetic_type> 
Geometry::minkowski_sum(const Rectangle<R>& r, const Orthotope<RC,RG>& z) 
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
  return Geometry::minkowski_sum(Orthotope<F>(r),z);
}


template<class R,class RC,class RG>
Geometry::Orthotope<typename Numeric::traits<R,RC>::arithmetic_type,typename Numeric::traits<R,RG>::arithmetic_type> 
Geometry::minkowski_sum(const Orthotope<RC,RG>& z, const Rectangle<R>& r) 
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
  return Geometry::minkowski_sum(z,Orthotope<F>(r));
}


template<class R,class RC,class RG>
Geometry::Orthotope<typename Numeric::traits<R,RC>::arithmetic_type,typename Numeric::traits<R,RG>::arithmetic_type> 
Geometry::minkowski_difference(const Rectangle<R>& r, const Orthotope<RC,RG>& z) 
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
  return Geometry::minkowski_difference(Orthotope<F>(r),z);
}


template<class R,class RC,class RG>
Geometry::Orthotope<typename Numeric::traits<R,RC>::arithmetic_type,typename Numeric::traits<R,RG>::arithmetic_type> 
Geometry::minkowski_difference(const Orthotope<RC,RG>& z, const Rectangle<R>& r) 
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
  return Geometry::minkowski_difference(z,Orthotope<F>(r));
}




template<class R> 
Geometry::Orthotope<Numeric::Interval<R>,R> 
Geometry::over_approximation(const Orthotope< Numeric::Interval<R>, Numeric::Interval<R> >& iz)
{
  // FIXME: This is incorrect; need over-approximations
  LinearAlgebra::Matrix< Numeric::Interval<R> > G(iz.generators());
  LinearAlgebra::Vector< Numeric::Interval<R> > e(iz.number_of_generators(),Numeric::Interval<R>(-1,+1));
  LinearAlgebra::Matrix<R> nG=midpoint(G);
  Geometry::Point< Numeric::Interval<R> > nc=iz.centre()+(G-nG)*e;
  
  return Orthotope< Numeric::Interval<R>, R >(nc,nG);
}

template<class R> 
Geometry::Orthotope<R,R> 
Geometry::over_approximation(const Orthotope<Numeric::Interval<R>, R>& ez)
{
  dimension_type d=ez.dimension();
  size_type ng=ez.number_of_generators();
  Point<R> c=midpoint(ez.centre());
  Matrix<R> g(d,ng+d);
  MatrixSlice<R>(d,ng,g.begin(),g.row_increment(),g.column_increment())=ez.generators();
  for(size_type i=0; i!=d; ++i) {
    g(i,i+ng)=ez.centre(i).radius();
  }
  return Orthotope<R,R>(c,g);
}



template<class R> 
Geometry::Orthotope<R,R> 
Geometry::over_approximation(const Orthotope<R,R>& z)
{
  return z;
}




template<class R> 
Geometry::Orthotope<R,R> 
Geometry::approximation(const Orthotope< Numeric::Interval<R>, Numeric::Interval<R> >& z)
{
  return Orthotope<R,R>(midpoint(z.centre()),
                       midpoint(z.generators()));
}


template<class R> 
Geometry::Orthotope<R,R> 
Geometry::approximation(const Orthotope<Numeric::Interval<R>,R>& z)
{
  return Orthotope<R,R>(midpoint(z.centre()),
                       z.generators());
}


template<class R> 
Geometry::Orthotope<R,R> 
Geometry::approximation(const Orthotope<R,R>& z)
{
  return z; 
}




template<class RC,class RG>
std::string
Geometry::Orthotope<RC,RG>::name()
{
  return std::string("Orthotope")+"<"+Numeric::name<RC>()+","+Numeric::name<RG>()+">";
}


template<class RC,class RG>
std::ostream&
Geometry::Orthotope<RC,RG>::write(std::ostream& os) const 
{
  const Orthotope<RC,RG>& z=*this;
  os << "["<<z.centre();
  for(size_type j=0; j!=z.number_of_generators(); ++j) {
    os << ";" << z.generator(j);
  }
  os << "]";
  /*
  if(z.dimension() > 0) {
    if (z.empty()) {
      os << "Orthotope( )" << std::endl;
    } else {
      os << "Orthotope( centre=" << z.centre();
      os << ", directions=";
      for(size_type i=0; i!=z.number_of_generators(); ++i) {
        os << (i==0 ? "[ " : ", ") << z.generator(i);
      }
      os << "] )";
    } 
  }
  */
  return os;
}


template<class RC,class RG>
std::istream& 
Geometry::Orthotope<RC,RG>::read(std::istream& is)
{
  Point<RC> centre;
  LinearAlgebra::Matrix<RG> generators;
  char c0,c1,c2;
  is >> c0 >> centre >> c1 >> generators >> c2;
  *this = Orthotope<RC,RG>(centre,generators);
  return is;
}








template<class RC,class RG>
std::ostream& 
Geometry::OrthotopeVerticesIterator<RC,RG>::write(std::ostream& os) const 
{
  return os << "OrthotopeVerticesIterator<" << Numeric::name<RC>() << "," << Numeric::name<RG>() << ">" 
            << "( &z=" << _z << ","
            << " i=" << _i << " p=" << _parity << ", v=" << **this << ")";
}


template<class RC,class RG>
void
Geometry::Orthotope<RC,RG>::_instantiate_geometry_operators() 
{
  typedef typename Numeric::traits<RC>::number_type R;
  ::instantiate_orthotope<R>();
}




}
                                                            
                                                            
                                                            

namespace {

using Geometry::verbosity;

template<class R>
void
adjoin_subdivision(ListSet< Orthotope<R,R> >& ls, const Orthotope<R,R>& z) 
{
  R two=2;
  dimension_type d=z.dimension();
  size_type m=z.number_of_generators();
  
  LinearAlgebra::Matrix<R> new_generators(d,m);
  for(size_type i=0; i!=d; ++i) {
    for(size_type j=0; j!=m; ++j) {
      new_generators(i,j)=div_up(z.generators()(i,j),two);
    }
  }
  
  Point<R> first_centre=z.centre();
  for(size_type i=0; i<m; i++) {
    first_centre=sub_approx(first_centre,LinearAlgebra::Vector<R>(new_generators.column(i)));
  }
  
  for(unsigned long k=0; k!=1u<<m; ++k) {
    Point<R> new_centre=first_centre;
    for(size_type i=0; i<m; i++) {
      if(k & 1u<<i) {
        new_centre=add_approx(new_centre,z.generator(i));
      }
    }
    ls.adjoin(Orthotope<R,R>(new_centre,new_generators));
  }
}


template<class R>
void
adjoin_subdivision(ListSet< Orthotope<Interval<R>,R> >& ls, const Orthotope<Interval<R>,R>& z) 
{
  // FIXME: Subdivide in zero order as well!
  static bool warn=true;
  if(warn) {
    std::cerr << std::endl << "WARNING: Orthotope<I,R>::subdivision()  does not subdivide the constant terms." << std::endl;
    warn=false;
  }

  typedef Numeric::Interval<R> I;
  
  R two=2;
  dimension_type d=z.dimension();
  size_type m=z.number_of_generators();
  
  LinearAlgebra::Matrix<R> new_generators(d,m);
  for(size_type i=0; i!=d; ++i) {
    for(size_type j=0; j!=m; ++j) {
      new_generators(i,j)=div_up(z.generators()(i,j),two);
    }
  }
  
  Vector<R> v(z.dimension());
  Point<I> first_centre=z.centre();
  for(size_type i=0; i<m; i++) {
    v=new_generators.column(i);
    first_centre=first_centre-v;
  }
  
  for(unsigned long k=0; k!=1u<<m; ++k) {
    Point<I> new_centre=first_centre;
    for(size_type i=0; i<m; ++i) {
      if(k & 1u<<i) {
        v=z.generator(i);
        new_centre=new_centre+v;
      }
    }
    ls.adjoin(Orthotope<I,R>(new_centre,new_generators));
  }
}


template<class R>
void
adjoin_subdivision(ListSet< Orthotope< Interval<R> > >& ls, const Orthotope< Interval<R>,Interval<R> >& z) 
{
  // FIXME: Subdivide in zero order as well!
  static bool warn=true;
  if(warn) {
    std::cerr << std::endl << "WARNING: Orthotope<I,I>::subdivision()  does not subdivide the constant terms." << std::endl;
    warn=false;
  }

  typedef Numeric::Interval<R> I;
  
  dimension_type d=z.dimension();
  size_type m=z.number_of_generators();
  
  LinearAlgebra::Matrix<I> new_generators(d,m);
  for(size_type i=0; i!=d; ++i) {
    for(size_type j=0; j!=m; ++j) {
      new_generators(i,j)=z.generators()(i,j)/2;
    }
  }
  
  Vector<I> v(z.dimension());
  Point<I> first_centre=z.centre();
  for(size_type i=0; i<m; i++) {
    v=new_generators.column(i);
    first_centre=first_centre-v;
  }
  
  Point<I> new_centre=first_centre;
  for(unsigned long k=0; k!=1u<<m; ++k) {
    Point<I> new_centre=first_centre;
    for(size_type i=0; i<m; ++i) {
      if(k & 1u<<i) {
        v=z.generator(i);
        new_centre=new_centre+v;
      }
    }
    ls.adjoin(Orthotope<I,I>(new_centre,new_generators));
  }

}



template<class R>
void
adjoin_division(ListSet< Orthotope<R,R> >& ls, const Orthotope<R,R>& z)
{
  size_type d=z.dimension();
  size_type m=z.number_of_generators();
  R two=2;
  
  LinearAlgebra::Matrix<R> new_generators=z.generators();
  
  R max_norm=0;
  size_type max_column=0;
  for(size_type j=0; j<m; j++) {
    R norm = LinearAlgebra::norm(LinearAlgebra::Vector<R>(new_generators.column(j)));
    if(norm>max_norm) {
      max_norm=norm;
      max_column=j;
    }
  }
  
  size_type j=max_column;
  for(size_type i=0; i!=d; ++i) {
    new_generators(i,j)=div_up(new_generators(i,j),two);
  }
  
  Point<R> new_centre=sub_approx(z.centre(),LinearAlgebra::Vector<R>(new_generators.column(j)));
  ls.adjoin(Orthotope<R,R>(new_centre,new_generators));
  new_centre=sub_approx(new_centre,LinearAlgebra::Vector<R>(new_generators.column(j)));
  ls.adjoin(Orthotope<R,R>(new_centre,new_generators));
}


template<class R>
void
adjoin_division(ListSet< Orthotope<Interval<R>,R> >& ls, const Orthotope<Numeric::Interval<R>,R>& z)
{
  typedef Numeric::Interval<R> I;
  using namespace LinearAlgebra;
  
  size_type d=z.dimension();
  size_type m=z.number_of_generators();
  R two=2;
  
  I max_radius=0;
  size_type max_value=0;
  for(size_type j=0; j<d; ++j) {
    R radius = z.centre()[j].radius();
    if(radius>max_radius) {
      max_radius=radius;
      max_value=j;
    }
  }
  
  I max_norm=0;
  size_type max_column=0;
  for(size_type j=0; j<m; ++j) {
    I norm = LinearAlgebra::norm(Vector<R>(z.generators().column(j)));
    if(norm>max_norm) {
      max_norm=norm;
      max_column=j;
    }
  }
  
  if(max_norm>max_radius) {
    LinearAlgebra::Matrix<R> new_generators=z.generators();
    size_type j=max_column;
    for(size_type i=0; i!=d; ++i) {
      new_generators(i,j)=div_up(new_generators(i,j),two);
    }
    
    Vector<R> v=new_generators.column(j);
    Point<I> new_centre=z.centre()-v;
    ls.adjoin(Orthotope<I,R>(new_centre,new_generators));
    new_centre=new_centre-v;
    ls.adjoin(Orthotope<I,R>(new_centre,new_generators));
  } else {
    Point<I> new_centre = z.centre();
    const LinearAlgebra::Matrix<R>& new_generators = z.generators();
    new_centre[max_value]=I(z.centre()[max_value].lower(),z.centre()[max_value].midpoint());
    ls.adjoin(Orthotope<I,R>(new_centre,new_generators));
    new_centre[max_value]=I(z.centre()[max_value].midpoint(),z.centre()[max_value].upper());
    ls.adjoin(Orthotope<I,R>(new_centre,new_generators));
  }
}


template<class R>
void
adjoin_division(ListSet< Orthotope< Interval<R>,Interval<R> > >& ls, const Orthotope< Interval<R>,Interval<R> >& z)
{
  typedef Numeric::Interval<R> I;
  using namespace LinearAlgebra;
  
  size_type d=z.dimension();
  size_type m=z.number_of_generators();
  
  I max_radius=0;
  size_type max_value=0;
  for(size_type j=0; j<d; ++j) {
    R radius = z.centre()[j].radius();
    if(radius>max_radius) {
      max_radius=radius;
      max_value=j;
    }
  }
  
  I max_norm=0;
  size_type max_column=0;
  for(size_type j=0; j<m; ++j) {
    I norm = LinearAlgebra::norm(Vector<I>(z.generators().column(j)));
    if(norm>max_norm) {
      max_norm=norm;
      max_column=j;
    }
  }
  
  if(max_norm>max_radius) {
    LinearAlgebra::Matrix<I> new_generators=z.generators();
    size_type j=max_column;
    for(size_type i=0; i!=d; ++i) {
      new_generators(i,j)=new_generators(i,j)/2;
    }
    
    Vector<I> v=new_generators.column(j);
    Point<I> new_centre=z.centre()-v;
    ls.adjoin(Orthotope<I,I>(new_centre,new_generators));
    new_centre=new_centre-v;
    ls.adjoin(Orthotope<I,I>(new_centre,new_generators));
  } else {
    Point<I> new_centre = z.centre();
    const LinearAlgebra::Matrix<I>& new_generators = z.generators();
    new_centre[max_value]=I(z.centre()[max_value].lower(),z.centre()[max_value].midpoint());
    ls.adjoin(Orthotope<I,I>(new_centre,new_generators));
    new_centre[max_value]=I(z.centre()[max_value].midpoint(),z.centre()[max_value].upper());
    ls.adjoin(Orthotope<I,I>(new_centre,new_generators));
  }
}







template<class R> inline
tribool 
contains_approx(const Orthotope<R,R>& z, const Point<R>& pt)
{
  Orthotope<Rational> qz(z);
  Point<Rational> qp(pt);
  return ::contains_exact(qz,qp);
}

template<class R> inline
tribool 
contains_approx(const Orthotope<Interval<R>,R>& z, const Point<R>& pt)
{
  Orthotope<Rational> qz;
  Point<Rational> qp(pt);
  convert(qz,z);
  return ::contains_exact(qz,qp);
}

template<class R> inline
tribool 
contains_approx(const Orthotope<Interval<R>,Interval<R> >& z, const Point<R>& pt)
{
  throw NotImplemented("tribool contains(Orthotope<I,I>,Point<R>): Undefined semantics\n");
}



template<class R> inline
tribool 
superset_approx(const Orthotope<R,R>& z, const Rectangle<R>& r)
{
  ARIADNE_LOG(6,"subset(Rectangle<R> r, Orthotope<R,R> z)\n");
  ARIADNE_LOG(7,"z="<<z<<", r="<<r<<"\n");
  Rectangle<Rational> qr(r);
  Orthotope<Rational> qz(z);
  return superset_exact(qz,qr);
}


template<class R> inline
tribool 
superset_approx(const Orthotope<Interval<R>,R>& z, const Rectangle<R>& r)
{
  ARIADNE_LOG(6,"subset(Rectangle<R> r, Orthotope<I,R> z)\n");
  ARIADNE_LOG(7,"z="<<z<<", r="<<r<<"\n");
  Rectangle<Rational> qr(r);
  Orthotope<Rational> qz;
  convert(qz,z);
  return superset_exact(qz,qr);
}


template<class R> inline
tribool 
superset_approx(const Orthotope<Interval<R>,Interval<R> >& z, const Rectangle<R>& r)
{
  throw NotImplemented("tribool subset(Rectangle<R>, Orthotope<I,I>): Undefined semantics\n");
}






template<class R>
tribool
disjoint_approx(const Orthotope<R,R>& z, const Rectangle<R>& r) {
  ARIADNE_LOG(6,"disjoint(Orthotope<R,R> z, Rectangle<R> r)\n");
  ARIADNE_LOG(7,"z="<<z<<", r="<<r<<"\n");
  return disjoint_exact(Orthotope<Rational,Rational>(z),Rectangle<Rational>(r));
}


template<class R>
tribool
disjoint_approx(const Orthotope<Interval<R>,R>& z, const Rectangle<R>& r) {
  ARIADNE_LOG(6,"disjoint(Orthotope<I,R> z, Rectangle<R> r)\n");
  ARIADNE_LOG(7,"z="<<z<<", r="<<r<<"\n");
  Point< Interval<Rational> > ic=z.centre();
  Matrix<Rational> g=z.generators();
  Point<Rational> c=midpoint(ic);
  Rectangle<Rational> qr=r;
  Rectangle<Rational> xr=qr+(ic-c);
  return disjoint_exact(Orthotope<Rational>(c,g),xr);
}


template<class R>
tribool
disjoint_approx(const Orthotope< Interval<R>, Interval<R> >& z, const Rectangle<R>& r)
{
  ARIADNE_LOG(6,"disjoint(Orthotope<I,I> z, Rectangle<R> r)\n");
  ARIADNE_LOG(7,"z="<<z<<", r="<<r<<"\n");
  Orthotope<Interval<R>,R> oaz=over_approximation(z);
  return disjoint_approx(oaz,r);
}





/* Test vertices individually. */
tribool 
superset_exact(const Orthotope<Rational>& z, const Rectangle<Rational>& r)
{
  tribool result=true;
  for(Rectangle<Rational>::vertices_const_iterator rv_iter=r.vertices_begin(); 
      rv_iter!=r.vertices_end(); ++rv_iter) 
  {
    const Point<Rational>& pt=*rv_iter;
    result=result && z.contains(pt);
    if(result==false) {
      break;
    }
  }
  return result;
}



/* Set up linear program to solve 
 *   \f[x=c+Ge;\ l\leq x\leq u;\ -1\leq e\leq1\f].
 *
 * Change variables to normalize \f$x\f$ and \f$e\f$
 *   \f[x'=x-l,\ e'=e+1;   x'-Ge' = c-G1-l;  0\leq x\leq u-l; \ 0\leq e\leq 2.\f] 
 * 
 * Introduce slack variables sx and se, and artificial variables ax. Problem in standard form
 *   \f[ \begin{matrix}I&0\\0&I\\\pm I&\mp G\end{matrix} \begin{matrix}x'\\e'\end{matrix}
 *        + \begin{matrix}I&&\\&I&\\&&I\end{matrix}\begin{matrix}sx\\se\\ax\end{matrix}
 *             = \begin{matrix}u-l\\2\\\pm(c-G1-l)\end{matrix} \f]
 * 
 */
tribool
disjoint_exact(const Orthotope<Rational,Rational>& z, const Rectangle<Rational>& r)
{
  ARIADNE_LOG(8,"disjoint(Orthotope<Rational> q, Rectangle<Rational> r)\n");
  ARIADNE_LOG(9,"z="<<z<<", r="<<r<<"\n");
  
  ARIADNE_CHECK_EQUAL_DIMENSIONS(z,r,"tribool disjoint(Orthotope<Rational> z, Rectangle<Rational> r)");
  dimension_type d=z.dimension();
  size_type m=z.number_of_generators();
  
  // Construct tableau for testing intersection of orthotope and rectangle
  // Rectangle  l<=x<=u
  // Orthotope  x==c+Ge,  -1<=e<=1
  // 
  // Translate x'=x-l,  e'=e+1
  //   0<=x'<=u-l      ->  x' +     + sx'               == u-l
  //   0<=e'<=2        ->     +  e' +     + se'         == 2
  //   x'+l==c+G(e'-1) ->  x' + Ge'             +/- ax' == c-l-G1
  //  
  // Change sign of RHS of first equality if necessary
  // Introduce slack variables for last two inequalities
  typedef Numeric::Rational Q;
  LinearAlgebra::Matrix<Q> T(2*d+m+1,d+m+1);
  
  const Geometry::Point<Q>& l=r.lower_corner();
  const Geometry::Point<Q>& u=r.upper_corner();
  const Geometry::Point<Q>& c=z.centre();
  const LinearAlgebra::Matrix<Q>& G=z.generators();
  
  const LinearAlgebra::Vector<Q> qo(m,Q(1));
  const LinearAlgebra::Vector<Q> ql=l.position_vector();
  const LinearAlgebra::Vector<Q> qu=u.position_vector();
  const LinearAlgebra::Vector<Q> qd=qu-ql;
  const LinearAlgebra::Vector<Q> qc=c.position_vector();
  const LinearAlgebra::Matrix<Q> qG=G;
  const LinearAlgebra::Vector<Q> qrhs=qc-ql-qG*qo;
  
  if(Geometry::verbosity>8) { std::clog << "ql=" << ql << ", qd=" << qd <<", qc=" << qc << ", qrhs=" << qrhs << std::endl; }
  
  // Set up constraints x+sx=u-l
  for(size_type i=0; i!=d; ++i) {
    T(i,i)=1;
    T(i,d+m)=qu(i)-ql(i);
  }
  
  // Set up constraints e+se=2
  for(size_type j=0; j!=m; ++j) {
    T(d+j,d+j)=1;
    T(d+j,d+m)=2;
  }
  
  // Set up constraints x-Ge \pm ax=c-l-G1 
  for(size_type i=0; i!=d; ++i) {
    if(qrhs(i)>=Q(0)) {
      T(i+d+m,i)=1;
      for(size_type j=0; j!=m; ++j) {
        T(i+d+m,d+j)=-qG(i,j);
      }
      T(i+d+m,d+m)=qrhs(i);
    }
    else {
      T(i+d+m,i)=-1;
      for(size_type j=0; j!=m; ++j) {
        T(i+d+m,d+j)=qG(i,j);
      }
      T(i+d+m,d+m)=-qrhs(i);
    }
  } 
  
  // Set up cost function ax^T1
  for(size_type i=0; i!=d; ++i) {
    T(2*d+m,i) -= T(i+d+m,i);
    for(size_type j=0; j!=m; ++j) {
      T(2*d+m,d+j) -= T(i+d+m,d+j);
    }
    T(2*d+m,d+m) -= T(i+d+m,d+m);
  }
  
  LinearProgramming::LinearProgram<Q> lp(T);
  tribool result=!lp.is_feasible();
  
  return result;
}



tribool
disjoint_exact(const Orthotope<Rational,Rational>& z1, const Orthotope<Rational,Rational>& z2)
{
  typedef Numeric::Rational Q;
  ARIADNE_CHECK_EQUAL_DIMENSIONS(z1,z2,"tribool disjoint(Orthotope<Rational> z1, Orthotope<Rational> z2)");
  
  dimension_type d=z1.dimension();
  Q one=1;
  size_type m1=z1.number_of_generators();
  size_type m2=z2.number_of_generators();
  
  LinearAlgebra::Matrix<Q> T(m1+m2+d+1,m1+m2+1);
  
  const LinearAlgebra::Vector<Q> qo1(m1,one);
  const LinearAlgebra::Vector<Q> qo2(m2,one);
  
  const Geometry::Point<Q>& qc1=z1.centre();
  const LinearAlgebra::Matrix<Q>& qG1=z1.generators();
  const Geometry::Point<Q>& qc2=z2.centre();
  const LinearAlgebra::Matrix<Q>& qG2=z2.generators();
  LinearAlgebra::Vector<Q> qrhs = qG1*qo1 - qG2*qo2 + (qc2 - qc1);
  
  // Set up constraints e1 + se1 = 2
  for(size_type j1=0; j1!=m1; ++j1) {
    T(j1,j1)=1;
    T(j1,m1+m2)=2;
  }
  
  // Set up constraints e2 + se2 = 2
  for(size_type j2=0; j2!=m2; ++j2) {
    T(m1+j2,m1+j2)=1;
    T(m1+j2,m1+m2)=2;
  }
  
  // Set up constraints G1*e1 - G2*e2 = (c2 - G2*1) - (c1 - G1*1)
  for(size_type i=0; i!=d; ++i) {
    if(qrhs(i)>=Q(0)) {
      for(size_type j1=0; j1!=m1; ++j1) {
        T(m1+m2+i,j1)=qG1(i,j1);
      }
      for(size_type j2=0; j2!=m2; ++j2) {
        T(m1+m2+i,m1+j2)=qG2(i,j2);
      }
      T(m1+m2+i,m1+m2)=qrhs(i);
    }
    else {
      for(size_type j1=0; j1!=m1; ++j1) {
        T(m1+m2+i,j1)=-qG1(i,j1);
      }
      for(size_type j2=0; j2!=m2; ++j2) {
        T(m1+m2+i,m1+j2)=-qG2(i,j2);
      }
      T(m1+m2+i,m1+m2)=-qrhs(i);
    }
  } 
  
  // Set up cost function ax^T1
  for(size_type i=0; i!=d; ++i) {
    for(size_type j1=0; j1!=m1; ++j1) {
      T(m1+m2+d,j1) -= T(m1+m2+i,j1);
    }
    for(size_type j2=0; j2!=m2; ++j2) {
      T(m1+m2+d,m1+j2) -= T(m1+m2+i,m1+j2);
    }
    T(m1+m2+d,m1+m2) -= T(m1+m2+i,m1+m2);
  }
  
  LinearProgramming::LinearProgram<Q> lp(T);
  
  tribool result=!lp.is_feasible();
  
  //std::clog << "disjoint(" << z1 << "," << z2 << ")=" << result << std::endl;
  return result;
}

/* Set up LP problem to solve \f$c+Ge=p\f$; \f$-1<=e<=1\f$.
 * Change variables so that the problem becomes \f$Ge=p-c-G1;\ 0\leq e\leq2\f$.
 * Change sign of \f$ Ge=p-c-G1\f$ to make right-hand side positive.
 */
tribool 
contains_exact(const Orthotope<Rational,Rational>& z, const Point<Rational>& pt)
{ 
  //std::clog << "Orthotope<RC,RG>::contains(const Point<R>& )" << std::endl;
  //typedef typename Numeric::traits<R,R>::arithmetic_type Q;
  typedef Rational Q;
  ARIADNE_CHECK_EQUAL_DIMENSIONS(z,pt,"tribool contains(Orthotope<Rational> z, Point<Rational> pt)");
  dimension_type d=z.dimension();
  dimension_type m=z.number_of_generators();
  
  Q zero=0;
  Q one=1;
  Q two=2;
  
  const Geometry::Point<Q>& qc=z.centre();
  const Geometry::Point<Q>& qp=pt;
  const LinearAlgebra::Matrix<Q>& qG=z.generators();
  const LinearAlgebra::Vector<Q> qo(m,one);
  const LinearAlgebra::Vector<Q> zv(m,zero);
  const LinearAlgebra::Vector<Q> tv(m,two);
  
  LinearAlgebra::Vector<Q> qrhs=qp-qc+qG*qo;
  
  LinearAlgebra::Matrix<Q> T(d+m+1,m+1);
  
  // Set up constraints e+se=2
  for(dimension_type j=0; j!=m; ++j) {
    T(j,j)=1;
    T(j,m)=2;
  }
  
  // Set up constraints Ge \pm ax = p-c+G1
  for(dimension_type i=0; i!=d; ++i) {
    if(qrhs(i)>=Q(0)) {
      for(dimension_type j=0; j!=m; ++j) {
        T(m+i,j)=qG(i,j); 
      }
      T(m+i,m)=qrhs(i);
    } else {
      for(dimension_type j=0; j!=m; ++j) {
        T(m+i,j)=-qG(i,j); 
      }
      T(m+i,m)=-qrhs(i);
    }
  }
  
  // Set up cost function ax^T 1
  for(dimension_type i=0; i!=d; ++i) {
    for(dimension_type j=0; j!=m; ++j) {
      T(m+d,j)-=T(m+i,j);
    }
    T(m+d,m)-=T(m+i,m);
  }
  
  LinearProgramming::LinearProgram<Q> lp(T);
  //std::clog << lp.tableau() << std::endl;
  tribool result=lp.is_feasible();
  //std::clog << lp.tableau() << std::endl;
  return result;
}

template<class R> inline
void
instantiate_orthotope()
{
  typedef Numeric::Interval<R> I;
  tribool tb;
  Geometry::Point<R> pt;
  Geometry::Rectangle<R> r;
  Geometry::Orthotope<I,R> ez;
  tb=Geometry::contains(ez,pt);
  tb=Geometry::disjoint(r,ez);
  tb=Geometry::disjoint(ez,r);
  tb=Geometry::disjoint(ez,ez);
  tb=Geometry::subset(r,ez);
  tb=Geometry::subset(ez,r);
  tb=Geometry::subset(ez,ez);
  
  Geometry::subdivide(ez);
  Geometry::divide(ez);
  
  Geometry::over_approximation(ez);
  Geometry::approximation(ez);
}

template<> inline
void
instantiate_orthotope<Rational>()
{
  typedef Rational R;
  tribool* tb=0;
  Geometry::Point<R>* pt=0;
  Geometry::Rectangle<R>* r=0;
  Geometry::Orthotope<R>* z=0;
  Geometry::ListSet< Geometry::Orthotope<R> >* zls=0;
  
  *tb=Geometry::contains(*z,*pt);
  *tb=Geometry::disjoint(*r,*z);
  *tb=Geometry::disjoint(*z,*r);
  *tb=Geometry::disjoint(*z,*z);
  *tb=Geometry::subset(*r,*z);
  *tb=Geometry::subset(*z,*r);
  *tb=Geometry::subset(*z,*z);

  *z=Geometry::minkowski_sum(*z,*z);
  *z=Geometry::minkowski_difference(*z,*z);
  
  *zls=Geometry::subdivide(*z);
  *zls=Geometry::divide(*z);
  
  *z=Geometry::approximation(*z);
  *z=Geometry::over_approximation(*z);
}



} // namespace



