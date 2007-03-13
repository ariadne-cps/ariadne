/***************************************************************************
 *            zonotope.code.h
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
 
#include <iostream>
#include <vector>
#include <algorithm>

#include "zonotope.h"

#include "../base/array.h"
#include "../exceptions.h"
#include "../numeric/conversion.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/lu_matrix.h"
#include "../linear_algebra/linear_program.h"

#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"
#include "../geometry/polytope.h"
#include "../geometry/polyhedron.h"
#include "../geometry/parallelotope.h"
#include "../geometry/list_set.h"

#include "../linear_algebra/vector.code.h"
#include "../linear_algebra/matrix.code.h"
#include "../linear_algebra/linear_program.code.h"
#include "../geometry/point.code.h"
#include "../geometry/rectangle.code.h"

namespace Ariadne {
  namespace Geometry {
    
    extern int verbosity;
    
    template<class R>
    tribool
    disjoint_approx(const Zonotope<R>& z, const Rectangle<R>& r);
   
    template<class R>
    tribool
    disjoint_exact(const Zonotope<R>& z, const Rectangle<R>& r);


    
    
    template<class R>
    Zonotope<R>::Zonotope(const std::string& s)
      : _dimension(), _number_of_generators(), _data()
    {
      std::stringstream ss(s);
      ss >> *this;
    }
    
     
    template<class R>
    typename Zonotope<R>::vertices_const_iterator
    Zonotope<R>::vertices_begin() const 
    {
      return ZonotopeVerticesIterator<R>(*this,false);
    }
    
    
    template<class R>
    typename Zonotope<R>::vertices_const_iterator
    Zonotope<R>::vertices_end() const 
    {
      return ZonotopeVerticesIterator<R>(*this,true);
    }
    
    
    template<class R>
    PointList<typename Numeric::traits<R>::arithmetic_type>
    Zonotope<R>::vertices() const
    {
      //std::cerr << "Zonotope<R>::vertices()" << std::endl;
      PointList<typename Numeric::traits<R>::arithmetic_type> v(this->dimension());
      for(typename Zonotope<R>::vertices_const_iterator vi=this->vertices_begin();
          vi!=this->vertices_end(); ++vi)
      {
        v.push_back(*vi);
      }
      return v;
    }
    
  
    
    template<class R>
    Zonotope<R>::operator Polytope<typename Numeric::traits<R>::arithmetic_type> () const
    {
      return Geometry::polytope(*this);
    }
    
    template<class R>
    Zonotope<R>::operator Polyhedron<typename Numeric::traits<R>::arithmetic_type> () const
    {
      return Geometry::polyhedron(*this);
    }
    
    
    template<class R>
    Polytope<typename Numeric::traits<R>::arithmetic_type> 
    polytope(const Zonotope<R>& z) 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      return Polytope<F>(z.vertices());
    }
  
    template<class R>
    Polyhedron<typename Numeric::traits<R>::arithmetic_type> 
    polyhedron(const Zonotope<R>& z) 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      return Polyhedron<F>(Polytope<F>(z.vertices()));
    }
  
    template<class R>
    Polyhedron<Rational> polyhedron(const Zonotope< Interval<R> >& z) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
  

    
    
    template<class R>
    Zonotope<typename Zonotope<R>::F> 
    Zonotope<R>::scale(const Zonotope<R>& z, const R& scale_factor) 
    {
      Geometry::Point<F> new_centre(z.centre().position_vector()*F(scale_factor));
      LinearAlgebra::Matrix<F> new_generators=z.generators()*F(scale_factor);
      return Geometry::Zonotope<F>(new_centre, new_generators);
    }
    
    
    template<class R>
    tribool 
    equal(const Zonotope<R>& A, const Zonotope<R>& B)
    {
      return subset(A,B) && subset(B,A);
    }
    
    
    template<class R> inline 
    Rectangle<R>
    bounding_box(const Zonotope<R>& z)
    {
      LinearAlgebra::Vector< Interval<R> > v(z.number_of_generators(),Interval<R>(-1,1));
      LinearAlgebra::Vector< Interval<R> > b=z.centre().position_vector()+z.generators()*v;
      return Rectangle<R>(b);
    }
    
    template<class R> inline
    Rectangle<R>
    bounding_box(const Zonotope< Interval<R> >& z)
    {
      LinearAlgebra::Vector< Interval<R> > e(z.number_of_generators(),Interval<R>(-1,1));
      LinearAlgebra::Vector< Interval<R> > b=z.centre().position_vector()+z.generators()*e;
      return Rectangle<R>(b);
    }

    template<class R>
    Rectangle<typename Numeric::traits<R>::number_type>
    Zonotope<R>::bounding_box() const
    {
      return Geometry::bounding_box(*this);
    }
    
    

    template<class R>
    ListSet< Zonotope<R> >
    subdivide(const Zonotope<R>& z) 
    {
      ListSet< Zonotope<R> > result(z.dimension());
      
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
        result.adjoin(Zonotope<R>(new_centre,new_generators));
      }
      return result;
    }
    
    
    template<class R>
    ListSet< Zonotope< Interval<R> > >
    subdivide(const Zonotope< Interval<R> >&) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    

    
    template<class R>
    ListSet< Zonotope<R> >
    divide(const Zonotope<R>& z)
    {
      size_type d=z.dimension();
      size_type m=z.number_of_generators();
      R two=2;
      ListSet< Zonotope<R> > result(d);
      
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
      result.adjoin(Zonotope<R>(new_centre,new_generators));
      new_centre=sub_approx(new_centre,LinearAlgebra::Vector<R>(new_generators.column(j)));
      result.adjoin(Zonotope<R>(new_centre,new_generators));

      return result;
    }
    
    
    template<class R>
    ListSet< Zonotope< Interval<R> > >
    divide(const Zonotope< Interval<R> >&)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    
    
    template<class R>
    void 
    Zonotope<R>::minimize_generators(void) 
    {
      return;
    }
    

    template<class R>
    inline
    tribool 
    norm_grtr(const LinearAlgebra::Vector<R>& v1, const LinearAlgebra::Vector<R>& v2) 
    {
      return LinearAlgebra::norm(v1)>LinearAlgebra::norm(v2);
    }
    
    
    template<class R>
    void 
    Zonotope<R>::sort_generators(void)
    {
      std::vector< LinearAlgebra::Vector<R> > generator_vectors;
      for(size_type j=0; j!=this->number_of_generators(); ++j) {
        generator_vectors.push_back(this->generators().column(j));
      }
      
      std::stable_sort(generator_vectors.begin(),generator_vectors.end(),&norm_grtr<R>);
      
      for(size_type j=0; j!=this->number_of_generators(); ++j) {
        for(size_type i=0; i!=this->dimension(); ++i) {
          this->_generators()(i,j)=generator_vectors[j](i);
        }
      }
    }
    
    
    /* Set up LP problem to solve \f$c+Ge=p\f$; \f$-1<=e<=1\f$.
     * Change variables so that the problem becomes \f$Ge=p-c-G1;\ 0\leq e\leq2\f$.
     * Change sign of \f$ Ge=p-c-G1\f$ to make right-hand side positive.
     */
    template<class R> inline
    tribool 
    contains(const Zonotope<R>& z, const Point<R>& pt)
    { 
      //std::cerr << "Zonotope<R>::contains(const Point<R>& )" << std::endl;
      //typedef typename Numeric::traits<R,R>::arithmetic_type F;
      typedef Rational F;
      check_equal_dimensions(z,pt,__PRETTY_FUNCTION__);
      dimension_type d=z.dimension();
      dimension_type m=z.number_of_generators();
      
      F zero=0;
      F one=1;
      F two=2;
      
      const Geometry::Point<F> qc=z.centre();
      const Geometry::Point<F> qp=pt;
      const LinearAlgebra::Matrix<F> qG=z.generators();
      const LinearAlgebra::Vector<F> qo(m,one);
      const LinearAlgebra::Vector<F> zv(m,zero);
      const LinearAlgebra::Vector<F> tv(m,two);

      LinearAlgebra::Vector<F> qrhs=qp-qc+qG*qo;
      
      LinearAlgebra::Matrix<F> T(d+m+1,m+1);
      
      // Set up constraints e+se=2
      for(dimension_type j=0; j!=m; ++j) {
        T(j,j)=1;
        T(j,m)=2;
      }
      
      // Set up constraints Ge \pm ax = p-c+G1
      for(dimension_type i=0; i!=d; ++i) {
        if(qrhs(i)>=F(0)) {
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
      
      LinearAlgebra::LinearProgram<F> lp(T);
      //std::cerr << lp.tableau() << std::endl;
      tribool result=lp.is_feasible();
      //std::cerr << lp.tableau() << std::endl;
      return result;
    }

    template<class R> inline
    tribool 
    contains(const Zonotope< Interval<R> >& z, const Point< Interval<R> >& pt)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool 
    Zonotope<R>::contains(const Point<R>& pt) const 
    {
      return Geometry::contains(*this,pt);
    }

  
    
    template<class R>
    tribool
    disjoint(const Zonotope<R>& z, const Rectangle<R>& r)
    {
      return disjoint_approx(z,r);
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
    template<class R>
    tribool
    disjoint_exact(const Zonotope<R>& z, const Rectangle<R>& r)
    {
      if(verbosity>7) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
      if(verbosity>8) { std::cerr << z << " " << r << std::endl; }
      
      check_equal_dimensions(z,r,__PRETTY_FUNCTION__);
      dimension_type d=z.dimension();
      size_type m=z.number_of_generators();
    
      //This method has some problems. It fails when the parameters are
      //[1,17/16]x[19/16,5/4] and 
      //  Zonotope(
      //    centre=(1/2, 1/10)
      //    directions=[ 1,1/2; 1/2,3/5 ]
      //  )
      
      // Construct tableau for testing intersection of zonotope and rectangle
      // Rectangle  l<=x<=u
      // Zonotope  x==c+Ge,  -1<=e<=1
      // 
      // Translate x'=x-l,  e'=e+1
      //   0<=x'<=u-l      ->  x' +     + sx'               == u-l
      //   0<=e'<=2        ->     +  e' +     + se'         == 2
      //   x'+l==c+G(e'-1) ->  x' + Ge'             +/- ax' == c-l-G1
      //  
      // Change sign of RHS of first equality if necessary
      // Introduce slack variables for last two inequalities
      typedef Rational F;
      LinearAlgebra::Matrix<F> T(2*d+m+1,d+m+1);

      const Geometry::Point<R>& l=r.lower_corner();
      const Geometry::Point<R>& u=r.upper_corner();
      const Geometry::Point<R> c=z.centre();
      const LinearAlgebra::Matrix<R> G=z.generators();

      const LinearAlgebra::Vector<F> qo(m,F(1));
      const LinearAlgebra::Vector<F> ql=l.position_vector();
      const LinearAlgebra::Vector<F> qu=u.position_vector();
      const LinearAlgebra::Vector<F> qd=qu-ql;
      const LinearAlgebra::Vector<F> qc=c.position_vector();
      const LinearAlgebra::Matrix<F> qG=G;
      const LinearAlgebra::Vector<F> qrhs=qc-ql-qG*qo;
      
      if(verbosity>8) { std::cerr << "ql=" << ql << ", qd=" << qd <<", qc=" << qc << ", qrhs=" << qrhs << std::endl; }
      
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
        if(qrhs(i)>=F(0)) {
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
      
      LinearAlgebra::LinearProgram<F> lp(T);
      tribool result=!lp.is_feasible();

      return result;
    }

    template<class R>
    tribool
    disjoint_approx(const Zonotope<R>& z, const Rectangle<R>& r)
    {
      //verbosity=9;
      
      if(verbosity>7) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
      if(verbosity>8) { std::cerr << z << " " << r << std::endl; }
      
      check_equal_dimensions(z,r,__PRETTY_FUNCTION__);
      dimension_type d=z.dimension();
      size_type m=z.number_of_generators();
    
      //This method has some problems. It fails when the parameters are
      //[1,17/16]x[19/16,5/4] and 
      //  Zonotope(
      //    centre=(1/2, 1/10)
      //    directions=[ 1,1/2; 1/2,3/5 ]
      //  )
      
      // Construct tableau for testing intersection of zonotope and rectangle
      // Rectangle  l<=x<=u
      // Zonotope  x==c+Ge,  -1<=e<=1
      // 
      // Translate x'=x-l,  e'=e+1
      //   0<=x'<=u-l      ->  x' +     + sx'               == u-l
      //   0<=e'<=2        ->     +  e' +     + se'         == 2
      //   x'+l==c+G(e'-1) ->  x' + Ge'             +/- ax' == c-l-G1
      //  
      // Change sign of RHS of first equality if necessary
      // Introduce slack variables for last two inequalities
      typedef typename Numeric::traits<R>::approximate_arithmetic_type B;
      B zero=0;
      B one=1;

      LinearAlgebra::Matrix<B> T(2*d+m+1,d+m+1);

      LinearAlgebra::Vector<B> l(d);
      LinearAlgebra::Vector<B> u(d);
      LinearAlgebra::Vector<B> c(d);
      LinearAlgebra::Matrix<B> G(d,m);
      for(dimension_type i=0; i!=d; ++i) {
        l(i)=r.lower_bound(i).get_base();
        u(i)=r.upper_bound(i).get_base();
        c(i)=z.centre(i).get_base();
        for(size_type j=0; j!=m; ++j) {
          G(i,j)=z.generators(i,j).get_base();
        }
      }

      if(verbosity>8) { std::cerr << l << u << c << G << std::endl; }

      const LinearAlgebra::Vector<B> o(m,one);
      const LinearAlgebra::Vector<B> rhs=c-l-G*o;

      if(verbosity>8) { std::cerr << "l=" << l << ", d=" << d <<", c=" << c << ", G=" << G << std::endl; }


      
      // Set up constraints x+sx=u-l
      for(size_type i=0; i!=d; ++i) {
        T(i,i)=1;
        T(i,d+m)=u(i)-l(i);
      }
      
      // Set up constraints e+se=2
      for(size_type j=0; j!=m; ++j) {
        T(d+j,d+j)=1;
        T(d+j,d+m)=2;
      }
      
      // Set up constraints x-Ge \pm ax=c-l-G1 
      for(size_type i=0; i!=d; ++i) {
        if(rhs(i)>=zero) {
          T(i+d+m,i)=1;
          for(size_type j=0; j!=m; ++j) {
            T(i+d+m,d+j)=-G(i,j);
          }
          T(i+d+m,d+m)=rhs(i);
        }
        else {
          T(i+d+m,i)=-1;
          for(size_type j=0; j!=m; ++j) {
            T(i+d+m,d+j)=G(i,j);
          }
          T(i+d+m,d+m)=-rhs(i);
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
      
      LinearAlgebra::LinearProgram<B> lp(T);
      tribool result=!lp.is_feasible();

      return result;
    }

   
    template<class R>
    tribool
    disjoint(const Rectangle<R>& r, const Zonotope<R>& z)
    {
      if(verbosity>7) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
      return disjoint(z,r);
    }
    
    
    template<class R>
    tribool
    disjoint(const Zonotope<R>& z1, const Zonotope<R>& z2)
    {
      typedef Rational F;
      check_equal_dimensions(z1,z2,__PRETTY_FUNCTION__);
      
      dimension_type d=z1.dimension();
      F one=1;
      size_type m1=z1.number_of_generators();
      size_type m2=z2.number_of_generators();
      
      LinearAlgebra::Matrix<F> T(m1+m2+d+1,m1+m2+1);

      const LinearAlgebra::Vector<F> qo1(m1,one);
      const LinearAlgebra::Vector<F> qo2(m2,one);

      Geometry::Point<F> qc1=z1.centre();
      LinearAlgebra::Matrix<F> qG1=z1.generators();
      Geometry::Point<F> qc2=z2.centre();
      LinearAlgebra::Matrix<F> qG2=z2.generators();
      LinearAlgebra::Vector<F> qrhs = qG1*qo1 - qG2*qo2 + (qc2 - qc1);
            
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
        if(qrhs(i)>=F(0)) {
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
      
      LinearAlgebra::LinearProgram<F> lp(T);
      
      tribool result=!lp.is_feasible();
      
      //std::cerr << "disjoint(" << z1 << "," << z2 << ")=" << result << std::endl;
      return result;
    }
    
    
    template<class R>
    tribool 
    subset(const Rectangle<R>& r, const Zonotope<R>& z) 
    {
      typedef typename Rectangle<R>::vertices_const_iterator RVIter;
      tribool result=true;
      for(RVIter rv_iter=r.vertices_begin(); rv_iter!=r.vertices_end(); ++rv_iter) {
        const Point<R>& pt=*rv_iter;
        result=result && z.contains(pt);
        if(result==false) {
          break;
        }
      }
      return result;
    }
    

    template<class R>
    tribool 
    subset(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Geometry::subset(A.bounding_box(),B);
    }
    

    template<class R>
    tribool 
    subset(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      return Geometry::subset(A.operator Polytope<F>(),B.operator Polyhedron<F>());
    } 
    
    
    template<class R> 
    Zonotope<typename Numeric::traits<R>::arithmetic_type>
    minkowski_sum(const Zonotope<R>& A, const Zonotope<R>& B)
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      
      Geometry::Point<F> new_centre=Geometry::minkowski_sum(A.centre(),B.centre());
      LinearAlgebra::Matrix<R> new_generators=LinearAlgebra::concatenate_columns(A.generators(),B.generators());
      return Zonotope<F>(new_centre,new_generators);
    }
     
    
    template<class R> 
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B)
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
     
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      
      return Zonotope<F>(Geometry::minkowski_difference(A.centre(),B.centre()),
                         LinearAlgebra::concatenate_columns(A.generators(),B.generators()));
    }
    
    
    template<class R> 
    Zonotope<R> 
    over_approximation(const Zonotope< Interval<R> >& iz)
    {
      if(iz.number_of_generators()==iz.dimension()) {
        Parallelotope< Interval<R> > ip=iz;
        Parallelotope<R> p=over_approximation(ip);
        return p;
      } else {
        return orthogonal_over_approximation(iz);
      }
    }
    
    
    template<class R> 
    Zonotope<R> 
    approximation(const Zonotope<R>& z)
    {
      return z; 
    }

    
    template<class R> 
    Zonotope<R> 
    approximation(const Zonotope< Interval<R> >& iz)
    {
      LinearAlgebra::Matrix< Interval<R> > G=iz.generators();

      return Zonotope<R>(approximate_value(iz.centre()),
                         approximate_value(G));
    }
    
    
    template<class R> 
    Zonotope< Interval<R> > 
    interval_over_approximation(const Zonotope< Interval<R> >& iz)
    {
      // FIXME: This is incorrect; need over-approximations
      LinearAlgebra::Matrix< Interval<R> > G(iz.generators());
      LinearAlgebra::Vector< Interval<R> > e(iz.number_of_generators(),Interval<R>(-1,+1));
      LinearAlgebra::Matrix<R> nG=approximate_value(G);
      Geometry::Point< Interval<R> > nc=iz.centre()+(G-nG)*e;
      
      return Zonotope<Interval<R> >(nc,nG);
    }

    template<class R>
    std::ostream&
    Zonotope<R>::write(std::ostream& os) const 
    {
      const Zonotope<R>& z=*this;
      if(z.dimension() > 0) {
        if (z.empty()) {
          os << "Zonotope( )" << std::endl;
        } else {
          os << "Zonotope( centre=" << z.centre();
          os << ", directions=" << z.generators();
          os << ")";
        } 
      }
      return os;
    }
    
    
    template<class R>
    std::istream& 
    Zonotope<R>::read(std::istream& is)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    
    template<class R>
    tribool
    disjoint(const Zonotope< Interval<R> >& z, const Rectangle<R>& r)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool
    disjoint(const Zonotope< Interval<R> >& z, const Rectangle< Interval<R> >& r)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool
    disjoint(const Zonotope< Interval<R> >& z, const Zonotope< Interval<R> >& r)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool
    subset(const Rectangle< R >& z, const Zonotope< Interval<R> >& r)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool
    subset(const Rectangle< Interval<R> >& z, const Zonotope< Interval<R> >& r)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool
    subset(const Zonotope< Interval<R> >& z, const Zonotope< Interval<R> >& r)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool
    subset(const Zonotope< Interval<R> >& z, const Rectangle<R>& r)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool
    subset(const Zonotope< Interval<R> >& z, const Rectangle< Interval<R> >& r)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }


    
       
    
    

    template<class R>
    std::ostream& 
    ZonotopeVerticesIterator<R>::write(std::ostream& os) const 
    {
      return os << "ZonotopeVerticesIterator<" << name<R>() << ">( &z=" << _z << ","
                << " i=" << _i << " p=" << _parity << ", v=" << **this << ")";
    }
  

    template<class R>
    void
    Zonotope<R>::_instantiate_geometry_operators() 
    {
      Rectangle<R> r;
      Zonotope<R> z;
      Zonotope<I> iz;
      Geometry::disjoint(r,z);
      Geometry::disjoint(z,r);
      Geometry::disjoint(z,z);
      Geometry::subset(r,z);
      Geometry::subset(z,r);
      Geometry::subset(z,z);
      Geometry::minkowski_sum(z,z);
      Geometry::minkowski_difference(z,z);
      
      Geometry::subset(iz,r);
      Geometry::disjoint(iz,r);

      typedef typename Numeric::traits<R>::arithmetic_type A;
      Zonotope<A> az;
      Geometry::over_approximation(az);
      Geometry::approximation(az);

      //Can't instantiate for Rational type
      //Geometry::interval_over_approximation(iz);

   }
    
    
  }
}
