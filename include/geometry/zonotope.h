/***************************************************************************
 *            zonotope.h
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
 
/*! \file zonotope.h
 *  \brief Zonotopes (affine images of cuboids).
 */

#ifndef ARIADNE_ZONOTOPE_H
#define ARIADNE_ZONOTOPE_H

#include <iosfwd>

#include "base/iterator.h"
#include "base/tribool.h"

#include "numeric/declarations.h"

#include "linear_algebra/declarations.h"
#include "linear_algebra/matrix.h"


namespace Ariadne {
  namespace Geometry {
    
    class basic_set_tag;
    template<class X> class Point;
    template<class E> class RectangleExpression;
    template<class X> class Rectangle;
    template<class X> class Polyhedron;
    template<class BS> class ListSet;
  
    class ExactTag;    
    class UniformErrorTag;    
    class IntervalTag;    
    class DomainTag;    
    template<class R,class Type=ExactTag> class Zonotope;
  

    template<class S>
    class OverApproximation 
    {
     public:
      const S& _s;
      OverApproximation(const S& s) : _s(s) { }
      template<class R> operator R() const { R r; over_approximate(r,_s); return r; }
    };
  
    template<class S> inline
    OverApproximation<S> over_approximation(const S& s) {
      return OverApproximation<S>(s); 
    }
        
    template<class R, class S> inline
    R over_approximation(const S& s) {
      R r; over_approximate(r,s); return r;
    }
     

#ifdef DOXYGEN
    /*!\ingroup BasicSet
     * \brief A zonotope of arbitrary dimension.
     * 
     * A zonotope is a set of the form \f$c+Ge\f$, where \f$||e||_{\infty}\leq1\f$.
     * The columns of the matrix \f$G\f$ are the <em>generators</em> of the 
     * zonotope. 
     *
     * Zonotopes are always bounded.
     * A zonotope always contains its centre point, so can never be empty.
     * However, it may not be regular.
     *
     *
     * The intersection and membership tests may be performed using algorithms from: <br>
     * Guibas, Leonidas J.; Nguyen, An; Zhang, Li, "Zonotopes as bounding volumes."  
     * <i>Proceedings of the Fourteenth Annual ACM-SIAM Symposium on Discrete Algorithms</i> 
     * (Baltimore, MD, 2003),  803--812, ACM, New York, 2003.
     *
     * \b Storage: A %Zonotope in dimension d with n generators is described by
     * d(n+1) real numbers. The ith component of the centre is given by the
     * ith element, and the ith component of the kth generator is given by the
     * (kd+i)th element.
     */
    template<class R, class Tag>
    class Zonotope {
     public:
      //@{
      //! \name Typedefs 
      /*! \brief The type of real number used to describe the zonotope. */
      typedef R real_type;
      //@}

      //@{
      //! \name Constructors
      /*! \brief Default constructor yields a zonotope with dimension zero and no generators. */
      explicit Zonotope();
      /*! \brief Construct a zonotope of dimension \a d with no generators. */
      explicit Zonotope(dimension_type d);
      /*! \brief Construct a zonotope of dimension \a n with centre at the origin and \a m generators. */
      explicit Zonotope(dimension_type d, size_type m);
     
      /*! \brief Construct from centre and directions. */
      template<class XC, class XG> 
      explicit Zonotope(const Point<XC>& c, const LinearAlgebra::Matrix<XG>& g);
      /*! \brief Copy constructor. */
      Zonotope(const Zonotope<R,Tag>& z);

      /*! \brief Copy assignment operator. */
      Zonotope<R>& operator=(const Zonotope<R,Tag>& z);
      /*! \brief Assign from a box. */
      Zonotope<R,Tag>& operator=(const Rectangle<R>& r);
      //@}
      
      //@{ 
      //! \name Data access
      /*! \brief The dimension of the zonotope. */
      dimension_type dimension() const;

      /*! \brief The number of generators of the zonotope. */
      size_type number_of_generators() const;

      /*! \brief The domain. */
      Rectangle<R> domain() const;

      /*! \brief The centre. */
      Point<R> centre() const;

      /*! \brief The matrix of principle directions. */
      LinearAlgebra::Matrix<R> generators() const;
     
      //@}
      
      
      //@{
      //! \name Geometric binary predicates
      /*! \brief Tests disjointness of \a z and \a r. */
      friend tribool disjoint(const Zonotope<R>& z, const Rectangle<R>& r);
      /*! \brief Tests inclusion of \a r in \a z. */
      friend tribool superset(const Zonotope<R>& z, const Rectangle<R>& r);
      /*! \brief Tests inclusion of \a r in \a z. */
      friend tribool subset(const Rectangle<R>& z, const Zonotope<R>& r);
      /*! \brief Tests inclusion of \a z in \a r. */
      friend tribool subset(const Zonotope<R>& z, const Rectangle<R>& r);
      /*! \brief Tests inclusion of \a z in \a p. */
      friend tribool subset(const Zonotope<R>& z, const Polyhedron<R>& p);
      //@}
    };
#endif
  

    template<class R, class Tag>
    tribool empty(const Zonotope<R,Tag>& z);

    template<class R, class Tag>
    tribool bounded(const Zonotope<R,Tag>& z);

    template<class R, class Tag> 
    R radius(const Zonotope<R,Tag>&);
    
    template<class R, class Tag> 
    Rectangle<R> bounding_box(const Zonotope<R,Tag>& z);
    
    template<class R, class Tag> 
    tribool subset(const Rectangle<R>& z, const Zonotope<R,Tag>& r);
    
    template<class R, class Tag> 
    tribool disjoint(const Rectangle<R>& r, const Zonotope<R,Tag>& z);
    


    template<class R, class Tag> 
    tribool subset(const Zonotope<R,Tag>& z, const Rectangle<R>& r);
    
    template<class R, class Tag> 
    tribool subset(const Zonotope<R,Tag>& z, const Polyhedron<R>& p);



    template<class R> 
    tribool contains(const Zonotope<R,ExactTag>& z, const Point<R>& pt);

    template<class R> 
    tribool disjoint(const Zonotope<R,ExactTag>& z, const Rectangle<R>& r);

    template<class R> 
    tribool superset(const Zonotope<R,ExactTag>& r, const Rectangle<R>& z);

    template<class R, class ExactTag> 
    ListSet< Zonotope<R,ExactTag> >
    subdivide(const Zonotope<R,ExactTag>& z);


    template<class R> 
    tribool contains(const Zonotope<R,UniformErrorTag>& z, const Point<R>& pt);

    template<class R> 
    tribool disjoint(const Zonotope<R,UniformErrorTag>& z, const Rectangle<R>& r);

    template<class R> 
    tribool superset(const Zonotope<R,UniformErrorTag>& r, const Rectangle<R>& z);

    template<class R, class UniformErrorTag> 
    ListSet< Zonotope<R,UniformErrorTag> >
    subdivide(const Zonotope<R,UniformErrorTag>& z);



    

    template<class R, class Tag> 
    void over_approximate(Zonotope<R,Tag>&, const Rectangle<R>&);


    template<class R, class Tag> 
    void approximate(Zonotope<R,ExactTag>&, const Zonotope<R,Tag>&);


    template<class R> 
    void over_approximate(Zonotope<R,ExactTag>&, const Zonotope<R,IntervalTag>&);

    template<class R> 
    void over_approximate(Zonotope<R,ExactTag>&, const Zonotope<R,UniformErrorTag>&);

    template<class R> 
    void over_approximate(Zonotope<R,ExactTag>&, const Zonotope<R,ExactTag>&);

    template<class R> 
    void over_approximate(Zonotope<R,UniformErrorTag>&, const Zonotope<R,UniformErrorTag>&);

    template<class R> 
    void over_approximate(Zonotope<R,UniformErrorTag>&, const Zonotope<R,IntervalTag>&);


    template<class R> 
    void orthogonal_over_approximate(Zonotope<R,UniformErrorTag>&, const Zonotope<R,UniformErrorTag>&);
    
    template<class R> 
    void orthogonal_over_approximate(Zonotope<R,UniformErrorTag>&, const Zonotope<R,IntervalTag>&);
    
    template<class R> 
    void orthogonal_over_approximate(Zonotope<R,IntervalTag>&, const Zonotope<R,IntervalTag>&);

    
    template<class R> 
    void
    nonsingular_over_approximate(Zonotope<R,ExactTag>&, const Zonotope<R,IntervalTag>&);


        
    template<class R, class Tag> inline
    Zonotope<R> approximation(const Zonotope<R,Tag>& z) {
      Zonotope<R> r; approximate(r,z); return r;
    }
        
    template<class R, class Tag> inline
    Zonotope<R,Tag> orthogonal_over_approximation(const Zonotope<R,Tag>& s) {
      Zonotope<R,Tag> r; orthogonal_over_approximate(r,s); return r;
    }


    template<class R, class Tag>  
    std::ostream& operator<<(std::ostream& os, const Zonotope<R,Tag>& z);
    
    template<class R, class Tag> 
    std::istream& operator>>(std::istream& is, Zonotope<R,Tag>& z);

 







   
    

    template<class R>
    class Zonotope<R,ExactTag>
    {
      typedef Numeric::Interval<R> I;
     public:
      typedef basic_set_tag set_category;
      typedef R real_type;
      Zonotope(const dimension_type& d=0)
        : _centre(d), _generators(d,0) { }
      Zonotope(const Point<R>& c, const LinearAlgebra::Matrix<R>& G)
        : _centre(c), _generators(G) { }
      Zonotope(const Rectangle<R>& r);
      Zonotope(const Zonotope<R,UniformErrorTag>& z);
      template<class RR> Zonotope(const Zonotope<RR>& z);
      dimension_type dimension() const { 
        return _centre.dimension(); }
      size_type number_of_generators() const { 
        return _generators.number_of_columns(); }
      Rectangle<R> domain() const { 
        return Rectangle<R>::unit_box(this->number_of_generators()); }
      const Point<R>& centre() const { 
        return this->_centre; }
      const LinearAlgebra::Matrix<R>& generators() const { 
        return this->_generators; }
      R radius() const {
        return Geometry::radius(*this); }
      tribool empty() const {
        return false; }
      tribool bounded() const {
        return true; }
      tribool contains(const Point<R>& pt) const {
        return Geometry::contains(*this,pt); }
      Rectangle<R> bounding_box() const {
        return Geometry::bounding_box(*this); }
     private:
      static void _instantiate();
     private:
      Point<R> _centre;
      LinearAlgebra::Matrix<R> _generators;
    };

    template<class R>
    class Zonotope<R,UniformErrorTag>
    {
      typedef Numeric::Interval<R> I;
     public:
      typedef basic_set_tag set_category;
      typedef R real_type;
      Zonotope(const dimension_type& d=0)
        : _centre(d), _generators(d,0) { }
      Zonotope(const Point<I>& c, const LinearAlgebra::Matrix<R>& G)
        : _centre(c), _generators(G) { }
      Zonotope(const Zonotope<R,ExactTag>& z)
        : _centre(z.centre()), _generators(z.generators()) { }
      Zonotope(const Rectangle<R>& r);
      Zonotope<R,UniformErrorTag>& operator=(const Rectangle<R>& r);
      dimension_type dimension() const { 
        return _centre.dimension(); }
      size_type number_of_generators() const { 
        return _generators.number_of_columns(); }
      Rectangle<R> domain() const { 
        return Rectangle<R>::unit_box(this->number_of_generators()); }
      const Point<I>& centre() const { 
        return this->_centre; }
      const LinearAlgebra::Matrix<R>& generators() const { 
        return this->_generators; }
      R radius() const {
        return Geometry::radius(*this); }
      tribool empty() const {
        return false; }
      tribool bounded() const {
        return true; }
      tribool contains(const Point<R>& pt) const {
        return Geometry::contains(*this,pt); }
      Rectangle<R> bounding_box() const {
        return Geometry::bounding_box(*this); }
     private:
      static void _instantiate();
     private:
      Point<I> _centre;
      LinearAlgebra::Matrix<R> _generators;
    };

    template<class R>
    class Zonotope<R,IntervalTag>
    {
      typedef Numeric::Interval<R> I;
     public:
      typedef basic_set_tag set_category;
      typedef R real_type;
      Zonotope(const dimension_type& d=0)
        : _centre(d), _generators(d,0) { }
      Zonotope(const Point<I>& c, const LinearAlgebra::Matrix<I>& G)
        : _centre(c), _generators(G) { }
      Zonotope(const Zonotope<R,ExactTag>& z) 
        : _centre(z.centre()), _generators(z.generators()) { }
      Zonotope(const Zonotope<R,UniformErrorTag>& z) 
        : _centre(z.centre()), _generators(z.generators()) { }
      dimension_type dimension() const { 
        return _centre.dimension(); }
      size_type number_of_generators() const { 
        return _generators.number_of_columns(); }
      Rectangle<R> domain() const { 
        return Rectangle<R>::unit_box(this->number_of_generators()); }
      const Point<I>& centre() const { 
        return this->_centre; }
      const LinearAlgebra::Matrix<I>& generators() const { 
        return this->_generators; }
      R radius() const {
        return Geometry::radius(*this); }
      tribool empty() const {
        return false; }
      tribool bounded() const {
        return true; }
      Rectangle<R> bounding_box() const {
        return Geometry::bounding_box(*this); }
     private:
      Point<I> _centre;
      LinearAlgebra::Matrix<I> _generators;
    };


  }
}

#include "zonotope.inline.h"
#include "zonotope.template.h"

#endif /* ARIADNE_ZONOTOPE_H */
